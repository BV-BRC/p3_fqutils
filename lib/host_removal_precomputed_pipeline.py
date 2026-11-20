#!/usr/bin/env python3
"""
Host removal pipeline.

Input (one of):
  -1 / --read1 [-2 / --read2]  exact read file(s): read1 alone = single-end/long reads,
                               read1 + read2 = paired-end
  -i / --src-dir               source directory to search for fastq files (legacy)

Required arguments:
  -s / --species-common-name   host common name for the reference genome
  -a / --aligner            aligner selection (bowtie2 for short reads; minimap2 for long reads)
  -o / --out-dir            output directory for filtered fastq files

This script will:
  1) search GCF accessions in a pre-computed index database, based on the aligner selection
     (human uses Hostile's built-in index, which Hostile downloads and caches itself;
     "Human - Complete" runs Hostile's human index first, then GRCh38 on the cleaned reads)
  2) remotely download reference genome and construct a custom index (if needed)
  3) filter host reads from input fastq files and output filtered fastq file(s)
"""

import argparse
import re
import subprocess
import sys
import zipfile
import csv
import os
import shutil
import threading
from contextlib import contextmanager
from pathlib import Path

# Base folder for all reference data; override with HOST_REMOVAL_REF_DIR
REF_DIR = Path(os.environ.get("HOST_REMOVAL_REF_DIR", 
"/vol/bvbrc/production/application-backend/fastqutils/reference_indices"))

BOWTIE2_DB = Path(os.environ.get("HOST_REMOVAL_BOWTIE2_DB", REF_DIR / "bowtie2"))
MINIMAP2_DB = Path(os.environ.get("HOST_REMOVAL_MINIMAP2_DB", REF_DIR / "minimap2"))
HOST_GCF_DICT = Path(os.environ.get("HOST_REMOVAL_GCF_DICT",REF_DIR / "reference_index_dict.csv"))

# Keep an existing HOSTILE_CACHE_DIR if one is set; otherwise use the default
os.environ.setdefault("HOSTILE_CACHE_DIR", str(REF_DIR / "hostile_cache"))

BOWTIE2_THREADS = 8

# Hostile's built-in standard indexes (keys are lower-case common names).
# Hostile downloads and caches these itself, for both Bowtie2 and Minimap2,
# so they bypass the GCF lookup and the local index database.
STANDARD_INDICES = {
    "human": "human-t2t-hla",
    "human - default": "human-t2t-hla"
    }

# Two-pass hosts (keys are lower-case common names): [Hostile standard index, GCF accession].
# Reads are first cleaned against the standard index, then the cleaned reads are
# cleaned again against the precomputed index for the GCF accession (GRCh38 here).
COMPLETE_INDICES = {
    "human - complete": ["human-t2t-hla", "GCF_000001405.40"]
}

# Optional suffixes added by fastq_utils' trimming step / Trim Galore
_TRIM_SUFFIX = r"(?:_ptrim|_strim|_val_[12]|_trimmed)?"
# FASTQ extensions: .fastq, .fq, .fastq.gz or .fq.gz
_FASTQ_EXT = r"\.f(?:ast)?q(?:\.gz)?$"

# Matches sample_1 / sample_R1 (optionally trimmed) with any FASTQ extension
FWD_READ_PATTERN = re.compile(r"_R?1" + _TRIM_SUFFIX + _FASTQ_EXT)
# Matches sample_2 / sample_R2 (optionally trimmed) with any FASTQ extension
RVS_READ_PATTERN = re.compile(r"_R?2" + _TRIM_SUFFIX + _FASTQ_EXT)
# Any FASTQ file
FASTQ_PATTERN = re.compile(_FASTQ_EXT)
# Trailing _1/_2/_R1/_R2, optionally followed by a trim suffix (captured)
READ_NUMBER_SUFFIX = re.compile(r"_R?[12](_ptrim|_strim|_val_[12]|_trimmed)?$")


def sample_name(read_path: Path) -> str:
    """
    Sample base name: file name without FASTQ extension or read number.
    A trim suffix is kept, so trimmed and untrimmed outputs stay distinguishable:
      sample_R1.fastq.gz       -> sample
      sample_R1_ptrim.fq.gz    -> sample_ptrim
      sample_R1_val_1.fq.gz    -> sample_val
      sample_strim.fq.gz       -> sample_strim   (single-end, no read number)
    """
    stem = FASTQ_PATTERN.sub("", Path(read_path).name)
    match = READ_NUMBER_SUFFIX.search(stem)
    if not match:
        return stem
    trim_suffix = match.group(1) or ""
    if trim_suffix.startswith("_val_"):
        trim_suffix = "_val"  # drop the read number Trim Galore puts in _val_1/_val_2
    return stem[:match.start()] + trim_suffix


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the host removal pipeline.")
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "-1", "--read1", dest="read1",
        help="Read file 1 (or the only file for single-end/long reads)",
    )
    inputs.add_argument(
        "-i", "--src-dir", dest="src_dir",
        help="Source directory to search for fastq files (legacy)",
    )
    parser.add_argument(
        "-2", "--read2", dest="read2",
        help="Read file 2 for paired-end reads (use with --read1)",
    )
    parser.add_argument(
        "-s", "--species-common-name", required=True, dest="species_common_name",
        help="Host common name for the reference genome (e.g. Human, Cow)",
    )
    parser.add_argument(
        "-a", "--aligner", required=True, choices=["bowtie2", "minimap2"],
        help="Aligner selection (bowtie2 for short reads; minimap2 for long reads)",
    )
    parser.add_argument(
        "-o", "--out-dir", required=True, dest="out_dir",
        help="Output directory for filtered fastq files",
    )
    return parser.parse_args()


def run(cmd: list, **kwargs) -> None:
    """Run a command, echoing it first, and raise on non-zero exit (mirrors `set -e`)."""
    print("+ " + " ".join(str(part) for part in cmd))
    subprocess.run(cmd, check=True, **kwargs)


def find_fastq_files(src_dir: Path):
    """
    Classifies each FASTQ file in src_dir (.fastq, .fq, .fastq.gz, .fq.gz) as
    forward (_1/_R1), reverse (_2/_R2), or single/unpaired (no read-number suffix),
    and derives the sample base name (extension and read-number suffix stripped).
    If there are multiple samples, the last matching file in each category "wins".
    """
    fwd_reads = rvs_reads = single_reads = base_name = None

    for file in sorted(src_dir.iterdir()):
        if not file.is_file() or not FASTQ_PATTERN.search(file.name):
            continue
        # Skip host-removal outputs (ours and Hostile's) left by earlier runs,
        # since src_dir can be the same folder as out_dir
        if file.name.startswith("cleaned_") or ".clean" in file.name:
            continue
        name = file.name
        base_name = sample_name(file)

        if FWD_READ_PATTERN.search(name):
            fwd_reads = file
        elif RVS_READ_PATTERN.search(name):
            rvs_reads = file
        else:
            single_reads = file

    return fwd_reads, rvs_reads, single_reads, base_name


def find_precomputed_index(db_dir: Path, species_common_name: str, suffix: str = ""):
    """
    Returns the list of indexes to pass to `hostile clean --index`, in the order
    they should be run, or None if no index is found:
      - [standard index name]            for built-in Hostile hosts (e.g. Human)
      - [standard name, precomputed path] for two-pass hosts (e.g. Human - Complete:
                                          Hostile's human index, then GRCh38)
      - [precomputed path]               for all other hosts: the index prefix
                                          (dir/dir_name + optional suffix) of the
                                          subdirectory of db_dir matching the GCF accession
    """
    key = species_common_name.strip().lower()
    # print("Current key is:", key)

    # Two-pass hosts: Hostile's standard index first, then the precomputed GCF index
    if key in COMPLETE_INDICES:
        standard_name, gcf_accession = COMPLETE_INDICES[key]
        print(f"'{species_common_name}' selected for host removal. This pipeline will filter using "
              f"Hostile's '{standard_name}' index first, followed by '{gcf_accession}'.")
        precomputed = find_index_dir(db_dir, gcf_accession, suffix)
        if precomputed is None:
            return None
        return [standard_name, precomputed]

    # Built-in Hostile indexes (e.g. human) skip the GCF lookup entirely
    if key in STANDARD_INDICES:
        index_name = STANDARD_INDICES[key]
        print(f"'{species_common_name}' uses Hostile's built-in index '{index_name}'. Using this index to filter host sequences...")
        return [index_name]

    print(f"'{species_common_name}' is not a Hostile built-in index...")
    gcf_accession = lookup_gcf_accession(species_common_name)
    if gcf_accession is None:
        return None
    precomputed = find_index_dir(db_dir, gcf_accession, suffix)
    return [precomputed] if precomputed is not None else None


def lookup_gcf_accession(species_common_name: str):
    """Returns the species' GCF accession from HOST_GCF_DICT, or None."""
    host_gcf_dict = {}

    with open(HOST_GCF_DICT, mode='r', encoding='utf-8') as file:
        csv_reader = csv.DictReader(file)

        for row in csv_reader:
            # Use the 'Common_Name' column as the main key
            key = row.pop('Common_Name')
            # The remaining row data becomes the value
            host_gcf_dict[key] = row

    gcf_accession = None
    for key, row in host_gcf_dict.items():
        if species_common_name.strip().lower() == key.strip().lower():
            gcf_accession = row['GCF_Accession']  # confirm this matches your actual CSV header exactly
            print(f"'{species_common_name}' found within reference species list. Paired GCF Accession is '{gcf_accession}'. Using this reference index to filter host sequences...")
            break

    if gcf_accession is None:
        print(f"'{species_common_name}' not found in reference species list.")
    return gcf_accession


def find_index_dir(db_dir: Path, gcf_accession: str, suffix: str = ""):
    """
    Returns the index prefix path (dir/dir_name + optional suffix) of the
    subdirectory under db_dir whose name starts with gcf_accession, or None.
    """
    for match in sorted(db_dir.glob(f"{gcf_accession}*")):
        if match.is_dir():
            return match / f"{match.name}{suffix}"
    print(f"No index directory for '{gcf_accession}' found in {db_dir}.")
    return None


def find_hostile_outputs(out_dir: Path, input_reads: list) -> list:
    """Return Hostile's <input_stem>.clean*.fastq.gz output for each input read file."""
    outputs = []
    for src in input_reads:
        stem = FASTQ_PATTERN.sub("", src.name)
        matches = sorted(out_dir.glob(f"{stem}.clean*.fastq.gz"))
        if len(matches) != 1:
            raise FileNotFoundError(
                f"Expected one Hostile output for {src.name} in {out_dir}, found {len(matches)}"
            )
        outputs.append(matches[0])
    return outputs


def hostile_clean(aligner: str, index, input_reads: list, out_dir: Path) -> list:
    """
    Runs one `hostile clean` pass against a single index and returns the cleaned
    read files (same order as input_reads: [fwd, rvs] for paired, [reads] otherwise).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["hostile", "clean", "--force", "--fastq1", str(input_reads[0])]
    if len(input_reads) == 2:
        cmd += ["--fastq2", str(input_reads[1])]
    # --aligner is always given: for bowtie2 it is required for unpaired short
    # reads, otherwise Hostile treats a single input file as long reads (Minimap2)
    cmd += [
        "--aligner", aligner,
        "--index", str(index),
        "--out-dir", str(out_dir),
    ]
    run(cmd)
    return find_hostile_outputs(out_dir, input_reads)


def run_hostile_passes(aligner: str, indexes: list, input_reads: list,
                       out_dir: Path, final_names: list) -> None:
    """
    Cleans input_reads against each index in turn, feeding each pass's cleaned
    reads into the next, and saves the last pass's output in out_dir as final_names.
    A single index runs straight into out_dir. With several indexes, each pass
    writes to its own temporary subfolder, removed once the final files are saved.
    """
    current = list(input_reads)
    pass_dirs = []
    for i, index in enumerate(indexes, start=1):
        if len(indexes) == 1:
            pass_dir = out_dir
        else:
            pass_dir = out_dir / f"hostile_pass{i}"
            pass_dirs.append(pass_dir)
            print(f"Host removal pass {i} of {len(indexes)}: index {index}")
        current = hostile_clean(aligner, index, current, pass_dir)

    for src, new_name in zip(current, final_names):
        target = out_dir / new_name
        src.rename(target)
        print(f"Renamed {src.name} -> {target.name}")

    for pass_dir in pass_dirs:
        shutil.rmtree(pass_dir, ignore_errors=True)


@contextmanager
def tee_to_file(path):
    """Show stdout/stderr on screen (including subprocess output) and save it to a file."""
    sys.stdout.flush()
    sys.stderr.flush()
    saved_out, saved_err = os.dup(1), os.dup(2)
    read_fd, write_fd = os.pipe()
    log_file = open(path, "wb")

    def pump():
        while True:
            chunk = os.read(read_fd, 4096)
            if not chunk:
                break
            os.write(saved_out, chunk)   # to the screen
            log_file.write(chunk)        # to the file
            log_file.flush()

    thread = threading.Thread(target=pump, daemon=True)
    thread.start()

    os.dup2(write_fd, 1)
    os.dup2(write_fd, 2)
    os.close(write_fd)
    try:
        yield
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(saved_out, 1)
        os.dup2(saved_err, 2)
        thread.join()
        log_file.close()
        os.close(read_fd)
        os.close(saved_out)
        os.close(saved_err)


def run_bowtie2(species_common_name: str, out_dir: Path, fwd_reads: Path, rvs_reads, base_name: str) -> None:
    """
    Runs Hostile with Bowtie2 on short reads.
    Paired mode when rvs_reads is given; unpaired (single-end) mode when rvs_reads is None.
    """
    paired = rvs_reads is not None
    mode = "paired" if paired else "unpaired"
    print(f"bowtie2 was selected as the aligner ({mode} short reads); searching for the host species within precomputed index database...")

    indexes = find_precomputed_index(BOWTIE2_DB, species_common_name)

    if indexes is not None:
        if paired:
            run_hostile_passes(
                "bowtie2", indexes, [fwd_reads, rvs_reads], out_dir,
                [f"cleaned_{base_name}_1.fastq.gz", f"cleaned_{base_name}_2.fastq.gz"],
            )
        else:
            run_hostile_passes(
                "bowtie2", indexes, [fwd_reads], out_dir,
                [f"cleaned_{base_name}.fastq.gz"],
            )
        return

    raise FileNotFoundError(
        f"No bowtie2 index found for '{species_common_name}' in {BOWTIE2_DB}; host removal not run"
    )


def run_minimap2(species_common_name: str, out_dir: Path, single_reads: Path, base_name: str) -> None:
    print("minimap2 was selected as the aligner; searching for the GCF Accession within precomputed index database...")

    indexes = find_precomputed_index(MINIMAP2_DB, species_common_name, suffix=".mmi")

    if indexes is not None:
        run_hostile_passes(
            "minimap2", indexes, [single_reads], out_dir,
            [f"cleaned_{base_name}.fastq.gz"],
        )
        return

    raise FileNotFoundError(
        f"No minimap2 index found for '{species_common_name}' in {MINIMAP2_DB}; host removal not run"
    )


def run_hostile(species_common_name: str, aligner: str, out_dir: str,
                read1: str = None, read2: str = None, src_dir: str = None) -> None:
    """
    Runs the full host-removal pipeline:
      1) use the given read1 (+ read2) files; or, if only src_dir is given,
         classify the fastq files in src_dir as forward/reverse/single-end
      2) search for a precomputed index matching gcf_accession, based on aligner
         (human uses Hostile's built-in index instead)
      3) if none exists, download the reference genome and build a custom index
      4) run `hostile clean` against the resolved index

    Can be called directly (e.g. imported into another script) or via the
    command-line entry point below.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "hostile_log.txt"

    if read1:
        # Exact files given: read1 + read2 = paired, read1 alone = single-end/long reads
        for path in (read1, read2):
            if path and not Path(path).is_file():
                raise FileNotFoundError(f"Read file not found: {path}")
        if read2:
            fwd_reads, rvs_reads, single_reads = Path(read1), Path(read2), None
        else:
            fwd_reads, rvs_reads, single_reads = None, None, Path(read1)
        base_name = sample_name(read1)
    elif src_dir:
        src_dir = Path(src_dir)
        fwd_reads, rvs_reads, single_reads, base_name = find_fastq_files(src_dir)
    else:
        raise ValueError("run_hostile needs read1 (and optionally read2) or src_dir")

    # A lone file (no matching mate, or no _1/_2 suffix) is treated as unpaired
    unpaired_reads = single_reads or fwd_reads or rvs_reads
    if unpaired_reads is None:
        raise FileNotFoundError(f"No FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz) found in {src_dir}")

    with tee_to_file(log_path):
        if aligner == "bowtie2":
            if fwd_reads is not None and rvs_reads is not None:
                run_bowtie2(species_common_name, out_dir, fwd_reads, rvs_reads, base_name)
            else:
                run_bowtie2(species_common_name, out_dir, unpaired_reads, None, base_name)
        elif aligner == "minimap2":
            run_minimap2(species_common_name, out_dir, unpaired_reads, base_name)
        else:
            raise ValueError(f"Unrecognized aligner '{aligner}' (expected bowtie2 or minimap2)")

    print("Host removal: complete")


def main() -> None:
    args = parse_args()
    if args.read2 and not args.read1:
        sys.exit("--read2 requires --read1")
    run_hostile(
        read1=args.read1,
        read2=args.read2,
        src_dir=args.src_dir,
        species_common_name=args.species_common_name,
        aligner=args.aligner,
        out_dir=args.out_dir,
    )


if __name__ == "__main__":
    main()