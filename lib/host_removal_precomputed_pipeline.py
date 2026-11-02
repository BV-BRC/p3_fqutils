#!/usr/bin/env python3
"""
Host removal pipeline.

Required arguments:
  -i / --src-dir            source directory containing fastq files
  -s / --species-common-name   host common name for the reference genome
  -a / --aligner            aligner selection (bowtie2 for short reads; minimap2 for long reads)
  -o / --out-dir            output directory for filtered fastq files

This script will:
  1) search GCF accessions in a pre-computed index database, based on the aligner selection
     (human uses Hostile's built-in index, which Hostile downloads and caches itself)
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
import threading
from contextlib import contextmanager
from pathlib import Path

# Base folder for all reference data; override with HOST_REMOVAL_REF_DIR
REF_DIR = Path(os.environ.get("HOST_REMOVAL_REF_DIR", 
"/vol/bvbrc/production/application-backend/fastqutils/reference_indices"))

BOWTIE2_DB = Path(os.environ.get("HOST_REMOVAL_BOWTIE2_DB", REF_DIR / "bowtie2"))
MINIMAP2_DB = Path(os.environ.get("HOST_REMOVAL_MINIMAP2_DB", REF_DIR / "minimap2"))
HOST_GCF_DICT = Path(os.environ.get("HOST_REMOVAL_GCF_DICT",REF_DIR / "reference_index_dict_31Aug2026.csv"))

# Keep an existing HOSTILE_CACHE_DIR if one is set; otherwise use the default
os.environ.setdefault("HOSTILE_CACHE_DIR", str(REF_DIR / "hostile_cache"))

BOWTIE2_THREADS = 2

# Hostile's built-in standard indexes (keys are lower-case common names).
# Hostile downloads and caches these itself, for both Bowtie2 and Minimap2,
# so they bypass the GCF lookup and the local index database.
STANDARD_INDEXES = {
    "human": "human-t2t-hla"
    }

# Matches sample_1.fastq, sample_R1.fastq, and .fastq.gz variants of both
FWD_READ_PATTERN = re.compile(r"_R?1\.fastq(\.gz)?$")
# Matches sample_2.fastq, sample_R2.fastq, and .fastq.gz variants of both
RVS_READ_PATTERN = re.compile(r"_R?2\.fastq(\.gz)?$")
FASTQ_PATTERN = re.compile(r"\.fastq(\.gz)?$")
# Strips a trailing _1/_2/_R1/_R2 from a file stem to get the sample name
READ_NUMBER_SUFFIX = re.compile(r"_R?[12]$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the host removal pipeline.")
    parser.add_argument(
        "-i", "--src-dir", required=True, dest="src_dir",
        help="Source directory containing fastq files",
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
    Classifies each .fastq file in src_dir as forward (_1/_R1), reverse (_2/_R2),
    or single/unpaired (no read-number suffix), and derives the sample base name
    (read-number suffix stripped).
    If there are multiple samples, the last matching file in each category "wins".
    """
    fwd_reads = rvs_reads = single_reads = base_name = None

    for file in sorted(src_dir.glob("*.fastq")):
        name = file.name
        base_name = READ_NUMBER_SUFFIX.sub("", file.stem)

        if FWD_READ_PATTERN.search(name):
            fwd_reads = file
        elif RVS_READ_PATTERN.search(name):
            rvs_reads = file
        else:
            single_reads = file

    return fwd_reads, rvs_reads, single_reads, base_name


def find_precomputed_index(db_dir: Path, species_common_name: str, suffix: str = ""):
    """
    Returns the index to pass to `hostile clean --index`:
      - a Hostile standard index name for built-in hosts (e.g. human), or
      - the index prefix path (dir/dir_name + optional suffix) of the subdirectory
        under db_dir matching the species' GCF accession, or
      - None if no index is found.
    """

    # Built-in Hostile indexes (e.g. human) skip the GCF lookup entirely
    key = species_common_name.strip().lower()

    if key in STANDARD_INDEXES:
        index_name = STANDARD_INDEXES[key]
        print(f"'{species_common_name}' uses Hostile's built-in index '{index_name}'. Using this index to filter host sequences...")
        return index_name
    else:
        print(f"'{species_common_name}' is not a Hostile built-in index...")

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
        return None

    for match in sorted(db_dir.glob(f"{gcf_accession}*")):
        if match.is_dir():
            return match / f"{match.name}{suffix}"
    return None


def rename_hostile_outputs(out_dir: Path, input_reads: list, new_names: list) -> None:
    """Rename Hostile's <input_stem>.clean*.fastq.gz outputs to the given names."""
    for src, new_name in zip(input_reads, new_names):
        stem = FASTQ_PATTERN.sub("", src.name)
        matches = sorted(out_dir.glob(f"{stem}.clean*.fastq.gz"))
        if len(matches) != 1:
            raise FileNotFoundError(
                f"Expected one Hostile output for {src.name} in {out_dir}, found {len(matches)}"
            )
        target = out_dir / new_name
        matches[0].rename(target)
        print(f"Renamed {matches[0].name} -> {target.name}")


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

    precomputed_index = find_precomputed_index(BOWTIE2_DB, species_common_name)

    if precomputed_index is not None:

        cmd = ["hostile", "clean", "--force", "--fastq1", str(fwd_reads)]
        if paired:
            cmd += ["--fastq2", str(rvs_reads)]
        # --aligner bowtie2 is required for unpaired short reads; otherwise
        # Hostile treats a single input file as long reads (Minimap2)
        cmd += [
            "--aligner", "bowtie2",
            "--index", str(precomputed_index),
            "--out-dir", str(out_dir),
        ]
        run(cmd)

        if paired:
            rename_hostile_outputs(
                out_dir,
                [fwd_reads, rvs_reads],
                [f"cleaned_{base_name}_1.fastq.gz", f"cleaned_{base_name}_2.fastq.gz"],
            )
        else:
            rename_hostile_outputs(
                out_dir,
                [fwd_reads],
                [f"cleaned_{base_name}.fastq.gz"],
            )
        return


def run_minimap2(species_common_name: str, out_dir: Path, single_reads: Path, base_name: str) -> None:
    print("minimap2 was selected as the aligner; searching for the GCF Accession within precomputed index database...")

    precomputed_index = find_precomputed_index(MINIMAP2_DB, species_common_name, suffix=".mmi")

    if precomputed_index is not None:

        run([
            "hostile", "clean", "--force",
            "--fastq1", str(single_reads),
            "--aligner", "minimap2",
            "--index", str(precomputed_index),
            "--out-dir", str(out_dir),
        ])
        rename_hostile_outputs(
            out_dir,
            [single_reads],
            [f"cleaned_{base_name}.fastq.gz"],
        )
        return


def run_hostile(src_dir: str, species_common_name: str, aligner: str, out_dir: str) -> None:
    """
    Runs the full host-removal pipeline:
      1) classify fastq files in src_dir as forward/reverse/single-end
      2) search for a precomputed index matching gcf_accession, based on aligner
         (human uses Hostile's built-in index instead)
      3) if none exists, download the reference genome and build a custom index
      4) run `hostile clean` against the resolved index

    Can be called directly (e.g. imported into another script) or via the
    command-line entry point below.
    """
    src_dir = Path(src_dir)
    out_dir = Path(out_dir)
    log_path = out_dir / "hostile_log.txt"

    fwd_reads, rvs_reads, single_reads, base_name = find_fastq_files(src_dir)

    # A lone file (no matching mate, or no _1/_2 suffix) is treated as unpaired
    unpaired_reads = single_reads or fwd_reads or rvs_reads
    if unpaired_reads is None:
        raise FileNotFoundError(f"No .fastq files found in {src_dir}")

    with tee_to_file(log_path):
        if aligner == "bowtie2":
            if fwd_reads is not None and rvs_reads is not None:
                run_bowtie2(species_common_name, out_dir, fwd_reads, rvs_reads, base_name)
            else:
                run_bowtie2(species_common_name, out_dir, unpaired_reads, None, base_name)
        elif aligner == "minimap2":
            run_minimap2(species_common_name, out_dir, unpaired_reads, base_name)
        else:
            print(f"Error: unrecognized aligner '{aligner}'", file=sys.stderr)
            sys.exit(1)

    print("Host removal: complete")


def main() -> None:
    args = parse_args()
    run_hostile(
        src_dir=args.src_dir,
        species_common_name=args.species_common_name,
        aligner=args.aligner,
        out_dir=args.out_dir,
    )


if __name__ == "__main__":
    main()