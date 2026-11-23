#!/usr/bin/env python3
"""
reference_index_pipeline.py

One pipeline that keeps the fastqutils reference genomes and their bowtie2,
minimap2 and HISAT2 indices in step with dict.csv.

Stages (run in this order):
  1. UPDATE  - Check every GCF accession in dict.csv against NCBI. Outdated
               accessions are replaced in dict.csv (original kept as .bak) and
               two lists are written, named with the count and today's date:
                   <N>_GCF_Accessions_Updated_<DDMonthYYYY>.txt   (new accessions)
                   <N>_GCF_Accessions_Removed_<DDMonthYYYY>.txt   (replaced accessions)
               e.g. 10_GCF_Accessions_Updated_21November2026.txt
  2. BUILD   - Download the new genomes with NCBI `datasets`, move the FASTA
               files into reference_genomes/, build bowtie2, minimap2 and
               HISAT2 indices (with timing logs), then move each FASTA to
               reference_genomes/done/.
  3. REMOVE  - Delete the bowtie2/minimap2/hisat2 index directories and the
               done/ FASTA files of the replaced accessions. An old accession
               is only removed once its replacement has all three indices, and
               never while it is still listed in dict.csv.
  4. VERIFY  - For every accession in dict.csv, confirm a complete bowtie2,
               minimap2 and HISAT2 index exists; build whatever is missing
               (downloading the genome again if its FASTA is not on disk).
  FINAL      - Write host_summary_manifesto.json (in reference_indices/), in the
               same format as BV-BRC's patric_host_summary.json, with
               "patric_ftp" replaced by "url_path": the path of the genome's
               .fna under reference_genomes/. Metadata comes from NCBI's RefSeq
               assembly_summary files (Datasets API as a fallback).

An accession counts as outdated when:
  1. A newer version of the same assembly exists
     (e.g. GCF_000001635.26 -> GCF_000001635.27), or
  2. The assembly has been suppressed/retired with no newer version. It is
     then replaced by the species' current RefSeq reference genome.
  3. (Only with --prefer-reference) The species' current RefSeq reference
     genome is a different assembly from the one in the CSV.

The CSV defaults to
    /vol/bvbrc/production/application-backend/fastqutils/reference_indices/reference_index_dict.csv
and every indexer runs with 8 threads unless --threads is given.

Usage:
    python reference_update_index_pipeline.py                     # default CSV, 8 threads
    python reference_update_index_pipeline.py --dry-run           # report only
    python reference_update_index_pipeline.py --yes               # no delete prompt
    python reference_update_index_pipeline.py --threads 16
    python reference_update_index_pipeline.py --skip-update       # build/verify only
    python reference_update_index_pipeline.py other.csv           # use a different CSV
    NCBI_API_KEY=xxxx python reference_update_index_pipeline.py   # faster rate limit
"""

import argparse
import csv
import datetime
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
DEFAULT_BASE = "/vol/bvbrc/production/application-backend/fastqutils"
DEFAULT_CSV = f"{DEFAULT_BASE}/reference_indices/reference_index_dict.csv"
DEFAULT_THREADS = 8
TOOLS = ("bowtie2", "minimap2", "hisat2")

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ACC_RE = re.compile(r"^(GC[FA]_\d{9})(?:\.(\d+))?$")
API_KEY = os.environ.get("NCBI_API_KEY")
# NCBI allows 5 requests/s without a key, 10/s with one.
MIN_INTERVAL = 0.11 if API_KEY else 0.22
_last_call = 0.0


class Paths:
    def __init__(self, base, datasets=None):
        base = Path(base)
        self.genomes = base / "reference_genomes"
        self.done = self.genomes / "done"
        self.index_base = base / "reference_indices"
        self.datasets = Path(datasets) if datasets else self.index_base / "datasets"

    def tool_root(self, tool):
        return self.index_base / tool

    def timing_log(self, tool):
        return self.index_base / f"{tool}indexing_timing_log.csv"


class PipelineError(Exception):
    pass


def log(msg=""):
    print(msg, flush=True)


def banner(title):
    log("\n" + "=" * 78)
    log(f"  {title}   [{datetime.datetime.now():%Y-%m-%d %H:%M:%S}]")
    log("=" * 78)


def date_tag(today=None):
    """21November2026 style date used in the list file names."""
    d = today or datetime.date.today()
    return f"{d.day}{d:%B%Y}"


# --------------------------------------------------------------------------- #
# NCBI helpers
# --------------------------------------------------------------------------- #
def ncbi_get(path, params=None, retries=4):
    """GET a Datasets API endpoint and return parsed JSON (or None on 404)."""
    global _last_call
    url = f"{API}/{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params, doseq=True)
    headers = {"Accept": "application/json"}
    if API_KEY:
        headers["api-key"] = API_KEY

    for attempt in range(retries):
        wait = MIN_INTERVAL - (time.time() - _last_call)
        if wait > 0:
            time.sleep(wait)
        _last_call = time.time()
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if e.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise
        except urllib.error.URLError:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise


def split_acc(acc):
    m = ACC_RE.match(acc.strip())
    if not m:
        return None, None
    return m.group(1), int(m.group(2)) if m.group(2) else None


def assembly_report(acc):
    """Dataset report for one accession, including non-current versions."""
    data = ncbi_get(
        f"genome/accession/{urllib.parse.quote(acc)}/dataset_report",
        {"filters.assembly_version": "all_assemblies"},
    )
    reports = (data or {}).get("reports", [])
    for rep in reports:
        if rep.get("accession") == acc:
            return rep
    return reports[0] if reports else None


def latest_version(acc):
    """Return (newest version of the same assembly or None, report)."""
    base, _ = split_acc(acc)
    candidates = set()

    rep = assembly_report(acc)
    if rep:
        for key in ("current_accession", "accession", "paired_accession"):
            val = rep.get(key)
            if val and val.startswith(base):
                candidates.add(val)

    hist = ncbi_get(f"genome/accession/{urllib.parse.quote(acc)}/revision_history")
    for rev in (hist or {}).get("assembly_revisions", []):
        for val in rev.values():
            if isinstance(val, str) and val.startswith(base):
                candidates.add(val)

    versions = [(split_acc(c)[1] or 0, c) for c in candidates if split_acc(c)[0] == base]
    return (max(versions)[1] if versions else None), rep


def species_reference(taxon):
    """Current RefSeq reference genome accession for a taxon (name or tax ID)."""
    data = ncbi_get(
        f"genome/taxon/{urllib.parse.quote(str(taxon))}/dataset_report",
        {
            "filters.reference_only": "true",
            "filters.assembly_source": "refseq",
            "filters.assembly_version": "current",
            "tax_exact_match": "true",
            "page_size": 5,
        },
    )
    for rep in (data or {}).get("reports", []):
        acc = rep.get("accession", "")
        if acc.startswith("GCF_"):
            return acc
    return None


def check_accession(acc, species, prefer_reference=False):
    """Return (new_accession_or_None, reason)."""
    base, version = split_acc(acc)
    if not base:
        return None, "not a valid GCF/GCA accession - skipped"

    newest, rep = latest_version(acc)
    status = ((rep or {}).get("assembly_info") or {}).get("assembly_status", "")
    taxid = ((rep or {}).get("organism") or {}).get("tax_id")

    # 1. Newer version of the same assembly
    if newest and newest != acc and (version is None or split_acc(newest)[1] > version):
        return newest, f"newer version of same assembly ({acc} -> {newest})"

    is_current = rep is not None and status in ("", "current")

    # 2. Suppressed / retired with nothing newer -> species reference genome
    # 3. Optionally, swap to the species' reference genome even if still current
    if not is_current or prefer_reference:
        ref = species_reference(taxid or species) if (taxid or species) else None
        if ref and split_acc(ref)[0] != base:
            why = (f"assembly status '{status or 'not found'}'" if not is_current
                   else "species reference genome is a different assembly")
            return ref, f"{why}; using current RefSeq reference {ref}"
        if not is_current:
            return None, (f"assembly status '{status or 'not found'}' and no "
                          f"RefSeq reference genome found - review manually")

    return None, "up to date"


# --------------------------------------------------------------------------- #
# CSV helpers
# --------------------------------------------------------------------------- #
def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(fh, dialect=dialect)
        return reader.fieldnames, list(reader), dialect


def guess_columns(header, rows):
    """Find the accession and species columns if not given explicitly."""
    acc_col = species_col = None
    for col in header:
        vals = [r[col] for r in rows if r.get(col)]
        if vals and sum(bool(ACC_RE.match(v.strip())) for v in vals) / len(vals) > 0.5:
            acc_col = col
            break
    for col in header:
        if col != acc_col and re.search(r"species|organism|name|taxon", col, re.I):
            species_col = col
            break
    return acc_col, species_col


def resolve_columns(header, rows, acc_col_arg, species_col_arg):
    acc_col, species_col = guess_columns(header, rows)
    acc_col = acc_col_arg or acc_col
    species_col = species_col_arg or species_col
    if not acc_col or acc_col not in header:
        raise PipelineError(f"Could not find the accession column. Columns: {header}. "
                            f"Use --acc-col.")
    return acc_col, species_col


def csv_accessions(rows, acc_col):
    """Unique, valid accessions from the CSV rows, in file order."""
    accs = [(r.get(acc_col) or "").strip() for r in rows]
    return list(dict.fromkeys(a for a in accs if split_acc(a)[0]))


# --------------------------------------------------------------------------- #
# Stage 1: update accessions
# --------------------------------------------------------------------------- #
def stage_update(args, header, rows, dialect, acc_col, species_col):
    """Check/replace accessions. Returns list of (old, new) pairs."""
    banner("STAGE 1/4: Checking dict.csv accessions against NCBI")
    log(f"Accession column: {acc_col!r}   Species column: {species_col!r}\n")

    pairs, problems = [], 0
    for i, row in enumerate(rows, 1):
        acc = (row.get(acc_col) or "").strip()
        species = (row.get(species_col) or "").strip() if species_col else ""
        if not acc:
            continue
        try:
            new_acc, reason = check_accession(acc, species, args.prefer_reference)
        except Exception as e:  # network or API error - keep going
            new_acc, reason = None, f"ERROR: {e}"
        label = species or f"row {i}"
        flag = "UPDATED" if new_acc else ("OK" if reason == "up to date" else "CHECK")
        log(f"[{flag:7}] {label}: {acc} - {reason}")
        if new_acc:
            row[acc_col] = new_acc
            pairs.append((acc, new_acc))
        elif flag == "CHECK":
            problems += 1

    log(f"\n{len(rows)} rows checked, {len(pairs)} updated, {problems} need review.")

    if not pairs:
        return pairs
    if args.dry_run:
        log("Dry run - dict.csv and accession lists not written.")
        return pairs

    backup = args.csv + ".bak"
    shutil.copy2(args.csv, backup)
    with open(args.csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, dialect=dialect)
        writer.writeheader()
        writer.writerows(rows)
    log(f"Updated {args.csv} (original saved as {backup})")

    tag = date_tag()
    list_dir = Path(args.list_dir)
    list_dir.mkdir(parents=True, exist_ok=True)
    new_accs = list(dict.fromkeys(n for _, n in pairs))
    old_accs = list(dict.fromkeys(o for o, _ in pairs))
    for accs, kind in ((new_accs, "Updated"), (old_accs, "Removed")):
        path = list_dir / f"{len(accs)}_GCF_Accessions_{kind}_{tag}.txt"
        path.write_text("".join(a + "\n" for a in accs))
        log(f"Wrote {len(accs)} accession(s) to {path}")
    return pairs


# --------------------------------------------------------------------------- #
# Download / index helpers
# --------------------------------------------------------------------------- #
def name_matches(name, acc):
    """True if name is acc, or acc followed by a non-digit (so .6 never matches .61)."""
    return name == acc or (name.startswith(acc) and not name[len(acc)].isdigit())


def download_genome(acc, staging, paths):
    """Download and unzip one genome into staging/. Returns True on success."""
    staging.mkdir(parents=True, exist_ok=True)
    zip_path = staging / f"{acc}_dataset.zip"
    cmd = [str(paths.datasets), "download", "genome", "accession", acc,
           "--include", "genome", "--filename", str(zip_path)]
    if API_KEY:
        cmd += ["--api-key", API_KEY]
    try:
        subprocess.run(cmd, check=True)
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(staging)
    except (subprocess.CalledProcessError, OSError, zipfile.BadZipFile) as e:
        log(f"DOWNLOAD FAILED: {acc} - {e}")
        return False
    finally:
        if zip_path.exists():
            zip_path.unlink()
    return True


def staged_fastas(staging):
    return sorted((staging / "ncbi_dataset" / "data").glob("GC[AF]_*/*.fna"))


def index_prefix(tool, base_name, paths):
    return paths.tool_root(tool) / f"{base_name}_index" / f"{base_name}_index"


def index_complete(tool, prefix):
    """Check that the index files a tool needs exist next to this prefix."""
    p = str(prefix)
    if tool == "minimap2":
        mmi = Path(p + ".mmi")
        return mmi.is_file() and mmi.stat().st_size > 0
    if tool == "bowtie2":
        parts, exts = ("1", "2", "3", "4", "rev.1", "rev.2"), ("bt2", "bt2l")
    else:  # hisat2
        parts, exts = tuple(str(i) for i in range(1, 9)), ("ht2", "ht2l")
    return any(all(Path(f"{p}.{part}.{ext}").is_file() for part in parts) for ext in exts)


def has_complete_index(tool, acc, paths):
    """Is there a complete <tool> index directory for this accession?"""
    root = paths.tool_root(tool)
    if not root.is_dir():
        return False
    for d in root.glob(f"{acc}*_index"):
        if d.is_dir() and name_matches(d.name, acc) and index_complete(tool, d / d.name):
            return True
    return False


def build_index(tool, fna, paths, threads=None):
    """Build one index for one FASTA, logging runtime and CPU%. Returns True on success."""
    base_name = fna.name[:-len(".fna")]
    prefix = index_prefix(tool, base_name, paths)
    outdir = prefix.parent
    outdir.mkdir(parents=True, exist_ok=True)

    if tool == "bowtie2":
        cmd = ["bowtie2-build"] + (["--threads", str(threads)] if threads else []) \
              + [str(fna), str(prefix)]
    elif tool == "minimap2":
        cmd = ["minimap2"] + (["-t", str(threads)] if threads else []) \
              + ["-d", f"{prefix}.mmi", str(fna)]
    else:  # hisat2
        cmd = ["hisat2-build"] + (["-p", str(threads)] if threads else []) \
              + [str(fna), str(prefix)]

    log(f"--- {tool}: {base_name}")
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.time()
    try:
        subprocess.run(cmd, check=True)
        ok = index_complete(tool, prefix)
        if not ok:
            log(f"INDEX INCOMPLETE: {tool} finished but expected files are missing in {outdir}")
    except (subprocess.CalledProcessError, OSError) as e:
        log(f"INDEX FAILED: {tool} {base_name} - {e}")
        ok = False
    runtime = int(time.time() - start)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)

    if not ok:
        shutil.rmtree(outdir, ignore_errors=True)  # so VERIFY rebuilds it cleanly
        return False

    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    cpu_pct = f"{cpu / max(time.time() - start, 1e-9) * 100:.0f}%"
    h, m, s = runtime // 3600, (runtime % 3600) // 60, runtime % 60
    with open(paths.timing_log(tool), "a") as fh:
        fh.write(f"{base_name},{h}h:{m}m:{s}s,{cpu_pct}\n")
    return True


# --------------------------------------------------------------------------- #
# Stage 2: download new genomes and build indices
# --------------------------------------------------------------------------- #
def stage_build(accessions, paths, args):
    """Returns True if every download and index build succeeded."""
    banner("STAGE 2/4: Downloading genomes and building indices")
    paths.genomes.mkdir(parents=True, exist_ok=True)
    paths.done.mkdir(parents=True, exist_ok=True)
    staging = Path(args.work_dir) / "ncbi_download_tmp"
    all_ok = True

    if args.dry_run:
        log(f"Dry run - would download {len(accessions)} genome(s) and index every "
            f".fna in {paths.genomes}.")
        return True

    # Download
    for acc in accessions:
        all_ok &= download_genome(acc, staging, paths)
    log("Downloading genome files: complete")

    # Move FASTA files to the reference genome directory
    for fasta in staged_fastas(staging):
        shutil.move(str(fasta), str(paths.genomes / fasta.name))
    log("All fasta files have been moved")

    # Check that nothing else was left behind before deleting the staging area
    data_dir = staging / "ncbi_dataset" / "data"
    dirs = sorted(d for d in data_dir.glob("GC[AF]_*") if d.is_dir()) if data_dir.is_dir() else []
    non_empty = [d for d in dirs if any(d.iterdir())]
    for d in non_empty:
        log(f"NOT EMPTY: {d}/")
        for f in sorted(d.iterdir()):
            log(f"    - {f.name}")
    log(f"\nChecked {len(dirs)} GC*_ directories.")
    if non_empty:
        raise PipelineError(f"{len(non_empty)} directory(ies) are NOT empty (see above). "
                            f"Staging left in place at {staging}")
    log("All GC*_ directories are empty.")
    if staging.exists():
        shutil.rmtree(staging)
        log("Download staging directory removed")

    # Build indices for every FASTA in reference_genomes/
    fastas = sorted(paths.genomes.glob("*.fna"))
    log(f"\nIndexing {len(fastas)} genome(s) with {', '.join(TOOLS)}")
    for fna in fastas:
        results = {tool: build_index(tool, fna, paths, args.threads) for tool in TOOLS}
        if all(results.values()):
            shutil.move(str(fna), str(paths.done / fna.name))
        else:
            failed = [t for t, ok in results.items() if not ok]
            log(f"KEEPING {fna.name} in {paths.genomes} - failed: {', '.join(failed)}")
            all_ok = False
    log("Indexing pipeline complete")
    return all_ok


# --------------------------------------------------------------------------- #
# Stage 3: remove outdated indices and genomes
# --------------------------------------------------------------------------- #
def stage_remove(pairs, current_accs, paths, args):
    banner("STAGE 3/4: Removing outdated indices and genomes")
    current = set(current_accs)
    targets = []
    for old, new in dict.fromkeys(pairs):
        if old in current:
            log(f"SKIP:      {old} is still listed in dict.csv")
            continue
        if not args.dry_run and not all(has_complete_index(t, new, paths) for t in TOOLS):
            log(f"SKIP:      {old} - replacement {new} is not fully indexed yet")
            continue
        targets.append(old)
    targets = list(dict.fromkeys(targets))

    dirs, files = [], []
    for acc in targets:
        found = False
        for tool in TOOLS:
            root = paths.tool_root(tool)
            if not root.is_dir():
                continue
            for d in sorted(root.glob(f"{acc}*")):
                if d.is_dir() and name_matches(d.name, acc):
                    dirs.append(d)
                    log(f"FOUND:     {acc} -> {d}/")
                    found = True
        if paths.done.is_dir():
            for f in sorted(paths.done.glob(f"{acc}*.fna")):
                if f.is_file() and name_matches(f.name, acc):
                    files.append(f)
                    log(f"FOUND:     {acc} -> {f}")
                    found = True
        if not found:
            log(f"NOT FOUND: {acc} (no index directory or .fna file)")

    for root in [paths.tool_root(t) for t in TOOLS] + [paths.done]:
        if not root.is_dir():
            log(f"WARNING: directory not found: {root}")

    total = len(dirs) + len(files)
    log(f"\n{len(dirs)} index director(ies) and {len(files)} genome .fna file(s) found")
    if total == 0:
        return
    if args.dry_run:
        log("Dry run - nothing deleted.")
        return
    if not args.yes:
        if not sys.stdin.isatty():
            log("No terminal to confirm deletion (use --yes) - nothing deleted.")
            return
        answer = input(f"Permanently delete these {total} item(s)? [y/N] ")
        if answer.strip().lower() not in ("y", "yes"):
            log("Nothing deleted.")
            return

    for d in dirs:
        try:
            shutil.rmtree(d)
            log(f"Deleted: {d}/")
        except OSError as e:
            log(f"FAILED to delete: {d}/ - {e}")
    for f in files:
        try:
            f.unlink()
            log(f"Deleted: {f}")
        except OSError as e:
            log(f"FAILED to delete: {f} - {e}")
    log("Removing old indices and genomes: complete")


# --------------------------------------------------------------------------- #
# Stage 4: verify every dict.csv accession has all three indices
# --------------------------------------------------------------------------- #
def find_fasta(acc, paths):
    for d in (paths.done, paths.genomes):
        if d.is_dir():
            for f in sorted(d.glob(f"{acc}*.fna")):
                if name_matches(f.name, acc):
                    return f
    return None


def stage_verify(accessions, paths, args):
    """Returns True if every accession ends up with all three indices."""
    banner("STAGE 4/4: Verifying indices for every accession in dict.csv")
    missing_any, still_missing = 0, []
    staging = Path(args.work_dir) / "ncbi_verify_tmp"

    for acc in accessions:
        missing = [t for t in TOOLS if not has_complete_index(t, acc, paths)]
        if not missing:
            log(f"[OK     ] {acc}")
            continue
        missing_any += 1
        log(f"[MISSING] {acc}: {', '.join(missing)}")
        if args.dry_run:
            continue

        fna = find_fasta(acc, paths)
        if fna is None:
            log(f"          no FASTA on disk - downloading {acc}")
            if download_genome(acc, staging, paths):
                fastas = staged_fastas(staging)
                if fastas:
                    paths.done.mkdir(parents=True, exist_ok=True)
                    fna = paths.done / fastas[0].name
                    shutil.move(str(fastas[0]), str(fna))
            shutil.rmtree(staging, ignore_errors=True)
        if fna is None:
            log(f"          could not get a FASTA for {acc}")
            still_missing.append(acc)
            continue

        for tool in missing:
            build_index(tool, fna, paths, args.threads)
        left = [t for t in TOOLS if not has_complete_index(t, acc, paths)]
        if left:
            still_missing.append(acc)
            log(f"          still missing after rebuild: {', '.join(left)}")
        elif fna.parent == paths.genomes:
            paths.done.mkdir(parents=True, exist_ok=True)
            shutil.move(str(fna), str(paths.done / fna.name))

    log(f"\n{len(accessions)} accession(s) checked, {missing_any} had missing indices"
        + ("" if args.dry_run else f", {len(still_missing)} still incomplete."))
    for acc in still_missing:
        log(f"  INCOMPLETE: {acc}")
    return not still_missing


# --------------------------------------------------------------------------- #
# Final step: host_summary_manifesto.json
# --------------------------------------------------------------------------- #
# Same fields and order as https://www.bv-brc.org/api/content/host/patric_host_summary.json,
# with "patric_ftp" replaced by "url_path" (the genome .fna under reference_genomes/).
MANIFEST_FIELDS = (
    "seq_rel_date", "asm_name", "gbrs_paired_asm", "species_taxid", "version_status",
    "assembly_accession", "isolate", "taxid", "url_path", "assembly_level",
    "refseq_category", "submitter", "infraspecific_name", "biosample", "bioproject",
    "wgs_master", "paired_asm_comp", "organism_name", "species_common_name",
    "release_type", "genome_rep",
)
ASSEMBLY_SUMMARIES = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt",
    "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq_historical.txt",
)


def iter_assembly_summary(url):
    """Yield one dict per row of an NCBI assembly_summary file, streamed."""
    req = urllib.request.Request(url, headers={"User-Agent": "reference_index_pipeline"})
    with urllib.request.urlopen(req, timeout=300) as r:
        header = None
        for raw in r:
            line = raw.decode("utf-8", "replace").rstrip("\n")
            if line.startswith("#"):
                if "assembly_accession" in line:
                    header = line.lstrip("#").strip().split("\t")
                continue
            if header and line:
                yield dict(zip(header, line.split("\t")))


def summary_from_ncbi_files(accessions):
    """Look the accessions up in NCBI's RefSeq assembly_summary files."""
    wanted, found = set(accessions), {}
    for url in ASSEMBLY_SUMMARIES:
        if not wanted - set(found):
            break
        log(f"Reading {url}")
        try:
            for row in iter_assembly_summary(url):
                acc = row.get("assembly_accession")
                if acc in wanted and acc not in found:
                    found[acc] = row
        except (urllib.error.URLError, OSError) as e:
            log(f"WARNING: could not read {url} - {e}")
    return found


def summary_from_datasets_api(acc):
    """Best-effort fallback from the Datasets API (some fields may be blank)."""
    rep = assembly_report(acc) or {}
    info = rep.get("assembly_info") or {}
    org = rep.get("organism") or {}
    names = dict(org.get("infraspecific_names") or {})
    isolate = names.pop("isolate", "")
    status = info.get("assembly_status", "")
    return {
        "seq_rel_date": (info.get("release_date") or "").replace("-", "/"),
        "asm_name": info.get("assembly_name", ""),
        "gbrs_paired_asm": rep.get("paired_accession", ""),
        "species_taxid": str(org.get("tax_id", "")),
        "version_status": "latest" if status == "current" else status,
        "assembly_accession": acc,
        "isolate": isolate,
        "taxid": str(org.get("tax_id", "")),
        "assembly_level": info.get("assembly_level", ""),
        "refseq_category": info.get("refseq_category") or "na",
        "submitter": info.get("submitter", ""),
        "infraspecific_name": next((f"{k}={v}" for k, v in names.items()), ""),
        "biosample": (info.get("biosample") or {}).get("accession", ""),
        "bioproject": info.get("bioproject_accession", ""),
        "organism_name": org.get("organism_name", ""),
    }


_common_names = {}


def common_name(acc, taxid):
    """Common name (e.g. 'house mouse') from NCBI, cached by tax ID."""
    key = taxid or acc
    if key in _common_names:
        return _common_names[key]
    name = ""
    try:
        if taxid:
            data = ncbi_get(f"taxonomy/taxon/{urllib.parse.quote(str(taxid))}/dataset_report")
            for rep in (data or {}).get("reports", []):
                tax = rep.get("taxonomy") or {}
                name = (tax.get("curator_common_name") or tax.get("genbank_common_name")
                        or tax.get("common_name") or "")
                if name:
                    break
        if not name:
            name = ((assembly_report(acc) or {}).get("organism") or {}).get("common_name", "")
    except Exception as e:
        log(f"WARNING: {acc} - common name lookup failed ({e})")
    _common_names[key] = name
    return name


def write_host_manifest(accessions, paths, args):
    banner("FINAL: Writing host_summary_manifesto.json")
    out = Path(args.manifest)
    if args.dry_run:
        log(f"Dry run - would write {len(accessions)} genome(s) to {out.resolve()}")
        return True

    rows = summary_from_ncbi_files(accessions)
    genomes, ok = [], True
    for acc in accessions:
        row = rows.get(acc)
        if row is None:
            log(f"{acc}: not in RefSeq assembly summaries - using Datasets API")
            try:
                row = summary_from_datasets_api(acc)
            except Exception as e:
                log(f"WARNING: {acc} - no metadata found ({e})")
                row, ok = {"assembly_accession": acc}, False
        fna = find_fasta(acc, paths)
        if fna is None:
            log(f"WARNING: {acc} - no genome .fna found under {paths.genomes}")
            ok = False
        entry = {f: str(row.get(f, "") or "") for f in MANIFEST_FIELDS}
        entry["assembly_accession"] = acc
        entry["url_path"] = str(fna.resolve()) if fna else ""
        entry["species_common_name"] = common_name(acc, entry["species_taxid"] or entry["taxid"])
        if not entry["species_common_name"]:
            log(f"NOTE: {acc} - no common name found at NCBI")
        genomes.append(entry)

    tmp = out.with_name(out.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump({"genomes": genomes}, fh, ensure_ascii=False, indent=4)
        fh.write("\n")
    os.replace(tmp, out)
    log(f"Wrote {len(genomes)} genome(s) to {out.resolve()}")
    return ok


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?", default=DEFAULT_CSV,
                    help=f"CSV of species and GCF accessions (default: {DEFAULT_CSV})")
    ap.add_argument("--acc-col", help="name of the accession column (auto-detected)")
    ap.add_argument("--species-col", help="name of the species column (auto-detected)")
    ap.add_argument("--prefer-reference", action="store_true",
                    help="also switch to the species' current RefSeq reference genome "
                         "when it is a different assembly")
    ap.add_argument("--base", default=DEFAULT_BASE,
                    help=f"fastqutils base directory (default: {DEFAULT_BASE})")
    ap.add_argument("--datasets", help="path to the NCBI datasets binary "
                                       "(default: <base>/reference_indices/datasets)")
    ap.add_argument("--list-dir", default=".",
                    help="where to write the Updated/Removed accession lists (default: .)")
    ap.add_argument("--work-dir", default=".",
                    help="where to stage downloads before moving them (default: .)")
    ap.add_argument("--threads", type=int, default=DEFAULT_THREADS,
                    help="threads for bowtie2-build, minimap2 and hisat2-build "
                         f"(default: {DEFAULT_THREADS})")
    ap.add_argument("--yes", "-y", action="store_true",
                    help="delete outdated indices without asking")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would happen without writing, building or deleting")
    ap.add_argument("--manifest", default="host_summary_manifesto.json",
                    help="host summary JSON to write (default: host_summary_manifesto.json "
                         "in reference_indices/)")
    for stage in ("update", "build", "remove", "verify", "manifest"):
        ap.add_argument(f"--skip-{stage}", action="store_true", help=f"skip the {stage} stage")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    # Work from the reference_indices directory
    # (/vol/bvbrc/production/application-backend/fastqutils/reference_indices/).
    # Relative paths given on the command line are relative to this directory.
    os.chdir(Path(args.base) / "reference_indices")
    log(f"Working directory: {os.getcwd()}")
    paths = Paths(args.base, args.datasets)
    ok = True
    try:
        header, rows, dialect = read_csv(args.csv)
        acc_col, species_col = resolve_columns(header, rows, args.acc_col, args.species_col)

        pairs = []
        if not args.skip_update:
            pairs = stage_update(args, header, rows, dialect, acc_col, species_col)
        current = csv_accessions(rows, acc_col)
        new_accs = list(dict.fromkeys(n for _, n in pairs))

        if not args.skip_build:
            if new_accs or list(paths.genomes.glob("*.fna")):
                ok &= stage_build(new_accs, paths, args)
            else:
                log("\nSTAGE 2/4: nothing new to download or index - skipped")

        if not args.skip_remove:
            if pairs:
                stage_remove(pairs, current, paths, args)
            else:
                log("\nSTAGE 3/4: no replaced accessions - nothing to remove")

        if not args.skip_verify:
            ok &= stage_verify(current, paths, args)

        if not args.skip_manifest:
            ok &= write_host_manifest(current, paths, args)
    except PipelineError as e:
        log(f"\nERROR: {e}")
        return 1
    except KeyboardInterrupt:
        log("\nInterrupted.")
        return 130

    banner("Pipeline finished " + ("successfully" if ok else "with problems (see above)"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())