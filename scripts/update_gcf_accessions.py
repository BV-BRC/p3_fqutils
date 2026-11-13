#!/usr/bin/env python3
"""
update_gcf_accessions.py

Check whether the NCBI RefSeq (GCF_) genome assembly accessions in a CSV are
up to date. Any outdated accession is replaced in the CSV, and the new
accession is appended to GCF_Accessions.txt.

An accession counts as outdated when:
  1. A newer version of the same assembly exists
     (e.g. GCF_000001635.26 -> GCF_000001635.27), or
  2. The assembly has been suppressed/retired with no newer version. It is
     then replaced by the species' current RefSeq reference genome.
  3. (Only with --prefer-reference) The species' current RefSeq reference
     genome is a different assembly from the one in the CSV.

Uses the NCBI Datasets v2 REST API and only the Python standard library.

Usage:
    python update_gcf_accessions.py species.csv
    python update_gcf_accessions.py species.csv --dry-run
    python update_gcf_accessions.py species.csv --acc-col Accession --species-col Species
    NCBI_API_KEY=xxxx python update_gcf_accessions.py species.csv   # faster rate limit
"""

import argparse
import csv
import json
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ACC_RE = re.compile(r"^(GC[FA]_\d{9})(?:\.(\d+))?$")
API_KEY = os.environ.get("NCBI_API_KEY")
# NCBI allows 5 requests/s without a key, 10/s with one.
MIN_INTERVAL = 0.11 if API_KEY else 0.22
_last_call = 0.0


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
    # Prefer the exact accession if several came back
    for rep in reports:
        if rep.get("accession") == acc:
            return rep
    return reports[0] if reports else None


def latest_version(acc):
    """Return the newest GCF_ version of the same assembly, or None."""
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
    return max(versions)[1] if versions else None, rep


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


# --------------------------------------------------------------------------- #
# Core logic
# --------------------------------------------------------------------------- #
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="CSV of species and GCF accessions")
    ap.add_argument("--acc-col", help="name of the accession column (auto-detected)")
    ap.add_argument("--species-col", help="name of the species column (auto-detected)")
    ap.add_argument("--out-list", default="GCF_Accessions.txt",
                    help="file to append new accessions to (default: GCF_Accessions.txt)")
    ap.add_argument("--remove-list", default="Remove_GCF_Indices.txt",
                    help="file to append replaced (old) accessions to "
                         "(default: Remove_GCF_Indices.txt)")
    ap.add_argument("--prefer-reference", action="store_true",
                    help="also switch to the species' current RefSeq reference genome "
                         "when it is a different assembly")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would change without writing anything")
    args = ap.parse_args()

    with open(args.csv, newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(fh, dialect=dialect)
        header = reader.fieldnames
        rows = list(reader)

    acc_col, species_col = guess_columns(header, rows)
    acc_col = args.acc_col or acc_col
    species_col = args.species_col or species_col
    if not acc_col or acc_col not in header:
        sys.exit(f"Could not find the accession column. Columns: {header}. "
                 f"Use --acc-col.")
    print(f"Accession column: {acc_col!r}   Species column: {species_col!r}\n")

    new_accessions, old_accessions, changes, problems = [], [], 0, 0
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
        print(f"[{flag:7}] {label}: {acc} - {reason}")
        if new_acc:
            row[acc_col] = new_acc
            new_accessions.append(new_acc)
            old_accessions.append(acc)
            changes += 1
        elif flag == "CHECK":
            problems += 1

    print(f"\n{len(rows)} rows checked, {changes} updated, {problems} need review.")

    if args.dry_run or not changes:
        if args.dry_run:
            print("Dry run - no files written.")
        return

    # Back up and rewrite the CSV
    backup = args.csv + ".bak"
    shutil.copy2(args.csv, backup)
    with open(args.csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, dialect=dialect)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Updated {args.csv} (original saved as {backup})")

    # New accessions -> GCF_Accessions.txt; replaced ones -> Remove_GCF_Indices.txt
    for path, accs, what in ((args.out_list, new_accessions, "new"),
                             (args.remove_list, old_accessions, "old")):
        n = append_unique(path, accs)
        print(f"Added {n} {what} accession(s) to {path}")


def append_unique(path, accessions):
    """Append accessions to a text file, skipping any already listed."""
    existing = set()
    if os.path.exists(path):
        with open(path) as fh:
            existing = {line.strip() for line in fh if line.strip()}
    to_add = [a for a in dict.fromkeys(accessions) if a not in existing]
    with open(path, "a") as fh:
        for a in to_add:
            fh.write(a + "\n")
    return len(to_add)


if __name__ == "__main__":
    main()