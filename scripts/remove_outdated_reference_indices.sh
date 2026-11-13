#!/usr/bin/env bash
#
# remove_old_indices.sh
#
# Delete the bowtie2 and minimap2 index directories, and the genome .fna
# files in reference_genomes/done, whose names begin with an accession
# listed in a file (e.g. Remove_GCF_Indices.txt).
#
# Usage:
#   ./remove_old_indices.sh Remove_GCF_Indices.txt             # find, confirm, delete
#   ./remove_old_indices.sh Remove_GCF_Indices.txt --dry-run   # only list matches
#   ./remove_old_indices.sh Remove_GCF_Indices.txt --yes       # delete without asking
#
# A directory or file matches an accession if its name is exactly the
# accession or starts with it followed by a non-digit
# (GCF_000003025.6_Sscrofa11.1, GCF_000003025.6_genomic.fna).
# So GCF_000003025.6 never matches .61.

index_base=/vol/bvbrc/production/application-backend/fastqutils/reference_indices
index_roots=("$index_base/bowtie2" "$index_base/minimap2")
genome_dir=/vol/bvbrc/production/application-backend/fastqutils/reference_genomes/done

if [[ $# -lt 1 || ! -f "$1" ]]; then
    echo "Usage: $0 <accession_list.txt> [--dry-run | --yes]" >&2
    exit 1
fi
accession_file="$1"
dry_run=false
assume_yes=false
case "$2" in
    --dry-run) dry_run=true ;;
    --yes|-y)  assume_yes=true ;;
    "") ;;
    *) echo "Unknown option: $2" >&2; exit 1 ;;
esac

name_matches() {  # name_matches <name> <accession>
    [[ "$1" == "$2" || "$1" == "$2"[!0-9]* ]]
}

dirs=()
files=()
while IFS= read -r acc || [[ -n "$acc" ]]; do
    acc="${acc//[[:space:]]/}"
    [[ -z "$acc" ]] && continue

    found=false
    # Index directories (bowtie2, minimap2)
    for root in "${index_roots[@]}"; do
        [[ -d "$root" ]] || continue
        for dir in "$root"/"$acc"*/; do
            [[ -d "$dir" ]] || continue
            if name_matches "$(basename "$dir")" "$acc"; then
                dirs+=("${dir%/}")
                echo "FOUND:     $acc -> ${dir%/}/"
                found=true
            fi
        done
    done
    # Genome .fna files
    if [[ -d "$genome_dir" ]]; then
        for f in "$genome_dir"/"$acc"*.fna; do
            [[ -f "$f" ]] || continue
            if name_matches "$(basename "$f")" "$acc"; then
                files+=("$f")
                echo "FOUND:     $acc -> $f"
                found=true
            fi
        done
    fi
    $found || echo "NOT FOUND: $acc (no index directory or .fna file)"
done < "$accession_file"

for root in "${index_roots[@]}" "$genome_dir"; do
    [[ -d "$root" ]] || echo "WARNING: directory not found: $root" >&2
done

total=$(( ${#dirs[@]} + ${#files[@]} ))
echo
echo "${#dirs[@]} index director(ies) and ${#files[@]} genome .fna file(s) found"
[[ $total -eq 0 ]] && exit 0

if $dry_run; then
    echo "Dry run - nothing deleted."
    exit 0
fi

if ! $assume_yes; then
    read -r -p "Permanently delete these $total item(s)? [y/N] " answer
    if [[ ! "$answer" =~ ^[Yy]$ ]]; then
        echo "Nothing deleted."
        exit 0
    fi
fi

for dir in "${dirs[@]}"; do
    if rm -rf -- "$dir"; then
        echo "Deleted: $dir/"
    else
        echo "FAILED to delete: $dir/" >&2
    fi
done
for f in "${files[@]}"; do
    if rm -f -- "$f"; then
        echo "Deleted: $f"
    else
        echo "FAILED to delete: $f" >&2
    fi
done
echo "Removing old indices and genomes: complete"