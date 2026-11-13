#!/usr/bin/env bash
#
# remove_old_indices.sh
#
# Delete the bowtie2 and minimap2 index directories whose names begin with
# an accession listed in a file (e.g. Remove_GCF_Indices.txt).
#
# Usage:
#   ./remove_old_indices.sh Remove_GCF_Indices.txt             # find, confirm, delete
#   ./remove_old_indices.sh Remove_GCF_Indices.txt --dry-run   # only list matches
#   ./remove_old_indices.sh Remove_GCF_Indices.txt --yes       # delete without asking
#
# A directory matches an accession if its name is exactly the accession
# (GCF_000003025.6) or starts with it followed by a non-digit
# (GCF_000003025.6_Sscrofa11.1). So GCF_000003025.6 never matches .61.

index_base=/vol/bvbrc/production/application-backend/fastqutils/reference_indices
index_roots=("$index_base/bowtie2" "$index_base/minimap2")

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

matches=()
while IFS= read -r acc || [[ -n "$acc" ]]; do
    acc="${acc//[[:space:]]/}"
    [[ -z "$acc" ]] && continue

    found=false
    for root in "${index_roots[@]}"; do
        [[ -d "$root" ]] || continue
        for dir in "$root"/"$acc"*/; do
            [[ -d "$dir" ]] || continue
            name=$(basename "$dir")
            if [[ "$name" == "$acc" || "$name" == "$acc"[!0-9]* ]]; then
                matches+=("${dir%/}")
                echo "FOUND:     $acc -> ${dir%/}"
                found=true
            fi
        done
    done
    $found || echo "NOT FOUND: $acc (no bowtie2 or minimap2 directory)"
done < "$accession_file"

for root in "${index_roots[@]}"; do
    [[ -d "$root" ]] || echo "WARNING: index directory not found: $root" >&2
done

echo
echo "${#matches[@]} matching director(ies) found"
[[ ${#matches[@]} -eq 0 ]] && exit 0

if $dry_run; then
    echo "Dry run - nothing deleted."
    exit 0
fi

if ! $assume_yes; then
    read -r -p "Permanently delete these ${#matches[@]} director(ies)? [y/N] " answer
    if [[ ! "$answer" =~ ^[Yy]$ ]]; then
        echo "Nothing deleted."
        exit 0
    fi
fi

for dir in "${matches[@]}"; do
    if rm -rf -- "$dir"; then
        echo "Deleted: $dir"
    else
        echo "FAILED to delete: $dir" >&2
    fi
done
echo "Removing old indices: complete"