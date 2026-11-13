set -euo pipefail

# Download reference genomes
if [[ $# -lt 1 || ! -f "$1" ]]; then
    echo "Usage: $0 <accession_list.txt>" >&2
    exit 1
fi

accession_file="$1"
outdir=./

while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "$line" ]] && continue
    /vol/bvbrc/production/application-backend/fastqutils/reference_indices/datasets download genome accession "$line" --include genome --filename "${outdir}/${line}_dataset.zip"
done < "$accession_file"

echo "Downloading genome files: complete"

for f in ./*.zip; do 
    unzip -o "$f" && rm "$f"; 
done

echo "Unzipping genome files: complete"

# Move FASTA to referemce genome directory
for fasta in ./ncbi_dataset/data/GCF_*/*.fna; do 
    mv "$fasta" /vol/bvbrc/production/application-backend/fastqutils/reference_genomes;
done 

echo "All fasta files have been moved"

SRC_DIR="./ncbi_dataset/data"

shopt -s nullglob

non_empty_found=0
checked=0

for dir in "$SRC_DIR"/GCF_*/; do
    [ -d "$dir" ] || continue
    checked=$((checked + 1))

    if [ -n "$(ls -A "$dir" 2>/dev/null)" ]; then
        non_empty_found=$((non_empty_found + 1))
        echo "NOT EMPTY: $dir"
        ls -A "$dir" | sed 's/^/    - /'
    fi
done

echo ""
echo "Checked $checked GCF_* directories."
if [ "$non_empty_found" -eq 0 ]; then
    echo "All GCF_* directories are empty."
else
    echo "$non_empty_found directory(ies) are NOT empty (see above)."
    exit 1
fi

rm -r ./ncbi_dataset

echo "ncbi_dataset directory removed"

# Create bowtie2 reference index
for file in /vol/bvbrc/production/application-backend/fastqutils/reference_genomes/*.fna; do 
    base_name=$(basename "$file" .fna)
    mkdir -p "/vol/bvbrc/production/application-backend/fastqutils/reference_indices/bowtie2/${base_name}_index/"
    SRC_DIR="/vol/bvbrc/production/application-backend/fastqutils/reference_indices/bowtie2/${base_name}_index"

    start=$(date +%s)
    /usr/bin/time -f "%P" -o /tmp/cpu_stats.txt \
        bowtie2-build "$file" "$SRC_DIR"/${base_name}_index
    end=$(date +%s)

    cpu_pct=$(cat /tmp/cpu_stats.txt)

    runtime=$((end - start))
    hours=$((runtime / 3600))
    minutes=$(( (runtime % 3600) / 60 ))
    seconds=$((runtime % 60))

    echo "${base_name},${hours}h:${minutes}m:${seconds}s,${cpu_pct}" >> /vol/bvbrc/production/application-backend/fastqutils/reference_indices/bowtie2indexing_timing_log.csv

done 

echo "Bowtie2 indexing complete"


# Create mnimap2 reference index
for file in /vol/bvbrc/production/application-backend/fastqutils/reference_genomes/*.fna; do 
    base_name=$(basename "$file" .fna)
    mkdir -p "/vol/bvbrc/production/application-backend/fastqutils/reference_indices/minimap2/${base_name}_index/"
    SRC_DIR="/vol/bvbrc/production/application-backend/fastqutils/reference_indices/minimap2/${base_name}_index"

    start=$(date +%s)
    /usr/bin/time -f "%P" -o /tmp/cpu_stats.txt \
       minimap2 -d "$SRC_DIR"/${base_name}_index.mmi "$file"
    end=$(date +%s)

    cpu_pct=$(cat /tmp/cpu_stats.txt)

    runtime=$((end - start))
    hours=$((runtime / 3600))
    minutes=$(( (runtime % 3600) / 60 ))
    seconds=$((runtime % 60))

    echo "${base_name},${hours}h:${minutes}m:${seconds}s,${cpu_pct}" >> /vol/bvbrc/production/application-backend/fastqutils/reference_indices/minimap2indexing_timing_log.csv

done 

echo "Indexing complete"

for file in /vol/bvbrc/production/application-backend/fastqutils/reference_genomes/*.fna; do 
    mv "$file" /vol/bvbrc/production/application-backend/fastqutils/reference_genomes/done/

done

echo "Indexing pipeline complete"
