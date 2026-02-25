#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" # We assume that the scripts are located in a subdirectory "scripts/"
# Works regardless of the current working directory.
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"  # Go up only 2 levels, not 3
RAW_DATA_DIR="${PROJECT_ROOT}/data/raw/total_RNA"

cd ${RAW_DATA_DIR}
ls -la ${RAW_DATA_DIR}
output_dir="${RAW_DATA_DIR}/cat_files"

mkdir -p "${output_dir}"

# Gunzip files and keep originals.
for file in *.gz; do
    echo "Extracting: $file"
    gunzip -k "$file"  # Keeps original .gz files
done

# For each sample, concatenate the lanes (only L001 and L002)
    for sample in $(ls *.fastq | grep -E '_L00[12]_' | cut -d'_' -f1-2 | sort -u); do
 
 # Concatenate R1 files (L001 + L002)
    cat "${sample}_L001_R1_001.fastq" "${sample}_L002_R1_001.fastq" > "${output_dir}/$(basename ${sample})_R1.fastq"
    
# Concatenate R2 files (L001 + L002)
    cat "${sample}_L001_R2_001.fastq" "${sample}_L002_R2_001.fastq" > "${output_dir}/$(basename ${sample})_R2.fastq"
    
    echo "Created: ${output_dir}/$(basename ${sample})_R1.fastq and ${output_dir}/$(basename ${sample})_R2.fastq"
done


# Remove the original lane files (optional)
rm *_L001_*.fastq *_L002_*.fastq