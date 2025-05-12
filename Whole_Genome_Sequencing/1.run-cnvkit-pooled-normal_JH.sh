#!/bin/bash

# Setting variables
anno='../99.Reference/refFlat_chr_removed.txt'
access='../99.Reference/access-5k-mappable.hg19_chr_removed.bed'
csv_path='./sarek_JH/csv/recalibrated.csv'
target='../99.Reference/SureSelect-V4-r2_chr_removed.bed'
outdir='result-pooled-normal-JH'
rscript_path='/opt/conda/bin/Rscript'

# Create output directory if it doesn't exist
mkdir -p "$outdir"

# Initialize arrays to hold BAM file paths
tumor_bams=()
normal_bams=()

# Read the CSV and prepare tumor and normal samples
while IFS=',' read -r patient sex status sample bam bai; do
    if [[ $status -eq 1 ]]; then # Tumor sample
        tumor_bams+=("$bam")
    else # Normal sample
        normal_bams+=("$bam")
    fi
done < <(tail -n +2 "$csv_path") # Skip the header line

# Convert arrays to space-separated strings for batch process
tumor_bam_paths=$(IFS=" " ; echo "${tumor_bams[*]}")
normal_bam_paths=$(IFS=" " ; echo "${normal_bams[*]}")

# Prepare output prefix for batch analysis
batch_output_prefix="${outdir}/pooled-analysis"

# Print cnvkit batch command for pooled analysis
cnvkit.py batch $tumor_bam_paths --normal $normal_bam_paths -t "$target" --annotate "$anno" -g "$access" -d "$batch_output_prefix" --rscript-path "$rscript_path" -p 100