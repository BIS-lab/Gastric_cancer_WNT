#!/bin/bash

folders=(
  "./result-pooled-Daniel/pooled-analysis"
  "./result-pooled-norma-HK/pooled-analysis"
  "./result-pooled-normal-JH/pooled-analysis"
)

output_dir="./CNV_visualizaiton"
mkdir -p "$output_dir"

cnvkit_cmd="cnvkit.py"

genes=("WNT2" "KRAS")

for folder in "${folders[@]}"; do
    echo "Searching in folder: $folder"
    
    find "$folder" -type f -name "*.cns" ! -name "*call.cns" ! -name "*bintest.cns" | while read cns_file; do
        base_name=$(basename "$cns_file" .cns)
        cns_dir=$(dirname "$cns_file")
        cnr_file="${cns_dir}/${base_name}.cnr"

        output_file="${output_dir}/${base_name}_gene.txt"

        echo "Processing $base_name -> $output_file"

        # gene metrics
        $cnvkit_cmd genemetrics "$cnr_file" -s "$cns_file" -t 0 -o "$output_file"
    done
done