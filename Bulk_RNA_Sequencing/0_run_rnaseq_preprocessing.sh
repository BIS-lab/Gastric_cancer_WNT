nextflow run nf-core/rnaseq -r 3.12.0 \
--input 0_RNA_seq_Sample_Sheet.csv \
-profile docker \
--genome GRCh37 \
--save_reference \
-c 0_rnaseq.config \
--outdir './0_result' \
--aligner 'star_salmon' \
-resume