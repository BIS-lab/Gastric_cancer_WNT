nextflow run nf-core/sarek -r 3.4.0 \
--input samplesheet_JH.csv \
-c sarek.config \
-profile docker \
--genome GATK.GRCh37 \
--save_reference \
--outdir sarek_JH \
-resume