## Get or build reference
mm10_ref="/w5home/bmoore/cut_and_tag/mouse_mm10/refdata-cellranger-arc-mm10-2020-A-2.0.0"

## Demultiplex
# cellranger-arc mkfastq --run ${bcl2_dir} --csv ${script_dir}/SampleSheet.atac.csv --output-dir ${fastq_dir}

## cellranger	
# cellranger-atac count --id=${s} --project=${bcl2##*_} --reference=${cellranger_ref} --fastqs=${fastq_dir}/ --sample=${s} --chemistry=ARC-v1