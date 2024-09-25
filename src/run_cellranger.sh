## Get or build reference
mm10_ref="/w5home/bmoore/cut_and_tag/mouse_mm10/refdata-cellranger-arc-mm10-2020-A-2.0.0"
fastq_dir="/w5home/bmoore/cut_and_tag/GSE224560/fastqs/"
s="FC_H3K27ac"
## Demultiplex
# cellranger-arc mkfastq --run ${bcl2_dir} --csv ${script_dir}/SampleSheet.atac.csv --output-dir ${fastq_dir}

## cellranger	
cellranger-atac count --id=${s} --reference=${mm10_ref} --fastqs=${fastq_dir}/ --sample=${s} --localcores 8 --localmem=64 --chemistry=ARC-v1