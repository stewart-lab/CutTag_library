# download bedops
wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2
# unpack
tar jxvf bedops_linux_x86_64-v2.4.41.tar.bz2
# copy to binary and export path
cp bin/* ~/bin/
export PATH=~/bin:$PATH
# get only genes from gtf file
awk '$3 == "gene"' human_genome_38/Homo_sapiens.GRCh38.108.filtered.gtf > human_genome_38/Homo_sapiens.GRCh38.108.filtered_genes.gtf 
# convert gtf to bed using bedops
convert2bed -i gtf --attribute-key="gene_name" < human_genome_38/Homo_sapiens.GRCh38.108.filtered_genes.gtf > human_genome_38/Homo_sapiens.GRCh38.108.filtered_genes.bed
# get 300bp past TSS
python parse_gene_bed.py /w5home/bmoore/genomes/human_genome_38/Homo_sapiens.GRCh38.108.filtered_genes.bed
# get flanking region for 400 bp upstream gene
bedtools flank -i /w5home/bmoore/genomes/human_genome_38/Homo_sapiens.GRCh38.108.filtered_genes.bed_300bpTSS.bed -g /w5home/bmoore/genomes/human_genome_38/chromSize_humanGR38.txt -l 400 -r 0 -s > /w5home/bmoore/genomes/human_genome_38/Homo_sapiens.GRCh38.108.filtered_genes.400bp.promoters.bed