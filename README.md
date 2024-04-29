# CutTag_library
Code for processing cut and tag data

## Make environment
```
conda env create -f environment_cut_tag.yml
```
## Activate environment
```
conda activate cut_tag
```
## Build Bowtie2 reference index if needed
```
nohup bowtie2-build "path to genome fasta file" "path to output with index name" &
# example:
nohup bowtie2-build AmexG_v6.0-DD_axolotl-omics_dataset/AmexG_v6.0-DD.fa AmexG_v6.0-DD_axolotl-omics_dataset/bowtie2_index/AmexG_v6.0-DD &
```
## get chromosome size file
```
samtools faidx "path to genome fasta file"
cut -f1,2 "path to genome fai file that was just created" > chromSize.txt
```
## Set variables in cut_tag.sh script
```
projPath="directory where data will be processes"
fastq_PE1="where fastq read 1 is located relative to projPath"
fastq_PE2="where fastq read 2 is located relative to projPath"
histName1="name of histone mark analyzing"
histName2="second histone mark (if needed)"
ref="path to bowtie2 index"
spikeInRef="path to spike-in reference (ie. E-coli)
binLen=500 # length of window for fragments
cores=16 # number of cpus
chromSize="path to chromSize.txt file"
seacr="path to SEACR .sh file"
minQualityScore=2 # minimum quality score for filtering reads
```

## Running cut and tag
```
source cut_tag.sh
```
