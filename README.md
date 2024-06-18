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
## Download SEACR 
(https://github.com/FredHutch/SEACR)
```
git clone git@github.com:FredHutch/SEACR.git
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
## Set variables in config file (found in src)
```
projPath = "top level directory where data will be processed"
Number_of_samples = "number of samples you have"
SAMPLE_1 = "each samples information. Copy and repeat for every sample"
  fastq_PE1 = "where fastq read 1 is located relative to projPath"
  fastq_PE2 = "where fastq read 2 is located relative to projPath"
  histName = "name of histone mark analyzing"
  rep = "replicate id- right now use a number, like 1 if it is rep 1"
runfastp = "if you want to run fastp (default), then true, else false"
ref = "path to bowtie2 index"
spikeInRef = "path to spike-in reference (ie. E-coli)"
spikeIn = "if you are using a spike-in like E.coli, then true (default), else false"
cores = "number of cpus to use"
chromSize = "path to chromSize.txt file"
seacr = "path to SEACR .sh file"
minQualityScore = "minimum quality score for filtering reads, typically 2"
binLen = "length of window for fragments, typically 500"
```

## Running cut and tag

### First we run a python wrapper that runs Bowtie2
```
# from the src folder run:
python run_fastqc_bowtie2.py <your output directory>
```
**Note** that this script automatically nohups your bowtie run. You must wait for this to finish before you go on to step 2.

### After Bowtie has finished, run the next python script to compute peaks and visualize
```
# from the src folder run:
python get_peaks-visualize.py <your output dir (can be same as bowtie2 dir)> <your bowtie2 dir from step 1> 
```
**Note** You may want to nohup this script, as this can take >10 hours to complete

## Outputs
* In alignment folder:
  * Alignment_summary.txt file and pdf summarizing Bowtie2 alignment
  * bam, bed, bedgraph, bigwig, and sam files. Bigiwg files can be uploaded to the genome browser for visualization of specific genes.
  * scale_boxplot.pdf summarizing scaling factor and normalization to spike-in
  * in sam/fragmentLen folder
    * fragment_length.pdf summarizing lengths of mapped fragments
* In peakCalling folder:
  * peak_plots.pdf summarizes peak results
  * In SEACR folder:
    * heatmaps centered around gene (_gene_cpm_smooth) and peak (_SEACR_heatmap) and their matrices and bed files
* log files from both bowtie run and peak run

