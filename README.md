# CutTag_library
Code for processing cut and tag data. Data can be bulk or single cell.

## For bulk cut&tag:
### Make environment
```
conda env create -f environment_cut_tag.yml
```
### Activate environment
```
conda activate cut_tag
```
### Download SEACR 
(https://github.com/FredHutch/SEACR)
```
git clone git@github.com:FredHutch/SEACR.git
```
### Build Bowtie2 reference index if needed
```
nohup bowtie2-build "path to genome fasta file" "path to output with index name" &
# example:
nohup bowtie2-build AmexG_v6.0-DD_axolotl-omics_dataset/AmexG_v6.0-DD.fa AmexG_v6.0-DD_axolotl-omics_dataset/bowtie2_index/AmexG_v6.0-DD &
```
### get chromosome size file
```
samtools faidx "path to genome fasta file"
cut -f1,2 "path to genome fai file that was just created" > chromSize.txt
```
### Set variables in config file (found in src)
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

### Running cut and tag

#### First we run a python wrapper that runs Bowtie2
```
# from the src folder run:
python run_fastqc_bowtie2.py <your output directory>
```
**Note** that this script automatically nohups your bowtie run. You must wait for this to finish before you go on to step 2.

#### After Bowtie has finished, run the next python script to compute peaks and visualize
```
# from the src folder run:
python get_peaks-visualize.py <your output dir (can be same as bowtie2 dir)> <your bowtie2 dir from step 1> 
```
**Note** You may want to nohup this script, as this can take >10 hours to complete

### Outputs
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

## For single-cell CUT&tag:
### Download cellranger-atac
* GO to: https://support.10xgenomics.com/single-cell-atac/software/downloads/latest
* Follow installation instructions

### Download SRAs if data is public
**NOTE** SRA-toolkit is needed for this: https://hpc.nih.gov/apps/sratoolkit.html
```
<sratoolkit bin folder/fasterq-dump SRR########>
# example:
~/programs/sratoolkit.3.0.5-centos_linux64/bin/fasterq-dump SRR23343778
```

### Download or build reference using cellranger arc:
*  If you have human or mouse data, you can download cellranger's prebuilt reference genomes here: https://www.10xgenomics.com/support/software/cell-ranger-arc/downloads
*  If you are building a reference, check https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/inputs/mkref
*  Then modify the following template and run
```
# input ensemble genome and gtf files, as well as motif file in build_cellranger_ref.sh
source src/build_cellranger_ref.sh
```

### Run cellranger-atac
```
# edit variables:
# cell_ref="" # location of cellranger reference genome
# fastq_dir="" # location of fastq files
# s="" # sample ID

source src/run_cellranger.sh
```

### process cut&tag counts with Signac and Seurat
#### Make environment
```
conda env create -f cut_tag_sc_environment.yml
```
#### Activate environment
```
conda activate cut_tag_sc
```
#### Set variables in config file (found in src, modify "sc_cut_tag" part)
```
projPath = "folder where count data is from cell ranger"
H5 = "name of H5 file" # if NA looks for matrix.mtx, barcodes.tsv, and peaks.bed files
metadata = "metadata file name"
fragments = "fragment file name"
genes = "genes to visualize for gene activity feature plots"
processed_sc_seurat_file = "processed scRNAseq object to help annotate cut&tag object"
```
#### Run Rmd file
```
Rscript -e "rmarkdown::render('process_scCutTag.rmd')"
```
Outputs:
* chromatin-assay_obj.rds --> object with processed peaks
* QC_plots.pdf --> plots for TSS enrichment, fragment length, blacklist ratio, nucleosome signal, pct reads in peaks
* depth_correlation.pdf --> Correlation between sequencing depth and reduced dimension components
* cluster_umap.pdf --> umap of clusters
* gene_feature_plots.pdf --> gene activity plots
* seurat_mapping_plot.pdf --> plot of scRNA annotations and their mapping for scCut&Tag annotations
* DEpeak_plot.pdf --> differentially expressed peaks between 2 celltypes
* GOterm_barplots.pdf --> Go enrichment plots for given celltype
* CoveragePlot.pdf --> For given genes, showing coverage for each cell type

