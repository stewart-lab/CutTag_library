# following cut and tag tutorial: https://yezhengstat.github.io/CUTTag_tutorial/
# make conda/mamba environment
# get bowtie2 via mamba
# mamba create --name cut_tag bowtie2
# mamba activate cut_tag
# conda install bioconda::fastqc
# mamba install fastp
# conda install bioconda::bedtools
# conda install bioconda::samtools
# conda install bioconda::seacr
# conda install bioconda::deeptools

mamba activate cut_tag

# set project path directory
projPath="/w5home/bmoore/cut_and_tag"

# make QC dir
mkdir -p ${projPath}/fastqFileQC/H3K27me3
mkdir -p ${projPath}/fastqFileQC/H3K36me3

# QC fastqs
fastqc -o ${projPath}/fastqFileQC/H3K27me3 -f fastq ${projPath}/fastqs/240226_EO_Ax1_K27me3_1_S26_L005_R1_001.fastq.gz
fastqc -o ${projPath}/fastqFileQC/H3K27me3 -f fastq ${projPath}/fastqs/240226_EO_Ax1_K27me3_1_S26_L005_R2_001.fastq.gz
fastqc -o ${projPath}/fastqFileQC/H3K36me3 -f fastq ${projPath}/fastqs/240226_EO_Ax1_K36me3_1_S27_L005_R1_001.fastq.gz
fastqc -o ${projPath}/fastqFileQC/H3K36me3 -f fastq ${projPath}/fastqs/240226_EO_Ax1_K36me3_1_S27_L005_R2_001.fastq.gz

# Trim fastqs with fastp
# get fastp:
#wget http://opengene.org/fastp/fastp
#chmod a+x ./fastp

fastp -i fastqs/240226_EO_Ax1_K27me3_1_S26_L005_R1_001.fastq.gz -I fastqs/240226_EO_Ax1_K27me3_1_S26_L005_R2_001.fastq.gz -o \
fastqs/240226_EO_Ax1_K27me3_1_S26_L005_R1_001_trimmed.fastq.gz -O fastqs/240226_EO_Ax1_K27me3_1_S26_L005_R2_001_trimmed.fastq.gz
fastp -i fastqs/240226_EO_Ax1_K36me3_1_S27_L005_R1_001.fastq.gz -I fastqs/240226_EO_Ax1_K36me3_1_S27_L005_R2_001.fastq.gz -o \
fastqs/240226_EO_Ax1_K36me3_1_S27_L005_R1_001_trimmed.fastq.gz -O fastqs/240226_EO_Ax1_K36me3_1_S27_L005_R2_001_trimmed.fastq.gz

# QC trimmed file again
fastqc -o ${projPath}/fastqFileQC/H3K27me3 -f fastq ${projPath}/fastqs/240226_EO_Ax1_K27me3_1_S26_L005_R1_001_trimmed.fastq.gz
astqc -o ${projPath}/fastqFileQC/H3K36me3 -f fastq ${projPath}/fastqs/240226_EO_Ax1_K36me3_1_S27_L005_R1_001_trimmed.fastq.gz


## Build the bowtie2 reference genome index if needed:
## bowtie2-build AmbMex60DD_axolotl_ncbi_dataset/data/GCA_002915635.3/GCA_002915635.3_AmbMex60DD_genomic.fna AmbMex60DD_axolotl_ncbi_dataset/data/GCA_002915635.3/bowtie2_index/AmbMex60DD
## for axolotl-omics genome:
## nohup bowtie2-build ${projPath}/AmexG_v6.0-DD_axolotl-omics_dataset/AmexG_v6.0-DD.fa ${projPath}/AmexG_v6.0-DD_axolotl-omics_dataset/bowtie2_index/AmexG_v6.0-DD &

# run bowtie2
cores=8
ref=${projPath}/AmexG_v6.0-DD_axolotl-omics_dataset/bowtie2_index/AmexG_v6.0-DD

mkdir -p ${projPath}/alignment/sam/bowtie2_summary
mkdir -p ${projPath}/alignment/bam
mkdir -p ${projPath}/alignment/bed
mkdir -p ${projPath}/alignment/bedgraph

# NOTE There is no need to trim reads from out standard 25x25 PE sequencing (--end-to-end), as adapter sequences will not be included in reads of inserts >25 bp. 
# However, for users performing longer sequencing, reads will need to be trimmed by Cutadapt and mapped by:
#  --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 to ignore any remaining adapter sequence at the 3’ ends of reads during mapping.

# align to axolotl genome
histName="H3K27me3"

nohup bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 ${projPath}/fastqs/240226_EO_Ax1_K27me3_1_S26_L005_R1_001_trimmed.fastq.gz \
-2 ${projPath}/fastqs/240226_EO_Ax1_K27me3_1_S26_L005_R2_001_trimmed.fastq.gz -S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt & > ${histName}.out

histName="H3K36me3"

nohup bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 ${projPath}/fastqs/240226_EO_Ax1_K36me3_1_S27_L005_R1_001_trimmed.fastq.gz \
-2 ${projPath}/fastqs/240226_EO_Ax1_K36me3_1_S27_L005_R2_001_trimmed.fastq.gz -S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt & > ${histName}.out

# align to the spike-in genome
# build ecoli refeernce
## bowtie2-build ${projPath}/ASM584v2_ecoli_ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna ${projPath}/ASM584v2_ecoli_ncbi_dataset/data/GCF_000005845.2/bowtie2Index/Ecoli

# set up
spikeInRef=${projPath}/ASM584v2_ecoli_ncbi_dataset/data/GCF_000005845.2/bowtie2Index/Ecoli
histName="H3K27me3"

# run bowtie2
bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} \
-1 ${projPath}/fastqs/240226_EO_Ax1_K27me3_1_S26_L005_R1_001_trimmed.fastq.gz -2 ${projPath}/fastqs/240226_EO_Ax1_K27me3_1_S26_L005_R2_001_trimmed.fastq.gz \
-S $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam &> $projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.txt
# For spike-in normalization, reads are aligned to the E. coli genome U00096.3 with two more parameters --no-overlap and --no-dovetail to avoid possible cross-mapping of the experimental genome to that of the carry-over E. coli DNA that is used for calibration.

seqDepthDouble=`samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth

#next histone
histName="H3K36me3"

bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} \
-1 ${projPath}/fastqs/240226_EO_Ax1_K36me3_1_S27_L005_R1_001_trimmed.fastq.gz -2 ${projPath}/fastqs/240226_EO_Ax1_K36me3_1_S27_L005_R2_001_trimmed.fastq.gz \
-S $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam &> $projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.txt

seqDepthDouble=`samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth

# see alignment_summary.rmd for visualization
# install r packages within R

# no need to remove duplicates: apparent duplication rate is low for high quality CUT&Tag datasets, 
# and even the apparent ‘duplicate’ fragments are likely to be true fragments

# assess mapped fragment size distribution
## CUT&Tag reactions targeting a histone modification predominantly results in fragments that are 
## nucleosomal lengths (~180 bp), or multiples of that length. CUT&Tag targeting transcription factors 
## predominantly produce nucleosome-sized fragments and variable amounts of shorter fragments, 
## from neighboring nucleosomes and the factor-bound site, respectively. Tagmentation of DNA on the 
## surface of nucleosomes also occurs, and plotting fragment lengths with single-basepair resolution 
## reveal a 10-bp sawtooth periodicity, which is typical of successful CUT&Tag experiments.
# conda install bioconda::samtools 
# 
mkdir -p $projPath/alignment/sam/fragmentLen
## Extract the 9th column from the alignment sam file which is the fragment length
histName="H3K27me3"
samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.sam | awk -F'\t' \
'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' \
>$projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.txt
histName="H3K36me3"
samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.sam | awk -F'\t' \
'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' \
>$projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.txt

# visualize in R
## fragment_length.rmd

# filter mapped reads by minimal quality score
histName="H3K27me3"
minQualityScore=2
samtools view -h -q $minQualityScore ${projPath}/alignment/sam/${histName}_bowtie2.sam \
>${projPath}/alignment/sam/${histName}_bowtie2.qualityScore$minQualityScore.sam
histName="H3K36me3"
minQualityScore=2
samtools view -h -q $minQualityScore ${projPath}/alignment/sam/${histName}_bowtie2.sam \
>${projPath}/alignment/sam/${histName}_bowtie2.qualityScore$minQualityScore.sam

# File format conversion to prep for peak calling
## conda install bioconda::bedtools

histName="H3K27me3"
## Filter and keep the mapped read pairs
samtools view -bS -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.qualityScore$minQualityScore.sam \
>$projPath/alignment/bam/${histName}_bowtie2.mapped.bam
## Convert into bed file format
bedtools bamtobed -i $projPath/alignment/bam/${histName}_bowtie2.mapped.bam -bedpe \
>$projPath/alignment/bed/${histName}_bowtie2.bed
## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
awk '$1==$4 && $6-$2 < 1000 {print $0}' $projPath/alignment/bed/${histName}_bowtie2.bed \
>$projPath/alignment/bed/${histName}_bowtie2.clean.bed
## Only extract the fragment related columns
cut -f 1,2,6 $projPath/alignment/bed/${histName}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  \
>$projPath/alignment/bed/${histName}_bowtie2.fragments.bed

# next histone mark
histName="H3K36me3"
samtools view -bS -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.qualityScore$minQualityScore.sam \
>$projPath/alignment/bam/${histName}_bowtie2.mapped.bam
bedtools bamtobed -i $projPath/alignment/bam/${histName}_bowtie2.mapped.bam -bedpe \
>$projPath/alignment/bed/${histName}_bowtie2.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' $projPath/alignment/bed/${histName}_bowtie2.bed \
>$projPath/alignment/bed/${histName}_bowtie2.clean.bed
cut -f 1,2,6 $projPath/alignment/bed/${histName}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  \
>$projPath/alignment/bed/${histName}_bowtie2.fragments.bed

# Assess replicate reproducibility
## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
binLen=500
histName="H3K27me3"
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $projPath/alignment/bed/${histName}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | \
awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >$projPath/alignment/bed/${histName}_bowtie2.fragmentsCount.bin$binLen.bed
histName="H3K36me3"
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $projPath/alignment/bed/${histName}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | \
awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >$projPath/alignment/bed/${histName}_bowtie2.fragmentsCount.bin$binLen.bed

# visulaize replicates
# replicate_corr_plot.rmd

# spike-in calibration
# The underlying assumption is that the ratio of fragments mapped to the primary genome to the E. coli genome is the same for a series of samples, 
# each using the same number of cells. Because of this assumption, we do not normalize between experiments or between batches of purified pATn5, which 
# can have very different amounts of carry-over E. coli DNA. Using a constant C to avoid small fractions in normalized data, we define a scaling factor S as
# S = C / (fragments mapped to E. coli genome)
# Normalized coverage is then calculated as:
# Normalized coverage = (primary_genome_coverage) * S
# The Constant is an arbitrary multiplier, typically 10,000. The resulting file will be comparatively small as a genomic coverage bedgraph file.

# to get chromosome size:
samtools faidx $projPath/AmexG_v6.0-DD_axolotl-omics_dataset/AmexG_v6.0-DD.fa
cut -f1,2 $projPath/AmexG_v6.0-DD_axolotl-omics_dataset/AmexG_v6.0-DD.fa.fai > chromSize_axl.txt

chromSize=${projPath}/AmexG_v6.0-DD_axolotl-omics_dataset/chromSize_axl.txt
histName="H3K27me3"

# assess via spike-in
seqDepth=`cat $projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth`
if [[ "$seqDepth" -gt "1" ]]; then
    
    mkdir -p $projPath/alignment/bedgraph

    scale_factor=`echo "10000 / $seqDepth" | bc -l`
    echo "Scaling factor for $histName is: $scale_factor!"
    bedtools genomecov -bg -scale $scale_factor -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize > \
    $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph
    
fi

# visualize scaling factor in R: scale.rmd

# if depth is 0, cannot scale, so making bed graph files without scaling
bedtools genomecov -bg -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize > \
$projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph
histName="H3K36me3"
bedtools genomecov -bg -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize > \
$projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph


# Peak calling using SEACR
# conda install bioconda::seacr

# Sparse Enrichment Analysis for CUT&RUN, SEACR, package is designed to call peaks and enriched regions 
# from chromatin profiling data with very low backgrounds (i.e., regions with no read coverage) that are 
# typical for CUT&Tag chromatin profiling experiments. SEACR requires bedGraph files from paired-end 
# sequencing as input and defines peaks as contiguous blocks of basepair coverage that do not overlap with 
# blocks of background signal. Since we have normalized fragment counts with the E. coli read count, we 
# set the normalization option of SEACR to “non”. Otherwise, the “norm” is recommended.

seacr="/w5home/bmoore/miniconda3/envs/cut_tag/bin/SEACR_1.3.sh"
mkdir -p $projPath/peakCalling/SEACR

# if using an IgG control, use the following command:
# histControl=$2
# $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
 #    $projPath/alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
 #   non stringent $projPath/peakCalling/SEACR/${histName}_seacr_control.peaks

# because we don't have a control file, just select top 1% of regions by AUC
histName="H3K27me3"
$seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph 0.01 non stringent \
$projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks
histName="H3K36me3"
$seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph 0.01 non stringent \
$projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks

# peak summary
# peak_sum.rmd

# visualize peaks in genome browser: 
# axolotl: https://genome.axolotl-omics.org/
# human/mouse etc: https://igv.org/app/

# for deeptools part, need annotation file

# make bigwig files normalized to genome
mkdir -p $projPath/alignment/bigwig
histName="H3K27me3" 
# sort                                                                                                                                       
samtools sort -o $projPath/alignment/bam/${histName}.sorted.bam \
$projPath/alignment/bam/${histName}_bowtie2.mapped.bam  
# index, -c for csi for genomes with chromosomes longer than 512Mb                                                   
samtools index -c $projPath/alignment/bam/${histName}.sorted.bam
# convert to bw and normalize                                                                                                              
bamCoverage -b $projPath/alignment/bam/${histName}.sorted.bam \
-o $projPath/alignment/bigwig/${histName}_norm.smooth.bw  --normalizeUsing CPM \
--binSize 500 --smoothLength 2000 --minMappingQuality 1 --numberOfProcessors $cores
# no smoothing
bamCoverage -b $projPath/alignment/bam/${histName}.sorted.bam \
-o $projPath/alignment/bigwig/${histName}_norm.bw  --normalizeUsing CPM \
--minMappingQuality 1 --numberOfProcessors $cores

# next mark
histName="H3K36me3"
# sort                                                                                                                                       
samtools sort -o $projPath/alignment/bam/${histName}.sorted.bam \
$projPath/alignment/bam/${histName}_bowtie2.mapped.bam  
# index, -c for csi for genomes with chromosomes longer than 512Mb                                                   
samtools index -c $projPath/alignment/bam/${histName}.sorted.bam
# convert to bw and normalize                                                                                                              
bamCoverage -b $projPath/alignment/bam/${histName}.sorted.bam \
-o $projPath/alignment/bigwig/${histName}_norm.smooth.bw  --normalizeUsing CPM \
--binSize 500 --smoothLength 2000 --minMappingQuality 1 --numberOfProcessors $cores
# no smoothing
bamCoverage -b $projPath/alignment/bam/${histName}.sorted.bam \
-o $projPath/alignment/bigwig/${histName}_norm.bw  --normalizeUsing CPM \
--minMappingQuality 1 --numberOfProcessors $cores
# no smoothing
bamCoverage -b $projPath/alignment/bam/${histName}.sorted.bam \
-o $projPath/alignment/bigwig/${histName}_norm.bw  --normalizeUsing CPM \
--minMappingQuality 1 --numberOfProcessors $cores

# use deeptools to compute matrix
cores=16
computeMatrix scale-regions -S $projPath/alignment/bigwig/H3K36me3_norm.smooth.bw \
                               $projPath/alignment/bigwig/H3K27me3_norm.smooth.bw \
                              -R $projPath/AmexG_v6.0-DD_axolotl-omics_dataset/AmexT_v47-AmexG_v6.0-DD.gtf \
                              --beforeRegionStartLength 5000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 10000 \
                              --skipZeros -o $projPath/output_plots/AmexT_v47-AmexG_v6.0-DD_H3K36me3-H3K27me3_gene_cpm_smooth10k.mat.gz \
                              -p $cores

# heat map on genes
stat=mean
plotHeatmap -m $projPath/output_plots/AmexT_v47-AmexG_v6.0-DD_H3K27me3-H3K36me3_gene_cpm_smooth10k.mat.gz \
-o $projPath/output_plots/AmexT_v47-AmexG_v6.0-DD_H3K27me3-H3K36me3_gene_cpm-smooth10k.pdf --averageTypeSummaryPlot $stat \
--sortUsing sum --heatmapHeight 16 --heatmapWidth 8 --outFileSortedRegions \
$projPath/output_plots/AmexT_v47-AmexG_v6.0-DD_H3K27me3-H3K36me3_gene.histone.cpm.smooth10k.bed

# sort by sample H3K36me3
stat=mean
plotHeatmap -m $projPath/output_plots/AmexT_v47-AmexG_v6.0-DD_H3K36me3-H3K27me3_gene_cpm_smooth10k.mat.gz \
-o $projPath/output_plots/AmexT_v47-AmexG_v6.0-DD_H3K36me3-H3K27me3_gene_cpm-smooth10k.pdf --averageTypeSummaryPlot $stat \
--sortUsing sum --heatmapHeight 16 --heatmapWidth 8 --outFileSortedRegions \
$projPath/output_plots/AmexT_v47-AmexG_v6.0-DD_H3K36me3-H3K27me3_gene.histone.cpm.smooth10k.bed

# Heatmap on CUT&Tag peaks
histName="H3K27me3"
# get summit region
awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' \
$projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks.stringent.bed >\
$projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks.summitRegion.bed
# compute matrix
computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_norm.bw \
              -R $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks.summitRegion.bed \
              --skipZeros -o $projPath/peakCalling/SEACR/${histName}_SEACR.mat.gz -p $cores -a 3000 -b 3000 \
              --referencePoint center

plotHeatmap -m $projPath/peakCalling/SEACR/${histName}_SEACR.mat.gz -out $projPath/peakCalling/SEACR/${histName}_SEACR_heatmap.pdf \
            --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel \
            "${histName}" --averageTypeSummaryPlot $stat 

# next mark
histName="H3K36me3"

awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' \
$projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks.stringent.bed >\
$projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks.summitRegion.bed
# compute matrix
computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_norm.bw \
              -R $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks.summitRegion.bed \
              --skipZeros -o $projPath/peakCalling/SEACR/${histName}_SEACR.mat.gz -p $cores -a 3000 -b 3000 \
              --referencePoint center

plotHeatmap -m $projPath/peakCalling/SEACR/${histName}_SEACR.mat.gz -out $projPath/peakCalling/SEACR/${histName}_SEACR_heatmap.pdf \
            --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel \
            "${histName}" --averageTypeSummaryPlot $stat 