# load libraries
library(reticulate)
use_condaenv("/w5home/bmoore/miniconda3/envs/cut_tag3")
#.libPaths(c("/w5home/bmoore/miniconda3/envs/cut_tag3/lib/R/library", .libPaths()))
library(rmarkdown)
library(ggplot2)
library(tidyverse)
library(stringr)
library(viridis)
library(ggpubr)
library(corrplot)
library(GenomicRanges)
library(chromVAR)
library(DESeq2)

# get arguments
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
projPath = args[1]
subDir <- "DESeq2_out"
# for testing
# projPath <- "/w5home/bmoore/kratz_cut_tag/"

# set working directory
setwd(projPath)
# create output directory
outDir = dir.create(file.path(projPath, subDir), showWarnings = FALSE)


# read in config.json
config = jsonlite::fromJSON(paste0(projPath, "/config.json"), simplifyDataFrame = FALSE)

# Extract keys that start with "SAMPLE"
sample_keys <- names(config)[grepl("^SAMPLE", names(config))]
samples <- lapply(sample_keys, function(key) config[[key]])
names(samples) <- sample_keys

# create master peak list
mPeak = GRanges()
histList = c()
sampleList = c()
## overlap with bam file to get count
for (sample_name in names(samples)) {
  sample <- samples[[sample_name]]
  # append histName to histList
  histList = c(histList, sample$histName)
  # get sample name
  samplename <- paste0(sample$histName, "Rep", as.character(sample$rep))
  sampleList = c(sampleList, samplename)
    peakRes = read.table(paste0(projPath, "/peakCalling/SEACR/S2/", samplename, "_seacr_top0.01.peaks.stringent.bed"), header = FALSE, fill = TRUE)
    mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
  }
masterPeak = reduce(mPeak)
## get counts
histL=c()
repL=c()
for (sample_name in names(samples)) {
  sample <- samples[[sample_name]]
    histL = c(histL, sample$histName)
    repL = c(repL, sample$rep)
}
histL <- unique(histL)
repL <- unique(repL)
bamDir = paste0(projPath, "/alignment/bam")
countMat = matrix(NA, length(masterPeak), length(histL)*length(repL))
## overlap with bam file to get count
i = 1
for(hist in histL){
  for(rep in repL){
    
    bamFile = paste0(bamDir, "/", hist, "Rep", rep, "_bowtie2.mapped.bam")
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
  }
}
colnames(countMat) = paste(rep(histL, 3), rep(repL, each = 2), sep = "_")
# get condition dataframe
samplelist <- c(paste(rep(histL, 3), rep(repL, each = 2), sep = "_"))
condition_df <- data.frame(sample = samplelist, condition = histL)
# DEseq
selectR = which(rowSums(countMat) > 5) ## remove low count genes
dataS = countMat[selectR,]
dds = DESeqDataSetFromMatrix(countData = dataS,
                              colData = condition_df,
                              design = ~ condition)
DDS = DESeq(dds)
# normalize with respect to seq depth
normDDS = counts(DDS, normalized = TRUE) 
colnames(normDDS) = paste0(colnames(normDDS), "_norm")
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
masterPeak <- as.data.frame(masterPeak)
masterPeakS <- masterPeak[selectR,]
countMatDiff = cbind(masterPeakS, dataS, normDDS, res)

# write out
write.csv(countMatDiff, paste0(projPath, "/", subDir,"/",histL[1],"vs.",histL[2],".csv"))
# make bedgraph file
bedgraphdf <- subset(countMatDiff, select=c(seqnames,start,end,log2FoldChange))
bedgraphdf$log2FoldChange <- as.numeric(bedgraphdf$log2FoldChange)
write.table(bedgraphdf, paste0(projPath, "/", subDir,"/",histL[1],"vs.",histL[2],".bedgraph"),
            row.names = FALSE, sep="\t", quote = FALSE, col.names = FALSE)
# write out packages
writeLines(capture.output(sessionInfo()), con = paste0(projPath, "/", subDir,"/DESeq2_sessionInfo.txt"))
