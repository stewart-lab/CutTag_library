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
#args = commandArgs(trailingOnly=TRUE)
#projPath = args[1]
subDir <- "DESeq2_out"
# for testing
projPath <- "/w5home/bmoore/kratz_cut_tag/"

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

# get bed files for genes and promoter
gene_bedfile <- "/w5home/bmoore/genomes/human_genome_38/Homo_sapiens.GRCh38.108.filtered_genes.bed"
genebed = read.table(gene_bedfile, header = FALSE, fill = TRUE)
promoter_bedfile <- "/w5home/bmoore/genomes/human_genome_38/Homo_sapiens.GRCh38.108.filtered_genes.1kb.promoters.bed"
promoterbed = read.table(promoter_bedfile, header = FALSE, fill = TRUE)

# read into Granges
# gene
geneRanges <- GRanges(seqnames = genebed$V1,
    ranges = IRanges(start = genebed$V2, end = genebed$V3,
                  names = genebed$V11),strand= genebed$V6,
                  genenames=genebed$V4)
# promoter
promoterRanges <- GRanges(seqnames = promoterbed$V1,
    ranges = IRanges(start = promoterbed$V2, end = promoterbed$V3,
                  names = promoterbed$V11),strand= promoterbed$V6,
                  genenames=promoterbed$V4)
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
# make empty dataframe for foldchange
FC_df <- data.frame(log2fc=numeric(),hist=character())

## overlap with bam file to get count
## loop through each hist/ mark
for(hist in histL){
  countMat = matrix(NA, length(geneRanges), length(repL)*2)
  i = 1
  namelist=c()
  conditionL=c()
  for(rep in repL){
    bamFile = paste0(bamDir, "/", hist, "Rep", rep, "_bowtie2.mapped.bam")
    fragment_counts <- getCounts(bamFile, geneRanges, paired = TRUE,
                                by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    name <- paste0(hist,"_",rep,"_gene")
    namelist <- append(namelist, name)
    conditionL <- append(conditionL, "gene")
    i = i + 1
    fragment_counts2 <- getCounts(bamFile, promoterRanges, paired = TRUE,
                                  by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts2)[,1]
    i = i + 1
    name2 <- paste0(hist,"_",rep,"_promoter")
    namelist <- append(namelist, name2)
    conditionL <- append(conditionL, "promoter")
  }
  colnames(countMat) = paste(namelist)
  # get condition dataframe
  condition_df <- data.frame(sample = namelist, condition = conditionL)
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
  geneRangesDF <- annoGR2DF(geneRanges)
  geneRangesDF2 <- geneRangesDF[selectR,]
  countMatDiff = cbind(geneRangesDF2, dataS, normDDS, res)

  # write out
  write.csv(countMatDiff, paste0(projPath, "/", subDir,"/",hist,"_genevs.promoter_logFC.csv"))
  # make bedgraph file
  bedgraphdf <- subset(countMatDiff, select=c(chr,start,end,log2FoldChange))
  bedgraphdf$log2FoldChange <- as.numeric(bedgraphdf$log2FoldChange)
  write.table(bedgraphdf, paste0(projPath, "/", subDir,"/",hist,"_genevs.promoter_logFC.bedgraph"),
            row.names = FALSE, sep="\t", quote = FALSE, col.names = FALSE)
  # add logFC to dataframe
  hist =c(hist)
  hist_vec= rep(hist, each=dim(countMatDiff)[1])
  FC_hist_df <- as.data.frame(countMatDiff$log2FoldChange)
  FC_hist_df <- cbind(FC_hist_df,hist_vec)
  colnames(FC_hist_df) <- c("log2fc","hist")
  # bind to big df
  FC_df <- rbind(FC_df,FC_hist_df)
}
write.csv(FC_df,paste0(projPath, "/", subDir,"/","all_FC_promotervs.gene.csv"))
# calculate cumalitive distribution (CDF)
#CDF <- ecdf(countMatDiff$log2FoldChange)
pdf(file="sample_CDF.pdf")
ggplot(FC_df, aes(log2fc, colour = hist)) + stat_ecdf()
dev.off()
# subset
FC_df_1 <- subset(FC_df, hist %in% c("ComboIgG","CtrlIgG","VipIgG","ZenIgG"))
FC_df_2 <- subset(FC_df, hist %in% c("ComboS2","CtrlS2","VipS2","ZenS2"))
# plot separately
pdf(file="IgG_CDF.pdf")
ggplot(FC_df_1, aes(log2fc, colour = hist)) + stat_ecdf()
dev.off()

pdf(file="S2_CDF.pdf")
ggplot(FC_df_2, aes(log2fc, colour = hist)) + stat_ecdf()
dev.off()


# write out packages
writeLines(capture.output(sessionInfo()), con = paste0(projPath, "/", subDir,"/DESeq2_2_sessionInfo.txt"))
