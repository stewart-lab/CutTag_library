# load libraries
library(reticulate)
use_condaenv("/w5home/bmoore/miniconda3/envs/cut_tag")
library(rmarkdown)
library(ggplot2)
library(tidyverse)
library(stringr)
library(viridis)
library(ggpubr)
library(corrplot)
library(GenomicRanges)
library(chromVAR)

# get arguments
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

projPath = args[1]
outDir = args[2]

# set working directory
setwd(projPath)

# read in config.json
config = jsonlite::fromJSON(paste0(projPath, "/config.json"), simplifyDataFrame = FALSE)

# read in alignment result
alignResult= read.table(paste0(projPath, "/alignment/alignment_summary.txt"), header = TRUE, fill = TRUE)

# Extract keys that start with "SAMPLE"
sample_keys <- names(config)[grepl("^SAMPLE", names(config))]
samples <- lapply(sample_keys, function(key) config[[key]])
names(samples) <- sample_keys

peakN = c()
peakWidth = c()
histList = c()
sampleList = c()

# get number of peaks called and peak width
for (sample_name in names(samples)) {
  sample <- samples[[sample_name]]
  # append histName to histList
  histList = c(histList, sample$histName)
  # get sample name
  samplename <- paste0(sample$histName, "_rep", as.character(sample$rep))
  sampleList = c(sampleList, samplename)
  if(as.character(sample$rep) != "IgG"){
    #for(type in peakType){
      peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/", samplename, "_seacr_top0.01.peaks.stringent.bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peakN = data.frame(peakN = nrow(peakInfo), peakType = "top0.01", Histone = sample$histName, Replicate = sample$rep) %>% rbind(peakN, .)
      peakWidth = data.frame(width = peakInfo$width, peakType = "top0.01", Histone = sample$histName, Replicate = sample$rep)  %>% rbind(peakWidth, .)
    #}
  }
}
peakN %>% select(Histone, Replicate, peakType, peakN)

# Reproducibility of the peak across biological replicates
#repL = paste0("rep", 1:2)
peakType = "top0.01" #c("control", "top0.01")
peakOverlap = c()
for(type in peakType){
  for (sample_name in names(samples)) {
    sample <- samples[[sample_name]]
    # get sample name
    samplename <- paste0(sample$histName, "_rep", as.character(sample$rep))
    overlap.gr = GRanges()
      peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/", samplename,"_seacr_", type, ".peaks.stringent.bed"), header = FALSE, fill = TRUE)
      peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
      if(length(overlap.gr) >0){
        overlap.gr = overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
      }else{
        overlap.gr = peakInfo.gr
      }
    peakOverlap = data.frame(peakReprod = length(overlap.gr), Histone = sample$histName, Replicate = sample$rep, peakType = type) %>% rbind(peakOverlap, .)
  }
}
peakReprod = left_join(peakN, peakOverlap, by = c("Histone", "peakType","Replicate")) %>% mutate(peakReprodRate = peakReprod/peakN * 100)
peakReprod %>% select(Histone, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)

# Calculate FRIPs: FRagment proportion in Peaks regions, as a measure of signal-to-noise
bamDir = paste0(projPath, "/alignment/bam")
inPeakData = c()
## overlap with bam file to get count
for (sample_name in names(samples)) {
    sample <- samples[[sample_name]]
    # get sample name
    samplename <- paste0(sample$histName, "_rep", as.character(sample$rep))
    peakRes = read.table(paste0(projPath, "/peakCalling/SEACR/", samplename, "_seacr_top0.01.peaks.stringent.bed"), header = FALSE, fill = TRUE)
    peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
    bamFile = paste0(bamDir, "/", samplename, "_bowtie2.mapped.bam")
    fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
    inPeakN = counts(fragment_counts)[,1] %>% sum
    inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN, Histone = sample$histName, Replicate = sample$rep))
}
frip = left_join(inPeakData, alignResult, by = c("Histone", "Replicate")) %>% mutate(frip = inPeakN/MappedFragNum_axl * 100)
frip %>% select(Histone, Replicate, SequencingDepth, MappedFragNum_axl, AlignmentRate_axl, FragInPeakNum = inPeakN, FRiPs = frip)

# visualize
figA = peakN %>% ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate)) + # , position = position_jitter(0.05)
    facet_grid(~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis() + #discrete = TRUE, begin = 0.1, end = 0.9
    theme_bw(base_size = 18) +
    ylab("Number of Peaks") +
    xlab("")

figB = peakWidth %>% ggplot(aes(x = Histone, y = width, fill = Histone)) +
    geom_violin() +
    facet_grid(Replicate~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis() + #discrete = TRUE, begin = 0.1, end = 0.9
    scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
    theme_bw(base_size = 18) +
    ylab("Width of Peaks") +
    xlab("")

figC = peakReprod %>% ggplot(aes(x = Histone, y = peakReprodRate, fill = Histone, label = round(peakReprodRate, 2))) +
    geom_bar(stat = "identity") +
    geom_text(vjust = 0.1) +
    facet_grid(Replicate~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis() + #discrete = TRUE, begin = 0.1, end = 0.9
    theme_bw(base_size = 18) +
    ylab("% of Peaks Reproduced") +
    xlab("")

figD = frip %>% ggplot(aes(x = Histone, y = frip, fill = Histone, label = round(frip, 2))) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate)) + # , position = position_jitter(0.05)
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis() + # discrete = TRUE, begin = 0.1, end = 0.9
    theme_bw(base_size = 18) +
    ylab("% of Fragments in Peaks") +
    xlab("")

pdf(file= paste0(projPath, "/peakCalling/peak_plots.pdf"), width=7, height=7)
ggarrange(figA, figB, figC, figD, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

writeLines(capture.output(sessionInfo()), paste0(outDir,"/pkgs_log/peaksum_sessionInfo.txt"))