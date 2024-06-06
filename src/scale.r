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
library(purrr)
library(jsonlite)
library(tidyjson)

# get arguments
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

projPath = args[1]
outDir = args[2]

# set working directory
setwd(projPath)

# read in config.json
config = jsonlite::fromJSON(paste0(projPath, "/config.json"), simplifyDataFrame = FALSE)

# Extract keys that start with "SAMPLE"
sample_keys <- names(config)[grepl("^SAMPLE", names(config))]
samples <- lapply(sample_keys, function(key) config[[key]])
names(samples) <- sample_keys
# multiplier
multiplier = 10000

# read in alignment summary
alignSummary = read.table(paste0(projPath, "/alignment/alignment_summary.txt"), header = TRUE, fill = TRUE)

# scale factor
scaleFactor = c()
histList = c()
sampleList = c()

for (sample_name in names(samples)) {
  sample <- samples[[sample_name]]
  # append histName to histList
  histList = c(histList, sample$histName)
  # get sample name
  samplename <- paste0(sample$histName, "_rep", as.character(sample$rep))
  sampleList = c(sampleList, samplename)
  spikeDepth = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", samplename, "_bowtie2_spikeIn.seqDepth"), header = FALSE, fill = TRUE)$V1[1]
  scaleFactor = data.frame(scaleFactor = multiplier/spikeDepth, Histone = sample$histName, Replicate = sample$rep)  %>% rbind(scaleFactor, .)
}
scaleFactor$Histone = factor(scaleFactor$Histone, levels = histList)
#alignSummary["Replicate"]= as.character(alignSummary["Replicate"])
left_join(alignSummary, scaleFactor, by = c("Histone", "Replicate"))
# write out alignment summary
write.table(alignSummary, file = paste0(projPath, "/alignment/alignment_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


## Generate sequencing depth boxplot
figA = scaleFactor %>% ggplot(aes(x = Histone, y = scaleFactor, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis() + #discrete = TRUE, begin = 0.1, end = 0.9
    theme_bw(base_size = 20) +
    ylab("Spike-in Scalling Factor") +
    xlab("")

normDepth = inner_join(scaleFactor, alignSummary, by = c("Histone", "Replicate")) %>% mutate(normDepth = MappedFragNum_axl * scaleFactor)

figB = normDepth %>% ggplot(aes(x = Histone, y = normDepth, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis() + # discrete = TRUE, begin = 0.1, end = 0.9
    theme_bw(base_size = 20) +
    ylab("Normalization Fragment Count") +
    xlab("") + 
    coord_cartesian(ylim = c(1000000, 130000000))

pdf(file= paste0(projPath, "/alignment/scale_boxplot.pdf"), width=7, height=7)
ggarrange(figA, figB, ncol = 2, common.legend = TRUE, legend="bottom")
dev.off()

writeLines(capture.output(sessionInfo()), paste0(outDir,"/pkgs_log/scale_sessionInfo.txt"))