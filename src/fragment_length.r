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

# set working directory
setwd(projPath)

# read in config.json
config = jsonlite::fromJSON(paste0(projPath, "/config.json"), simplifyDataFrame = FALSE)

# Extract keys that start with "SAMPLE"
sample_keys <- names(config)[grepl("^SAMPLE", names(config))]
samples <- lapply(sample_keys, function(key) config[[key]])
names(samples) <- sample_keys

## Collect the fragment size information
fragLen = c()
histList = c()
sampleList = c()
# Use the extracted sample information
for (sample_name in names(samples)) {
  sample <- samples[[sample_name]]
  cat(paste("Information for", sample_name, ":\n"))
  cat(paste("  FASTQ PE1:", sample$fastq_PE1, "\n"))
  cat(paste("  FASTQ PE2:", sample$fastq_PE2, "\n"))
  cat(paste("  Histone Name:", sample$histName, "\n"))
  cat(paste("  Replicate:", sample$rep, "\n\n"))
  # append histName to histList
  histList = c(histList, sample$histName)
  samplename <- paste0(sample$histName, "_rep", as.character(sample$rep))
  sampleList = c(sampleList, samplename)
  fragLen = read.table(paste0(projPath, "/alignment/sam/fragmentLen/", samplename,"_fragmentLen.txt"), header = FALSE) %>% 
  mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Histone = sample$histName, 
  Replicate = sample$rep, sampleInfo = samplename) %>% rbind(fragLen, .) 
}
fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone = factor(fragLen$Histone, levels = histList)

## Generate the fragment size density plot (violin plot)
figA = fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
    geom_violin(bw = 5) +
    scale_y_continuous(breaks = seq(0, 800, 50)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis() + # discrete = TRUE, begin = 0.1, end = 0.9
    theme_bw(base_size = 14) +
    ggpubr::rotate_x_text(angle = 20) +
    ylab("Fragment Length") +
    xlab("")

figB = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo)) + #linetype = Replicate)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 14) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

pdf(file= paste0(projPath, "/alignment/sam/fragmentLen/fragment_length.pdf"), width=10, height=7)
ggarrange(figA, figB, ncol = 2)
dev.off()

