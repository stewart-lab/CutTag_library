# load packages
library(reticulate)
use_condaenv("/w5home/bmoore/miniconda3/envs/cut_tag2")
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
print(names(samples))
# get alignment result vector
alignResult = c()
histList = c()
replist = c()
# Use the extracted sample information
for (sample_name in names(samples)) {
  sample <- samples[[sample_name]]
  cat(paste("Information for", sample_name, ":\n"))
  cat(paste("  FASTQ PE1:", sample$fastq_PE1, "\n"))
  cat(paste("  FASTQ PE2:", sample$fastq_PE2, "\n"))
  cat(paste("  Histone Name:", sample$histName, "\n"))
  # append histName to histList
  histList = c(histList, sample$histName)
  replist = c(replist, paste0(sample$histName,sample$rep))
  cat(paste("  Replicate:", sample$rep, "\n\n"))
  alignRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", 
            sample$histName, "Rep", as.character(sample$rep),"_bowtie2.txt"), 
            header = FALSE, fill = TRUE) # skip 1 to skip no hup line
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  print(alignRes)
  print(alignRate)
  alignResult = data.frame(Histone = sample$histName, Replicate = sample$rep, 
                  SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
                  MappedFragNum_axl = alignRes$V1[4] %>% as.character %>% as.numeric + 
                  alignRes$V1[5] %>% as.character %>% as.numeric, 
                  AlignmentRate_axl = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}
print(alignResult)
#alignResult$Histone = factor(alignResult$Histone, levels = histList)
alignResult %>% mutate(AlignmentRate_axl = paste0(AlignmentRate_axl, "%"))
#print(alignResult)
# get spike in alignment
spikeAlign = c()
# for (sample_name in names(samples)) {
#   sample <- samples[[sample_name]]
#   spikeRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", sample$histName, "_rep", 
#                as.character(sample$rep), "_bowtie2_spikeIn.txt"), skip=1, header = FALSE, fill = TRUE)
#   alignRate = substr(spikeRes$V1[6], 1, nchar(as.character(spikeRes$V1[6]))-1)
#   spikeAlign = data.frame(Histone = sample$histName, Replicate = sample$rep, 
#                     SequencingDepth = spikeRes$V1[1] %>% as.character %>% as.numeric, 
#                     MappedFragNum_spikeIn = spikeRes$V1[4] %>% as.character %>% as.numeric + 
#                     spikeRes$V1[5] %>% as.character %>% as.numeric, 
#                     AlignmentRate_spikeIn = alignRate %>% as.numeric)  %>% rbind(spikeAlign, .)
# }
# spikeAlign$Histone = factor(spikeAlign$Histone, levels = histList)
# spikeAlign %>% mutate(AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))
# # join the two alignment results
# alignSummary = left_join(alignResult, spikeAlign, by = c("Histone", "Replicate", "SequencingDepth")) %>%
#   mutate(AlignmentRate_axl = paste0(AlignmentRate_axl, "%"), 
#          AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))
alignSummary <- alignResult
# write alignment summary
write.table(alignSummary, file = paste0(projPath, "/alignment/alignment_summary.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

## Generate sequencing depth boxplot
figA = alignResult %>% ggplot(aes(x = Histone, y = SequencingDepth/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis() + #discrete = TRUE, begin = 0.1, end = 0.9
    theme_bw(base_size = 18) +
    ylab("Sequencing Depth per Million") +
    xlab("") + 
    ggtitle("A. Sequencing Depth")

figB = alignResult %>% ggplot(aes(x = Histone, y = MappedFragNum_axl/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis() + #discrete = TRUE, begin = 0.1, end = 0.9
    theme_bw(base_size = 18) +
    ylab("Mapped Fragments per Million") +
    xlab("") +
    ggtitle("B. Alignable Fragment (Axl)")

figC = alignResult %>% ggplot(aes(x = Histone, y = AlignmentRate_axl, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis() + #discrete = TRUE, begin = 0.1, end = 0.9
    theme_bw(base_size = 18) +
    ylab("% of Mapped Fragments") +
    xlab("") +
    ggtitle("C. Alignment Rate (Axl)")

# figD = spikeAlign %>% ggplot(aes(x = Histone, y = AlignmentRate_spikeIn, fill = Histone)) +
#     geom_boxplot() +
#     geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
#     scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
#     scale_color_viridis() + # discrete = TRUE, begin = 0.1, end = 0.9
#     theme_bw(base_size = 18) +
#     ylab("Spike-in Alignment Rate") +
#     xlab("") +
#     ggtitle("D. Alignment Rate (E.coli)")

pdf(file= paste0(projPath, "/alignment/alignment_summary.pdf"))
ggarrange(figA, figB, figC, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

writeLines(capture.output(sessionInfo()), paste0(outDir,"/pkgs_log/alignsum_sessionInfo.txt"))