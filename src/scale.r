# load libraries
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
# multiplier
multiplier = 10000

# read in alignment summary
alignSummary = read.table(paste0(projPath, "/alignment/alignment_summary.txt"), header = TRUE, fill = TRUE)

# scale factor
scaleFactor = c()
histList = c()
sampleList = c()
# initialize dataframe
normDepthDF <- data.frame(normDepth = numeric(), Histone = character(), Replicate = numeric(), stringsAsFactors = FALSE)
scaleFactorDF <- data.frame(scaleFactor = numeric(), Histone = character(), Replicate = numeric(), stringsAsFactors = FALSE)

for (sample_name in names(samples)) {
  sample <- samples[[sample_name]]
  # append histName to histList
  histList = c(histList, sample$histName)
  subdf <- alignSummary[alignSummary$Histone == sample$histName, ]
  if (nrow(subdf) == 0) next
  # get sample name
  samplename <- paste0(sample$histName, "Rep", as.character(sample$rep))
  sampleList = c(sampleList, samplename)
   # Read spikeDepth
  spikeDepthFile <- paste0(projPath, "/alignment/sam/bowtie2_summary/", samplename, "_bowtie2_spikeIn.seqDepth")
  if (file.exists(spikeDepthFile)) {
    spikeDepth <- read.table(spikeDepthFile, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)$V1[1]
  } else {
    spikeDepth <- NA
    warning(paste("File not found:", spikeDepthFile))
  }
  # Skip if spikeDepth is NA
  if (is.na(spikeDepth)) next
  # Calculate normalized depth
  sF <- multiplier/spikeDepth
  normDepth <- subdf$MappedFragNum_axl * sF
  print(normDepth)
  # Skip if normDepth is empty
  if (length(normDepth) == 0) next
  #scaleFactor = data.frame(scaleFactor = multiplier/spikeDepth, Histone = sample$histName, Replicate = sample$rep)  %>% rbind(scaleFactor, .)
  #normDepth = data.frame(normDepth = normDepth,  Histone = sample$histName, Replicate = sample$rep)  %>% rbind(normDepth, .)
  # Create a temporary data frame to store current iteration results
  tempScaleFactor <- data.frame(scaleFactor = multiplier / spikeDepth, Histone = sample$histName, Replicate = sample$rep, stringsAsFactors = FALSE)
  tempNormDepth <- data.frame(normDepth = normDepth, Histone = sample$histName, Replicate = sample$rep, stringsAsFactors = FALSE)
  
  # Append the temporary data frame to the main data frame
  scaleFactorDF <- rbind(scaleFactorDF, tempScaleFactor)
  normDepthDF <- rbind(normDepthDF, tempNormDepth)
}
# Print final results
print("Final NormDepth DataFrame:")
print(normDepthDF)

print("Final ScaleFactor DataFrame:")
print(scaleFactorDF)

#scaleFactorDF$Histone = factor(scaleFactorDF$Histone, levels = histList)
#normDepthDF$Histone = factor(normDepthDF$Histone, levels = histList)
#alignSummary["Replicate"]= as.character(alignSummary["Replicate"])
alignSummary <- left_join(alignSummary, scaleFactorDF, by = c("Histone", "Replicate"))
alignSummary <- left_join(alignSummary, normDepthDF, by = c("Histone", "Replicate"))
#normDepth = inner_join(scaleFactor, alignSummary, by = c("Histone", "Replicate")) %>% mutate(normDepth = MappedFragNum_axl * scaleFactor)
# write out alignment summary
write.table(alignSummary, file = paste0(projPath, "/alignment/alignment_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


## Generate sequencing depth boxplot
figA = scaleFactorDF %>% ggplot(aes(x = Histone, y = scaleFactor, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis() + #discrete = TRUE, begin = 0.1, end = 0.9
    theme_bw(base_size = 20) +
    ylab("Spike-in Scalling Factor") +
    xlab("")



figB = normDepthDF %>% ggplot(aes(x = Histone, y = normDepth, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis() + # discrete = TRUE, begin = 0.1, end = 0.9
    theme_bw(base_size = 20) +
    ylab("Normalization Fragment Count") +
    xlab("")
    #+ coord_cartesian(ylim = c(1000000, 130000000))

pdf(file= paste0(projPath, "/alignment/scale_boxplot.pdf"), width=7, height=7)
ggarrange(figA, figB, ncol = 2, common.legend = TRUE, legend="bottom")
dev.off()

writeLines(capture.output(sessionInfo()), paste0(outDir,"scale_sessionInfo.txt"))