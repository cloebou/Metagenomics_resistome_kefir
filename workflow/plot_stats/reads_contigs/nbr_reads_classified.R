# 17/05/2026
# Barplot of the numbers of classified reads by Kraken generated
# with the 'taxonomy workflow'

# Part 1: Importing the necessary packages and data -----------------------

library(readr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(tibble)

reads <- read.csv("read_classified_kraken.csv",dec = ".", sep = ",", header = FALSE)
reads.env <- read.csv("sample_collection.csv",dec = ",", sep = ";", header = TRUE)

# Part 2: rearranging data ------------------------------------------------

reads$bioproject <- reads$V2
reads$bioproject <- reads.env$bioproject[match(reads$V1,reads.env$sample)]
reads[reads$V1=="SRR13801805", "bioproject"] <- "PRJNA704713"

# Part 3:  coloring by Bioproject ------------------------------------------

color_map <- c(
  "PRJEB65292" = "red",
  "PRJEB77409" = "steelblue",
  "PRJNA388572" = "green",
  "PRJNA704713" = "purple"
)

bar_colors <- color_map[reads$bioproject]
n = as.character(nrow(reads))

# Part 4:  Plotting --------------------------------------------------------
barplot(
  reads$V2,
  names.arg = reads$V1,
  col = bar_colors,
  las = 3, # labels verticaux
  main = paste("Number of classified reads by Kraken2 on the PlusPF-16 database,\n", n,"samples"),
  sub = "Metagenome samples from kefir products",
  ylab = "Number of reads",
  cex.names = 0.5
)

legend(
  "topleft",
  legend = names(color_map),
  fill = color_map,
  title = "Bioproject",
  cex = 0.5
)
dev.off()

# Saving
pdf("barplot_reads_classified.pdf", width = 12, height = 4)
layout(matrix(c(1,2), ncol = 1), heights = c(3,1.5))


