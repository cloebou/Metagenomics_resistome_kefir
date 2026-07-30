# 17/05/2026
# Barplot of the number of reads per sample 

# Part 1: Importing the necessary packages and data -----------------------
library(readr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(tibble)

reads <- read.csv("read.csv",dec = ".", sep = ",", header = FALSE)
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

par(mar = c(4, 4, 2, 2)) 
barplot(
  reads$V2,
  names.arg = reads$V1,
  col = bar_colors,
  las = 3, # labels verticaux
  main = paste("Number of total reads per samples\n", n,"samples"),
  sub = "Metagenome samples from kefir products",
  ylab = "Number of reads",
  cex.names = 0.5
)

legend(
  "bottom",
  inset = c(0, -0.6),
  legend = names(color_map),
  fill = color_map,
  title = "Bioproject",
  cex = 0.4,
  ncol = 3,        # 🔑 MULTI colonnes au lieu de 1 ligne infinie
  bty = "n",
  xpd = TRUE
)

# Saving
pdf("barplot_reads.pdf", width = 12, height = 4)
layout(matrix(c(1,2), ncol = 1), heights = c(3,1.5))
