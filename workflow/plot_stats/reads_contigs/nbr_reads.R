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

# Information on data
taille <-read.csv("sample_nbreads_length.csv",dec = ".", sep = ",", header = FALSE)
taille$bioproject <- taille$V2
taille$bioproject <- reads.env$bioproject[match(taille$V1,reads.env$sample)]
taille[taille$V1=="SRR13801805", "bioproject"] <- "PRJNA704713"
grand=taille[taille$bioproject=="PRJEB65292",]
mean(grand$V3)
petit=taille[taille$bioproject=="PRJEB77409",]
mean(petit$V3)
tpetit=taille[taille$bioproject=="PRJNA704713",]
mean(tpetit$V3)
ttpetit=taille[taille$bioproject=="PRJNA388572",]
mean(ttpetit$V3)
sum(reads$V2)

# Part 3:  coloring by Bioproject ------------------------------------------
color_map <- c(
  "PRJEB65292" = "red",
  "PRJEB77409" = "steelblue",
  "PRJNA388572" = "green",
  "PRJNA704713" = "purple"
)

# Assign colors to each bar
bar_colors <- color_map[reads$bioproject]
n = as.character(nrow(reads))

# save plot
pdf("barplot_reads.pdf", width = 12, height = 4)
layout(matrix(c(1,2), ncol = 1), heights = c(3,1.5))

# Part 4:  Plotting --------------------------------------------------------
par(mar = c(1, 4, 6, 0.05))
barplot(
  reads$V2,
  names.arg = reads$V1,
  col = bar_colors,
  las = 3,
  main = paste("Number of total reads per samples\n", n, "samples"),
  sub = "Metagenome samples from kefir products",
  ylab = "Number of reads",
  cex.names = 0.2
)

# Legende
par(mar = c(0, 2, 0, 2))
plot.new()

legend(
  "center",
  legend = names(color_map),
  fill = color_map,
  title = "Bioproject",
  cex = 0.5,
  ncol = 4, 
  bty = "n"
)

