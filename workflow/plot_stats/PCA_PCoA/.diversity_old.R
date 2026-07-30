# 11/05/2026
# Analysis of the microbial diversity of 212 Kefirs via α and β indice. 
# The α and β diversity were calculated via a Python script :
# (https://github.com/jenniferlu717/KrakenTools/tree/master/DiversityTools)
# Out : Barplots and PCoA Plots

# Part 1: Importing the necessary packages and data -----------------------
library(dplyr)
library(readr)
library(tidyr)
library(mixOmics)
library(vegan) 
library(ggplot2)
library(ggforce)
library(ggrepel)

alpha <- read.csv("alpha.csv.csv",dec = ",", sep = ",", header = TRUE)

beta <- read.csv(
  "betadiv.csv",dec = ",", sep = ",", header = TRUE)

# 1 line = number of ARGs classified as 'Only for Humans' by samples.
beta.env <- read.csv(
  "count_arg_dfs.csv",dec = ",", sep = ",", header = TRUE)

# renamed file 'species_relative_abundance.csv', obtained in taxonomy/abundancy_csv
tab <- read.csv(
  "intrasamples_abund_relativ.csv",dec = ",", sep = ",", header = TRUE) 

rownames(tab) <- tab[,1]
tab <- tab[,-1]
tab <- tab[, colnames(tab) != "Total"]

# Part 2: α diversity in terms of Kefir microbial composition --------------

# Formatting data 
alpha$Shannon <- gsub(",", ".", alpha$Shannon)     
alpha$Shannon[alpha$Shannon %in% c("NA", "N/A")] <- NA 
alpha$Shannon <- as.numeric(alpha$Shannon)
alpha_clean <- alpha[!is.na(alpha$Shannon), ]
nrow(alpha_clean)
alpha_clean$Colocalization <- as.factor(alpha_clean$Colocalization)
alpha_clean$Bioproject <- beta.env$bioproject[match(alpha_clean$sample,beta.env$sample)]

# Plot
par(mar = c(4, 4, 1, 4))
boxplot(Shannon ~ Colocalization, data = alpha_clean)
title("Comparing Alpha Diversity by abundance of ARG colocated with at least 1 MGE.
      \n n = 212, All types of ARG are included")
dev.off()

unique(alpha_clean$Bioproject)
color_biop <- c("PRJEB65292" = "red", 
                "PRJEB77409" = "steelblue",
                "PRJNA388572" = "green",
                "PRJNA704713" = "purple")

par(mar = c(4, 4, 1, 4))
barplot(Shannon ~ sample, data = alpha_clean, col = color_biop[alpha_clean$Bioproject])
title("\n Alpha Diversity by samples\n n = 212, All types of ARG are included")
dev.off()

alpha_clean$Colocalization <- as.numeric(alpha_clean$Colocalization)
par(mar = c(4, 4, 1, 4))
barplot(Colocalization ~ sample, data = alpha_clean, 
        col = color_biop[alpha_clean$Bioproject])
title("\n Count of ARGs colocated with an MGE by samples\n n = 212, 
      All types of ARG are included")
dev.off()

# Part 2:  β diversity in terms of Kefir microbial composition --------------

# Formatting dissimilarity matrix
rownames(beta) <- beta$x
beta <- beta[, colnames(beta) != "x"]
beta <- beta[, order(colnames(beta))]
beta <- beta[order(rownames(beta)), ]

dist_mat <- as.matrix(beta)
mode(dist_mat) <- "numeric"
n = nrow(dist_mat)

# PCoA
pcoa_result <- cmdscale(dist_mat, eig = TRUE, k = 2)
points <- as.data.frame(pcoa_result$points)
colnames(points) <- c("PCoA1", "PCoA2")
variance <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 2)

# Eigenvalues fall
barplot(pcoa_result$eig)
title(main = "Eigenvalues fall for β diversity \n in terms of Kefir microbial composition")

# Variables projection
# Fits environmental vectors (tab) onto an ordination(pcoa_result$points)
# To check "strong" and "weak" predictors : View(fit$vectors$r)
fit <- envfit(pcoa_result$points, tab, perm = 0) 
species_coords <- as.data.frame(scores(fit, display = "vectors"))
species_coords$species <- rownames(species_coords)
# Add transparency to vectors with weak quality representation
species_coords$length <- sqrt(species_coords$Dim1^2 + species_coords$Dim2^2)
species_coords$quality <- species_coords$length / max(species_coords$length)

# Grouped by Bioproject
points$bioproject <- beta.env$bioproject[match(rownames(points),beta.env$sample)]

# Plot with variables
ggplot(points, aes(x = PCoA1, y = PCoA2, color = bioproject)) +
  geom_point(size = 3, alpha = 0.8) +
  # Arrows for species
  geom_segment(data = species_coords,
               aes(x = 0, y = 0, xend = Dim1, yend = Dim2, alpha = quality),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black",
               inherit.aes = FALSE) +
  # Species labels
  geom_text_repel(data = species_coords,
                  aes(x = Dim1, y = Dim2, label = species),
                  size = 3,
                  color = "black",
                  inherit.aes = FALSE) +
  scale_color_manual(values = color_biop) +
  labs(
    title = paste("PCoA (Bray-Curtis) – Kefir microbial composition"),
    subtitle = paste(n, "samples"),
    x = paste0("PCoA1 (", variance[1], "%)"),
    y = paste0("PCoA2 (", variance[2], "%)")
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

ggplot(points, aes(x = PCoA1, y = PCoA2, color = bioproject))+
 geom_point(size = 2.7) +
 scale_color_manual(values = color_biop) +
 geom_text_repel(aes(label = rownames(points)), size = 1.2, color = "black") +
 labs(
 title = paste("PCoA Analysis (Bray-Curtis)– Kefir microbial composition" ),
 subtitle = paste(n, "samples"),
 x = paste0("PCoA1 (", variance[1], "%)"),
 y = paste0("PCoA2 (", variance[2], "%)")
 ) +
 theme_minimal()

# Grouped by tFermentation
points$tFermentation<- beta.env$tFermentation[match(rownames(points),beta.env$sample)]

ggplot(points, aes(x = PCoA1, y = PCoA2, color = tFermentation))+
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linetype = 2, alpha = 0.3) +
  labs(
    title = paste("PCoA Analysis (Bray-Curtis)\n", n, " samples"),
    x = paste0("PCoA1 (", variance[1], "%)"),
    y = paste0("PCoA2 (", variance[2], "%)")
  ) +
  theme_minimal()

# Grouped by ARGs Count
points$count<- beta.env$count[match(rownames(points),beta.env$sample)]
c = sum(!is.na(points$count))
points$count <- as.character(points$count)
points$count[is.na(points$count)] <- "NA"

ggplot(points, aes(x = PCoA1, y = PCoA2, shape = count))+
  geom_point(size = 3) +
  labs(
    title = paste("PCoA Analysis (Bray-Curtis), Count of ARGs colocated with at least 1 MGE by Samples"),
    subtitle = paste(n, " samples", c, " transferable ARGs"),
    x = paste0("PCoA1 (", variance[1], "%)"),
    y = paste0("PCoA2 (", variance[2], "%)")
  ) +
  theme_minimal()

# Plot with colors by Bioproject, names of the samples, shape of points by count of transferable ARG

ggplot(points, aes(x = PCoA1, y = PCoA2, color = bioproject, shape = count)) +
  geom_point(size = 2.7, alpha = 0.8) +
  # Arrows for species
  geom_segment(data = species_coords,
               aes(x = 0, y = 0, xend = Dim1, yend = Dim2, alpha = quality),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black",
               inherit.aes = FALSE) +
  # Species labels
  geom_text_repel(data = species_coords,
                  aes(x = Dim1, y = Dim2, label = species),
                  size = 2.7,
                  color = "black",
                  inherit.aes = FALSE) +
  scale_color_manual(values = color_biop) +
  labs(
    title = paste("PCoA Bray-Curtis of Kefir microbial composition,\nCount of ARGs colocated with at least 1 MGE by samples"),
    subtitle = paste(n, " samples", c, " transferable ARGs"),
    x = paste0("PCoA1 (", variance[1], "%)"),
    y = paste0("PCoA2 (", variance[2], "%)")
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

