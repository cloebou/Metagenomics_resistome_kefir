# 23/07/2026
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
library(vegan)

# Metadata import
taxo <-  read.csv("intrasamples_abund_relativ.csv",dec = ",", sep = ",", header = TRUE)
taxo_norm <-  read.csv("intrasample_abund.csv",dec = ",", sep = ",", header = TRUE)
colnames(taxo)
taxo = c("Lactococcus.cremoris", "Lactobacillus.kefiranofaciens","Lactobacillus.helveticus","Lactococcus.lactis","Lentilactobacillus.kefiri",
         "Leuconostoc.mesenteroides","Pseudolactococcus.raffinolactis","Acinetobacter.sp..TTH0.4","Pseudomonas.helleri","Acetobacter.orientalis",
         "Other..0.005", "Streptococcus.thermophilus","Pseudomonas.lundensis", "Macrococcoides.caseolyticum","Hafnia.paralvei",
         "Pseudomonas.brenneri","Lactococcus.sp","Streptococcus.parauberis","Thermus.thermophilus","Leuconostoc.falkenbergense",
         "Pseudolactococcus.piscium","Acetobacter.ghanensis","Pseudomonas.saxonica", "Carnobacterium.maltaromaticum", "Kluyvera.intermedia",
         "Saccharomyces.cerevisiae","Kluyveromyces.marxianus","Lelliottia.amnigena","Pseudomonas.putida","Klebsiella.sp..HC6","Bacillus.mycoides",
         "Pantoea.agglomerans","Clostridium.saccharobutylicum","Pseudolactococcus.carnosus","Pseudomonas.lurida","Paenibacillus.sp..FSL.H8.0332",
         "Pseudomonas.sp..B21128","Propionibacterium.freudenreichii","Klebsiella.grimontii","Yersinia.alsatica","Duffyella.gerundensis","Clostridium.beijerinckii",
         "Pseudolactococcus.paracarnosus","Yarrowia.lipolytica", "Pseudomonas.azotoformans","Clostridium.tyrobutyricum", "Pseudomonas.sp..FP453",
         "Pseudomonas.moraviensis", "Acinetobacter.towneri","Leclercia.adecarboxylata","Rahnella.aceris", "Pseudomonas.sp..S3_B08","Kluyvera.cryocrescens",
         "Pseudomonas.rustica","Superficieibacter.sp..BNK.5","Pseudomonas.sp..SD17.1","Janthinobacterium.sp..J1.1", "Maaswegvirus.Kp24","Staphylococcus.sp..IVB6181",
         "Citrobacter.freundii.complex.sp..CFNIH3", "Macrococcus.psychrotolerans", "Lactococcus.phage.63301","Bacillus.sp..WC2507","Enterobacter.sp..RHBSTW.00901",
         "Pseudomonas.sp..SC3.2021","Bordetella.bronchiseptica")

# Importing data

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

range(alpha_clean$Shannon)

#---Stats Alpha---
# 2 categories : bioproject and resistance

# Category 1 : bioproject 
# 3 groups big: "PRJEB65292"  middle:"PRJEB77409"  low:"PRJNA704713", rm low 1 sample
df <- alpha_clean[,c("sample","Colocalization","Shannon","Bioproject")]
unique(df$Bioproject)
big=df[df$Bioproject=="PRJEB65292",]
middle=df[df$Bioproject=="PRJEB77409",]
low=df[df$Bioproject=="PRJNA704713",]

color_map <- c(
  "PRJEB65292" = "red",
  "PRJEB77409" = "steelblue"
)
boxplot(Shannon ~ Bioproject, data = alpha_clean[alpha_clean$Bioproject!="PRJNA704713",],col=color_map)
# outlier : 	ERR11865364, Shannon 2.2092480

# Test Gaussian Variable 
shapiro.test(big$Shannon) #p-value 0.29>0.05 => do NOT reject H0 (gaussienne)
qqnorm(big$Shannon)
hist(big$Shannon)
shapiro.test(middle$Shannon) #p-value 0.88>0.05 => do NOT reject H0 
qqnorm(middle$Shannon)
hist(middle$Shannon)

# Equal variance test
var.test(big$Shannon,middle$Shannon) #p-value  5.13e-08 => reject H0 (so do not have equal variances)

# But they have more then 30 samples and Student test is robust
length(big$Shannon) # > 30
length(middle$Shannon) # >30

student <- t.test(big$Shannon,middle$Shannon,var.equal=TRUE) #p-value = 7.043e-13 => Reject H0  (so do not have equality of the same expected value)

stat.test <- data.frame(
  group1 = "PRJEB65292",
  group2 = "PRJEB77409",
  p = student$p.value,
  y.position = max(alpha_clean$Shannon, na.rm = TRUE) + 0.2
)
stat.test$label <- paste0("p = ", signif(stat.test$p, 3)) # Keep only 3 decimals

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  )

# Boxplot with Stats
p <- ggplot(alpha_clean[alpha_clean$Bioproject!="PRJNA704713",],aes(x = Bioproject, y = Shannon, color = Bioproject)) +
  geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2) +
  scale_color_manual(values = color_map) +
  stat_pvalue_manual(stat.test, label = "label") +
  labs(
    x = "Bioproject",
    y = "Shannon Diversity Index",
    title = "Alpha Diversity (Shannon Index)"
  ) +
  My_Theme+
  theme(legend.position="non")

p

# Category 2 : resistance 
dfr_uncl <- alpha_clean[,c("sample","Colocalization","Shannon","Bioproject")]
dfr <- dfr_uncl %>%
  mutate(Colocalization = ifelse(Colocalization != 0, "resistant", "not_resistant"))
res=dfr[dfr$Colocalization=="resistant",]
notres=dfr[dfr$Colocalization=="not_resistant",]

color_map <- c(
  "resistant" = "orange",
  "not_resistant" = "darkgreen"
)
boxplot(Shannon ~ Colocalization, data = dfr,col=color_map)
# outlier : 	ERR11865364, Shannon 2.2092480

# Test Gaussian Variable 
shapiro.test(res$Shannon) #p-value 0.0006<0.05 => reject H0 (gaussienne)
qqnorm(res$Shannon)
hist(res$Shannon)
shapiro.test(notres$Shannon) #p-value 0.01176<0.05 => reject H0 
qqnorm(notres$Shannon)
hist(notres$Shannon)

# Equal variance test
var.test(res$Shannon,notres$Shannon) #p-value  0.4759=> reject H0 (so do not have equal variances)

# But they have more then 30 samples and Student test is robust
length(res$Shannon) # 193 > 30
length(notres$Shannon) # 19 < 30

will <- wilcox.test(res$Shannon,notres$Shannon) #p-value = 0.1492>0.05 => No NOT Reject H0  (so have equality of the distribution functions of these two random variables.)

stat.test <- data.frame(
  group1 = "resistant",
  group2 = "not_resistant",
  p = will$p.value,
  y.position = max(dfr$Shannon, na.rm = TRUE) + 0.2
)
stat.test$label <- paste0("p = ", signif(stat.test$p, 3)) # Keep only 3 decimels

p <- ggplot(dfr,aes(x = Colocalization, y = Shannon, color = Colocalization)) +
  geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2) +
  scale_color_manual(values = color_map) +
  stat_pvalue_manual(stat.test, label = "label") +
  labs(
    x = "Resistance",
    y = "Shannon Diversity Index",
    title = "Alpha Diversity (Shannon Index)"
  ) +
  theme_classic()+
  My_Theme+
  theme(legend.position="non")

p

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


# Grouped by Bioproject, adapted for rapport
ggplot(points, aes(x = PCoA1, y = PCoA2, color = bioproject))+
 geom_point(size = 2.7) +
 scale_color_manual(values = color_biop) +
# geom_text_repel(aes(label = rownames(points)), size = 1.2, color = "black") +
 labs(
 title = paste("PCoA of Bray-Curtis distance" ),
 x = paste0("PCoA1 (", variance[1], "%)"),
 y = paste0("PCoA2 (", variance[2], "%)")
 ) +
# Arrows for species
 #geom_segment(data = species_coords,
             # aes(x = 0, y = 0, xend = Dim1, yend = Dim2, alpha = quality),
             # arrow = arrow(length = unit(0.2, "cm")),
             # color = "black",
             # inherit.aes = FALSE) +
#geom_text_repel(data = species_coords,
              #  aes(x = Dim1, y = Dim2, label = species),
              #  size = 2.7,
              #  color = "black",
              #  inherit.aes = FALSE) +
My_Theme

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

ggplot(points, aes(x = PCoA1, y = PCoA2, color = count))+
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

##---stats BETA---
#2 categories : bioproject and resistance
# Category 1 : bioproject : 
# 3 groups big: "PRJEB65292"  middle:"PRJEB77409"  low:"PRJNA704713", rm low 1 sample

meta <- beta.env[match(rownames(dist_mat), beta.env$sample), ]
nrow(meta)
all(meta$sample == rownames(dist_mat)) # Check same order for 2 tables

keep <- meta$bioproject != "PRJNA704713"
dist_sub <- as.dist(dist_mat[keep, keep])
meta_sub <- meta[keep, ]
meta_sub$bioproject <- factor(meta_sub$bioproject)

bd <- betadisper(dist_sub, meta_sub$bioproject)
anova(bd) # mean sq : 1.1, F-value : 122.3, p-val : 2.2e-16
permutest(bd, permutations = 9999) # mean sq : 1.12, F-value : 122.3, p-val : 1e-04

# PERMANOVA, H0 : 
# Assumption : homogeneous dispersion (similar spread within groups), same variability => HERE the Anova is not possible
#R2: 0.085, F-val 19.3, p-val 1e-04<0.05 The effect of the bioproject on the eta div in unlikely due to random chance.
perm <- adonis2(
  dist_sub ~ bioproject,
  data = meta_sub,
  permutations = 9999,
  method="bray"
)

perm

# Category 2 : resistance: 
meta_sub_res <- meta_sub %>%
  mutate(count = ifelse(is.na(count), "not_resistant", "resistant"))
meta_sub_res

# Dispersion (Leven's test equivalent but multivariate)
bd <- betadisper(dist_sub, meta_sub_res$count)
anova(bd) # mean sq : 0.01, F-value : 0.79, p-val : 0.3752
permutest(bd, permutations = 9999) # mean sq : 0.01, F-value : 0.79, p-val : 0.3737

# PERMANOVA, 
# Assumption : homogeneous dispersion (similar spread within groups), same variability => HERE the Anova is not possible
#R2: 0.014, F-val : 3.0, p-val 0.0044<0.05 The resistance on the eta div in unlikely due to random chance.
perm <- adonis2(
  dist_sub ~ count,
  data = meta_sub_res,
  permutations = 9999,
  method="bray"
)

perm

resistance <- c(
  "resistant" = "orange",
  "not_resistant" = "darkgreen"
)
points$count<- beta.env$count[match(rownames(points),beta.env$sample)]
c = sum(!is.na(points$count))
points$count <- as.character(points$count)
points$count[is.na(points$count)] <- "NA"

points_res <- points%>%
  mutate(count = ifelse(points$count=="NA", "not_resistant", "resistant"))


My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
)

ggplot(points_res, aes(x = PCoA1, y = PCoA2, color = count))+
  geom_point(size = 3) +
  scale_color_manual(values = resistance) +
  labs(
    title = paste("PCoA of Bray-Curtis distance"),
    subtitle = paste(n, " samples", c, " transferable ARGs"),
    x = paste0("PCoA1 (", variance[1], "%)"),
    y = paste0("PCoA2 (", variance[2], "%)")
  ) +
  My_Theme

