# 19/05/2026 
# Multivariate analysis of ARGs targeting antibiotics classified as 'Only for Humans'. 
# RQ : Why a PCA or PCoA ? 
  #PCA : use of FactomineR or MixOmics package and  PCoA use of Vegan(cmdscale) package
  #For a PCA we suppose the data is an Euclidienne distance not for a PCoA
# Part 1: Importing and loading the necessary packages and data
# Part 2: Data cleaning and formatting
# Part 3: PCA
# Part 4: Dissimilarity matrix with Bray-Curtis index
# Part 5: PCoA


#  Part 1: Importing the necessary packages and data ----------------------

library(dplyr)
library(tidyr)
library(mixOmics) # PCA
library(vegan) # Dissimilarity matrix
library(ggplot2)
library(ggrepel)

# Obtained in merging/sankey/exposition.R :

beta.env <- read.csv(
  "count_arg_dfs.csv",dec = ",", sep = ",", header = TRUE)

df <- read.csv(
  "amrsample_met.csv",dec = ",", sep = ",", header = TRUE)

# Part 2: Data cleaning and formatting ------------------------------------

df <- df[df$"WHO_classif"=="Authorized for use in humans only" | df$"WHO_classif"=="Authorized for use in humans only ?", c("sample","contig","tFermentation","bioproject","sequenceur","origin","species","short_name")]
df <- df[!is.na(df$short_name),]
species_list <- unique(df$species)
short_names <- unique(df$short_name)
df_clean =  df[c("sample","short_name")]

# Abundance matrix
tab_name <- df_clean %>%
  mutate(value = 1) %>%
  distinct(sample, short_name, .keep_all = TRUE) %>%  # éviter doublons
  pivot_wider(
    names_from = short_name,
    values_from = value,
    values_fill = 0
  )

tab_name_biop <- tab_name
tab_name_biop <- column_to_rownames(tab_name_biop, var = "sample")
tab_name_biop$bioproject <- beta.env$bioproject[match(rownames(tab_name_biop),beta.env$sample)]


# Part 3: sPCA ------------------------------------------------------------

res.acp = mixOmics::spca(tab_name_biop[,-44])
plotIndiv(res.acp, ind.names = TRUE, group = tab_name_biop$bioproject, legend = TRUE, title = "sPCA on ARG abundance targeting antibiotics classified as 'Only for Humans', n = 105")
plotVar(res.acp, title = "sPCA on ARG abundance targeting antibiotics classified as 'Only for Humans', n = 105")

# Extract the variables used to construct the first PC
selectVar(res.acp, comp = 1)$name 
# Depict weight assigned to each of these variables
plotLoadings(res.acp, method = 'mean', contrib = 'max', title = "Loading on comp1 : \n sPCA on ARG abundance, n = 105") 
plotLoadings(res.acp,comp = 2 ,method = 'mean', contrib = 'max', title = "Loading on comp2 : \n sPCA on ARG abundance, n = 105")


#  Part 4: Dissimilarity matrix with Bray-Curtis index --------------------

tab <- tab_name[,-1]
beta <- vegdist(tab, method = "bray")
n = nrow(tab)


# Part 5: PCoA ------------------------------------------------------------

pcoa_result <- cmdscale(beta, eig = TRUE, k = 2)
barplot(pcoa_result$eig)

# Grouped by Bioproject
points <- as.data.frame(pcoa_result$points,tab_name$sample)
colnames(points) <- c("PCoA1", "PCoA2")
points$bioproject <- beta.env$bioproject[match(rownames(points),beta.env$sample)]
variance <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 2)

ggplot(points, aes(x = PCoA1, y = PCoA2, color = bioproject)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_text_repel(
    aes(label = tab_name$sample),
    size = 2.5,
    max.overlaps = 20, # Change "Inf" by for eg 20 to have more clarity
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.2
  ) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.8, alpha = 0.2) +
  labs(
    title = "PCoA of ARG Diversity (Bray-Curtis)",
    subtitle = paste(n, "samples – Human-only classification"),
    x = paste0("PCoA1 (", variance[1], "%)"),
    y = paste0("PCoA2 (", variance[2], "%)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


# Grouped by sequenceur 
points$sequenceur <- beta.env$sequenceur[match(rownames(points),beta.env$sample)]

variance <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 2)

ggplot(points, aes(x = PCoA1, y = PCoA2, color = sequenceur))+
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linetype = 2, alpha = 0.3) +
  labs(
    title = paste("PCoA Analysis (Bray-Curtis)\n", n, " samples"),
    x = paste0("PCoA1 (", variance[1], "%)"),
    y = paste0("PCoA2 (", variance[2], "%)")
  ) +
  theme_minimal()  

