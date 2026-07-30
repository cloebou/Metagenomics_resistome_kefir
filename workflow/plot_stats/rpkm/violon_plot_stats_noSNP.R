# 15/06/2026
# Generating violon plot : RPKM for every ARGs in colocalisation with an 
# MGE, classified as OfH, no mutations.
# Plotting every samples containing this ARG.

# Part 1: Importing the necessary packages and data -----------------------
library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)
library(rstatix)

df <- read_csv("rpkm_biop.csv")
df <- df %>%
  mutate(species = ifelse(species == "other", "Other", species)) %>%
  mutate(bioproject = ifelse(bioproject == "PRJEB65292", "Traditional", bioproject)) %>%
  mutate(bioproject = ifelse(bioproject == "PRJEB77409", "Industrial", bioproject)) %>%
  mutate(bioproject = ifelse(bioproject == "PRJNA388572", "Traditional", bioproject)) %>%
  mutate(bioproject = ifelse(bioproject == "PRJNA704713", "Traditional", bioproject))
  

# Part 2: Filtering data ---------------------------------------------------

# ARG with mutation
snp <- c("A16S",
         "FABG",
         "FUSA",
         "GLPT",
         "GYRA",
         "GYRB",
         "MEXT",
         "MLS23S",
         "MURA",
         "O23S",
         "OMP36",
         "OMPK",
         "OMPK36",
         "OMPFB",
         "PARE",
         "PHOB",
         "PHOP",
         "PTSL",
         "RPOB",
         "TUFAB",
         "UL3",
         "UHPT",
         "ACRR",
         "AMPCR",
         "PMRAR",
         "RPSL",
         "PARC")


df1 <- df %>%
  filter(
    WHO_classif %in% c(
    "Authorized for use in humans only",
    "Authorized for use in humans only ?"),
    !is.na(short_name),
    !short_name %in% snp
    )
df2 <- df1 %>%
  group_by(short_name) %>%
  mutate(n_samples = n_distinct(sample)) %>%  # compter nb de samples uniques
  ungroup() %>%
  filter(n_samples > 5)

df_filtered <- df2 %>%
  group_by(species) %>%
  mutate(n_samples = n_distinct(sample)) %>%  # compter nb de samples uniques
  ungroup() %>%
  filter(n_samples > 5)

#df_filtered[1:20,]

color_map <- c(
  "PRJEB65292" = "red",
  "PRJEB77409" = "steelblue",
  "PRJNA388572" = "green",
  "PRJNA704713" = "purple"
)

# Part 4: Stats------------------------------------------------------------

# Normal distribution (H0 : Follows a Gaussian law)
for (i in unique(df_filtered$species)) {
        print(i)
        print(shapiro.test(df_filtered$rpkm[df_filtered$species==i]))
}
# => Do not accept H0 for every species and few samples so no ANOVA

# Kruskal-Wallis test : (H0 : there is no statistically significant difference
# between the median rpkm ratings in all species)
kruskal.test(rpkm ~ species, data = df_filtered)
# => globale p-value 0.0001< 0.05 : H0 rejected, possible to do a Dunn's Test

# Paired comparison -> Dunn's Test with BH correction for p-values
# H0: there is no difference between groups 
dunn_res <- df_filtered %>%
  dunn_test(rpkm ~ species, p.adjust.method = "BH")
dunn_res

# Part 3: Ploting data ------------------------------------------------------
# Log-transform
df_filtered <- df_filtered %>%
  mutate(log_rpkm = log10(rpkm))

# Dunn test
dunn_res <- df_filtered %>%
  dunn_test(log_rpkm ~ species, p.adjust.method = "BH") %>%  # BH > Bonferroni
  add_y_position()

# Keep only significatifs (< 0.05)
dunn_res <- dunn_res %>%
  filter(p.adj < 0.05)

# Plot 1 - shapes related to grain origin
My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
)

p <- ggplot(df_filtered, aes(x = species, y = log_rpkm)) +
  geom_violin(fill = "grey90", color = NA, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +

  # points 
  geom_jitter(aes(color = short_name, shape = bioproject),
              width = 0.15, size = 1.5, alpha = 0.5) +
  scale_color_brewer(
    palette = "Paired") +

  # stats
  stat_pvalue_manual(
    dunn_res,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE
  ) +

  # Axis and title
  labs(
    x = NULL,
    y = expression(log[10]*"(RPKM)"),
    title = "ARG abundance across species",
    subtitle = "Only ARGs in >5 samples and WHO-restricted antibiotics and no mutations\nSpecies in >5 samples with these ARG. All kefir are handmade.",
   # 'Only' means <0.005% diversity according to Kraken2.",
    color = "Study",
    shape = "Grain origin"
   
  ) +

  theme_classic(base_size = 12) +
  theme(
    title = element_text(size = 18,face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1,size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "right",
  )
  

p


# Plot 2 - no particular shapes
My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
)

p <- ggplot(df_filtered, aes(x = species, y = log_rpkm)) +
  geom_violin(fill = "grey90", color = NA, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  
  # points 
  geom_jitter(aes(color = short_name),
              width = 0.15, size = 1.5) +
  scale_color_brewer(
    palette = "Paired") +
  
  # stats
  stat_pvalue_manual(
    dunn_res,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE
  ) +
  
  # Axis and title
  labs(
    x = NULL,
    y = expression(log[10]*"(RPKM)"),
    title = "ARG abundance across species",
    #subtitle = "Only ARGs in >5 samples and WHO-restricted antibiotics and no mutations\nSpecies in >5 samples with these ARG. All kefir are handmade.",
    # 'Only' means <0.005% diversity according to Kraken2.",
    color = "ARG"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    title = element_text(size = 18,face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1,size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "right",
    legend.text = element_text(size = 16)
  )

p

