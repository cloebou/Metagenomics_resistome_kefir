# 15/06/2026
# Generating violon plot : RPKM for every ARGs in colocalisation with an 
# MGE and classified as OfH. Plotting every samples containg this ARG.

# Part 1: Importing the necessary packages and data -----------------------
library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)
library(rstatix)

df <- read_csv("rpkm_biop.csv")

# Part 2: Filtering data and basic stats -----------------------------------
count_shortname <- df %>%
  group_by(species) %>%
  summarise(
    nb_short_name = n_distinct(short_name)
  ) %>%
  arrange(desc(nb_short_name))

count_shortname
length(unique(df$sample)) # 207
length(unique(df$species))
unique(df$species)

df1 <- df %>%
  filter(WHO_classif %in% c(
    "Authorized for use in humans only",
    "Authorized for use in humans only ?"
  ))

length(unique(df1$sample))
length(unique(df1$species))

# Mean RPKM
rpkm_summary <- df1 %>%
  group_by(species) %>%
  summarise(
    mean_rpkm = mean(rpkm, na.rm = TRUE),
    sd_rpkm = sd(rpkm, na.rm = TRUE),
    .groups = "drop"
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
# => globale p-value < 0.05 : H0 rejected, possible to do a Dunn's Test

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

# Plot
p <- ggplot(df_filtered, aes(x = species, y = log_rpkm)) +
  geom_violin(fill = "grey90", color = NA, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +

  # points 
  geom_jitter(aes(color = short_name),
              width = 0.15, size = 1.5) +
  scale_color_brewer(palette = "Paired") +
  
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
   # subtitle = "Only AMRs in >5 samples and WHO-restricted antibiotics\nSpecies in >5 samples with theses AMR.",
   # 'Only' means <0.005% diversity according to Kraken2."
    color = "AMR"
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

