# 15/06/2026
# Generating violon plot : RPKM for every ARGs in colocalisation with an 
# MGE and classified as OfH. Plotting every samples containg this ARG.

# Part 1: Importing the necessary packages and data -----------------------
library(ggplot2)
library(dplyr)
library(readr)

df <- read_csv("rpkm_biop.csv")

# Part 2: Filtering data ---------------------------------------------------
df1 <- df %>%
  filter(WHO_classif %in% c(
    "Authorized for use in humans only",
    "Authorized for use in humans only ?"
  ))
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

# Part 3: Plotting---------------------------------------------------------
ggplot(df_filtered, aes(x = species, y = rpkm)) +
  geom_violin(fill = "lightblue", trim = FALSE) +
  geom_jitter(width = 0.5, size = 2, alpha = 0.6, aes(color = short_name)) + # handling overplotting
  scale_y_log10() + 
  labs(
    x = "ARG (short_name)",
    y = "RPKM (log10 scale)",
    title = "Distribution of ARG abundance (RPKM)", 
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
