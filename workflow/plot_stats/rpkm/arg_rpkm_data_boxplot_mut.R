# 15/06/2026
# Distribution of the most prevalent AMR determinants co-localized 
# with mobile genetic elements (minimum occurrence: five samples). 
# Boxes are coloured according to resistance mechanism and points to its BioProject.

# Part 1: Importing the necessary packages and data -----------------------
library(ggplot2)
library(dplyr)
library(readr)

df <- read_csv("rpkm_biop.csv")
length(unique(df$sample))

# Part 2: Filtering data ---------------------------------------------------
df1 <- df %>%
  filter(WHO_classif %in% c(
    "Authorized for use in humans only",
    "Authorized for use in humans only ?"
  ))

length(unique(df1$sample))
length(unique(df1$short_name))

#Mean RPKM
rpkm_summary <- df1 %>%
  group_by(short_name) %>%
  summarise(
    mean_rpkm = mean(rpkm, na.rm = TRUE),
    sd_rpkm = sd(rpkm, na.rm = TRUE),
    .groups = "drop"
  )

rpkm_summary

color_map <- c(
  "PRJEB65292" = "red",
  "PRJEB77409" = "steelblue",
  "PRJNA388572" = "green",
  "PRJNA704713" = "purple"
)

biopcol <- color_map[df_filtered$bioproject]

# Part 3: Plotting---------------------------------------------------------
ggplot(df1, aes(x = short_name, y = rpkm)) +
  geom_violin(fill = "lightblue", trim = FALSE) +
  geom_jitter(width = 0.5, size = 2, alpha = 0.6, aes(color = bioproject)) + # handling overplotting
  scale_y_log10() + 
  labs(
    x = "ARG (short_name)",
    y = "RPKM (log10 scale)",
    title = "Distribution of ARG abundance (RPKM)", 
    subtitle = "ARG in colocalization with MGE, classified as 'Only for Humans'"
  ) +
  scale_color_manual(values = color_map) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
# Part 4 - Mutation-------------------
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

df2 <- df1%>%
  mutate(mutation_status = ifelse(short_name %in% snp,
                                  "mutation",
                                 "not_mutation"))
df3 <- df2[df2$mutation_status=="mutation",]
length(unique(df3$short_name))
mutation <- c(
  "mutation" = "lightgreen",
  "not_mutation" = "darkgrey"
)

ggplot(df2, aes(x = short_name, y = rpkm, fill = mutation_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.5,
    size = 2,
    alpha = 0.6,
    aes(color = bioproject)
  ) +
  scale_fill_manual(values = mutation) +
  scale_color_manual(values = color_map) +
  scale_y_log10() +
  labs(
    x = "ARG (short_name)",
    y = "RPKM (log10 scale)",
    title = "Distribution of ARG abundance (RPKM)",
    subtitle = "ARG in colocalization with MGE, classified as 'Only for Humans'",
    fill = "Mutation status"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

