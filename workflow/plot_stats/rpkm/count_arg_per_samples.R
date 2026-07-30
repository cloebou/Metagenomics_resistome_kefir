# 15/06/2026
# Count of samples containing each ARG in coloc with MGE classif OfH

# Part 1: Importing the necessary packages and data -----------------------
library(dplyr)
library(readr)

df <- read_csv("rpkm.csv")
mean(df$rpkm)
nrow(df)
length(unique(df$sample))
unique(df$AMR)
unique(df$short_name)

# Part 2: Filtering data ---------------------------------------------------
df_filtered <- df %>%
  filter(WHO_classif %in% c(
    "Authorized for use in humans only",
    "Authorized for use in humans only ?"
  ))
unique(df_filtered $AMR)
unique(df_filtered$short_name)

# Part 3: Extracting data ---------------------------------------------------
result <- df_filtered %>%
  group_by(short_name) %>%
  summarise(n_samples = n_distinct(sample)) %>%
  arrange(desc(n_samples))

result
