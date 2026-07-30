# 30/07/2026
# Getting some insight into the data

# Part 1: Importing the necessary packages and data -----------------------
library(readr)
library(dplyr)

df <- read_csv("~/Desktop/out/amrsample_met.csv")
View(amrsample_met)

all <- read_delim("all_samples_one_ARG_per_line.csv", 
                  delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(all)

# Part 2: Manipulation ----------------------------------------------------
samples = unique(amrsample_met$sample)
length(samples)

df %>%
  group_by(bioproject) %>%
  summarise(n_samples = n_distinct(sample))

count_class = all %>%
  group_by(AMR) %>%
  summarise(n_samples = n_distinct(sample))
sum(count_class$n_samples)
