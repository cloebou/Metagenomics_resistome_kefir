#!/usr/bin/env Rscript
# Getting the species that represent more then 0.005% of the sample. Not normalized. 
# Out : "species_relative_abundance.csv" 

# Importing the necessary packages and config parameters -----------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))

library(data.table)
library(dplyr)
library(tidyr)

kraken_dir <- "results/kraken"
out_dir    <- "results/"
threshold  <- 0.005   # relative abundance (%)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Reading Kraken reports---------------------------------------------------------
files <- list.files(
  kraken_dir,
  pattern = "\\.kraken.report$",
  recursive = TRUE,
  full.names = TRUE
)

read_kraken <- function(file) {
  fread(
    file,
    col.names  = c("percent", "reads", "reads_clade",
                   "rank", "taxid", "name"),
    colClasses = list(character = "rank")
  ) %>%
    mutate(sample = basename(dirname(file)))
}
kraken <- bind_rows(lapply(files, read_kraken))

# Building species abundance table-----------------------------------------------
species_table <- kraken %>%
  filter(rank == "S") %>% # Species level only
  group_by(sample, name) %>% 
  summarise(percent = sum(percent), .groups = "drop") %>% # Sum abundance per sample × species
  group_by(name) %>% 
  mutate(mean_p = mean(percent)) %>% # Global mean abundance per species
  ungroup() %>%
  # Group rare species
  mutate(
    name = ifelse(
      mean_p < threshold,
      paste0("Other <", threshold, "%"),
      name
    )
  ) %>%
  # Re-sum after grouping "Other"
  group_by(sample, name) %>%
  summarise(percent = sum(percent), .groups = "drop") %>%
  # Wide format (one column per species)
  pivot_wider(
    names_from  = name,
    values_from = percent,
    values_fill = 0
  ) %>%

  arrange(sample)

# Order species columns by global abundance--------------------------------------
species_order <- species_table %>%
  select(-sample) %>%
  colMeans() %>%
  sort(decreasing = TRUE) %>%
  names()

species_table <- species_table %>%
  select(sample, all_of(species_order))

# Write CSV----------------------------------------------------------------------
write.csv(
  species_table,
  file = file.path(out_dir, "species_relative_abundance.csv"),
  row.names = FALSE,
  quote = FALSE
)
