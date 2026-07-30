# 25/04/2026
# Joining the Blast Telcomb output and Taxonomy output :
# for every species with ARGs found wih Blast, the abundance of each species 
# within every sample is mention. 
# The abundance intrasample, i.e  mean abundance conditional on the presence, 
# and the abundance intersample, i.e the frequency of presence of the species
# between the samples, is calculated.
# Input : 
# 'intrasamples_abund_relativ.csv' : 
# renamed file 'species_relative_abundance.csv', obtained via taxonomy/abundancy_csv
# 'intrasample_abund.csv' : 
# file 'species_relative_abundance.csv' normalized in Excel
# 'contig_arg_classif.csv' : 
# 1 line by contig with ARGs by samples with WHO classification
# Out : 'joint_relativ.csv'
# Out : 'joint_norm.csv'

# Part 1: Importing the necessary packages and data -----------------------
library(dplyr)
library(readr)
library(tidyr)
library(networkD3)
library(mixOmics)
library(tibble)

# Warnings : names of the species should NOT include coma ',' and should be replace with a dot '.'. Remove also the quotes as they should not be necessary.
df_abund <- read_csv(
  "intrasample_abund.csv",
  col_types = cols(.default = "c")
)

df_abund_relat <-read_csv(
  "intrasamples_abund_relativ.csv",
  col_types = cols(.default = "c")
)

df_amr <- read_csv(
  "contig_arg_classif.csv",
  col_types = cols(.default = "c")
)

amrsample <- read_csv(
  "contig_arg_classif.csv",
  col_types = cols(.default = "c")
)
#-------------------- NOT NORMALIZED RELATIVE ABUNDANCE--------------------

# Part 2: Raw count matrix species/ARG ------------------------------------
# Handling of the table with the AMR abundances for all contigs.
# The species are those found by Blast for the sequence of the first ARG 
# found by contig. Transformation to have a raw count matrix with 1 row, 
# 1 species, and in column the short names of the ARGs.

# Selection of ARG "only for humans"
df_amr <- df_amr[4:9]
df_amrclean <- df_amr[df_amr$"WHO_classif"=="Authorized for use in humans only" | df_amr$"WHO_classif"=="Authorized for use in humans only ?", c("species", "short_name")]

df_amrclean <- df_amrclean[
  !is.na(df_amrclean$species) &
  !is.na(df_amrclean$short_name),
]
species_list <- unique(df_amrclean$species)
short_names <- unique(df_amrclean$short_name)

# Raw count matrix
result <- matrix(0,
                 nrow = length(species_list),
                 ncol = length(short_names))

rownames(result) <- species_list
colnames(result) <- short_names

for (i in seq_along(species_list)) {
  for (j in seq_along(short_names)) {
    result[i, j] <- sum(
      df_amrclean$species == species_list[i] &
        df_amrclean$short_name == short_names[j],
      na.rm = TRUE   # ← CRUCIAL
    )
  }
}
result
result_df <- as.data.frame(result)
result_df$species <- rownames(result_df)
result_df <- result_df[, c("species", setdiff(colnames(result_df), "species"))]

# Part 3: Intra and intersamples abundance---------------------------------
# Handling of the table with the species abundances by standardized samples.
# The species present are those with a cutoff of 0.005 in kraken.
# Transformation to have on line the species, in column the samples, 
# the average intrasample abundance and the intersample abundances.

tdf_abund <- as.data.frame(t(df_abund[ , -1]), stringsAsFactors = FALSE)
colnames(tdf_abund) <- df_abund[[1]]

tdf_abund$species <- rownames(tdf_abund)
tdf_abund <- tdf_abund[, c("species", setdiff(colnames(tdf_abund), "species"))]
tdf_abund <- tdf_abund[tdf_abund$species != "Total", ]
                                              
# Abondance intersample (= frequency of presence of the species between samples)

sample_cols <- setdiff(colnames(tdf_abund), "species")
n_samples <- length(sample_cols)

tdf_abund$abondance_intersamples <- apply(
  tdf_abund[, sample_cols],
  1,
  function(x) (sum(!is.na(x) & x != 0) / n_samples)*100
)

tdf_abund <- tdf_abund %>%
  relocate(abondance_intersamples, .after = species)

# Abondance intrasample. The average is only taken from samples containing ARGs 
# (= average abundance conditional on presence)

sample_cols <- setdiff(
  colnames(tdf_abund),
  c("species", "abondance_intersamples","abondance_intrasample")
)

tdf_abund[ , sample_cols] <- lapply(
  tdf_abund[ , sample_cols],
  function(x) as.numeric(gsub(",", ".", x))
)

tdf_abund$abondance_intrasample <- apply(
  tdf_abund[, sample_cols],
  1,
  function(x) 
  {
    x <- as.numeric(x)
    x_non_nuls <- x[!is.na(x) & x != 0]
    
    if (length(x_non_nuls) == 0) {
      return(NA_real_)   # ou NA si tu préfères
    } else {
      return(mean(x_non_nuls))
    }
  }
  
)

tdf_abund <- tdf_abund %>%
  relocate(abondance_intrasample, .after = species)

# Part 4: Join between the 2 tables ---------------------------------------

join <- left_join(result_df, tdf_abund, by="species")
summary(join)

write.csv2(join,
          file = "joint_norm.csv",
          row.names = FALSE)
# 'Na' in the colomn 'sample' means the species is absent in Kraken taxonomy

#---------------------- NORMALIZED RELATIVE ABUNDANCY----------------------

# Handling of the table with the species abundances by standardized samples.
# The species present are those with a cutoff of 0.005 in kraken. Transformation
# to have on line the species, in column the samples, the average intrasample
# abundance and the intersample abundances.

tdf_abund_relat <- as.data.frame(t(df_abund_relat[ , -1]), stringsAsFactors = FALSE)
colnames(tdf_abund_relat) <- df_abund_relat[[1]]

tdf_abund_relat$species <- rownames(tdf_abund_relat)
tdf_abund_relat <- tdf_abund_relat[, c("species", setdiff(colnames(tdf_abund_relat), "species"))]
tdf_abund_relat <- tdf_abund_relat[tdf_abund_relat$species != "Total", ]

# Part 3 bis : Intra and intersamples normalized abundance----------------
# Intersample abundnace

sample_cols <- setdiff(colnames(tdf_abund_relat), "species")
n_samples <- length(sample_cols)

tdf_abund_relat$abondance_intersamples <- apply(
  tdf_abund_relat[, sample_cols],
  1,
  function(x) (sum(!is.na(x) & x != 0) / n_samples)*100
)

tdf_abund_relat <- tdf_abund_relat %>%
  relocate(abondance_intersamples, .after = species)

# Intrasample abundance

sample_cols <- setdiff(
  colnames(tdf_abund_relat),
  c("species", "abondance_intersamples","abondance_intrasample")
)

tdf_abund_relat[ , sample_cols] <- lapply(
  tdf_abund_relat[ , sample_cols],
  function(x) as.numeric(gsub(",", ".", x))
)

tdf_abund_relat$abondance_intrasample <- apply(
  tdf_abund_relat[, sample_cols],
  1,
  function(x) 
  {
    x <- as.numeric(x)
    x_non_nuls <- x[!is.na(x) & x != 0]
    
    if (length(x_non_nuls) == 0) {
      return(NA_real_)   # ou NA si tu préfères
    } else {
      return(mean(x_non_nuls))
    }
  }
  
)

tdf_abund_relat <- tdf_abund_relat %>%
  relocate(abondance_intrasample, .after = species)

# Part 4 bis : Join between the 2 tables ----------------------------------

join_relativ <- left_join(result_df, tdf_abund_relat, by="species")
summary(join_relativ)

write.csv2(join_relativ,
          file = "joint_relativ.csv",
          row.names = FALSE)
