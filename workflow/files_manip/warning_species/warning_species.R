# 6/05/2026
# Creating a csv for species with arg and their abundance: 
# Out : 'warning_species.csv'

# Part 1: Importing the necessary packages and data -----------------------

library(dplyr)
library(readr)
library(tidyr)


db1 <- read.csv2(
  "joint_relativ.csv"
)
db <- db1[,1:46]
antibio <- read.csv2(
  "WHO_AMR.csv",
)
arg <- read.csv2(
  "all_samples_one_ARG_per_line.csv",
)
arginsam <- arg %>%
  distinct(sample, species, short_name)

spe <- read_csv(
  "species_relative_abundance_0.005.csv",
  col_types = cols(.default = "c")
  )

# Part 2: Extracting Species, ARG name, abondance intra and intersamples
col = c("species", "abondance_intrasample","abondance_intersamples")

df2 <- db %>%
  pivot_longer(
    cols = -all_of(col),
    names_to = "Short_name",
    values_to = "count"
  ) %>%
  filter(count > 0) %>%
  dplyr::select(species, Short_name, abondance_intrasample, abondance_intersamples)

# Part 3: adding the name of the antibiotitic(s) targeted

df3 <- left_join(df2, antibio, by="Short_name")
df3 <- df3[,1:5]

# Part 4: calculating if the species is an automatic resistance carrier

nbspe <- spe %>%
  pivot_longer(
    cols = -sample,
    names_to = "species",
    values_to = "abundance"
  ) %>%
  filter(!is.na(abundance) & abundance > 0)%>%
  distinct(sample, species) %>% 
  count(species, name = "species_count")


counts <- arginsam %>%
  count(species, short_name, name = "count")

df4 <- df3 %>%
  left_join(counts, by = c("species", "Short_name" = "short_name")) %>%
  mutate(
    count = ifelse(is.na(count), 0, count),
    # Percentage of samples in which species Y carries ARG Z over all samples analyzed with TELcOMB:
    stat_resist_samp = (count / 226) * 100  
  ) %>%
  left_join(nbspe, by = "species") %>%
  mutate(
    # Percentage of samples in which species Y carries ARG Z over samples carrying species y
    stat_resist_spe = (count / species_count) * 100
    # Interpretation should be done with caution, as species detection depends on Kraken results.
    # Although this metric considers presence/absence rather than abundance within each sample,
    # low-abundance species may lead to unreliable estimates (NA values may occur when species are too rare).
  )

col4 = c("species", "Short_name", "antibiotic_class","abondance_intrasample", "abondance_intersamples", "stat_resist_samp","stat_resist_spe")
df4 <- df4[,col4]
col5 = c("species", "arg", "antibiotic_class","species_mean_abund_intrasample", "species_abund_intersample", "count_samples_with_spec&arg/226","count_species_with_arg/count_species")
colnames(df4) <- col5
write.csv(df4, "warning_species.csv", row.names = FALSE)
  
samples <- unique(arginsam$sample)
samples_with_argmge = length(samples) - 1 


