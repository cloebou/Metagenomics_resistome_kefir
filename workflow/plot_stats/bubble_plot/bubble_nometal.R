# 02/05/2026
# Generating bubble plot : in x ARGs, y species, circle for abundance

# Part 1: Importing the necessary packages and data -----------------------

library(dplyr)
library(ggplot2)
library(readr)
library(forcats)

df <- read_csv(
	         "all_samples_one_ARG_per_line_v2.csv",
		   col_types = cols(.default = "c")
		 )

# Part 2: Cleaning and filtering data --------------------------------------

df_filtered <- df %>%
	  filter(
		     !is.na(type_resistance),
		         type_resistance != "Metals",
		         type_resistance != "Biocide",
		     !is.na(AMR),
		         AMR != "Biocide_and_metal_resistance",
		     !is.na(species),
			       species != "",
			   !is.na(short_name),
			       short_name != ""
			       )

# Part 3: Calculation of Species x ARG frequencies-------------------------
freq_df <- df_filtered %>%
	  count(species, short_name, name = "frequency")

  # Optional : Keep only most frequent AMR for clarity
  freq_df <- freq_df %>%
	    filter(frequency >= 2)

   freq_df <- freq_df %>%
	     mutate(
		         species = fct_reorder(species, frequency, sum),
			     short_name = fct_reorder(short_name, frequency, sum)
			   )

# Part 4: Plotting---------------------------------------------------------
    p <- ggplot(freq_df, aes(
			       x = short_name,
			         y = species,
			         size = frequency
				 )) +
  geom_point(
	         alpha = 0.7,
		     color = "#2c7fb8"
		   ) +
  scale_size_continuous(range = c(1, 12)) +
    labs(
	     x = "AMR",
	         y = "Species",
	         size = "Frequency",
		     title = "AMR Distribution by species, only drug and Multi-compound target. ARG in colocalization with MGE(s)."
		   ) +
  theme_minimal(base_size = 8) +
    theme(
	      axis.text.x = element_text(angle = 90, hjust = 0.5),
	          panel.grid.major = element_line(color = "grey90")
	        )

  print(p)

  # Save
  ggsave(
	   "AMR_species_bubbleplot_no_Metals_long.png",
	     p,
	     width = 14,
	       height = 10,
	       dpi = 300
	       )


