# 02/05/2026
# Generating bubble plot : in x ARGs, y species, circle for abundance

# Part 1: Importing the necessary packages and data -----------------------

library(dplyr)
library(ggplot2)
library(readr)
library(forcats)
library(stringr)

df <- read_csv(
	         "all_samples_one_ARG_per_line_v2.csv",
		   col_types = cols(.default = "c")
		 )

# Part 2: Cleaning data ---------------------------------------------------

df_clean <- df %>%
	  filter(
		     !is.na(species),
		         species != "",
		         !is.na(AMR),
			     AMR != ""
			   )

# Part 3: Calculation of Species x ARG frequencies-------------------------

freq_df <- df_clean %>%
	  count(species, AMR, name = "frequency")

  # Optional : Only most frequent AMR for clarity
  freq_df <- freq_df %>%
	    filter(frequency >= 2)

    # Reorder factors according to total frequency
    freq_df <- freq_df %>%
	      mutate(
		         species = fct_reorder(species, frequency, .fun = sum),
			     arg = fct_reorder(AMR, frequency, .fun = sum)
			   )

# Part 4: Plotting---------------------------------------------------------

    p <- ggplot(freq_df, aes(
			       x = arg,
			         y = species,
			         size = frequency
				 )) +
  geom_point(alpha = 0.7, color = "#1f78b4") +
    scale_size_continuous(range = c(1, 12)) +
      labs(
	       x = "ARG",
	           y = "Species",
	           size = "Frequency",
		       title = "ARG Distribution by species, for all kinds of ARG in colocalization with MGE(s)."
		     ) +
  theme_minimal(base_size = 12) +
    theme(
	      axis.text.x = element_text(angle = 45, hjust = 1),
	          panel.grid.major = element_line(color = "grey90")
	        )

  print(p)

  # Save
  ggsave(
	   "ARG_all_bubble.png",
	     p,
	     width = 14,
	       height = 10,
	       dpi = 300
	       )
























