# 30/07/2026
# Generating bubble plot : in x ARGs, y species, circle for abundance
# Selecting only the mutation WITHOUT mutations and classified as "Only for Humans"

# Part 1: Importing the necessary packages and data -----------------------
library(dplyr)
library(ggplot2)
library(readr)
library(forcats)

df <- read_csv(
	         "all_samples_one_ARG_per_line-AMR_who.csv",
		   col_types = cols(.default = "c")
		 )

# Part 2: Cleaning data ---------------------------------------------------
df <- df %>%
  mutate(species = ifelse(species == "other", "Other", species))
         
# ARG with mutation
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

# Filter to avoid having the type of metal resistance and mutational resistance
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
			       short_name != "",
			       !short_name %in% snp
			       )
# Part 3: Calculation of Species x ARG frequencies-------------------------

freq_df <- df_filtered %>%
	  count(species, short_name, WHO_classif, name = "frequency")

# Optional : Only most frequent AMR for clarity
freq_df <- freq_df %>%
	  filter(frequency >= 2)

# Reorder factors by total frequency
freq_df <- freq_df %>%
	   mutate(
		     species = fct_reorder(species, frequency, sum),
			   short_name = fct_reorder(short_name, frequency, sum)
			   )
# Part 4: Plotting---------------------------------------------------------

who_colors <- c(
  "Authorized for use in humans only" = "#7f0000",
  "Authorized for use in animals only " = "#7f0000",
  "Authorized for use in humans only ?" = "#7f0000",
  "HPCIA" = "#d73027",                      
  "CIA" = "#fc8d59",
  "CIA ?" = "#fc8d59", 
  "HIA" = "#fee08b", 
  "HIA ?" = "#fee08b",  
  "IA" = "#d9ef8b",
  "Not medically important" = "#91cf60",
  "0" = "#bdbdbd",
  "NA" = "#bdbdbd",
  "NA " = "#cccccc"
  )
    
  
    # ---- BUBBLE PLOT ----
    p <- ggplot(freq_df, aes(
			       x = short_name,
			       y = species,
			       size = frequency,
			       color = WHO_classif
				 )) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(1, 12)) +
  scale_color_manual(
    values = who_colors,
    name = "WHO classification"
    ) +
  labs(
    x = "ARG",
    y = "Species",
    size = "Frequency",
    title = "AMR Distribution by species, only drug and Multi-compound target"
    ) +
  theme_minimal(base_size = 8) +
    theme(
	      axis.text.x = element_text(angle = 90, hjust = 0.5),
	      legend.position = "right",
	      panel.grid.major = element_line(color = "grey90")
	        )

  print(p)

  # Save
  ggsave(
	   "bubbleplot_short_who_nosnp.png",
	     p,
	     width = 14,
	       height = 10,
	       dpi = 300
	       )


