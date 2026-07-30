#!/usr/bin/env Rscript
# Plotting barplot at family level : 20 samples per graph for more visibility

# Importing the necessary packages and config parameters -----------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))

library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)

kraken_dir <- "results/kraken"
out_dir    <- "results/visu"
threshold  <- 0.05      
group_size <- 20        # samples per plot

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
		        col.names = c("percent","reads","reads_clade",
				                        "rank","taxid","name"),
		    colClasses = list(character = "rank")
		      ) %>%
    mutate(sample = basename(dirname(file)))
}

kraken <- bind_rows(lapply(files, read_kraken))

samples <- sort(unique(kraken$sample))
sample_groups <- split(samples, ceiling(seq_along(samples) / group_size))

# Building a global palette------------------------------------------------------
all_taxa <- kraken %>%
	  filter(rank %in% c("S")) %>%
	    group_by(name) %>%
	      summarise(mean_p = mean(percent), .groups = "drop") %>%
	        filter(mean_p >= threshold) %>%
		  pull(name) %>%
		    sort()

	    n_colors <- length(all_taxa)
	    base_palette <- colorRampPalette(brewer.pal(12, "Paired"))(n_colors)

	    palette <- setNames(base_palette, all_taxa)
	    palette["Other <0.05%"] <- "grey70"

# Function to plot stacked bars--------------------------------------------------
	    make_stacked <- function(rank_level, samples_subset, outfile, title) {

		      df <- kraken %>%
			          filter(rank == rank_level, sample %in% samples_subset) %>%
				      group_by(sample, name) %>%
				          summarise(percent = sum(percent), .groups = "drop")

				    # group rare taxa
				    df <- df %>%
					        group_by(name) %>%
						    mutate(mean_p = mean(percent)) %>%
						        ungroup() %>%
							    mutate(
								         name = ifelse(mean_p < threshold,
										                           paste0("Other <", threshold, "%"),
													                       name)
									     ) %>%
				        group_by(sample, name) %>%
					    summarise(percent = sum(percent), .groups = "drop")

				      ggplot(df,
					              aes(x = sample,
							               y = percent,
								                    fill = name)) +
    geom_col(position = "fill", width = 0.85) +
        scale_y_continuous(labels = scales::percent) +
	    scale_fill_manual(values = palette, drop = FALSE) +
	        labs(
		           x = NULL,
			         y = "Relative abundance (%)",
			         fill = NULL,
				       title = title
				     ) +
    theme_minimal(base_size = 14) +
        theme(
	            axis.text.x = element_text(angle = 45, hjust = 1),
		          panel.grid.major.x = element_blank()
		        )

      ggsave(outfile, width = 12, height = 7)
	    }

	    # Generate plots per group
	    c=0
	    for (i in seq_along(sample_groups)) {

		      group_samples <- sample_groups[[i]]
	      idx <- sprintf("%02d", i)

	        make_stacked(
			         rank_level = "S",
				     samples_subset = group_samples,
				     outfile = file.path(out_dir, paste0("species_stacked_", idx, ".png")),
				         title = paste("Species-level composition (group", idx, ")")
				       )
	      c=c+1
	    }

	    message(c," Stacked barplots generated")
