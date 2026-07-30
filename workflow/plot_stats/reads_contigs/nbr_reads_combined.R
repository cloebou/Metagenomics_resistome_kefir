# Barplot with Number of classified reads by Kraken2 on the PlusPF-16 database 
# and raw reads
library(readr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(tibble)
library(scales)

reads <- read.csv("read.csv",dec = ".", sep = ",", header = FALSE)
reads.env <- read.csv("sample_collection.csv",dec = ",", sep = ";", header = TRUE)
reads$bioproject <- reads$V2
reads$bioproject <- reads.env$bioproject[match(reads$V1,reads.env$sample)]
reads[reads$V1=="SRR13801805", "bioproject"] <- "PRJNA704713"

creads <- read.csv("read_classified_kraken.csv",dec = ".", sep = ",", header = FALSE)
creads.env <- read.csv("sample_collection.csv",dec = ",", sep = ";", header = TRUE)
creads$bioproject <- creads$V2
creads$bioproject <- creads.env$bioproject[match(creads$V1,creads.env$sample)]
creads[creads$V1=="SRR13801805", "bioproject"] <- "PRJNA704713"

# Assign colors to each bar by BioProjects
color_cr <- c(
  "PRJEB65292" = "red",
  "PRJEB77409" = "steelblue1",
  "PRJNA388572" = "green",
  "PRJNA704713" = "purple"
)
color_r<- c(
  "PRJEB65292" = "red4",
  "PRJEB77409" = "steelblue",
  "PRJNA388572" = "green3",
  "PRJNA704713" = "purple4"
)

# Merging 2 tables
df <- merge(
  reads[, c("V1", "V2", "bioproject")],
  creads[, c("V1", "V2")],
  by = "V1",
  suffixes = c("_reads", "_creads")
)

# Long format
df_long <- pivot_longer(
  df,
  cols = c(V2_reads, V2_creads),
  names_to = "Condition",
  values_to = "Value"
)

# Variable combinant condition et bioproject
df_long$Bioprojects <- paste(df_long$Condition,
                            df_long$bioproject,
                            sep = "_")

# Colors and names
Bioproject <- c()

for(bp in names(color_r)) {
  Bioproject[paste0("V2_reads_", bp)]  <- color_r[bp]
  Bioproject[paste0("V2_creads_", bp)] <- color_cr[bp]
}

labels_fill <- c()

for(bp in names(color_r)) {
  labels_fill[paste0("V2_reads_", bp)]  <- paste(bp, "- Raw reads")
  labels_fill[paste0("V2_creads_", bp)] <- paste(bp, "- Cleaned and classified reads")
}

ggplot(
  df_long,
  aes(
    x = V1,
    y = Value,
    fill = Bioprojects
  )
) +
  geom_bar(
    stat = "identity",
    position = "stack"
  ) +
  scale_fill_manual(
    values = Bioproject,
    labels = labels_fill
  ) +
  scale_y_log10(
    labels = label_scientific()
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 3
    )
  ) +
  labs(
    x = "Sample",
    y = "Number of reads (log10)",
    title = "Abundance of raw reads and clean and classified reads by Kraken2 on the PlusPF-16 database"
  )