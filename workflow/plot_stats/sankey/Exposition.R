# 6/05/2026
# Formatting and joining output from TELcOMB and taxonomy : 
# Out : 'amrsample_met.csv' file : 1 line = 1 ARG in 1 contig + species + WHO classif
# Out : 'count_arg_dfs.csv' file : 1 line = 1 sample and 1 ARG classified as 'Only for humans'
# Out : Sankey plot : ARG -> Samples 
#  Color by Bioproject type : 'sankey_arg_samples.html' 
#  Color by fermentation time : 'sankey_arg_samples_ferm.html' 
#  Color by sequencer type : 'sankey_arg_samplesseq.html' 
# And some other visual representation

# Part 1: Importing the necessary packages and data -----------------------

library(dplyr)
library(readr)
library(tidyr)
library(networkD3)
library(mixOmics)
library(tibble)

join <- read.csv2(
  "joint_norm.csv"
)
amrsample <- read_csv(
  "contig_arg_classif.csv",
  col_types = cols(.default = "c")
)

df_sample <- read_delim(
  "sample_collection.csv", 
  delim = ";", escape_double = FALSE, trim_ws = TRUE)
df_sample <- df_sample[1:5]


# Part 2: Data cleaning and formatting ------------------------------------

# Join samples metadata (Bioproject, Fermentation time, Type of sequencer, Origin)
# PS : 'NA' in the 'species' column means that their was no match in BLAST for 
# the contig's identification (e-value threshold = 1e-20)

amrsample_met <- left_join(df_sample, amrsample, by="sample")
amrsample_met <- amrsample_met %>%
  relocate(sample, .before = tFermentation)

write.csv(amrsample_met,
          file = "amrsample_met.csv",
          row.names = FALSE)

# Counts the number of times each ARGs classied as 'Only for Humans' 
# appear in 1 sample

amrsample =  amrsample_met[c(1,2,3,4,5,6,8,12,13)]
df <- amrsample[amrsample$"WHO_classif"=="Authorized for use in humans only" |
                  amrsample$"WHO_classif"=="Authorized for use in humans only ?", 
                c("sample","contig","tFermentation","bioproject","sequenceur",
                  "origin","species","short_name")]
df <- df[!is.na(df$species) & !is.na(df$short_name),]
species_list <- unique(df$species)
short_names <- unique(df$short_name)

div = join[,c("species",grep("^ERR|^SRR", colnames(join), value = TRUE))]
col_sums <- colSums(div[, -1], na.rm = TRUE)
samples_keep <- names(col_sums[col_sums > 0])

div_filt <- div[, c("species", samples_keep)]
argcount = join[,c("species",short_names)]

join$genus <- sapply(
  strsplit(join$species, " "),
  function(x) x[1]
)
join$genus[join$genus == "other"] <- "Other"
table(join$genus)

count_arg = table(df$short_name,df$sample)
count_sample = table(df$sample,df$short_name)
count_arg_df <- as.data.frame(count_arg)
colnames(count_arg_df) <- c("short_name", "sample", "count")
count_arg_df <- count_arg_df[count_arg_df$count > 0.5, ] 
count_arg_dfs <- left_join(df[1:6], count_arg_df, by="sample","short_name")
count_arg_dfs <- df %>%
  count(sample, tFermentation, bioproject, sequenceur, origin, short_name, name = "count")
summary(count_arg_dfs)

# Part 3: Sankey plot, ARG -> Samples --------------------------------------------

arg_samples <- count_arg_dfs %>%
  transmute(
    from  = short_name,
    to    = sample,
    value = count
  )

# Combine links
links <- bind_rows(arg_samples)

# Node table
nodes <- data.frame(
  name = unique(c(links$from, links$to))
)

# Colors by Bioproject
links$source <- match(links$from, nodes$name) - 1
links$target <- match(links$to,   nodes$name) - 1

links$group <- count_arg_dfs$bioproject
my_color <- 'd3.scaleOrdinal() .domain(["PRJEB65292", "PRJEB77409","PRJNA388572"]) 
.range(["red", "steelblue", "green"])'

sankey4= tagList(
  tags$h4("Sankey diagram: ARG → Samples"),
  tags$h6("Transferable ARG classified as 'Only for Humans'"),
  tags$h6(" red : PRJEB65292, blue : PRJEB77409, green : PRJNA388572"),
  
  sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value  = "value",
    NodeID = "name",
    fontSize = 10,
    nodeWidth = 1,
    nodePadding = 10,
    colourScale = my_color,
    LinkGroup = "group",
    iterations = 32
  )
)
save_html(sankey4, "sankey_arg_samples.html")

# Colors by fermentation 
links$groupferm <- trimws(count_arg_dfs$tFermentation)
links$groupferm <- gsub("\\s+", " ", links$groupferm)

links$groupferm <- recode(links$groupferm,
                          "24h every 1-3 day for 1 week" = "1w",
                          "24h every 1-3 day for 5 week" = "5w",
                          "24h every 1-3 day for 9 week" = "9w",
                          "24h every 1-3 day for 13 week" = "13w",
                          "24h every 1-3 day for 17 week" = "17w",
                          "24h every 1-3 day for 21 week" = "21w"
)


my_colorferm <- '
d3.scaleOrdinal()
  .domain([
    "8",
    "24",
    "48",
    "1w",
    "5w",
    "9w",
    "13w",
    "17w",
    "21w"
  ])
  .range([
    "#1f77b4",  // bleu
    "#ff7f0e",  // orange
    "#2ca02c",  // vert
    "#d62728",  // rouge
    "#9467bd",  // violet
    "#8c564b",  // marron
    "#e377c2",  // rose
    "#7f7f7f",  // gris
    "#17becf"   // cyan
  ])
'
sankey5= tagList(
  tags$h4("Sankey diagram: ARG → Samples"),
  tags$h6("Transferable ARG classified as 'Only for Humans'"),
  tags$div(
    style = "margin-bottom: 20px;",
    tags$h6("Fermentation time:"),
    
    tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 10px; margin-top: 10px;",
      
      tags$span(style = "color:#1f77b4;", "■ 8"),
      tags$span(style = "color:#ff7f0e;", "■ 24"),
      tags$span(style = "color:#2ca02c;", "■ 48"),
      tags$span(style = "color:#d62728;", "■ 1w"),
      tags$span(style = "color:#9467bd;", "■ 5w"),
      tags$span(style = "color:#8c564b;", "■ 9w"),
      tags$span(style = "color:#e377c2;", "■ 13w"),
      tags$span(style = "color:#7f7f7f;", "■ 17w"),
      tags$span(style = "color:#17becf;", "■ 21w")
    )
  ),
  
  sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value  = "value",
    NodeID = "name",
    fontSize = 10,
    nodeWidth = 1,
    nodePadding = 10,
    colourScale = my_colorferm,
    LinkGroup = "groupferm",
    iterations = 32
  )
)
save_html(sankey5, "sankey_arg_samples_ferm.html")

# Colors by sequencer
links$groupseq <- count_arg_dfs$sequenceur
links$groupseq <- recode(links$groupseq,
                          "ILLUMINA, NextSeq 500, paired-end" = "500",
                          "Illumina HiSeq 2500, paired-end" = "2500",
                          "Illumina HiSeq 4000, paired-end" = "4000",
                          "Illumina NovaSeq 6000, paired-end" = "6000"
)

my_colorseq <- 'd3.scaleOrdinal() .domain(["500", "2500", "4000","6000"]) .range(["red", "black", "green","steelblue"])'

sankey6= tagList(
  tags$h4("Sankey diagram: ARG → Samples"),
  tags$h6("Transferable ARG classified as 'Only for Humans'"),
  tags$h6(" red : Illumina NextSeq 500, blue : Illumina NovaSeq 600, green : Illumina HiSeq 4000, black : Illumina HiSeq 2500"),
  
  sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value  = "value",
    NodeID = "name",
    fontSize = 10,
    nodeWidth = 1,
    nodePadding = 10,
    colourScale = my_colorseq,
    LinkGroup = "groupseq",
    iterations = 32
  )
)
save_html(sankey6, "sankey_arg_samplesseq.html")

# Part 4: Descriptive statistics ------------------------------------------

# ARGs count
somme <- rowSums(count_arg)
df_sum <-as.data.frame(somme)
df_sum %>% arrange(desc(somme))

# Samples count per ARG
count_samples <- count_arg_df  %>%
  group_by(short_name) %>%
  summarise(n_samples = n_distinct(sample)) %>%
  arrange(desc(n_samples))

# Species count per samples
arg_species_counts <- colSums(argcount[,-1] > 0)
arg_species_counts <- data.frame(
  ARG = names(arg_species_counts),
  n_species = arg_species_counts
)
arg_species_counts <- arg_species_counts[order(-arg_species_counts$n_species), ]


# Some visual plot
boxplot(div_filt[-1])
hist(join$abondance_intersamples)
hist(join$abondance_intrasample)
plot(join$abondance_intrasample, join$O23S)
plot(join$abondance_intrasample, join$TEM)
plot(join$abondance_intrasample, join$OMP36)
plot(join$abondance_intrasample, join$ACRB)
plot(join$abondance_intrasample, join$ACRA)
plot(join$abondance_intersamples, join$O23S)
plot(join$abondance_intersamples, join$TEM)
plot(join$abondance_intersamples, join$OMP36)
plot(join$abondance_intersamples, join$ACRB)
plot(join$abondance_intersamples, join$ACRA)

barplot(count_arg, beside=TRUE, legend=FALSE)
barplot(count_sample, beside=TRUE, legend=FALSE)
chisq.test(count_sample) # 

