# 7/05/2026
# Generating Sankey plot : ARG -> Species, Species -> Samples, ARG -> Species -> Samples
# Out : sankey_arg_species.html
# Out : sankey_species_sample.html
# Out : sankey_ARG_species_samples.html


# Part 1: Importing the necessary packages and data -------------- --------

library(dplyr)
library(readr)
library(tidyr)
library(networkD3)

joint <- read_csv(
  "joint.csv",
  col_types = cols(.default = "c")
)

# Part 2: Sankey plot, ARG -> Species -----------------------------------

arg_cols <- c(
  "O23S","ACRB","TEM","RPOCL","BC","BLA1","BLAZ","BCII","CTX","CRP","ACRA",
  "ACRR","ROBA","ACRE","ACRF","ACRS","MARA","MARR","EVGA","EVGS","HNS",
  "OMPA","OMPFB","PBP4B","OMP36","BLAC","OQXA","CMY","OMP37","OMPK36",
  "LIAFSR","RAHN","OQXB","OMPF","GADW","GADX","MDTE","MDTF","MDEA",
  "OXY","RAMR","MIR","LAQ")

sample_cols <- grep("^ERR|^SRR", colnames(join), value = TRUE)

# Links ARG -> species
arg_species <- joint[,c("species", all_of(arg_cols))] %>%
  pivot_longer(
    cols = all_of(arg_cols),
    names_to  = "ARG",
    values_to = "value"
  ) %>%
  filter(value > 0) %>%
  transmute(
    from  = ARG,
    to    = species,
    value = value
  )

# Combine links
links <- bind_rows(arg_species)

# Node table
nodes <- data.frame(
  name = unique(c(links$from, links$to))
)

links$source <- match(links$from, nodes$name) - 1
links$target <- match(links$to,   nodes$name) - 1

# Plot
sankey1= tagList(
  tags$h4("Sankey diagram: ARG → Species"),
  tags$h6("Transferable ARG classified as 'Only for Humans'"),
  tags$h6("Species : identification with Blast (E-Threshold = 1e-20)"),
  
  sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value  = "value",
    NodeID = "name",
    fontSize = 7,
    nodeWidth = 40
  )
)
save_html(sankey1, "sankey_arg_species.html")


# Part 2: Sankey plot, Species -> Samples ---------------------------------

# Links species -> samples
species_samples_filt <- joint[,c("species", all_of(sample_cols))] %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to  = "sample",
    values_to = "value"
  ) %>%
  filter(!is.na(value), value > 10) %>%
  transmute(
    from  = species,
    to    = sample,
    value = value
  )

# Creation of an artificial knot for species without samples
species_with_sample <- unique(species_samples_filt$from)
all_species <- unique(arg_species$to)

species_without_sample <- setdiff(all_species, species_with_sample)
dummy_node <- "NO_SAMPLE_inf0.005"

dummy_links <- data.frame(
  from  = species_without_sample,
  to    = dummy_node,
  value = 1e-6
)

# Combine links
links <- bind_rows(species_samples_filt, dummy_links)

# Node table
nodes <- data.frame(
  name = unique(c(links$from, links$to))
)

links$source <- match(links$from, nodes$name) - 1
links$target <- match(links$to,   nodes$name) - 1

# Plot
sankey2= tagList(
  tags$h3("Sankey diagram: Species → Samples"),
  tags$h6("Species present are >= 10% normalized diversity"),
  
  sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value  = "value",
    NodeID = "name",
    fontSize = 7,
    nodeWidth = 40
  )
)
save_html(sankey2, "sankey_species_sample.html")


# Part 3: Sankey plot, ARG -> Species -------------------------------------

# Links ARG -> species
arg_species <- joint[,c("species", all_of(arg_cols))]%>%
  pivot_longer(
    cols = all_of(arg_cols),
    names_to  = "ARG",
    values_to = "value"
  ) %>%
  filter(value > 0) %>%
  transmute(
    from  = ARG,
    to    = species,
    value = value
  )

# Links species -> sample
species_sample <- joint[,c("species", all_of(sample_cols))] %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to  = "sample",
    values_to = "value"
  ) %>%
  filter(value > 0) %>%
  transmute(
    from  = species,
    to    = sample,
    value = value
  )

# Top ARGs
arg_species$value <- as.numeric(arg_species$value)
typeof(arg_species$value)
top_args <- arg_species %>%
  group_by(from) %>%
  summarise(total = sum(value)) %>%
  slice_max(total, n = 10) %>%
  pull(from)

arg_species_filt <- arg_species %>%
  filter(from %in% top_args)

# Top species
top_species <- arg_species_filt %>%
  group_by(to) %>%
  summarise(total = sum(value)) %>%
  slice_max(total, n = 10) %>%
  pull(to)

arg_species_filt <- arg_species_filt %>%
  filter(to %in% top_species)

# Top samples
species_sample$value <- as.numeric(species_sample$value)
typeof(arg_species$value)
top_samples <- species_sample %>%
  group_by(to) %>%
  summarise(total = sum(value)) %>%
  slice_max(total, n = 20) %>%
  pull(to)

species_sample_filt <- species_sample %>%
  filter(to %in% top_samples)

# Creation of an artificial knot for species without samples
species_with_sample <- unique(species_sample_filt$from)
all_species <- unique(arg_species_filt$to)

species_without_sample <- setdiff(all_species, species_with_sample)
dummy_node <- "NO_SAMPLE"

dummy_links <- data.frame(
  from  = species_without_sample,
  to    = dummy_node,
  value = 1e-6
)

# Combine links
links <- bind_rows(arg_species_filt, species_sample_filt, dummy_links)

# Knot selection
nodes <- data.frame(
  name = unique(c(links$from, links$to))
)

nodes$level <- ifelse(
  nodes$name %in% arg_species$from, 1,
  ifelse(nodes$name %in% arg_species$to, 2, 3)
)
# Node table
nodes <- data.frame(
  name = unique(c(links$from, links$to))
)

links$source <- match(links$from, nodes$name) - 1
links$target <- match(links$to,   nodes$name) - 1

# Plot
sankey3= tagList(
  tags$h3("Sankey diagram: ARG → Species → Samples"),
  tags$h6("10 top ARG, 10 top species, 20 top samples"),
  
  sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value  = "value",
    NodeID = "name",
    fontSize = 7,
    nodeWidth = 40
  )
)
save_html(sankey3, "sankey_ARG_species_samples.html")


