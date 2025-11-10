### Figure 4
### Author: Yunjeong So

# Figure4 a. Bubble plot ------------------------------------------------------
# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(tidyverse)

# Load data
metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered2.tsv')
taxa_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/meta-wholetaxa-converted.txt')

# Add 'prefix' and 'newID' to metadata
metadata <- metadata %>%
  mutate(prefix = paste(participant, timepoint, substrate, sep = "_")) %>% 
  mutate(newID = paste(participant, timepoint, sep = "_"))

# Define key classes and genus
key_classes <- c("Actinomycetia", "Bacilli", "Clostridia", "Negativicutes", "Bacteroidia", "Gammaproteobacteria", "Others", "Unclassified")

key_genus <- c(
  "Collinsella",
  "Enterococcus", "Erysipelatoclostridium", "Limosilactobacillus",  "Streptococcus", 
  "Blautia", "Butyribacter", "Clostridium", 
  "Dorea", "Enterocloster", "Eubacterium", "Faecalibacterium", "Lachnospira", "Paraclostridium", "Roseburia", 
  "Ruminococcus", "Megamonas", "Veillonella", "Bacteroides", "Phocaeicola"
)

# Separate the taxonomic levels from 'taxa'
taxa_abundance_long <- taxa_abundance %>%
  pivot_longer(-ID, names_to = "taxa", values_to = "taxa_abundance") %>%
  filter(taxa_abundance >= 0.01, ID %in% metadata$ID) %>% 
  separate(taxa, into = c("kingdom", "phyla", "classes", "orders", "family", "genus", "species"),
           sep = ";", fill = "right") %>%
  mutate(across(kingdom:species, ~str_remove(.x, "^[dpcofgs]__"))) %>%
  mutate(
    genus = str_remove(genus, "_.*"),
    species = str_remove(species, "_[A-Za-z]+"),
    genus = if_else(is.na(genus) | genus == "", "unassigned_genus", genus),
    species = if_else(is.na(species) | species == "", "unassigned_species", species)
  ) %>% 
  filter(taxa_abundance >= 0.01) %>% 
  mutate(
    classes = ifelse(classes %in% key_classes, classes, 'Others'),
    grouping = ifelse(
      genus == "Bifidobacterium", species,
      ifelse(
        genus == "Clostridium" & classes == "Bacilli", paste0(classes, '-others'),
        ifelse(genus %in% key_genus, genus, paste0(classes, '-others'))
      )))

# Extract unique Bifidobacterium species
bifido_species <- taxa_abundance_long %>% 
  filter(genus == "Bifidobacterium") %>% 
  pull(species) %>% 
  unique()

# Define orders
order_grouping <- c(
  "Bifidobacterium longum.infantis", "Bifidobacterium longum.longum", "Bifidobacterium bifidum", "Bifidobacterium breve", 
  "Bifidobacterium catenulatum", "Bifidobacterium dentium", "Bifidobacterium pseudocatenulatum",  
  "Collinsella", "Actinomycetia-others", 
  "Enterococcus", "Erysipelatoclostridium", "Limosilactobacillus", "Streptococcus", "Bacilli-others",
  "Agathobacter", "Blautia", "Butyribacter", "Clostridium", "Copromonas",
  "Dorea", "Enterocloster", "Eubacterium", "Faecalibacterium", 
  "Lachnospira", "Mediterraneibacter", "Paraclostridium", "Roseburia", "Ruminococcus", "Clostridia-others",
  "Megamonas", "Veillonella", "Negativicutes-others", 
  "Bacteroides", "Phocaeicola", "Bacteroidia-others", 
  "Gammaproteobacteria-others", "Others-others", "Unclassified-others"
)

# Merge metadata and reshaped taxa abundance data
merged_data <- merge(metadata, taxa_abundance_long, by = "ID")
merged_data$timepoint <- factor(merged_data$timepoint, levels = c('T1', 'T2', 'T3'))

# Filter faeces
F_data <- merged_data %>%
  filter(infant == "I", substrate == "F") %>%
  group_by(newID, grouping) %>%  
  summarise(F_abundance = sum(taxa_abundance, na.rm = TRUE)) %>%
  ungroup()

# Filter HMOs and mucin group
non_F_data <- merged_data %>%
  filter(infant == "I", substrate %in% c("HMOs", "MUC")) %>%
  group_by(prefix, newID, participant, timepoint, classes, grouping, substrate) %>%
  summarise(non_F_abundance = sum(taxa_abundance, na.rm = TRUE)) %>%
  ungroup()

# Merge faeces and non-faeces data for fold change calculation
merged_subset <- non_F_data %>% 
  left_join(F_data, by = c("newID", "grouping")) %>% 
  mutate(
    F_abundance = coalesce(F_abundance, 0),
    non_F_abundance = coalesce(non_F_abundance, 0)
  )

# Calculate fold change (log2 transformation)
filtered_subset <- merged_subset %>%
  mutate(
    non_F_abundance = round(as.numeric(non_F_abundance), 2),  
    F_abundance = round(as.numeric(F_abundance), 2),  
    fold_change = round(log2(non_F_abundance + 1) - log2(F_abundance + 1), 2)
  ) 

# Reorder substrates in the desired order
sub_list <-  c("HMOs", "MUC")
filtered_subset$substrate <- factor(filtered_subset$substrate, levels = sub_list)

# Set classes as a factor with levels defined by 'order_class'
filtered_subset <- filtered_subset %>%
  mutate(classes = factor(classes, levels = key_classes))

# Arrange by classes and grouping
sorted_subset <- filtered_subset %>%
  arrange(classes, grouping) %>%
  mutate(grouping = factor(grouping, levels = rev(order_grouping)))

# Define the color scale with grey80 at the midpoint
my_colours <- c("#003171", "white", "#ae3127")

# Min and max values for abundance and fold change
abundance_min <- 0.01
abundance_max <- 100
fold_change_min <- min(sorted_subset$fold_change, na.rm = TRUE)
fold_change_max <- max(sorted_subset$fold_change, na.rm = TRUE)
fold_change_range <- c(fold_change_min, 0, fold_change_max)

# Plot with substrate and participant facets
p <- ggplot(sorted_subset, aes(x = timepoint, y = grouping, size = non_F_abundance, fill = fold_change)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +  
  scale_size(range = c(0.5, 6), limits = c(abundance_min, abundance_max)) +  # Fix bubble size range
  scale_fill_gradientn(
    colours = my_colours, 
    limits = c(fold_change_min, fold_change_max), 
    values = scales::rescale(fold_change_range, to = c(0, 1))
  ) +
  #geom_hline(yintercept = seq(1.5, length(unique(sorted_subset$grouping)) - 0.5, by = 1), 
  #           color = "grey70", linewidth = 0.2) +
  theme_minimal() +
  labs(x = "Timepoint", y = "Genus/Class/Species", size = "Non-F Abundance", fill = "Fold Change") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add tile outline
    strip.background = element_rect(fill = "white", color = "black", size = 1)
  ) +
  facet_grid(substrate ~ participant)

# Display the plot
print(p)

# Save the bubble plot
filename <- paste0("figures/infants_genus_class_abundance_HMOs-MUC.svg")
#ggsave(plot = p, filename = filename, width = 8.5, height = 11.5, bg = "transparent")


# Figure4 c,d Bar plots of CAZymes------------------------------------------------------
# Load necessary libraries
library(data.table);library(tidyverse)

# Read data
metadata <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered2.tsv')
gene_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/gene-abundance-nonlog-table.txt')
mapping <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sample_mapping_all.txt')
grouped_cazymes <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/cazyme-subgroup-for-box.2025.tsv')

# Merge and process data
metadata <- metadata %>%
  mutate(prefix = paste(infant, participant, timepoint, sep = "_")) %>%
  left_join(mapping, by = "ID")

# Set cazymes and substrates
substrate_list <- c("F", "HMOs", "MUC")
HMOs_list <- c(
  "GH112" = "#ffa500",
  "GH136" = "dodgerblue",
  "GH20" = "pink",
  "GH42" = "#aa98a9",
  "GH29" = "#a12232",
  "GH95" = "#af593e",
  "GH33" = "purple4")

# Filter metadata
filtered_metadata <- metadata %>%
  filter(substrate %in% substrate_list & infant == "I") 

# Filter abundance
filtered_abundance <- gene_abundance %>%
  filter(mappedID %in% filtered_metadata$mappedID) %>%
  select(mappedID, any_of(names(HMOs_list))) 

# Use the extracted cazymes levels to order the factor
abundance_long <- filtered_abundance %>%
  pivot_longer(cols = -mappedID, names_to = "cazymes", values_to = "abundance") %>% 
  merge(filtered_metadata, by = "mappedID") %>% 
  merge(grouped_cazymes, by = "cazymes") %>%
  mutate(
    abundance = abundance,
    timepoint = factor(timepoint, levels = c('T1', 'T2', 'T3')),
    cazymes = factor(cazymes, levels = names(HMOs_list)),
    participant = factor(participant, levels = c("A", "B", "C", "D", "E", "F", "G")),
    timepoint_numeric = as.numeric(timepoint))

# Update the plot to use the custom colors
p <- ggplot(abundance_long, aes(x = prefix, y = abundance, fill = cazymes)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.85, color = "black", linewidth = 0.1) +
  labs(y = "GPM") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = "black"),  
    panel.grid = element_blank(), 
    plot.title = element_text(hjust = 0.5),
    axis.ticks = element_blank(), 
    axis.text.y = element_text(),
    axis.line.y = element_blank(),
    legend.position = "right"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000)) +
  scale_fill_manual(values = HMOs_list) +  
  facet_wrap(~ substrate, ncol = 1, scales = "free_y")

# Display the plot
print(p)

# Save the plot
filename <- paste0("figures/bar_infants_HMOs.svg")
#ggsave(plot = p, filename = filename, width = 8, height = 5, bg = "transparent")


