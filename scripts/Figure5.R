### Figure 5
### Author: Yunjeong So

# Figure5 a. Bubble plot------------------------------------------------------
# Load necessary libraries
library(tidyverse)
library(RColorBrewer)

# Load the metadata and genome abundance data
metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered2.tsv')
genome_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/meta-wholetaxa-converted.txt')

# Define substrate list
sub_list <- c("F", "AX", "PEC", "XG")

# Prepare metadata
metadata <- metadata %>%
  filter(substrate %in% sub_list, infant == "I") %>% 
  mutate(
    participant = as.factor(participant),
    prefix = paste(participant, timepoint, substrate, sep = "_"), 
    newID = paste(participant, timepoint, sep = "_")
  ) %>%
  arrange(timepoint, match(substrate, sub_list), participant) %>%  
  mutate(prefix = factor(prefix, levels = unique(paste(participant, timepoint, substrate, sep = "_")))) 

# Define class order and key genus
key_classes <- c(
  "Actinomycetia-others",  "Actinomycetia", "Coriobacteriia", "Bacilli", "Bacilli-others", "Clostridia", "Clostridia-others", "Negativicutes", "Negativicutes-others",
  "Bacteroidia", "Bacteroidia-others", "Gammaproteobacteria", "Gammaproteobacteria-others", "Others", "Others-others", "Unclassified", "Unclassified-others")

# Reshape genome abundance data to long format
genome_abundance_long <- genome_abundance %>%
  pivot_longer(-ID, names_to = "taxa", values_to = "genome_abundance") %>%
  filter(genome_abundance >= 0.01, ID %in% metadata$ID) %>% 
  separate(taxa, into = c("kingdom", "phyla", "classes", "orders", "family", "genus", "species"),
           sep = ";", fill = "right") %>%
  mutate(across(kingdom:species, ~str_remove(.x, "^[dpcofgs]__"))) %>%
  mutate(
    genus = str_remove(genus, "_.*"),
    species = str_remove(species, "_[A-Za-z]+"),
    genus = if_else(is.na(genus) | genus == "", "unassigned_genus", genus),
    species = if_else(is.na(species) | species == "", "unassigned_species", species)
  )

# Merge metadata and genome abundance data
merged_data <- metadata %>%
  left_join(genome_abundance_long, by = "ID") %>%
  mutate(
    timepoint = factor(timepoint, levels = c('T1', 'T2', 'T3'))
  )

# Filter faeces
F_data <- merged_data %>%
  filter(substrate == "F") %>%
  group_by(newID, classes, genus, species) %>%  
  summarise(F_abundance = sum(genome_abundance, na.rm = TRUE)) %>%
  ungroup()

# Filter HMOs and mucin group
non_F_data <- merged_data %>%
  filter(substrate %in% sub_list) %>%
  group_by(prefix, newID, participant, timepoint, classes, genus, species, substrate) %>%
  summarise(non_F_abundance = sum(genome_abundance, na.rm = TRUE)) %>%
  ungroup()

# Merge faeces and non-faeces data for fold change calculation
merged_subset <- non_F_data %>% 
  left_join(F_data, by = c("newID", "classes", "genus", "species")) %>% 
  mutate(
    F_abundance = coalesce(F_abundance, 0),
    non_F_abundance = coalesce(non_F_abundance, 0))

# Calculate fold change (log2 transformation)
filtered_subset <- merged_subset %>%
  mutate(
    non_F_abundance = round(as.numeric(non_F_abundance), 2),  
    F_abundance = round(as.numeric(F_abundance), 2),  
    fold_change = round(log2(non_F_abundance + 1) - log2(F_abundance + 1), 2),
    substrate = factor(substrate, levels = sub_list),
    classes = factor(classes, levels = key_classes)
  ) 

# Species list
filter_list <- filtered_subset %>%
  filter(fold_change >= 1, non_F_abundance >= 10) %>%
  distinct(species) %>%
  pull(species)

# Filtered subset
filtered_subset2 <- merged_subset %>%
  mutate(
    classes = if_else(classes %in% key_classes, classes, 'Others'),
    grouping = case_when(
      genus == "Enterococcus" ~ paste0(genus),
      classes == "Bacilli" & genus == "Clostridium" ~ paste0('z-', classes),
      classes == "Negativicutes"  ~ paste0(classes),
      !species %in% filter_list ~ paste0('z-', classes),
      classes == "Gammaproteobacteria" ~ paste0('z-', classes),
      TRUE ~ species)) %>%
  group_by(prefix, classes, substrate, grouping) %>%
  summarise(
    non_F_abundance = sum(non_F_abundance),
    F_abundance = sum(F_abundance),
    .groups = "drop") %>%
  mutate(
    fold_change = log2(non_F_abundance + 1) - log2(F_abundance + 1),
    classes = factor(classes, levels = key_classes),
    substrate = factor(substrate, levels = sub_list)) %>%
  arrange(classes, grouping) %>% 
  mutate(grouping = factor(grouping, levels = rev(unique(grouping))))


# Define the color scale with grey80 at the midpoint
my_colours <- c("#003171", "white", "#ae3127")
fold_change_min <- min(filtered_subset$fold_change, na.rm = TRUE)
fold_change_max <- max(filtered_subset$fold_change, na.rm = TRUE)
fold_change_range <- c(fold_change_min, 0, fold_change_max)

# Plot with substrate and participant facets
p <- ggplot(filtered_subset2, aes(x = prefix, y = grouping, size = non_F_abundance, fill = fold_change)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_size(range = c(0.5, 6), limits = c(0.5, max(merged_subset$non_F_abundance))) +
  scale_fill_gradientn(
    colours = my_colours,
    limits = c(fold_change_min, fold_change_max),
    values = scales::rescale(fold_change_range, to = c(0, 1))
  ) +
  theme_minimal() +
  labs(size = "Relative abundance (%)", fill = "Log2 fold change") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "white", color = "black", size = 1)
  ) 

# Display the plot
print(p)

# Save the bubble plot
#ggsave(plot = p, filename = "figures/plant-species-bubbles_FC1_RA10.svg", width = 12, height = 4.5, bg = "transparent")


# Figure5 b. CAZymes gene enrichment ------------------------------------------------------
library(data.table);library(tidyverse);library(ggpubr)

# Read data
metadata <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered2.tsv')
gene_abundance <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/gene-abundance-nonlog-table.txt')
mapping <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sample_mapping_all.txt')
grouped_cazymes <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/cazyme-subgroup-for-box.2025.tsv', header = TRUE)

# Define colors
AX_colors <- c(
  "GH10" = "red",
  "GH11" = "#87789c",
  "GH5_21" = "#ed9897",
  "GH5_35" = "#994D78",
  "Xylan-side-chains" = "#e4ddd4",
  "xylan/cellulose/b-glucans" = "#8b5d2f")
PEC_colors <- c(
  "GH28" = "#6f8520",
  "GH105" = "#989279",
  "PL" = "#006400",
  "Pectin-side-chains" = "#E7E178")
XG_colors <- c(
  "GH74" = "#1C0CA1",
  "GH5_4" = "#5599cc",
  "GH12" = "#6FE8FD",
  "xylan/cellulose/b-glucans" = "#8b5d2f")


# Function 
extract_substrate <- function(target) {
  target_list <- switch(target,
                        "AX"  = c("F", "AX"),
                        "PEC" = c("F", "PEC"),
                        "XG"  = c("F", "XG"),
                        stop("target must be one of: 'AX', 'PEC', or 'XG'"))
  
  target_colors <- switch(target,
                          "AX"  = AX_colors,
                          "PEC" = PEC_colors,
                          "XG"  = XG_colors)
  filtered_metadata <- metadata %>%
    left_join(mapping, by = "ID") %>% 
    filter(infant == "I", substrate %in% target_list) %>% 
    mutate(
      participant = as.factor(participant),
      substrate = factor(substrate, levels = target_list),
      prefix = paste(participant, timepoint, substrate, sep = "_"), 
      newID = paste(timepoint, participant, sep = "_")
    ) %>%
    arrange(timepoint, substrate, participant) %>%  
    mutate(prefix = factor(prefix, levels = unique(paste(participant, timepoint, substrate, sep = "_"))))
  
  # filter gene
  filtered_abundance <- gene_abundance %>%
    filter(mappedID %in% filtered_metadata$mappedID) %>%
    select(mappedID, starts_with(c("GH", "PL"))) 
  
  # Make long
  abundance_long <- filtered_abundance %>%
    pivot_longer(cols = -mappedID, names_to = "cazymes", values_to = "abundance") %>% 
    merge(filtered_metadata, by = "mappedID") %>% 
    merge(grouped_cazymes, by = "cazymes") %>%
    mutate(
      abundance = abundance,
      timepoint = factor(timepoint, levels = c("T1", "T2", "T3"))) %>% 
    filter(Figure_5 %in% names(target_colors))
  
  # Make wide
  abundance_wide <- abundance_long %>%
    select(cazymes, newID, Figure_5, substrate, abundance) %>%
    pivot_wider(names_from = substrate, values_from = abundance) %>%
    filter(newID %in% filtered_metadata$newID[filtered_metadata$substrate == target]) %>% 
    replace_na(list(F = 0, AX = 0, PEC = 0, XG = 0)) %>%
    mutate(change = (get(target) - F),
           prefix = paste0(
             target, "_", sub("_.*$", "", newID), "_",  sub(".*_", "", newID)))
  
  sum_data <- abundance_wide %>%
    group_by(prefix, newID, Figure_5) %>%
    filter(change > 0) %>%
    summarise(sum_change = sum(change, na.rm = TRUE), .groups = "drop") %>%
    mutate(Figure_5 = factor(Figure_5, levels = names(target_colors)))
  
  sum_data
}

merged_data <- bind_rows(extract_substrate("AX"), extract_substrate("PEC"), extract_substrate("XG"))

target_colors <- c(AX_colors, PEC_colors, XG_colors)

p <- ggplot(merged_data, aes(x = prefix, y = (sum_change), fill = Figure_5)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.85, color = NA) +
  geom_hline(yintercept = 0, color = "black", size = 0.25) +
  scale_fill_manual(values = target_colors) +
  labs(
    x = "Sample (newID)",
    fill = "CAZyme group"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = "black"),  
    panel.grid = element_blank(), 
    axis.ticks = element_blank(), 
    axis.text.y = element_text(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line.y = element_blank(),
    legend.position = "right"
  )

print(p)

#ggsave("figures/new-plants.svg", p, width = 10, height = 3, bg = "transparent")
