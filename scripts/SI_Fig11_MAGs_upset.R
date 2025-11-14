# SI Figure 11 MAGs Upset
# Load required libraries
library(tidyverse);library(paletteer);library(ComplexUpset)

# File paths
metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/selected-das-bin0.txt')
bin_metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/bin-metadata.txt')
gene_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/bins-cazymes.txt')
grouped_cazymes <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/cazyme-subgroup-for-MAGs.2025.tsv')

# Bin metadata
age_mapping <- bin_metadata %>% 
  select(mappedID, age, detail_age) %>% 
  distinct() %>%
  group_by(mappedID) %>%
  summarise(
    has_M = any(age == "M" | str_starts(detail_age, "M")),
    has_I = any(age == "I" | str_starts(detail_age, "I")),
    .groups = "drop"
  ) %>%
  mutate(age_category = case_when(
    has_M & has_I ~ "Transit",
    has_M ~ "M",
    has_I ~ "I",
    TRUE ~ "None"
  ))

# Define custom order
substrate_order <- c(
  "HMOs", "GH42", "GH20",  "Fucose", "Sialic acid", "Host glycans",
  "Inulin", "β-glucans", "Xylan-backbone", "Xylan-side-chains", "Pectin-backbone", "Pectin-side-chains")

# Define classes order
classes_order <-c ("Actinomycetes" = "#608cd8", 
                   "Coriobacteriia" = "#2941a8", 
                   "Bacilli" = "#c69c6d",
                   "Clostridia" = "#a26990", 
                   "Negativicutes" = "#a88e7c",
                   "Bacteroidia" = "#307b16", 
                   "Verrucomicrobiae" = "darkviolet",
                   "Gammaproteobacteria" = "coral4",
                   "Alphaproteobacteria" = "grey30",
                   "Brachyspirae" = "black")

# Clean and prepare metadata
split_metadata <- metadata %>%
  separate(taxa, into = c("kingdom", "phyla", "classes", "orders", "family", "genus", "species"), sep = ";", fill = "right") %>%
  mutate(across(c(kingdom, phyla, classes, orders, family, genus, species),
                ~ str_replace(., "^[a-z]__", "")),
         classes = str_replace(classes, "_[A-Z]$", ""),
         genus = ifelse(is.na(genus) | genus == "", "unassigned_genus", genus),
         species = ifelse(is.na(species) | species == "", "unassigned_species", species),
         prefix = paste(species, mappedID, sep = "_")) %>%
  left_join(age_mapping, by = "mappedID") %>% 
  distinct(mappedID, species, classes, prefix, age_category) %>% 
  mutate(classes = factor(classes, levels = names(classes_order)),
         age_category = factor(age_category, levels = c("I", "Transit", "M")))

# Filter and reshape gene abundance
for (age in unique(split_metadata$age_category)) {
  filtered_abundance <- gene_abundance %>%
    filter(mappedID %in% (split_metadata %>% filter(age_category == age) %>% pull(mappedID))) %>% 
    select(mappedID, any_of(grouped_cazymes$cazymes)) %>%
    pivot_longer(-mappedID, names_to = "CAZyme", values_to = "Abundance") %>%
    left_join(grouped_cazymes, by = c("CAZyme" = "cazymes")) %>% 
    filter(substrates %in% substrate_order) %>% 
    group_by(substrates, mappedID) %>% 
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Presence/absence matrix per species
  species_matrix <- filtered_abundance %>%
    left_join(split_metadata %>% filter(age_category == age), by = "mappedID") %>%
    mutate(present = Abundance > 0) %>% 
    select(-Abundance) %>% 
    group_by(prefix, substrates) %>%
    pivot_wider(names_from = substrates, values_from = present, values_fill = FALSE) %>%
    ungroup() 
  # Upset data
  upset_data <- species_matrix %>%
    select(prefix, all_of(substrate_order)) %>%
    mutate(across(all_of(substrate_order), ~ replace_na(., FALSE))) %>%
    group_by(prefix) %>%
    summarise(across(all_of(substrate_order), ~ as.integer(any(.))), .groups = "drop") %>%
    mutate(None = as.integer(rowSums(across(all_of(substrate_order))) == 0)) %>%  
    left_join(
      split_metadata %>% select(prefix, classes) %>% distinct(),
      by = "prefix"
    ) %>% 
    column_to_rownames("prefix")
  
  # Plot
  p <- upset(
    upset_data,
    intersect = rev(substrate_order),
    name = "Substrates",
    sort_sets = FALSE,
    sort_intersections_by = 'degree',
    base_annotations = list(
      'Intersection size' = intersection_size(
        aes(fill = classes)
      ) + scale_fill_manual(values = classes_order)
    ), 
    width_ratio = 0.2
  )  + ggtitle(paste("Age group:", age))
  
  print(p)
  
  # Save plot with age in filename
  #ggsave(filename = paste0("figures/complex_upset_plot2_", age, ".svg"),plot = p, width = 14, height = 10, bg = "transparent")
}

compiled <- gene_abundance %>%
  filter(mappedID %in% unique(split_metadata$mappedID)) %>% 
  select(mappedID, any_of(grouped_cazymes$cazymes)) %>%
  pivot_longer(-mappedID, names_to = "CAZyme", values_to = "Abundance") %>%
  left_join(grouped_cazymes, by = c("CAZyme" = "cazymes")) %>% 
  filter(substrates %in% substrate_order) %>% 
  select(-c(substrates, groups, subgroups, Figure_1c, Figure_1e)) %>% 
  pivot_wider(id_cols = c(mappedID), names_from = CAZyme, values_from = Abundance) %>% 
  left_join(split_metadata, by = "mappedID")

write.table(compiled, "MAGs.tsv", 
            sep = "\t",        
            quote = FALSE,     
            row.names = FALSE)
