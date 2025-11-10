### Figure 6
### Author: Yunjeong So

# Figure6 b. Bubble plot------------------------------------------------------
# Load necessary libraries
library(tidyverse)
library(RColorBrewer)

# Load the metadata and taxa abundance data
metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered2.tsv')
taxa_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/meta-wholetaxa-converted.txt')

# Define substrate list
sub_list <- c("F", "HMOs", "MUC")

# Prepare metadata
metadata <- metadata %>%
  filter(substrate %in% sub_list, infant == "M", timepoint == "T2") %>% 
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
  "Bacteroidia", "Bacteroidia-others", "Verrucomicrobiae", "Gammaproteobacteria", "Gammaproteobacteria-others", "Others", "Others-others", "Unclassified", "Unclassified-others")


# Reshape taxa abundance data to long format
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

# Merge metadata and taxa abundance data
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
  filter(fold_change >= 1, non_F_abundance >= 1) %>%
  distinct(species) %>%
  pull(species)

# Filter data
filtered_subset2 <- merged_subset %>%
  mutate(
    classes = if_else(classes %in% key_classes, classes, 'Others'),
    grouping = case_when(
      classes == "Bacilli" & genus == "Clostridium" ~ paste0('z-', classes),
      genus %in% c("Enterococcus", "Streptococcus") ~ paste0(genus),
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
  scale_size(range = c(0.5, 6), limits = c(0.5, max(filtered_subset2$non_F_abundance))) +
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
#ggsave(plot = p, filename = "figures/2025.mother-bubbles.0.01.svg", width = 6, height = 8, bg = "transparent")


# Figure6 c. CAZymes gene enrichment ------------------------------------------------------
library(data.table)
library(tidyverse)

# Read data
metadata <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered.tsv')
gene_abundance <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/gene-abundance-nonlog-table.txt')
mapping <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sample_mapping_all.txt')
grouped_cazymes <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/cazymes_groups_mah.txt', header = T)

# Define substrate list
sub_list <- c("F", "HMOs", "MUC")

# Prepare metadata
metadata <- metadata %>%
  left_join(mapping, by = "ID") %>% 
  filter(infant == "M") %>% 
  mutate(
    participant = as.factor(participant),
    substrate = factor(substrate, levels = sub_list),  
    prefix = paste(participant, timepoint, substrate, sep = "_"), 
    newID = paste(participant, timepoint, sep = "_")
  ) %>%
  filter(!(prefix %in% c("A_T2_MUC", "B_T2_MUC", "C_T2_MUC", "G_T2_MUC"))) %>% 
  arrange(participant, timepoint, substrate) %>% 
  mutate(prefix = factor(prefix, levels = unique(paste(participant, timepoint, substrate, sep = "_")))) 


filtered_abundance <- gene_abundance %>%
  filter(mappedID %in% metadata$mappedID) %>%
  select(mappedID, starts_with(c("GH", "PL"))) 

filtered_cazymes <- setdiff(names(filtered_abundance), "mappedID")
unmapped_cazymes <- setdiff(filtered_cazymes, grouped_cazymes$cazymes)
unmapped_df <- data.frame(
  cazymes = unmapped_cazymes,
  groups = "Others",
  subgroups = "Others"
)

grouped_cazymes <- grouped_cazymes %>%
  mutate(groups = ifelse(groups == "Mixed", "Others", groups)) %>%
  bind_rows(unmapped_df) 

HMOs_list <- c("GH112", "GH136", "GH20", "GH42", "GH29", "GH95", "GH33")

# Filter cazymes
cazyme_levels <- grouped_cazymes %>%
  filter(subgroups %in% HMOs_list) %>%
  pull(cazymes)  

# Arrange data
abundance_long <- filtered_abundance %>%
  pivot_longer(cols = -mappedID, names_to = "cazymes", values_to = "abundance") %>% 
  merge(metadata, by = "mappedID") %>% 
  filter(substrate %in% c("F", "HMOs", "MUC")) %>% 
  merge(grouped_cazymes, by = "cazymes") %>%
  mutate(
    abundance = (abundance),
    timepoint = factor(timepoint, levels = c('T1', 'T2', 'T3')),
    cazymes = factor(cazymes, levels = cazyme_levels)) %>%
  filter(subgroups %in% HMOs_list) 

sum_data <- abundance_long %>%
  group_by(prefix, subgroups, timepoint) %>%
  summarise(abundance = (sum(abundance, na.rm = TRUE)), .groups = 'drop') %>%
  mutate(timepoint_numeric = as.numeric(timepoint)) %>% 
  mutate(
    subgroups = factor(subgroups, levels = (HMOs_list)))


# Define a custom color
subgroup_colors <- c(
  "GH112" = "#ffa500",
  "GH136" = "dodgerblue",
  "GH20" = "pink",
  "GH42" = "#aa98a9",
  "GH29" = "#a12232",
  "GH95" = "#af593e",
  "GH33" = "purple4",
  "GAGs" = "green",
  "O-glycans" ="darkgreen",
  "N-glycans" = "black"
)

# plot
p <- ggplot(sum_data, aes(x = prefix, y = abundance, fill = subgroups)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.85, color = "black", linewidth = 0.1) +
  labs(y = "(GPM)") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = "black"),  
    panel.grid = element_blank(), 
    plot.title = element_text(hjust = 0.5),
    axis.ticks = element_blank(), 
    axis.text.y = element_text(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line.y = element_blank(),
    legend.position = "right"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2500)) +
  scale_fill_manual(values = subgroup_colors)   

# Display the plot
print(p)

# Figure6 d. MAGs CAZYme heatmap------------------------------------------------------
# # Load required libraries 
library(tidyverse)
library(paletteer)
library(pheatmap)

# File paths 
metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/selected-das-bin0.txt')
bin_metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/bin-metadata.txt')
gene_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/bins-cazymes.txt')
grouped_cazymes <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/cazyme-subgroup-for-MAGs.2025.tsv')
metageome_metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered.tsv')
mp_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/meta-species-converted.txt')
mapping <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sample_mapping_all.txt')

# Define orderings 
substrate_order <- c("HMOs", "GH42", "GH20", "Fucose", "Sialic acid", "Host glycans", "Starch",
                     "Inulin", "β-glucans", "Xylan-backbone", "Xylan-side-chains",
                     "Pectin-backbone", "Pectin-side-chains", "Food additives", "β-glucans, Xylan-backbone, Chitosanase")

classes_order <- c("Actinomycetes" = "#608cd8", 
                   "Coriobacteriia" = "#2941a8", 
                   "Bacilli" = "#c69c6d",
                   "Clostridia" = "#a26990", 
                   "Negativicutes" = "#a88e7c",
                   "Bacteroidia" = "#307b16", 
                   "Verrucomicrobiae" = "darkviolet",
                   "Gammaproteobacteria" = "coral4",
                   "Alphaproteobacteria" = "grey30",
                   "Brachyspirae" = "black")

# Prepare age category 
age_mapping <- bin_metadata %>%
  select(mappedID, age, detail_age) %>%
  distinct() %>%
  group_by(mappedID) %>%
  summarise(
    has_M = any(age == "M" | str_starts(detail_age, "M")),
    has_I = any(age == "I" | str_starts(detail_age, "I")),
    .groups = "drop") %>%
  mutate(age_category = case_when(
    has_M & has_I ~ "Transit",
    has_M ~ "M",
    has_I ~ "I",
    TRUE ~ "None"))

# Prepare metadata
metageome_metadata <- metageome_metadata %>%
  mutate(prefix = paste(infant, participant, sep = "_")) %>%
  left_join(mapping, by = "ID")

md_faeces <- metageome_metadata %>%
  filter(substrate == "F") %>%
  select(ID, mappedID, infant) %>%
  left_join(mapping, by = "ID") %>%
  distinct()

# Reshape abundance table and clean species names
long_abundance <- mp_abundance %>%
  filter(ID %in% md_faeces$ID) %>%
  pivot_longer(-ID, names_to = "species", values_to = "abundance") %>%
  mutate(
    species = str_remove(species, "^s__"),
    species = case_when(
      species == "Bifidobacterium longum.infantis" ~ "Bifidobacterium infantis",
      species == "Bifidobacterium longum.longum" ~ "Bifidobacterium longum",
      TRUE ~ species)
  ) %>%
  filter(species != "Unclassified") %>%
  left_join(select(metageome_metadata, ID, infant), by = "ID")

# Sum abundance across all samples for infants and mothers separately
top_species <- long_abundance %>%
  mutate(group = if_else(str_detect(infant, "^I"), "Infant", "Mother")) %>%
  group_by(group, species) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>% 
  arrange(group, desc(total_abundance)) %>%
  group_by(group) %>% 
  slice_head(n = 25) %>%
  ungroup()

# Prepare split_metadata
split_metadata <- metadata %>%
  separate(taxa, into = c("kingdom", "phyla", "classes", "orders", "family", "genus", "species"), sep = ";", fill = "right") %>%
  mutate(across(c(kingdom, phyla, classes, orders, family, genus, species), ~ str_replace(., "^[a-z]__", "")),
         classes = str_replace(classes, "_[A-Z]$", ""),
         genus = ifelse(is.na(genus) | genus == "", "unassigned_genus", genus),
         species = ifelse(is.na(species) | species == "", "unassigned_species", species),
         prefix = paste(species, mappedID, sep = "_")) %>%
  filter(species %in% top_species$species | species == "Clostridium butyricum" | species == "Mediterraneibacter faecis") %>%
  left_join(age_mapping, by = "mappedID") %>%
  distinct(mappedID, species, classes, prefix, age_category) %>%
  mutate(classes = factor(classes, levels = names(classes_order)),
         age_category = factor(age_category, levels = c("I", "Transit", "M")))

# Filter gene abundance for selected MAGs
filtered_abundance <- gene_abundance %>%
  filter(mappedID %in% split_metadata$mappedID) %>%
  select(mappedID, any_of(grouped_cazymes$cazymes)) %>%
  mutate(across(-mappedID, as.numeric)) %>%
  # keep mappedID plus numeric cols with colSums > 0
  { bind_cols(
    select(., mappedID),
    select(., -mappedID)[, colSums(select(., -mappedID)) > 4]
  ) } %>%
  pivot_longer(-mappedID, names_to = "CAZyme", values_to = "Abundance") %>%
  left_join(grouped_cazymes, by = c("CAZyme" = "cazymes"))

# Create abundance matrix (raw values, not log-transformed)
abundance_matrix <- filtered_abundance %>%
  left_join(split_metadata, by = "mappedID") %>%
  arrange(classes, species, prefix) %>%
  group_by(CAZyme, prefix) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = prefix, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("CAZyme") %>%
  as.matrix()

# Column annotation (age_category and classes)
column_annotation <- split_metadata %>%
  filter(prefix %in% colnames(abundance_matrix)) %>%
  select(prefix, age_category, classes) %>%
  distinct() %>%
  mutate(classes = factor(classes, levels = names(classes_order)),
         age_category = factor(age_category, levels = c("I", "Transit", "M"))) %>%
  arrange(age_category, classes, prefix) %>%
  column_to_rownames("prefix")

# Row annotation (substrate group)
row_annotation <- filtered_abundance %>%
  select(CAZyme, substrates) %>%
  distinct() %>%
  filter(CAZyme %in% rownames(abundance_matrix)) %>%
  mutate(substrates = factor(substrates, levels = substrate_order)) %>%
  arrange(substrates) %>%
  column_to_rownames("CAZyme")

# Log-transform for color
log_abundance_matrix <- log10(abundance_matrix + 1)

# Raw values as cell labels, hide zeros
label_matrix <- abundance_matrix
label_matrix[label_matrix == 0] <- "" 
label_matrix <- formatC(label_matrix, format = "f", digits = 0)

# Define log-scaled breaks
breaks <- seq(0, ceiling(max(log_abundance_matrix, na.rm = TRUE)), length.out = 100)
raw_ticks <- round(10^breaks - 1)

log_abundance_matrix <- log_abundance_matrix[
  rownames(row_annotation), 
  rownames(column_annotation)
]

label_matrix <- label_matrix[
  rownames(row_annotation), 
  rownames(column_annotation)
]

# Define color palette 
my_colours <- colorRampPalette(c("snow2", "red"))(length(breaks) - 1)

p <- pheatmap(
  log_abundance_matrix,
  col = my_colours,
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = label_matrix,
  annotation_col = column_annotation,
  annotation_row = row_annotation,
  annotation_names_row = TRUE,
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 90,
  show_colnames = TRUE,
  show_rownames = TRUE,
  annotation_legend = TRUE,
  legend_breaks = breaks[seq(1, length(breaks), length.out = 6)],
  legend_labels = raw_ticks[seq(1, length(breaks), length.out = 6)],
  treeheight_row = 10,
  treeheight_col = 20
)

# Save the plot 
#ggsave("figures/cazyme_heatmap_common_abundant_MAGs.svg", plot = p, width = 10, height = 8)
