### Figure 2

# Fig2 a. Phylum area ------------------------------------------------------
# Load required libraries
library(tidyverse)
library(ggplot2)

# Load data
metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered.tsv')
taxa_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/meta-phyla-converted.txt')

# Reshape data
taxa_abundance_long <- taxa_abundance %>%
  pivot_longer(cols = -ID, names_to = "taxa", values_to = "taxa_abundance")

# Split taxa
taxa_abundance_long <- taxa_abundance_long %>%
  separate(taxa, into = c("phylum", "species"), sep = ";") %>%
  mutate(
    phylum = str_replace(phylum, "^p__", ""),  
    species = str_replace(species, "^s__", "") 
  ) 

# Merge the metadata and reshaped taxa abundance data
merged_data <- merge(metadata, taxa_abundance_long, by = "ID")

# Convert timepoint and substrate to factors for ordered plotting
merged_data$timepoint <- factor(merged_data$timepoint, levels = c('T1', 'T2', 'T3'))
timepoint_mapping <- c("T1" = 1, "T2" = 2, "T3" = 3)

# Function non target phyla into "Others"
group_others <- function(phylum) {
  if (!(phylum %in% c("Actinobacteriota", "Firmicutes", "Bacteroidota",  "Proteobacteria", "Others",  "Unclassified"))) {
    return("Others")
  }
  return(phylum)
}

# Filter data
filtered_data <- merged_data %>%
  filter(infant == "I", substrate %in% c("F")) %>%
  mutate(phylum = sapply(phylum, group_others)) %>% 
  group_by(participant, phylum, timepoint) %>%      
  summarise(taxa_abundance = sum(taxa_abundance, na.rm = TRUE)) %>% 
  select(participant, timepoint, phylum, taxa_abundance) 

# Phyla level and colours
p_levels <- c("Actinobacteriota", "Firmicutes", "Bacteroidota",  "Proteobacteria", "Others",  "Unclassified")
p_colours <- c(
  "Actinobacteriota" = "#171b60", 
  "Bacteroidota" = "#065e06",  
  "Firmicutes" = "#a54846", 
  "Proteobacteria" = "wheat3", 
  "Verrucomicrobiota" = "#5b264b",
  "Others" = "#8a9eb4",
  "Unclassified" = "grey80"
)

# Reassign the factor levels for phylum
filtered_data$phylum <- factor(filtered_data$phylum, levels = p_levels)
filtered_data$timepoint_numeric <- timepoint_mapping[filtered_data$timepoint]

# Continue with the plotting as before
p <- ggplot(filtered_data, aes(x = timepoint_numeric, y = taxa_abundance, fill = phylum)) +
  geom_area(alpha=0.85 , linewidth=0.0, colour=NA) +
  scale_fill_manual(values = p_colours) +  
  labs( y = "Relative Abundance") +  
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = "black"),  
    panel.grid = element_blank(), 
    plot.title = element_text(hjust = 0.5),  
    axis.text.y = element_text(),  
    #legend.position = "none"
  )  +
  facet_wrap(~participant, nrow=1, scales ="free_y")

# Print and save plot
print(p)
filename <- paste0("figures/area-Infant-faeces_others.svg")
#ggsave(plot = p, filename = filename, width = 8, height = 5, bg = "transparent")


# Fig2 b. Bifidobacterium species heatmap ----------------------------------
# Load necessary libraries
library(RColorBrewer)
library(tidyverse)
library(paletteer)
library(pheatmap)
library(grDevices)

# Load data
metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered.tsv')
taxa_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/meta-wholetaxa-converted.txt')

# Reshape data
taxa_abundance_long <- taxa_abundance %>%
  pivot_longer(cols = -ID, names_to = "taxa", values_to = "taxa_abundance") %>%
  separate(taxa, into = c("kingdom", "phyla", "classes", "orders", "family", "genus", "species"), sep = ";", fill = "right") %>%
  mutate(
    kingdom = str_replace(kingdom, "^d__", ""),
    phyla = str_replace(phyla, "^p__", ""),  
    classes = str_replace(classes, "^c__", ""),  
    orders = str_replace(orders, "^o__", ""),  
    family = str_replace(family, "^f__", ""),
    genus = str_replace(genus, "^g__", ""),
    species = str_replace(species, "^s__", "")) %>%
  replace_na(list(genus = "unassigned_genus", species = "unassigned_species")) %>% 
  filter(taxa_abundance >= 0.01)

# Merge metadata and taxa abundance data
merged_data <- merge(metadata, taxa_abundance_long, by = "ID")
merged_data$timepoint <- factor(merged_data$timepoint, levels = c('T1', 'T2', 'T3'))

# Define a function
class_func <- function(classes, species) {
  if (!(classes %in% c("Actinomycetia", "Bacilli", "Bacteroidia", "Clostridia", "Gammaproteobacteria", "Negativicutes", "Verrucomicrobiae", "Unclassified"))) 
  { return("Others") 
  } else if (classes == "Actinomycetia" & startsWith(as.character(species), "Bifidobacterium" ) & !(species %in% c("Bifidobacterium animalis", "Bifidobacterium angulatum")) ) 
  { return(species)
  } else 
  { return(classes)}}

phylum_func <- function(phyla, fill_taxa) {
  if (fill_taxa%in% c("Others", "Unclassified")) {
    return("Others") 
  } else {
    return(phyla)}}

# Define the order of taxa
taxa_levels <- c("Bifidobacterium longum.infantis", "Bifidobacterium longum.longum", "Bifidobacterium bifidum", "Bifidobacterium breve",  
                 "Bifidobacterium adolescentis", "Bifidobacterium angulatum", "Bifidobacterium animalis", "Bifidobacterium catenulatum", "Bifidobacterium dentium", "Bifidobacterium pseudocatenulatum",
                 "Actinomycetia", "Bacilli", "Clostridia", "Negativicutes", "Bacteroidia", "Gammaproteobacteria", "Verrucomicrobiae", "Others", "Unclassified")

# format the abundance data
bif_phyla <- merged_data %>%
  filter(infant == "I", substrate %in% c("F")) %>% 
  mutate(
    fill_taxa = mapply(class_func, classes, species),
    newID = paste(participant, timepoint, sep = "_")) %>%
  mutate(anno_phyla = mapply(phylum_func, phyla, fill_taxa)) %>%
  select(newID, fill_taxa, anno_phyla, taxa_abundance) 

bif_phyla_log <- bif_phyla %>% 
  select(newID, fill_taxa, anno_phyla, taxa_abundance) %>%
  group_by(newID, fill_taxa) %>%
  summarise(fill_abundance = sum(taxa_abundance, na.rm = TRUE), .groups = 'drop')  %>% 
  mutate(fill_abundance = ((fill_abundance))) %>% 
  pivot_wider(names_from = fill_taxa, values_from = fill_abundance) %>%
  column_to_rownames(var = "newID") 

# Create a phyla annotation
phyla_anno <- bif_phyla %>%
  select(fill_taxa, anno_phyla) %>%               
  distinct(fill_taxa, .keep_all = TRUE) %>%  
  mutate(fill_taxa = factor(fill_taxa, levels = taxa_levels)) %>% 
  arrange(fill_taxa) %>% 
  column_to_rownames(var = "fill_taxa") 

# Set colours
my_colours <- colorRampPalette(c("snow2", "darkmagenta"))(100)
p_colours <- c(
  "Actinobacteriota" = "#171b60", 
  "Bacteroidota" = "#065e06",  
  "Firmicutes" = "#a54846", 
  "Proteobacteria" = "wheat3", 
  "Verrucomicrobiota" = "#4e193e",
  "Others" = "grey80")

# Generate heatmap with the correct phyla annotation
matrix <- as.matrix(t(bif_phyla_log))
species_names <- rownames(matrix)
matrix <- matrix[match(rownames(phyla_anno), rownames(matrix)), ]
matrix[is.na(matrix)] <- 0

# Define breaks and labels manually
breaks <- seq(0, max(matrix, na.rm = TRUE), length.out = 100)
label_matrix <- ifelse(matrix == 0, "0", "")

# Plot the heatmap
p <- pheatmap(matrix, 
              border_color = FALSE, 
              col = my_colours,
              breaks = breaks,
              cluster_rows = F, 
              cluster_cols = F, 
              annotation_row = phyla_anno, 
              annotation_colors = list(anno_phyla = p_colours), 
              gaps_col = c(3,6,9,11,14,15),
              gaps_row = c(9,12,13,14,15),
              fontsize_row = 6,         
              fontsize_col = 7,          
              angle_col = 45, 
              show_colnames = T, 
              show_rownames = T,
              annotation_legend = T,
              treeheight_row = 10,  
              treeheight_col = 20,
              display_numbers = label_matrix, 
              number_color = "grey50")

print(p)
filename <- paste0("figures/heatmap_bif_class_by_subject.svg")
#ggsave(plot = p, filename = filename, width = 5, height = 2.5, bg = "transparent")


# Fig2 c, d. Violin plot - Counts of CAZymes over weaning ----------------------------------
# Load necessary libraries
library(data.table);library(tidyverse);library(ggpubr)

# Read data (this remains unchanged)
metadata <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered.tsv')
gene_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/gene-abundance-nonlog-table.txt')
mapping <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sample_mapping_all.txt')
grouped_cazymes <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/cazyme-subgroup-for-box.2025.tsv')


# Merge metadata
md_infant_F <- metadata %>%
  mutate(prefix = paste(infant, participant, timepoint, sep = "_")) %>%
  left_join(mapping, by = "ID") %>%
  filter(substrate == "F" & infant == "I")


order_substrates <- c(
  "HMOs", "HMOs/host glycans", "Host glycans",
  "Starch", "Xylan-backbone", "Pectin-backbone", "β-glucans", "Food additives")


# Filter gene abundance data based on CAZyme_ID 
filtered_abundance <- gene_abundance %>%
  left_join(md_infant_F, by = "mappedID") %>% 
  filter(!is.na(prefix)) %>%
  select(mappedID, prefix, timepoint, participant, all_of(grouped_cazymes$cazymes)) %>% 
  pivot_longer(cols = -c(mappedID, prefix, timepoint, participant), names_to = "cazymes", values_to = "abundance") %>% 
  left_join(grouped_cazymes, by = "cazymes") %>% 
  filter(Figure_2c %in% order_substrates) %>% 
  group_by(prefix, timepoint, participant, mappedID, Figure_2c)  %>% 
  summarise(total_abundance = log10(sum(abundance , na.rm = TRUE) + 1), .groups = "drop") %>% 
  mutate(Figure_2c = factor(Figure_2c, levels = order_substrates))


participant_colours <- c("A" = "#4c89cd",
                         "B" = "#e99215",
                         "C" = "#478c2c",
                         "D" = "#d64a28",
                         "E" = "#713d91",
                         "F" = "#8c564b",
                         "G" = "#ce9ab9")


# Cazymes bar plots with y-axis break
p <- ggplot(filtered_abundance, aes(x = timepoint, y = total_abundance)) +
  geom_violin(fill = "#998675", alpha = 0.1) +
  stat_summary(fun = median,
               geom = "crossbar",
               colour = "black",
               fatten = 1) +  
  geom_point(aes(colour = participant), 
             position = position_dodge2(width = 0.8), size = 1.5, alpha = 0.8) +
  scale_colour_manual(values = participant_colours) +
  scale_y_continuous(
    breaks = log10(c(1, 6, 11, 101, 501, 1001, 5001)),
    labels = c("0", "5", "10", "100", 500, "1000", "5000")) +
  labs(y = "GPM") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  facet_wrap(~ Figure_2c, nrow = 1) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3")),
                     label.y = c(3.5, 3.6, 3.7),
                     format.args = list(digits = 5))

print(p)


# Save the plot# Save theposition_jitterdodge() plot
filename <- paste0("figures/infant-faeces-cazympes-plants-violin-plot.svg")
#ggsave(plot = p, filename = filename, width = 5, height = 2.5, bg = "transparent")

# Fig2 e. Heatmap of CAZymes ----------------------------------
# Load necessary libraries
library(data.table);library(tidyverse);library(paletteer);library(pheatmap);library(RColorBrewer)

# Read data
metadata <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered.tsv')
gene_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/gene-abundance-nonlog-table.txt')
mapping <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sample_mapping_all.txt')
grouped_cazymes <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/cazyme-subgroup-for-box.2025.tsv')

# Merge metadata
md_infant_F <- metadata %>%
  mutate(prefix = paste(infant, participant, timepoint, sep = "_")) %>%
  left_join(mapping, by = "ID") %>%
  filter(substrate == "F" & infant == "I")

# Filter gene abundance data based on CAZyme_ID 
log_abundance <- gene_abundance %>%
  filter(mappedID %in% md_infant_F$mappedID) %>%
  column_to_rownames(var = "mappedID") %>% 
  select(all_of(grouped_cazymes$cazymes)) %>% 
  select(which(colSums(.) > 10)) %>% 
  mutate(across(everything(), ~ log10(as.numeric(.) + 1))) 

# Convert to matrix and transpose
abundance_matrix <- as.matrix(t(log_abundance))

# Set group colours
cazy_colours <- c("HMOs/host glycans" = "yellow",
                  "Host glycans" = "skyblue3",
                  "HMOs/host glycans/plant" = "red",
                  "Starch" = "#ebc4c1", 
                  "Dietary fibres" = "darkgreen",
                  "Food additives" = "grey")

# Order of substrates
order_substrates <- c(
  "HMOs", "GH20", "HMOs/host glycans", "HMOs/plant", 
  "O-glycans", "Glycosaminoglycans", "Starch",
  "Inulin", "Xylan-backbone", 
  "Pectin-backbone", "Pectin-backbone-PL",
  "β-glucans", "Mannan", "Gum Arabic")

# Set substrate colours
substrate_colours <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(length(order_substrates)),
  order_substrates)

# Create a dataframe for row annotation
cazy_anno <- grouped_cazymes %>%
  as.data.frame() %>% 
  filter(cazymes %in% colnames(log_abundance)) %>% 
  column_to_rownames(var = "cazymes") %>% 
  filter(Figure_1e %in% order_substrates) %>% 
  mutate(groups = factor(groups, levels = names(cazy_colours)),
         Figure_1e = factor(Figure_1e, levels = order_substrates)) %>%
  arrange(groups, Figure_1e) %>%
  select(groups, Figure_1e)

# Reorder abundance_matrix rows based on cazy_anno
abundance_matrix <- abundance_matrix[match(rownames(cazy_anno), rownames(abundance_matrix)), ]
abundance_matrix[is.na(abundance_matrix)] <- 0

# Color palette and legends
my_colours <- colorRampPalette(c("snow2", "#b22222"))(100)
breaks <- seq(0, max(abundance_matrix, na.rm = TRUE), length.out = 100)
label_matrix <- ifelse(abundance_matrix == 0, "0", "")
legend_values <- c(0, 1, 10, 100, 500, 1000)
legend_breaks <- log10(legend_values + 1)
legend_labels <- as.character(legend_values)

# Annotation color mapping
annotation_colors <- list(
  groups = cazy_colours,
  substrates = substrate_colours
)

# Generate the heatmap
p <- pheatmap(
  abundance_matrix, 
  col = my_colours,
  border_color = FALSE, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE, 
  gaps_col = c(3, 6, 9, 11, 14, 15),
  gaps_row = c(2, 5, 8, 13, 14, 16),
  annotation_row = cazy_anno, 
  annotation_colors = annotation_colors, 
  annotation_names_row = TRUE, 
  fontsize_row = 7,         
  fontsize_col = 7,          
  angle_col = 45, 
  show_colnames = TRUE, 
  show_rownames = TRUE,
  annotation_legend = TRUE,
  display_numbers = label_matrix, 
  number_color = "grey50",
  legend_breaks = legend_breaks,
  legend_labels = legend_labels
)

# Save the heatmap
#ggsave("figures/selected_cazymes-heatmap.svg", plot = p, width = 10, height = 7.8)


# Fig2 f. Shannon diversity of microbiota and CAZymes ----------------------------------
# Load libraries
library(vegan); library(tidyverse); library(rstatix); library(dplyr)

# Read input files
metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered.tsv')
gene_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/gene-abundance-nonlog-table.txt')
mapping <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sample_mapping_all.txt')
taxa_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/meta-wholetaxa-converted.txt')
grouped_cazymes <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/cazyme-subgroup-for-box.2025.tsv')

# Define colors
participant_colours <- c("A" = "#4c89cd", "B" = "#e99215", "C" = "#478c2c",
                         "D" = "#d64a28", "E" = "#713d91", "F" = "#8c564b", "G" = "#ce9ab9")

# Prepare metadata
md_infant_F <- metadata %>%
  mutate(prefix = paste(infant, participant, timepoint, sep = "_")) %>%
  left_join(mapping, by = "ID") %>%
  filter(substrate == "F" & infant == "I")

all_ids <- md_infant_F %>% select(ID) %>% distinct()

# taxa abundance
taxa_abundance_long <- taxa_abundance %>%
  pivot_longer(cols = -ID, names_to = "taxa", values_to = "taxa_abundance") %>%
  separate(taxa, into = c("kingdom", "phyla", "classes", "orders", "family", "genus", "species"), sep = ";", fill = "right") %>%
  mutate(across(c(kingdom, phyla, classes, orders, family, genus, species),
                ~ str_replace(., "^[a-z]__", ""))) %>%
  mutate(
    genus = ifelse(is.na(genus) | genus == "", "unassigned_genus", genus),
    species = ifelse(is.na(species) | species == "", "unassigned_species", species)
  ) %>%
  filter(taxa_abundance >= 0.01) %>%
  filter(ID %in% md_infant_F$ID)

filtered_taxa <- taxa_abundance_long %>%
  select(ID, species, taxa_abundance) %>%
  filter(species != "Unclassified") %>%
  group_by(ID, species) %>%
  summarise(total_abundance = sum(taxa_abundance, na.rm = TRUE), .groups = "drop")

shannon_taxa <- filtered_taxa %>%
  group_by(ID) %>%
  summarise(shannon = diversity(total_abundance, index = "shannon"), .groups = "drop") %>%
  right_join(all_ids, by = "ID") %>%
  mutate(shannon = replace_na(shannon, 0),
         source = "taxa abundance")

# PLANT gene abundance (Starch + Dietary fibres)
plant_cazymes <- grouped_cazymes %>%
  filter(groups %in% c("Starch", "Dietary fibres"))

filtered_abundance_plant <- gene_abundance %>%
  left_join(md_infant_F, by = "mappedID") %>%
  filter(!is.na(prefix)) %>%
  select(ID, prefix, timepoint, participant, all_of(plant_cazymes$cazymes)) %>%
  pivot_longer(cols = -c(ID, prefix, timepoint, participant), 
               names_to = "cazymes", values_to = "abundance") %>%
  left_join(grouped_cazymes, by = "cazymes")

shannon_gene <- filtered_abundance_plant %>%
  group_by(ID) %>%
  summarise(shannon = diversity(abundance, index = "shannon"), .groups = "drop") %>%
  right_join(all_ids, by = "ID") %>%
  mutate(shannon = replace_na(shannon, 0),
         source = "Plant gene abundance")

# FOOD ADDITIVES gene abundance
additives_cazymes <- grouped_cazymes %>%
  filter(groups %in% c("Food additives"))

filtered_abundance_add <- gene_abundance %>%
  left_join(md_infant_F, by = "mappedID") %>%
  filter(!is.na(prefix)) %>%
  select(ID, prefix, timepoint, participant, all_of(additives_cazymes$cazymes)) %>%
  pivot_longer(cols = -c(ID, prefix, timepoint, participant), 
               names_to = "cazymes", values_to = "abundance") %>%
  left_join(grouped_cazymes, by = "cazymes")

shannon_additives <- filtered_abundance_add %>%
  group_by(ID) %>%
  summarise(shannon = diversity(abundance, index = "shannon"), .groups = "drop") %>%
  right_join(all_ids, by = "ID") %>%
  mutate(shannon = replace_na(shannon, 0),
         source = "Food additives abundance")

# COMBINE and plot
combined_shannon <- bind_rows(shannon_gene, shannon_additives, shannon_taxa)

# Add metadata
combined_shannon <- combined_shannon %>%
  left_join(md_infant_F %>% select(ID, timepoint, participant), by = "ID")

# Plot
p <- ggplot(combined_shannon, aes(x = timepoint, y = shannon)) +
  geom_violin(alpha = 0.1, lwd = 0.5, fill = "#998675") +
  stat_summary(fun = median,
               geom = "crossbar",
               colour = "black",
               fatten = 1) +  
  geom_point(aes(colour = participant),
             position = position_jitter(width = 0.25),
             size = 1.5, alpha = 0.8) +
  scale_colour_manual(values = participant_colours) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = "black"),
    plot.background = element_rect(fill = "transparent", colour = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(),
    axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(y = "Shannon Index") +
  facet_wrap(~ source, nrow = 1) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3")),
    label = "p.format",
    label.y = c(3.47, 3.59, 3.71))

print(p)

# save
#ggsave("figures/shannon_diversity_taxa_gene-with-additives.svg", plot = p, width = 6, height = 3)


# Fig2 g. NMDS with CAZymes ----------------------------------
# Load libraries
library(tidyverse);library(vegan);library(data.table);library(ggrepel)

# Read data
metadata <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sequencing-metadata-ordered.tsv')
taxa_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/meta-wholetaxa-converted.txt')
mapping <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/sample_mapping_all.txt')
gene_abundance <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/gene-abundance-nonlog-table.txt')
grouped_cazymes <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/cazyme-subgroup-for-box.2025.tsv')

# Merge and filter metadata
md_infant_F <- metadata %>%
  mutate(prefix = paste(infant, participant, sep = "_")) %>%
  left_join(mapping, by = "ID") %>%
  filter(substrate == "F" & infant == "I")

# Make data long
taxa_abundance_long <- taxa_abundance %>%
  pivot_longer(cols = -ID, names_to = "taxa", values_to = "taxa_abundance") %>%
  separate(taxa, into = c("kingdom", "phyla", "classes", "orders", "family", "genus", "species"), sep = ";", fill = "right") %>%
  mutate(across(c(kingdom, phyla, classes, orders, family, genus, species),
                ~ str_replace(., "^[a-z]__", ""))) %>%
  mutate(
    genus = ifelse(is.na(genus) | genus == "", "unassigned_genus", genus),
    species = ifelse(is.na(species) | species == "", "unassigned_species", species)) %>%
  filter(ID %in% md_infant_F$ID) %>%
  filter(taxa_abundance >= 0.01) %>%
  select(ID, species, taxa_abundance) %>% 
  filter(species != "Unclassified")  

# Convert to matrix
df_abundance <- taxa_abundance_long %>%
  pivot_wider(names_from = species, values_from = taxa_abundance, values_fill = 0) %>%
  column_to_rownames("ID")

# NMDS analysis
abundance_matrix <- as.matrix(df_abundance)
set.seed(123)
nmds <- metaMDS(abundance_matrix, distance = "bray")

# Extract NMDS scores
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores2 <- merge(data.scores, md_infant_F, by.x =  0, by.y = 'ID')

grouped_cazymes <- grouped_cazymes %>%
  filter(substrates %in% c("HMOs", "HMOs/host glycans", "HMOs/plant", "Host glycans",
                           "Starch", "Inulin", "β-glucans", "β-glucans, Xylan-backbone, Food additives",       
                           "Xylan-backbone", "Pectin-backbone", "Food additives"))

# Filter gene abundance data
gene_filtered_abundance <- gene_abundance %>%
  filter(mappedID %in% md_infant_F$mappedID) %>%
  left_join(mapping, by = "mappedID") %>%
  select(ID, any_of(grouped_cazymes$cazymes)) %>%
  mutate(across(-ID, as.numeric)) %>%
  column_to_rownames("ID") %>%
  select(where(~ sum(.) > 0)) 

# Perform envfit for species-environment relationships
en <- envfit(nmds, gene_filtered_abundance, permutations = 999, na.rm = TRUE)
en_coord_cont <- as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cont$pvals <- en$vectors$pvals

# Adjust p-value for multiple testing
en_coord_cont_plot <- en_coord_cont %>%
  filter(pvals <= 0.05) %>%
  rownames_to_column("cazymes") %>%
  left_join(grouped_cazymes, by = "cazymes")

# Top cazymes per substrate by p-values
selected_cazymes <- c(
  # HMOs/host-glycans
  "GH181", "GH29",
  # Host glycans
  "PL35", "PL33_1",
  # Starch
  "GH13_46",
  # Xylan
  "GH5_21", "GH10",
  # Pectin
  "PL1_2", "PL11",
  # b-glucans
  "GH5_2", "GH5_5",
  # Food additives
  "GH26", "GH130_1")

# Subset the en_coord_cont_plot to include only the selected CAZymes
en_coord_cont_plot <- en_coord_cont_plot %>%
  filter(cazymes %in% selected_cazymes)

# Set colours
participant_colours <- c("A" = "#4c89cd", "B" = "#e99215", "C" = "#478c2c", "D" = "#d64a28",
                         "E" = "#713d91", "F" = "#8c564b", "G" = "#ce9ab9")
timepoint_colours <- c(  "T1" = "#a6cee3",
                         "T2" = "pink",
                         "T3" = "#b2df8a")

# Plot with ggplot2
p <- ggplot(data = data.scores2, aes(x = NMDS1, y = NMDS2)) + 
  stat_ellipse(aes(group = timepoint), 
               type = "t", level = 0.68, alpha = 0.2) +
  geom_point(aes(colour = as.factor(participant), shape = timepoint), 
             size = 4, alpha = 0.8, stroke = 0.5) +
  scale_colour_manual(values = participant_colours) +
  geom_segment(data = en_coord_cont_plot, 
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, colour = substrates), 
               linewidth = 0.5, alpha = 0.8, show.legend = FALSE) +
  geom_text_repel(data = en_coord_cont_plot, 
                  aes(x = NMDS1, y = NMDS2, colour = substrates, label = cazymes),
                  size = 3, show.legend = FALSE) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "transparent", colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

print(p)

# PERMANOVA
set.seed(2010)
model <- adonis2(abundance_matrix ~ timepoint, data = md_infant_F, method = "bray", permutations = 9999, by = "terms")
print(model)

# Save the plot
#ggsave("figures/NMDS_taxas-with-cazylabel.svg", plot = p, width = 8, height = 6)

