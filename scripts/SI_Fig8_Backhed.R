### SI Figure 8 Backhed
### Author: Yunjeong So

# Extended Figure 5. Violin plot------------------------------------------------------
# Load libraries 
library(dplyr);library(DT);library(curatedMetagenomicData);library(data.table);library(tidyverse);library(paletteer);library(ggbreak);library(pheatmap)

# Load data
metadata <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/public-cohort/data/metadata.tsv')
gene_abundance <- fread('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/public-cohort/data/gene-abundance-nonlog-table.txt')
grouped_cazymes <- read_tsv('/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/sequencing/data/cazyme-subgroup-for-box.2025.tsv')
exc_bm_backhed <- sampleMetadata |> 
  filter(study_name == "BackhedF_2015" &
           feeding_practice == "exclusively_breastfeeding" & 
           body_site =="stool" & disease == "healthy" &  
           antibiotics_current_use == "no" &
           born_method != "c_section")

order_substrates <- c(
  "HMOs" = "#252e9a", 
  "HMOs/host glycans" = "#70a6d9",
  "Host glycans" = "#a0cac7",
  "Starch" = "#ffff99", 
  "Inulin" = "darkgreen", 
  "Xylan-backbone" = "#e99215",
  "Pectin-backbone" = "#6f9940",
  "β-glucans" = "#b19cd9", 
  "Food additives" = "#ddc5c1")

# Your desired order of months
month_order <- c("B", "4M")

new_id_levels <- c()  # initialize empty vector

for (substrate in names(order_substrates)) {
  for (month in month_order) {
    new_id_levels <- c(new_id_levels, paste(substrate, month, sep = "_"))
  }
}

# Filter gene abundance data based on CAZyme_ID 
filtered_abundance <- gene_abundance %>%
  filter(ID %in% exc_bm_backhed$NCBI_accession) %>%
  left_join(metadata, by = c("ID")) %>%
  pivot_longer(cols = -c(ID, Month, Sample_name, personalID), names_to = "cazymes", values_to = "abundance") %>%
  left_join(grouped_cazymes, by = "cazymes") %>%
  filter(Figure_1c %in% names(order_substrates)) %>%
  mutate(new_ID = paste(Figure_1c, Month, sep = "_")) %>%
  group_by(new_ID, ID, Figure_1c, Month) %>%
  summarise(total_abundance = log10(sum(abundance, na.rm = TRUE) + 1), .groups = "drop") %>%
  mutate(new_ID = factor(new_ID, levels = new_id_levels))

# Calculate Wilcoxon p-values per group (new_ID)
wilcox_results <- filtered_abundance %>%
  group_by(new_ID) %>%
  summarise(
    p_value = wilcox.test(total_abundance, mu = 0, alternative = "greater")$p.value,
    .groupd = "drop") %>%
  mutate(p_label = ifelse(p_value < 0.01, 
                          formatC(p_value, format = "e", digits = 2),  
                          round(p_value, 4)))

# Plot with p-value labels
violin_p <- ggplot(filtered_abundance, aes(x = new_ID, y = total_abundance, fill = Figure_1c)) +
  geom_violin(alpha = 0.5) +
  stat_summary(fun = median,
               geom = "crossbar",
               colour = "black",
               fatten = 1,
               width = 0.5) +
  geom_point(colour = "black", size = 1, alpha = 0.2, position = position_jitter(width = 0.1)) +
  scale_y_continuous(
    breaks = log10(c(0, 10, 100, 1000, 10000) + 1),
    labels = c("0", "10", "100", "1,000", "10,000")
  ) +
  scale_fill_manual(values = order_substrates) +
  theme_minimal() +
  labs(y = "Genes per million") +
  geom_text(data = wilcox_results,
            aes(x = new_ID, y = log10(10000), label = p_label),
            inherit.aes = FALSE,
            size = 4,
            color = "black")

print(violin_p)
ggsave(plot = violin_p, filename = "figures/backhed_fibres.svg", width = 15, height = 4, bg = "transparent")

# SI Figure 3. Heatmap------------------------------------------------------
# Heatmap 
heatmap_matrix <- filtered_abundance %>%
  select(ID, Figure_1c, total_abundance) %>%
  pivot_wider(names_from = Figure_1c, values_from = total_abundance, values_fill = 0) %>%
  column_to_rownames("ID") %>%
  as.matrix()

row_annotation <- filtered_abundance %>%
  distinct(ID, Month) %>%
  mutate(Month = factor(Month, levels = c("B", "4M", "12M"))) %>% 
  column_to_rownames("ID") %>% 
  arrange(Month)

heatmap_matrix <- heatmap_matrix[rownames(row_annotation), names(order_substrates)]

# Define desired breakpoints and labels
my_colours <- colorRampPalette(c("snow2",  "thistle3", "brown"))(100)
breaks <- seq(0, max(heatmap_matrix, na.rm = TRUE), length.out = 100)
label_matrix <- ifelse(heatmap_matrix == 0, "0", "")
legend_values <- c(0, 10, 100, 1000, 6000)
legend_breaks <- log10(legend_values + 1)
legend_labels <- as.character(legend_values)

# Plot
p <- pheatmap(
  heatmap_matrix, 
  col = my_colours,
  border_color = FALSE, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE, 
  annotation_row = row_annotation, 
  fontsize_row = 7,         
  fontsize_col = 7,          
  angle_col = 45, 
  gaps_row = c(53),
  gaps_col = c(3, 7),
  display_numbers = label_matrix, 
  number_color = "grey50",
  legend_breaks = legend_breaks,
  legend_labels = legend_labels
)

print(p)
#ggsave("figures/backhed-heatmap.svg", p, width = 5, height = 6, bg = "transparent")
