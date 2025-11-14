### Figure SI Fig 9
### Author: Yunjeong So

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
  "xylan/cellulose/b-glucans" = "#8b5d2f",
  "GH74" = "#1C0CA1",
  "GH5_4" = "#5599cc",
  "GH12" = "#6FE8FD")
PEC_colors <- c(
  "GH28" = "#6f8520",
  "GH105" = "#989279",
  "PL" = "#006400",
  "Pectin-side-chains" = "#E7E178")
GM_colors <- c(
  "GH5_7" = "#363700",
  "GH5_8" = "#4A8A4C",
  "GH5_41" = "#E5EAE7",
  "GH26" = "#9AC3AD",
  "GH113" = "#9C6769",
  "GH130_1" = "#2B4C40",
  "GH130_2" = "#D8C978",
  "GH130_3" = "#515B50",
  "GH130_4" = "#FDD4D3",
  "GH130_5" = "#1F2020")
GA_colors <- c(
  "GH39" = "#D5D254",
  "GH43_24" = "#A46F53",
  "GH154" = "#7AA1C3",
  "PL27" = "#6B7566",
  "PL42" = "#4656A2")

# Function 
extract_substrate <- function(target) {
  target_list <- switch(target,
                        "AX"  = c("F", "AX", "PEC", "XG"),
                        "PEC" = c("F", "AX", "PEC", "XG"),
                        "GM"  = c("F", "GM", "GA"),
                        "GA"  = c("F", "GM", "GA"),
                        stop("target must be one of: 'GM', or 'GA'"))
  
  target_colors <- switch(target,
                          "AX"  = AX_colors,
                          "PEC" = PEC_colors,
                          "GM"  = GM_colors,
                          "GA"  = GA_colors,)
  
  filtered_metadata <- metadata %>%
    left_join(mapping, by = "ID") %>% 
    filter(infant == "I", substrate %in% target_list) %>% 
    mutate(
      participant = as.factor(participant),
      substrate = factor(substrate, levels = target_list),
      prefix = paste(timepoint, substrate, participant, sep = "_"), 
      newID = paste(participant, timepoint, sep = "_")
    ) %>%
    arrange(timepoint, substrate, participant) %>%  
    mutate(prefix = factor(prefix, levels = unique(paste(timepoint, substrate, participant,sep = "_"))))
  
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

  abundance_long
}


sub_list <- c("AX", "PEC", "GM", "GA")

for (sub in sub_list) {
  target_colors <- get(paste0(sub, "_colors"))
  
  merged_data <- bind_rows(extract_substrate(sub)) %>% 
    mutate(Figure_5 = factor(Figure_5, levels = names(target_colors)))
  
  p <- ggplot(merged_data, aes(x = prefix, y = abundance, fill = Figure_5)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.85, color = NA) +
    geom_hline(yintercept = 0, color = "black", size = 0.25) +
    scale_fill_manual(values = target_colors, drop = FALSE) +
    labs(
      x = "Sample (newID)",
      fill = "CAZyme group",
      title = paste("Substrate:", sub)
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
  #ggsave(paste0("figures/with-faeses_", sub, "_gene-enrichment.svg"), p, width = 10, height = 3, bg = "transparent")
}


