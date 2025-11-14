### SI Figure 5 First enrichment
### Author: Yunjeong So

# Figure6 a. OD bar plot------------------------------------------------------

# load  packages
library(ggplot2);library(dplyr);library(tidyr);library(data.table);library(glmmTMB)

data <- fread("/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/first-enrichment/data/all-OD-split-ID.txt")

# Set factor levels
data$Condition <- factor(data$Condition, levels = c("NC", "Monos", "HMOs", "Mucin", "PX", "XG", "GG"))

# Pivot longer
data_long <- data %>%
  pivot_longer(
    cols = c("Pre-weaning", "Early-weaning", "Late-weaning", "Mother"),
    names_to = "timepoints",
    values_to = "OD"
  )

# Mixed models by timepoint
model_pre <- glmmTMB(`Pre-weaning` ~ Condition + (1|DonorID) + (1|Replicate), data = data, family = "gaussian")
model_early <- glmmTMB(`Early-weaning` ~ Condition + (1|DonorID) + (1|Replicate), data = data, family = "gaussian")
model_late <- glmmTMB(`Late-weaning` ~ Condition + (1|DonorID) + (1|Replicate), data = data, family = "gaussian")
model_mother <- glmmTMB(Mother ~ Condition + (1|DonorID) + (1|Replicate), data = data, family = "gaussian")

summary(model_pre)
summary(model_early)
summary(model_late)
summary(model_mother)

# Extract and format coefficients with 3 significant digits
coef_pre <- as.data.frame(summary(model_pre)$coefficients$cond)
coef_pre[] <- lapply(coef_pre, function(x) if (is.numeric(x)) signif(x, 3) else x)

coef_early <- as.data.frame(summary(model_early)$coefficients$cond)
coef_early[] <- lapply(coef_early, function(x) if (is.numeric(x)) signif(x, 3) else x)

coef_late <- as.data.frame(summary(model_late)$coefficients$cond)
coef_late[] <- lapply(coef_late, function(x) if (is.numeric(x)) signif(x, 3) else x)

coef_mother <- as.data.frame(summary(model_mother)$coefficients$cond)
coef_mother[] <- lapply(coef_mother, function(x) if (is.numeric(x)) signif(x, 3) else x)

# View the results
coef_pre
coef_early
coef_late
coef_mother


my_colours <- c("#5a86c8", "#a77aa7","#8bb956", "#e8b462" )

p <- ggplot(data_long, aes(Condition, OD, color = timepoints)) +
  geom_boxplot(position = position_dodge2(), colour = "grey30",  
               aes(fill = as.factor(timepoints)),  alpha = 0.1, , outlier.shape = NA) +
  geom_point(position = position_dodge2(width = .8), 
             size = 1, 
             aes(colour = as.factor(timepoints))) +
  ylab("OD") +
  scale_color_manual(values = my_colours) +
  scale_fill_manual(values = my_colours) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black")
  ) 
print(p)

# save the plot
#ggsave(paste0("/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/scripts/first-enrichment/figures/all_pH.svg"), width = 24, height = 10, dpi = 600, units = "cm", plot = p)

#+facet_wrap(. ~ variable, scales = "free_x", row_number(1)) 

# Figure6 b. pH bar plot------------------------------------------------------

# Load data
data <- fread("/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/first-enrichment/data/2025.all-pH.txt")

# Set factor levels
data$Condition <- factor(data$Condition, levels = c("NC", "Monos", "HMOs", "Mucin", "PX", "XG", "GG"))

# Pivot longer
data_long <- data %>%
  pivot_longer(
    cols = c("Pre-weaning", "Early-weaning", "Late-weaning", "Mother"),
    names_to = "timepoints",
    values_to = "pH"
  )

# Mixed models by timepoint
model_pre <- glmmTMB(`Pre-weaning` ~ Condition + (1|DonorID) + (1|Replicate), data = data, family = "gaussian")
model_early <- glmmTMB(`Early-weaning` ~ Condition + (1|DonorID) + (1|Replicate), data = data, family = "gaussian")
model_late <- glmmTMB(`Late-weaning` ~ Condition + (1|DonorID) + (1|Replicate), data = data, family = "gaussian")
model_mother <- glmmTMB(Mother ~ Condition + (1|DonorID) + (1|Replicate), data = data, family = "gaussian")

summary(model_pre)
summary(model_early)
summary(model_late)
summary(model_mother)

# Extract and format coefficients with 3 significant digits
coef_pre <- as.data.frame(summary(model_pre)$coefficients$cond)
coef_pre[] <- lapply(coef_pre, function(x) if (is.numeric(x)) signif(x, 3) else x)

coef_early <- as.data.frame(summary(model_early)$coefficients$cond)
coef_early[] <- lapply(coef_early, function(x) if (is.numeric(x)) signif(x, 3) else x)

coef_late <- as.data.frame(summary(model_late)$coefficients$cond)
coef_late[] <- lapply(coef_late, function(x) if (is.numeric(x)) signif(x, 3) else x)

coef_mother <- as.data.frame(summary(model_mother)$coefficients$cond)
coef_mother[] <- lapply(coef_mother, function(x) if (is.numeric(x)) signif(x, 3) else x)

# View the results
coef_pre
coef_early
coef_late
coef_mother


my_colours <- c("#5a86c8", "#a77aa7","#8bb956", "#e8b462" )

p <- ggplot(data_long, aes(Condition, pH, color = timepoints)) +
  geom_boxplot(position = position_dodge2(), colour = "grey30",  
               aes(fill = as.factor(timepoints)),  alpha = 0.1, , outlier.shape = NA) +
  geom_point(position = position_dodge2(width = .8), 
             size = 1, 
             aes(colour = as.factor(timepoints))) +
  ylab("pH") +
  scale_color_manual(values = my_colours) +
  scale_fill_manual(values = my_colours) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black")
  ) 
print(p)

# save the plot
#ggsave(paste0("/Users/yunjeongso/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD/milkome-cohort/scripts/first-enrichment/figures/all_pH.svg"), width = 24, height = 10, dpi = 600, units = "cm", plot = p)

#+facet_wrap(. ~ variable, scales = "free_x", row_number(1)) 

