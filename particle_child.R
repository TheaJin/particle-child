##############################################
# Particle Deposition Simulation Analysis
# Author: Theodora Jin
#
# Description:
#   This script processes particle deposition data from simulation experiments
#   on children (virtual and constricted conditions) and clustered data. It:
#     • Reads and processes virtual child data (with various constriction levels)
#       and cluster data.
#     • Performs paired comparisons (e.g., paired t-test of intra-thoracic deposition
#       between VC (1.0) and Cluster_12).
#     • Generates a series of ggplot2 visualizations for different deposition measures:
#         - Intra-thoracic (TDF)
#         - Bronchial (BDE)
#         - Alveolar (ADF)
#     • Merges adult simulation data (from TDF, BDF, ADF sheets) with virtual child data
#       for multi-measure comparisons.
#     • Runs pairwise comparisons using emmeans and computes compact letter displays (CLD).
#     • Optionally, reads external "lobe summary" files, summarizes by condition and lobe,
#       and plots deposition by lobe.
#
# Data Requirements:
#   - Excel files containing virtual child deposition data (e.g., "particle_deposition_results_virtual_old.xlsx",
#     "particle_deposition_results_constrict_90.xlsx", "particle_deposition_results_constrict_85.xlsx",
#     "particle_deposition_results_constrict_80.xlsx").
#   - An Excel file for cluster deposition data (e.g., "particle_deposition_results_cluster_12.xlsx").
#   - An adult simulation Excel file ("result_new.xlsx") with sheets "TDF_Lung", "BDF", "ADF".
#   - A demographic dataset ("demo") with subject IDs, gender, etc.
##############################################

# ------------------------------
# 1. Load Required Libraries
# ------------------------------
library(readxl)
library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(reshape2)
library(tidyr)
library(purrr)
library(broom)
library(emmeans)
library(multcomp)
library(pheatmap)
library(scales)

# ------------------------------
# 2. Import and Process Virtual Child and Cluster Data
# ------------------------------
# Virtual child data (VC Ref)
virtual_child_200 <- read_xlsx("/hpc/gjin250/refine_1d/particle_deposition_results_virtual_old.xlsx", col_names = TRUE)
colnames(virtual_child_200) <- c('sbj', 'dp', 'TDE', 'DE_ET', 'BDE', 'ADE', 'flow', 'intra')
virtual_child_200$group <- "VC (1.0)"
virtual_child_200_t <- virtual_child_200[, c("sbj", "group", "dp", "intra")]
virtual_child_200_t$variable <- "intra"

# Constricted conditions
constrict_90 <- read_xlsx("/hpc/gjin250/refine_1d/particle_deposition_results_constrict_90.xlsx", col_names = TRUE)
colnames(constrict_90) <- c('sbj', 'dp', 'TDE', 'DE_ET', 'BDE', 'ADE', 'flow', 'intra')
constrict_90$group <- "VC (0.9)"
constrict_90_t <- constrict_90[, c("sbj", "group", "dp", "intra")]
constrict_90_t$variable <- "intra"

constrict_85 <- read_xlsx("/hpc/gjin250/refine_1d/particle_deposition_results_constrict_85.xlsx", col_names = TRUE)
colnames(constrict_85) <- c('sbj', 'dp', 'TDE', 'DE_ET', 'BDE', 'ADE', 'flow', 'intra')
constrict_85$group <- "VC (0.85)"
constrict_85_t <- constrict_85[, c("sbj", "group", "dp", "intra")]
constrict_85_t$variable <- "intra"

constrict_80 <- read_xlsx("/hpc/gjin250/refine_1d/particle_deposition_results_constrict_80.xlsx", col_names = TRUE)
colnames(constrict_80) <- c('sbj', 'dp', 'TDE', 'DE_ET', 'BDE', 'ADE', 'flow', 'intra')
constrict_80$group <- "VC (0.8)"
constrict_80_t <- constrict_80[, c("sbj", "group", "dp", "intra")]
constrict_80_t$variable <- "intra"

# Cluster data
cluster_12 <- read_xlsx("/hpc/gjin250/refine_1d/particle_deposition_results_cluster_12.xlsx", col_names = TRUE)
colnames(cluster_12) <- c('sbj', 'dp', 'TDE', 'DE_ET', 'BDE', 'ADE', 'flow', 'intra')
cluster_12$group <- "cluster_12"
cluster_12_t <- cluster_12[, c("sbj", "group", "dp", "intra")]
cluster_12_t$variable <- "intra"

# ------------------------------
# 3. Paired Comparison of Intra-Thoracic Deposition: VC (1.0) vs. Cluster_12
# ------------------------------
vc <- virtual_child_200_t[, c("sbj", "dp", "intra")]
colnames(vc)[3] <- "intra_vc"
c12 <- cluster_12_t[, c("sbj", "dp", "intra")]
colnames(c12)[3] <- "intra_c12"
paired_data <- merge(vc, c12, by = c("sbj", "dp"))
t_test_intra <- t.test(paired_data$intra_vc, paired_data$intra_c12, paired = TRUE)
print(t_test_intra)
print(summary(paired_data[, c("intra_vc", "intra_c12")]))

# ------------------------------
# 4. Visualization: Intra-Thoracic Deposition Plot
# ------------------------------
# Reshape paired data to long format for plotting
long_data <- paired_data %>%
  pivot_longer(cols = c("intra_vc", "intra_c12"), names_to = "group", values_to = "intra") %>%
  mutate(group = recode(group, "intra_vc" = "VC (1.0)", "intra_c12" = "Cluster_12"))

ggplot(long_data, aes(x = factor(dp), y = intra, group = group, color = group)) +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 15) +
  stat_summary(fun = mean, geom = "line", linetype = "dotted") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "Intrathoracic Deposition Across Particle Sizes",
       subtitle = "VC (1.0) vs. Cluster_12",
       x = expression("Particle Diameter "~(d[p])~" (µm)"),
       y = "Deposition Fraction",
       color = "Condition") +
  scale_x_log10() +
  theme_classic(base_size = 16)

# ------------------------------
# 5. Processing and Plotting BDE and ADE Data (Virtual Child)
# ------------------------------
# Process BDE
virtual_child_200_b <- virtual_child_200[, c("sbj", "group", "dp", "BDE")]
colnames(virtual_child_200_b)[4] <- "value"
virtual_child_200_b$variable <- "Bronchial"

ggplot(virtual_child_200_b, aes(x = dp, y = value, color = group)) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 15) +
  stat_summary(fun = mean, geom = "line", linetype = "dotted") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  labs(x = "dp (log scale)", y = "Bronchial Deposition Fraction",
       title = "Bronchial Deposition by dp (Virtual Child)",
       color = "Condition") +
  scale_x_log10() +
  theme_classic()

# Process ADE
virtual_child_200_a <- virtual_child_200[, c("sbj", "group", "dp", "ADE")]
colnames(virtual_child_200_a)[4] <- "value"
virtual_child_200_a$variable <- "Alveolar"

ggplot(virtual_child_200_a, aes(x = dp, y = value, color = group)) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 15) +
  stat_summary(fun = mean, geom = "line", linetype = "dotted") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  labs(x = "dp (log scale)", y = "Alveolar Deposition Fraction",
       title = "Alveolar Deposition by dp (Virtual Child)",
       color = "Condition") +
  scale_x_log10() +
  theme_classic()

# ------------------------------
# 6. Merge Adult Simulation and Virtual Child Data for Multi-Measure Plotting
# ------------------------------
# Adult simulation data from TDF, BDF, ADF
tdf <- sbj_tdf_melt %>% mutate(Type = "TDF")
bdf <- sbj_bdf_melt %>% mutate(Type = "BDF")
adf <- sbj_adf_melt %>% mutate(Type = "ADF")
adult_df <- bind_rows(tdf, bdf, adf) %>%
  mutate(group = ifelse(Sex == "f", "Female Adult", ifelse(Sex == "m", "Male Adult", NA))) %>%
  select(dp, value, group, Type)

# Virtual child data (using virtual_child_200 data)
vc_tdf <- virtual_child_200_t %>% mutate(Type = "TDF", group = "Virtual Child", value = intra)
vc_bdf <- virtual_child_200 %>% mutate(Type = "BDF", group = "Virtual Child", value = BDE)
vc_adf <- virtual_child_200 %>% mutate(Type = "ADF", group = "Virtual Child", value = ADE)
vc_df <- bind_rows(vc_tdf[, c("dp", "value", "group", "Type")],
                   vc_bdf[, c("dp", "value", "group", "Type")],
                   vc_adf[, c("dp", "value", "group", "Type")])
plot_df <- bind_rows(adult_df, vc_df)
plot_df$Type <- factor(plot_df$Type, levels = c("TDF", "BDF", "ADF"))

ggplot(plot_df, aes(x = dp, y = value, color = group)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  scale_x_log10() +
  facet_wrap(~ Type, ncol = 1,
             labeller = labeller(Type = c("TDF" = "(a) Total deposition fraction",
                                          "BDF" = "(b) Bronchial deposition fraction",
                                          "ADF" = "(c) Alveolar deposition fraction"))) +
  labs(x = "dp (log scale)", y = "Deposition Fraction",
       title = "Deposition Comparison by Particle Size (Children vs Adults)",
       color = "Age Group") +
  scale_color_manual(values = c("Female Adult" = "plum", "Male Adult" = "forestgreen", "Virtual Child" = "#EFC000FF")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom", legend.box = "horizontal")

# ------------------------------
# 7. Pairwise Comparisons Using Emmeans and CLD (ADF Data)
# ------------------------------
anova_model_adf <- aov(value ~ Sex * factor(dp), data = sbj_adf_melt)
emm <- emmeans(anova_model_adf, pairwise ~ Sex | dp, adjust = "tukey")
cld_df <- cld(emm$emmeans, Letters = letters, alpha = 0.05)
head(cld_df)

plot_data_adf <- sbj_adf_melt %>% group_by(Sex, dp) %>% summarise(mean = mean(value), .groups = "drop")
plot_with_labels <- left_join(plot_data_adf, cld_df[, c("Sex", "dp", ".group")])
ggplot(sbj_adf_melt, aes(x = Sex, y = value, fill = Sex)) +
  geom_boxplot(width = 0.5) +
  geom_text(data = plot_with_labels, aes(label = .group, y = mean + 0.1), vjust = 0, size = 4, color = "black") +
  facet_wrap(~ dp, labeller = label_both) +
  theme_classic(base_size = 14)

# ------------------------------
# 8. Paired T-Tests Across dp and Variable
# ------------------------------
data_long <- combined_data %>% select(sbj, dp, group, variable, value)
paired_wide <- data_long %>%
  pivot_wider(names_from = group, values_from = value) %>%
  rename(vc = `VC (1.0)`, c12 = cluster_12)
ttest_results_all <- paired_wide %>%
  group_by(dp, variable) %>%
  summarise(ttest = list(t.test(vc, c12, paired = TRUE)), .groups = "drop") %>%
  mutate(tidy_res = map(ttest, tidy)) %>%
  unnest(tidy_res)
ttest_summary <- ttest_results_all %>% select(dp, variable, estimate, statistic, p.value, conf.low, conf.high, method)
print(ttest_summary)

# ------------------------------
# 9. Lobe Summary Processing and Plotting (Optional)
# ------------------------------
read_lobe_summary_merged <- function(dp_value) {
  file_list <- list(
    "VC (1.0)"   = paste0("/hpc/gjin250/refine_1d/particle_development/output/025-virtual/", dp_value, "/airway_result.txt"),
    "VC (0.90)"  = paste0("/hpc/gjin250/refine_1d/particle_development/output/025-constrict-0.9/", dp_value, "/airway_result.txt"),
    "VC (0.85)"  = paste0("/hpc/gjin250/refine_1d/particle_development/output/025-constrict-0.85/", dp_value, "/airway_result.txt"),
    "VC (0.80)"  = paste0("/hpc/gjin250/refine_1d/particle_development/output/025-constrict/", dp_value, "/airway_result.txt"),
    "Cluster 12" = paste0("/hpc/gjin250/refine_1d/particle_development/output/025-cluster-12/", dp_value, "/airway_result.txt")
  )
  
  column_names <- c(
    "ne", "x", "y", "z", "distance", "gen", "horsf_ord", "lobe", "ne_length",
    "ne_radius", "ne_Vdot", "ne_mass", "nj_loss_dif", "nj_loss_imp",
    "nj_loss_sed", "nj_conc1", "ne_part_vel"
  )
  
  load_and_process <- function(path, label) {
    df <- read.table(path, header = FALSE, sep = "", fill = FALSE, comment.char = "")
    colnames(df) <- column_names
    df <- df %>%
      mutate(across(c(ne_Vdot, ne_mass, nj_loss_dif, nj_loss_imp, nj_loss_sed, nj_conc1, ne_part_vel),
                    ~ as.numeric(gsub("D", "E", .)))) %>%
      mutate(total_loss = nj_loss_dif + nj_loss_imp + nj_loss_sed,
             condition  = label)
    # Optionally override 'lobe' if needed
    df
  }
  
  all_data <- bind_rows(lapply(names(file_list), function(name) {
    load_and_process(file_list[[name]], name)
  }))
  
  lobe_summary <- all_data %>%
    group_by(condition, lobe) %>%
    summarise(
      total_sed      = sum(nj_loss_sed, na.rm = TRUE),
      total_imp      = sum(nj_loss_imp, na.rm = TRUE),
      total_dif      = sum(nj_loss_dif, na.rm = TRUE),
      total_loss_sum = sum(total_loss, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(condition) %>%
    mutate(
      sed_percent  = 100 * total_sed      / sum(total_loss_sum),
      imp_percent  = 100 * total_imp      / sum(total_loss_sum),
      dif_percent  = 100 * total_dif      / sum(total_loss_sum),
      lobe_percent = 100 * total_loss_sum / sum(total_loss_sum)
    ) %>%
    ungroup() %>%
    mutate(dp = dp_value)
  
  return(lobe_summary)
}

dp_values <- c("0.80", "1.00", "3.00", "5.00")
df_all <- bind_rows(lapply(dp_values, read_lobe_summary_merged))

df_long <- df_all %>%
  select(condition, lobe, dp, total_sed, total_imp, total_dif) %>%
  pivot_longer(cols = c("total_sed", "total_imp", "total_dif"),
               names_to = "mechanism", values_to = "loss_value")

ggplot(df_long, aes(x = lobe, y = loss_value, fill = mechanism)) +
  geom_col(position = "stack") +
  facet_grid(dp ~ condition, scales = "free_y") +
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  labs(
    title = "Deposition by Lobe and Condition",
    x = "Lobe (Region)",
    y = "Total Deposition (mass)",
    fill = "Mechanism"
  ) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))

# ------------------------------
#  End of Script
##############################################
