# Script 3. Target gene longitudinal expression (Pre vs. On anti-PD-1 treatment)
# Dataset: BioKey Clinical Study (NCT03197389; Bassez et al., 2021)
# Author: Clayton Voidelo

# ---

#packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)

# Cohort parameters
# ER+ (E): c("BIOKEY_5", "BIOKEY_12", "BIOKEY_18") | genes: c("CXCL9", "SPP1", "MDK")
# TNBC (E): c("BIOKEY_10", "BIOKEY_11", "BIOKEY_15", "BIOKEY_16", "BIOKEY_31") | genes: c("CXCL10", "CXCL9", "IGF2", "LGALS9", "SPP1")
patients_of_interest <- c("BIOKEY_5", "BIOKEY_12", "BIOKEY_18")
genes_of_interest    <- c("CXCL9", "SPP1", "MDK")
tumor_subtype        <- "ER+"

# Load counts and metadata
counts_matrix <- readRDS("/home/labceca/DocumentosSSD/Clayton Voidelo/1863-counts_cells_cohort1.rds")
metadata      <- read.csv("/home/labceca/DocumentosSSD/Clayton Voidelo/1872-BIOKEY_metaData_cohort1_web.csv")
rownames(metadata) <- metadata$Cell
metadata           <- metadata[colnames(counts_matrix), ]

# Extract normalized expression per timepoint
extract_timepoint <- function(tp) {
  idx <- metadata$patient_id %in% patients_of_interest & metadata$timepoint == tp
  obj <- CreateSeuratObject(
    counts       = counts_matrix[, idx],
    meta.data    = metadata[idx, ],
    min.cells    = 3,
    min.features = 200
  )
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  FetchData(obj, vars = c(genes_of_interest, "patient_id", "timepoint", "expansion"))
}

plot_data_pre <- extract_timepoint("Pre")
plot_data_on  <- extract_timepoint("On")
rm(counts_matrix, metadata)
gc(verbose = FALSE)

long_data <- bind_rows(plot_data_pre, plot_data_on) %>%
  mutate(timepoint = factor(timepoint, levels = c("Pre", "On"))) %>%
  pivot_longer(
    cols      = all_of(genes_of_interest),
    names_to  = "gene",
    values_to = "expression"
  ) %>%
  mutate(gene = factor(gene, levels = genes_of_interest))

# Wilcoxon rank-sum test (Pre vs. On per patient x gene) with BH correction
stat_results <- long_data %>%
  group_by(gene, patient_id) %>%
  rstatix::wilcox_test(expression ~ timepoint, detailed = TRUE) %>%
  ungroup() %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

effect_sizes <- long_data %>%
  group_by(gene, patient_id) %>%
  rstatix::wilcox_effsize(expression ~ timepoint) %>%
  ungroup() %>%
  select(gene, patient_id, effsize, magnitude)

patient_summary <- long_data %>%
  group_by(gene, patient_id, expansion, timepoint) %>%
  summarise(
    n_cells     = n(),
    mean_expr   = mean(expression, na.rm = TRUE),
    median_expr = median(expression, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  pivot_wider(
    names_from  = timepoint,
    values_from = c(n_cells, mean_expr, median_expr)
  ) %>%
  mutate(
    delta_mean   = mean_expr_On - mean_expr_Pre,
    delta_median = median_expr_On - median_expr_Pre
  )

supp_table <- patient_summary %>%
  left_join(stat_results, by = c("gene", "patient_id")) %>%
  left_join(effect_sizes, by = c("gene", "patient_id")) %>%
  mutate(
    subtype          = tumor_subtype,
    statistical_unit = "Single cells stratified by patient (patient x gene)"
  ) %>%
  select(
    subtype, gene, patient_id, expansion, statistical_unit,
    n_cells_Pre, n_cells_On,
    mean_expr_Pre, mean_expr_On, delta_mean,
    median_expr_Pre, median_expr_On, delta_median,
    estimate, conf.low, conf.high,
    statistic, p, p.adj, p.adj.signif,
    effsize, magnitude
  )

write.csv(supp_table, paste0("Supplementary_Table_Stats_", tumor_subtype, ".csv"), row.names = FALSE)

# Single-cell violin plots (Fig. 8A / 8B)
expansion_colors <- c("E" = "#F8766D", "NE" = "#00BFC4")
stat_results_plot <- stat_results %>%
  add_xy_position(x = "timepoint", fun = "max", data = long_data, step.increase = 0.15)

violin_plot <- ggplot(long_data, aes(x = timepoint, y = expression)) +
  geom_violin(aes(fill = expansion), trim = FALSE, scale = "width", alpha = 0.8) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.2, size = 0.5) +
  facet_grid(gene ~ patient_id, scales = "free_y") +
  scale_fill_manual(name = "Clonal Expansion Status", values = expansion_colors) +
  stat_pvalue_manual(
    stat_results_plot,
    label       = "p.adj.signif",
    tip.length  = 0.01,
    size        = 3,
    hide.ns     = FALSE,
    inherit.aes = FALSE
  ) +
  labs(
    title    = paste("Longitudinal Single-Cell Gene Expression:", tumor_subtype, "Cohort"),
    subtitle = "Pre-treatment vs. On-treatment (anti-PD-1 pembrolizumab)",
    x        = "Treatment Timepoint",
    y        = "Log-Normalized Expression",
    caption  = "Two-sided Wilcoxon rank-sum test (BH-adjusted across gene x patient comparisons): ns: p.adj > 0.05, * p.adj <= 0.05, ** p.adj <= 0.01, *** p.adj <= 0.001, **** p.adj <= 0.0001"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    plot.caption    = element_text(size = 8, hjust = 1),
    strip.text.x    = element_text(size = 10, face = "bold"),
    strip.text.y    = element_text(size = 10, face = "bold", angle = 0),
    axis.text.x     = element_text(size = 9, face = "bold"),
    legend.position = "top"
  )

ggsave(paste0("Fig8_SingleCell_Violin_", tumor_subtype, ".pdf"), violin_plot, width = 12, height = 10, units = "in", dpi = 300)

# Patient-level paired trajectories (Supp.Fig. 13)
patient_long <- long_data %>%
  group_by(gene, patient_id, expansion, timepoint) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")

paired_plot <- ggplot(patient_long, aes(x = timepoint, y = mean_expr, group = patient_id)) +
  geom_line(aes(color = patient_id), linewidth = 1, alpha = 0.85) +
  geom_point(aes(color = patient_id), size = 3) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  labs(
    title    = paste("Patient-Level Paired Expression Trajectories:", tumor_subtype, "Cohort (E Group)"),
    subtitle = "Mean log-normalized expression per patient (Pre-treatment vs. On-treatment)",
    x        = "Treatment Timepoint",
    y        = "Mean Log-Normalized Expression",
    color    = "Patient ID"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    strip.text      = element_text(size = 10, face = "bold"),
    axis.text.x     = element_text(size = 9, face = "bold"),
    legend.position = "bottom"
  )

ggsave(paste0("SupFig13_PatientPaired_Trajectory_", tumor_subtype, ".pdf"), paired_plot, width = 10, height = 4.2, units = "in", dpi = 300)
