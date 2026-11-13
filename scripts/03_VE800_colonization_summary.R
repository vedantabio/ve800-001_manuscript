# =============================================================================
# 03_VE800_colonization_summary.R
# -----------------------------------------------------------------------------
# Summarise VE800 strain colonization in melanoma patients using the marker
# panel: per subject and strain, mean relative abundance and mean prevalence
# (fraction of visits with the strain detected) over the dosing window
# D7 - D30.
#
# Input:  Processed_data/VE303_abun.csv
# Output: Results/03_VE800_colonization_summary/
#           Timeline_Mic_SD_PD.pdf               sampling timeline (D7-D30 shaded)
#           VE800_Abundance_RA_Subject_Time.csv  per subject / visit / strain
#           VE800_Mean_RA_Subject.csv            per subject / strain means
#                                                (used by 07_colonization_predictors_ordinalNet.R)
#           VE800_Mean_RA_Subject.pdf            heatmap: mean relative abundance
#           VE800_Mean_Prev_Subject.pdf          heatmap: mean prevalence
# =============================================================================

source(here::here("R/utils.R"))
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggprism)
  library(lemon)
  library(ComplexHeatmap)
})

results_folder <- create_results_dir("Results/03_VE800_colonization_summary")

time_levels <- c("Baseline", "D1", "D3", "D7", "D14", "D21", "D30", "D58", "D86")


# ---- 1. Marker-panel data (melanoma only) -----------------------------------
m_panel_dt <- read.csv(here("Processed_data/VE303_abun.csv")) %>%
  select(organism, Subject_ID, Visit.Name.Norm, Cohort, Response_Wk8, detection_status,
         targeted_panel_relative_abundance, normalized_marker_depth) %>%
  filter(Cohort %in% "Melanoma")
m_panel_dt$Visit <- as.character(m_panel_dt$Visit.Name.Norm)


# ---- 2. Sampling timeline (dosing window D7-D30 shaded) ---------------------
timeline_dt <- m_panel_dt %>%
  select(Subject_ID, Visit.Name.Norm, Cohort, Response_Wk8) %>%
  filter(Visit.Name.Norm %in% time_levels) %>%
  unique()
timeline_dt$Subject_ID   <- as.character(timeline_dt$Subject_ID)
timeline_dt$Timepoint    <- factor(timeline_dt$Visit.Name.Norm, levels = time_levels)
timeline_dt$Response_Wk8 <- gsub("PD/Off TX", "PD", timeline_dt$Response_Wk8)

p_time <- ggplot(timeline_dt, aes(x = Timepoint, y = Subject_ID,
                                  group = Subject_ID, color = Response_Wk8)) +
  annotate("rect", xmin = 3.5, xmax = 7.5, ymin = -Inf, ymax = Inf, alpha = .3) +
  geom_point(size = 3, alpha = 0.7) +
  geom_line() +
  scale_color_manual(name = "Response", values = col_resp) +
  facet_rep_wrap(Cohort ~ Response_Wk8, scales = "free_y", ncol = 2, repeat.tick.labels = "x") +
  ylab("Subjects") +
  theme_prism() +
  theme_bold_axes()
save_pdf(p_time, file.path(results_folder, "Timeline_Mic_SD_PD.pdf"), height = 5, width = 8)


# ---- 3. Per subject x strain colonization over D7-D30 -----------------------
num_dt <- m_panel_dt %>%
  filter(Visit %in% c("D7", "D14", "D21", "D30")) %>%
  mutate(N_det = ifelse(detection_status == "Detected", 1, 0),
         # Abundance / depth are set to 0 when the strain is not detected
         est_rel_abundance_adj       = ifelse(N_det %in% 0, 0, targeted_panel_relative_abundance),
         normalized_marker_depth_adj = ifelse(N_det %in% 0, 0, normalized_marker_depth))

abun_dt <- num_dt %>%
  select(organism, Subject_ID, Visit.Name.Norm, Response_Wk8, N_det, est_rel_abundance_adj)
write.csv(abun_dt, file.path(results_folder, "VE800_Abundance_RA_Subject_Time.csv"))

sum_dt <- num_dt %>%
  group_by(Subject_ID, organism) %>%
  summarise(mean_abun = mean(est_rel_abundance_adj),
            mean_prev = sum(N_det) / n())
write.csv(sum_dt, file.path(results_folder, "VE800_Mean_RA_Subject.csv"))


# ---- 4. Heatmaps -------------------------------------------------------------
#' Subjects x strains matrix of one summary column from sum_dt
to_matrix <- function(value_col) {
  wide <- sum_dt %>%
    pivot_wider(id_cols = Subject_ID, names_from = organism, values_from = all_of(value_col)) %>%
    data.frame()
  mat <- as.matrix(wide[2:ncol(wide)])
  rownames(mat) <- wide$Subject_ID
  colnames(mat) <- gsub("\\.", "-", colnames(mat))   # undo make.names: VE800.01 -> VE800-01
  mat
}

#' Shared heatmap settings
make_heatmap <- function(mat, cluster_columns, legend_param) {
  Heatmap(mat,
          column_title_gp = gpar(fontsize = 14, fontface = "bold"),
          row_names_gp = gpar(fontsize = 10, fontface = "bold"),
          show_parent_dend_line = FALSE,
          column_names_rot = 60,
          border = TRUE,
          col = viridis::magma(10),
          cluster_columns = cluster_columns,
          cluster_rows = FALSE,
          column_names_side = "top",
          show_column_dend = FALSE,
          heatmap_legend_param = c(legend_param,
                                   list(legend_width = unit(5, "cm"),
                                        legend_direction = "horizontal",
                                        title_gp = gpar(fontsize = 14, fontface = "bold"),
                                        labels_gp = gpar(fontsize = 12, fontface = "bold"))))
}

# Mean relative abundance (log10); strains clustered
ht_abun <- make_heatmap(log10(to_matrix("mean_abun") + 1e-5), cluster_columns = TRUE,
                        legend_param = list(title = "RA", at = c(0, -1, -2, -3, -4, -5),
                                            labels = c(expression("10"^0), expression("10"^-1),
                                                       expression("10"^-2), expression("10"^-3),
                                                       expression("10"^-4), expression("10"^-5))))
pdf(file.path(results_folder, "VE800_Mean_RA_Subject.pdf"), height = 5, width = 4, useDingbats = FALSE)
draw(ht_abun, heatmap_legend_side = "top")
dev.off()

# Mean prevalence; strains in the same order as the abundance heatmap
mat_prev <- to_matrix("mean_prev")[, column_order(ht_abun)]
ht_prev <- make_heatmap(mat_prev, cluster_columns = FALSE,
                        legend_param = list(title = "Mean_Prev", at = c(0, 0.25, 0.5, 0.75, 1),
                                            labels = c(0, 0.25, 0.5, 0.75, 1)))
pdf(file.path(results_folder, "VE800_Mean_Prev_Subject.pdf"), height = 5, width = 4, useDingbats = FALSE)
draw(ht_prev, heatmap_legend_side = "top")
dev.off()
