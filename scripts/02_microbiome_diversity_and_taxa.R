# =============================================================================
# 02_microbiome_diversity_and_taxa.R
# -----------------------------------------------------------------------------
# Microbiome sampling timeline, alpha diversity over time and species that
# change relative to baseline, for every cohort.
#
# Input:  Processed_data/phy_mic.rds
# Output: Results/02_microbiome_diversity_and_taxa/
#           Timeline_Mic_SD_PD.pdf, Timeline_Mic_HE_LE.pdf
#           Shannon_diversity_all_cohs.pdf, SD_PD_HC_LC_Melanoma.pdf
#         Results/02_microbiome_diversity_and_taxa/<Cohort>/
#           Timeline_Mic_<Cohort>.pdf
#           Shannon_diversity_<Cohort>.pdf/.csv        (time vs baseline)
#           Shannon_diversity_PD_SD_<Cohort>.pdf/.csv  (Time x Response LME)
#           Shannon_diversity_HE_LE_<Cohort>.pdf/.csv  (Time x Engraftment LME)
#           Tab_LME_Species_<Cohort>.csv               (per-species LME)
#           Fig_Taxa_<Cohort>.pdf                      (heatmap of sig. species)
#
# Statistics: linear mixed models (nlme::lme) with a random intercept per
# subject; BH correction. Wilcoxon tests from rstatix are only used to get
# bracket positions for the plots - the p-values shown come from the LMEs.
# =============================================================================

source(here::here("R/utils.R"))
suppressPackageStartupMessages({
  library(phyloseq)
  library(purrr)
  library(vegan)
  library(tidyverse)
  library(yingtools2)    # get.samp(): sample data + diversity metrics
  library(ggprism)
  library(lemon)         # facet_rep_wrap()
  library(nlme)
  library(rstatix)       # wilcox_test(), add_xy_position(), add_significance()
  library(circlize)
  library(ComplexHeatmap)
})

results_folder <- create_results_dir("Results/02_microbiome_diversity_and_taxa")

time_levels     <- c("Baseline", "D1", "D3", "D7", "D14", "D21", "D30", "D58", "D86")
div_time_levels <- time_levels[time_levels != "D86"]   # D86 excluded from diversity


# =============================================================================
# 1. Data
# =============================================================================
phy_ph1 <- readRDS(here("Processed_data/phy_mic.rds"))

phy_sel <- subset_samples(phy_ph1, Visit.Name.Norm %in% time_levels)
phy_sel <- prune_taxa(taxa_sums(phy_sel) > 0, phy_sel)
taxa_names(phy_sel) <- make.names(taxa_names(phy_sel))
sm_dt <- data.frame(sample_data(phy_sel))

#' Shannon / InvSimpson per sample, restricted to the diversity time points.
#' @param cohort optional cohort name to subset to
get_diversity <- function(phy, cohort = NULL) {
  div_dt <- get.samp(phy, stats = TRUE, measures = c("Shannon", "InvSimpson")) %>%
    as.data.frame()
  if (!is.null(cohort)) div_dt <- div_dt %>% filter(Cohort %in% cohort)
  div_dt <- div_dt %>% filter(Visit.Name.Norm %in% div_time_levels)
  div_dt$Visit.Name.Norm <- factor(div_dt$Visit.Name.Norm, levels = div_time_levels)
  div_dt$Response_Wk8 <- gsub("PD/Off TX", "PD", div_dt$Response_Wk8)
  div_dt$Subject_ID <- factor(div_dt$Subject_ID)
  div_dt
}


# =============================================================================
# 2. Sampling timeline (all cohorts)
# =============================================================================
timeline_dt <- sm_dt %>%
  select(Subject_ID, Visit.Name.Norm, Response_Wk8, Cohort, Engraftment) %>%
  unique()
timeline_dt$Subject_ID   <- as.character(timeline_dt$Subject_ID)
timeline_dt$Timepoint    <- factor(timeline_dt$Visit.Name.Norm, levels = time_levels)
timeline_dt$Response_Wk8 <- gsub("PD/Off TX", "PD", timeline_dt$Response_Wk8)
timeline_dt$Engraftment  <- factor(timeline_dt$Engraftment)
levels(timeline_dt$Engraftment) <- c("HE", "LE")

#' Subject x time point dot/line plot, coloured and faceted by `group`
plot_timeline <- function(dt, x, group, colours, legend_name) {
  ggplot(dt, aes(x = .data[[x]], y = Subject_ID, group = Subject_ID,
                 color = .data[[group]])) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    scale_color_manual(name = legend_name, values = colours) +
    facet_rep_wrap(as.formula(paste("Cohort ~", group)), scales = "free_y",
                   ncol = 2, repeat.tick.labels = "x") +
    ylab("Subjects") +
    theme_prism() +
    theme_bold_axes()
}

save_pdf(plot_timeline(timeline_dt, "Timepoint", "Response_Wk8", col_resp, "Response"),
         file.path(results_folder, "Timeline_Mic_SD_PD.pdf"), height = 12, width = 8)
save_pdf(plot_timeline(timeline_dt, "Timepoint", "Engraftment", col_eng, "Engraftment"),
         file.path(results_folder, "Timeline_Mic_HE_LE.pdf"), height = 12, width = 8)


# =============================================================================
# 3. Per-cohort analysis
# =============================================================================

#' Fit Shannon ~ Visit * <group> LME, save the coefficient table, and return
#' a bracket table (interaction p < 0.05) for plotting.
fit_group_lme <- function(div_dt, group, csv_file, cohort) {
  model_lme <- lme(as.formula(paste("Shannon ~ Visit.Name.Norm *", group)),
                   random = ~ 1 | Subject_ID, data = div_dt)
  dt_sum <- data.frame(summary(model_lme)$tTable)
  dt_sum$Var <- rownames(dt_sum)
  write.csv(dt_sum, csv_file)
  
  # Interaction terms -> one row per visit
  dt_sum <- dt_sum[grep(":", dt_sum$Var), ]
  dt_sum$Visit.Name.Norm <- gsub(":.*", "", gsub("Visit.Name.Norm", "", dt_sum$Var))
  dt_sum$group2 <- dt_sum$Visit.Name.Norm
  
  div_dt %>%
    wilcox_test(Shannon ~ Visit.Name.Norm, ref.group = "Baseline") %>%
    add_xy_position(fun = "mean_sd") %>%
    select(!c(p, p.adj, p.adj.signif)) %>%
    left_join(dt_sum) %>%
    mutate(p_label = sprintf("%.4f", p.value)) %>%
    add_significance(p.col = "p.value") %>%
    mutate(!!group := "") %>%
    filter(p.value <= 0.05) %>%
    mutate(Cohort = cohort)
}

#' Shannon box plot split by `group`, with LME interaction p-values if any
plot_group_diversity <- function(div_dt, group, colours, legend_name, df_p_val) {
  p <- ggplot(div_dt, aes(x = Visit.Name.Norm, y = Shannon, fill = .data[[group]])) +
    geom_boxplot(outlier.colour = NA, alpha = 0.2) +
    geom_point(position = position_dodge(width = 0.75),
               aes(group = .data[[group]]),
               shape = 21, color = "black", size = 2, alpha = 0.5) +
    scale_fill_manual(name = legend_name, values = colours) +
    ylab("Shannon") + xlab("Time") +
    theme_prism() +
    theme_bold_axes()
  if (nrow(df_p_val) > 0) p <- p + add_pvalue(df_p_val, label = "p.value.signif", tip.length = 0.01)
  p
}

cohorts <- unique(timeline_dt$Cohort)

p_val_list_div <- list()   # Time vs baseline, per cohort
p_val_sd_pd    <- list()   # Time x Response interaction, per cohort
p_val_hc_lc    <- list()   # Time x Engraftment interaction, per cohort

for (cohs in cohorts) {
  
  coh_folder <- create_results_dir("Results/02_microbiome_diversity_and_taxa", cohs)
  div_dt <- get_diversity(phy_sel, cohort = cohs)
  
  # ---- 3a. Cohort timeline ---------------------------------------------------
  p_time <- ggplot(div_dt, aes(x = Visit.Name.Norm, y = Subject_ID,
                               group = Subject_ID, color = Response_Wk8)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    scale_color_manual(name = "Response", values = col_resp) +
    facet_rep_wrap(Cohort ~ Response_Wk8, scales = "free", ncol = 2, repeat.tick.labels = "x") +
    ylab("Subjects") +
    theme_prism() +
    theme_bold_axes()
  save_pdf(p_time, file.path(coh_folder, paste0("Timeline_Mic_", cohs, ".pdf")), height = 6, width = 8)
  
  # ---- 3b. Shannon over time (each visit vs baseline) -----------------------
  fit_div <- lme(Shannon ~ Visit.Name.Norm, random = ~ 1 | Subject_ID, data = div_dt)
  dt_sum <- data.frame(summary(fit_div)$tTable)
  dt_sum$Var <- rownames(dt_sum)
  dt_sum <- dt_sum[grep("Visit", dt_sum$Var), ] %>%
    mutate(padj = p.adjust(p.value, method = "BH"),
           Subject_ID = "", Response_Wk8 = "")
  dt_sum$Visit.Name.Norm <- gsub("Visit.Name.Norm", "", dt_sum$Var)
  write.csv(dt_sum, file.path(coh_folder, paste0("Shannon_diversity_", cohs, ".csv")))
  dt_sum$group2 <- dt_sum$Visit.Name.Norm
  
  df_p_val <- div_dt %>%
    wilcox_test(Shannon ~ Visit.Name.Norm, ref.group = "Baseline") %>%
    add_xy_position(fun = "mean_sd") %>%
    select(!c(p, p.adj, p.adj.signif)) %>%
    left_join(dt_sum) %>%
    mutate(p_label = sprintf("%.4f", padj)) %>%
    add_significance(p.col = "padj") %>%
    mutate(Cohort = cohs)
  p_val_list_div[[cohs]] <- df_p_val
  
  p_div <- ggplot(div_dt, aes(x = Visit.Name.Norm, y = Shannon)) +
    geom_boxplot(outlier.colour = NA, alpha = 0.2) +
    geom_point(size = 2, alpha = 0.7) +
    ylab("Shannon") + xlab("Time") +
    theme_prism() +
    theme_bold_axes()
  if (nrow(df_p_val) > 0) p_div <- p_div + add_pvalue(df_p_val, label = "padj.signif", tip.length = 0.01)
  save_pdf(p_div, file.path(coh_folder, paste0("Shannon_diversity_", cohs, ".pdf")), height = 5, width = 4)
  
  # ---- 3c. Shannon: Time x Response (SD vs PD) ------------------------------
  sel_div_dt <- div_dt
  sel_div_dt$Response_Wk8 <- factor(sel_div_dt$Response_Wk8, levels = c("PD", "SD"))
  df_p_val <- fit_group_lme(sel_div_dt, "Response_Wk8",
                            file.path(coh_folder, paste0("Shannon_diversity_PD_SD_", cohs, ".csv")), cohs)
  p_val_sd_pd[[cohs]] <- df_p_val
  save_pdf(plot_group_diversity(sel_div_dt, "Response_Wk8", col_resp, "Response", df_p_val),
           file.path(coh_folder, paste0("Shannon_diversity_PD_SD_", cohs, ".pdf")), height = 5, width = 5)
  
  # ---- 3d. Shannon: Time x Engraftment (HE vs LE) ---------------------------
  sel_div_dt$Engraftment <- factor(sel_div_dt$Engraftment)
  levels(sel_div_dt$Engraftment) <- c("HE", "LE")
  df_p_val <- fit_group_lme(sel_div_dt, "Engraftment",
                            file.path(coh_folder, paste0("Shannon_diversity_HE_LE_", cohs, ".csv")), cohs)
  p_val_hc_lc[[cohs]] <- df_p_val
  save_pdf(plot_group_diversity(sel_div_dt, "Engraftment", col_eng, "Engraftment", df_p_val),
           file.path(coh_folder, paste0("Shannon_diversity_HE_LE_", cohs, ".pdf")), height = 5, width = 6)
  
  # ---- 3e. Species changing over time (LME on arcsine-sqrt RA) --------------
  phy_coh <- prune_samples(div_dt$sample, phy_sel)
  phy_coh <- prune_taxa(taxa_sums(phy_coh) > 0, phy_coh)
  
  # Keep species present in > 5% of samples AND > 0.1% of total abundance
  prevalence <- apply(otu_table(phy_coh),
                      MARGIN = ifelse(taxa_are_rows(phy_coh), 1, 2),
                      FUN = function(x) sum(x > 0))
  prevdf <- data.frame(OTU = names(prevalence), Prevalence = prevalence,
                       TotalAbundance = taxa_sums(phy_coh))
  prevdf$Prop <- prevdf$TotalAbundance / sum(prevdf$TotalAbundance)
  prev_cut_off <- ceiling(0.05 * nsamples(phy_coh))
  keep_taxa <- prevdf$OTU[prevdf$Prevalence > prev_cut_off & prevdf$Prop > 1e-3]
  phy_fil <- prune_taxa(keep_taxa, phy_coh)
  
  melt_dt <- psmelt(phy_fil)
  result_lme <- list()
  for (taxa in taxa_names(phy_fil)) {
    phy_m <- melt_dt[melt_dt$OTU %in% taxa, ]
    phy_m$PID <- factor(phy_m$Subject_ID)
    phy_m$t_Abundance <- asin(sqrt(phy_m$Abundance))
    phy_m$Visit.Name.Norm <- factor(phy_m$Visit.Name.Norm, levels = div_time_levels)
    mod_bac <- lme(t_Abundance ~ Visit.Name.Norm, random = ~ 1 | PID, data = phy_m)
    sum_mod_dt <- data.frame(summary(mod_bac)$tTable)
    sum_mod_dt$Bac <- taxa
    sum_mod_dt$Var <- rownames(sum_mod_dt)
    result_lme[[taxa]] <- sum_mod_dt
  }
  
  result_lme_dt <- do.call("rbind", result_lme)
  names(result_lme_dt)[5] <- "pval"
  result_lme_dt <- result_lme_dt[grep("Visit", result_lme_dt$Var), ] %>%
    group_by(Var) %>%                                   # BH within each visit
    mutate(p_adj = p.adjust(pval, method = "BH"))
  result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj, decreasing = FALSE), ]
  write.csv(result_lme_dt, file.path(coh_folder, paste0("Tab_LME_Species_", cohs, ".csv")))
  
  # ---- 3f. Heatmap of significant species (p_adj < 0.05) --------------------
  result_sig <- result_lme_dt[result_lme_dt$p_adj < 0.05, ]
  ps_sig <- prune_taxa(unique(result_sig$Bac), phy_fil)
  
  mat <- data.frame(otu_table(ps_sig))
  met <- data.frame(sample_data(ps_sig))
  met <- met[order(met$Visit.Name.Norm, met$Subject_ID), ]
  mat <- mat[, match(rownames(met), colnames(mat))]
  
  # Row annotation: LME coefficient of each visit vs baseline
  coeff_dt_w <- result_sig[, c("Bac", "Var", "Value")] %>%
    pivot_wider(names_from = Var, values_from = Value) %>%
    data.frame()
  coeff_dt_w <- coeff_dt_w[match(rownames(mat), coeff_dt_w$Bac), ]
  
  col_coef <- colorRamp2(c(-0.5, -0.25, 0, 0.25, 0.5),
                         colorRampPalette(c("blue", "white", "red"))(5))
  ha_left <- rowAnnotation(D1  = coeff_dt_w$Visit.Name.NormD1,
                           D3  = coeff_dt_w$Visit.Name.NormD3,
                           D7  = coeff_dt_w$Visit.Name.NormD7,
                           D14 = coeff_dt_w$Visit.Name.NormD14,
                           D21 = coeff_dt_w$Visit.Name.NormD21,
                           D30 = coeff_dt_w$Visit.Name.NormD30,
                           D58 = coeff_dt_w$Visit.Name.NormD58,
                           col = list(D1 = col_coef, D3 = col_coef, D7 = col_coef, D14 = col_coef,
                                      D21 = col_coef, D30 = col_coef, D58 = col_coef),
                           show_legend = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                           show_annotation_name = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                           annotation_name_gp = gpar(fontface = "bold"),
                           annotation_name_rot = 45,
                           annotation_legend_param = list(title = "Coeff",
                                                          title_gp = gpar(fontsize = 14, fontface = "bold"),
                                                          labels_gp = gpar(fontsize = 12, fontface = "bold"),
                                                          legend_height = unit(10, "cm"),
                                                          at = c(-0.5, -0.25, 0, 0.25, 0.5)))
  
  tax_dt <- tax_table(ps_sig)@.Data %>% data.frame()
  tax_dt <- tax_dt[match(rownames(mat), rownames(tax_dt)), ]
  split_rows <- factor(tax_dt$Class)
  split_cols <- factor(met$Visit.Name.Norm, levels = div_time_levels)
  colnames(mat) <- met$Subject_ID
  
  # VE800 strains highlighted in red
  fontcolors <- ifelse(grepl("VE800", rownames(mat)), "red", "black")
  
  ht <- Heatmap(log10(mat + 0.00001), name = "Log10(RA)",
                column_split = split_cols,
                row_split = split_rows,
                left_annotation = ha_left,
                col = viridis::viridis(10),
                row_names_side = "left",
                row_gap = unit(2, "mm"),
                row_title_gp = gpar(fontsize = 14, fontface = "bold"),
                row_title_rot = 0,
                column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                column_names_gp = gpar(fontsize = 8, fontface = "bold"),
                row_names_gp = gpar(fontsize = 10, fontface = "bold", col = fontcolors),
                column_title_rot = 0,
                cluster_rows = TRUE,
                show_parent_dend_line = FALSE,
                show_row_dend = FALSE,
                cluster_columns = FALSE,
                border = TRUE,
                heatmap_legend_param = list(title = "RA", at = c(0, -1, -2, -3, -4, -5),
                                            legend_width = unit(10, "cm"),
                                            legend_direction = "horizontal",
                                            title_gp = gpar(fontsize = 14, fontface = "bold"),
                                            labels_gp = gpar(fontsize = 12, fontface = "bold"),
                                            labels = c(expression("10"^0), expression("10"^-1),
                                                       expression("10"^-2), expression("10"^-3),
                                                       expression("10"^-4), expression("10"^-5))))
  
  pdf(file.path(coh_folder, paste0("Fig_Taxa_", cohs, ".pdf")), height = 12, width = 20, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "top", annotation_legend_side = "right", legend_grouping = "original")
  dev.off()
}


# =============================================================================
# 4. Combined figures
# =============================================================================

# ---- 4a. Shannon over time, all cohorts side by side ------------------------
p_val_dt_comb <- do.call("rbind", p_val_list_div)
# Use the Gastric bracket heights for all panels.
# NOTE: assumes 3 cohorts with the same set of comparisons.
y_pos_gastric <- p_val_dt_comb$y.position[p_val_dt_comb$Cohort == "Gastric"]
p_val_dt_comb$y.position <- rep(y_pos_gastric, 3)

div_dt <- get_diversity(phy_sel)

# D58 is displayed as "D56" (nominal study day)
x_labels <- c("Baseline", "D1", "D3", "D7", "D14", "D21", "D30", "D56")

p_time <- ggplot(div_dt, aes(x = Visit.Name.Norm, y = Shannon)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.2) +
  geom_point(size = 2, alpha = 1) +
  facet_rep_wrap(~Cohort, repeat.tick.labels = TRUE) +
  ylab("Shannon") + xlab("Time") +
  theme_prism() +
  theme(axis.text.x  = element_text(size = 14, angle = 90, hjust = 1),
        axis.text.y  = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold")) +
  add_pvalue(p_val_dt_comb, label = "padj.signif", tip.length = 0.01) +
  scale_x_discrete(labels = x_labels)
save_pdf(p_time, file.path(results_folder, "Shannon_diversity_all_cohs.pdf"), height = 5, width = 9)

# ---- 4b. Melanoma: SD vs PD and HC vs LC in one figure ----------------------
comb_p_val_dt <- rbind(p_val_sd_pd$Melanoma, p_val_hc_lc$Melanoma)
comb_p_val_dt$Cat <- "COLONIZATION OF VE800"
comb_p_val_dt$Grp <- ""

div_dt <- get_diversity(phy_sel, cohort = "Melanoma")
div_dt_1 <- div_dt %>% mutate(Cat = "CLINICAL EFFICACY")
div_dt_2 <- div_dt %>% mutate(Cat = "COLONIZATION OF VE800")
div_dt_comb <- rbind(div_dt_1, div_dt_2)
div_dt_comb$Grp <- c(div_dt_1$Response_Wk8, div_dt_1$Engraftment)
div_dt_comb$Grp <- factor(div_dt_comb$Grp, levels = c("SD", "PD", "High Engraft", "Low Engraft"))
levels(div_dt_comb$Grp) <- c("SD", "PD", "HC", "LC")

col_grp <- c("PD" = "#ff5e1e", "SD" = "#0099b0", "HC" = "#009031", "LC" = "#7f1f98")

p_comb_lab <- ggplot(div_dt_comb, aes(x = Visit.Name.Norm, y = Shannon, fill = Grp)) +
  geom_boxplot(outlier.colour = NA, alpha = 1) +
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 1) +
  scale_fill_manual(name = "Response", values = col_grp) +
  facet_rep_wrap(~Cat, repeat.tick.labels = TRUE) +
  ylab("Shannon") + xlab("Time") +
  theme_prism() +
  theme(axis.text.x  = element_text(size = 14, angle = 90, hjust = 1),
        axis.text.y  = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold")) +
  add_pvalue(comb_p_val_dt, label = "p.value.signif", tip.length = 0.01) +
  scale_x_discrete(labels = x_labels)
save_pdf(p_comb_lab, file.path(results_folder, "SD_PD_HC_LC_Melanoma.pdf"), height = 5, width = 8)

