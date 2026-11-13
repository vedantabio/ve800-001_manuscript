# =============================================================================
# 07_colonization_predictors_ordinalNet.R
# -----------------------------------------------------------------------------
# Ordinal elastic-net (ordinalNet) models: baseline features vs VE800 colonization
# -----------------------------------------------------------------------------
# For each VE800 strain, a penalized cumulative-logit (proportional-odds) model
# predicts the strain's colonization category (mean prevalence across visits)
# from each baseline data set separately:
#   serum SCFA, serum bile acids, and microbiome (Phylum/Class/Order/Genus/Species).
# Strains are fitted in parallel (one strain per worker).
#
# Inputs: Processed_data/phy_mic.rds, SCFA_serum_melanoma.csv, BA_serum_melanoma.csv,
#         VE800_taxonomy.csv (strain lineage; reference file, not patient data)
#         Results/03_VE800_colonization_summary/VE800_Mean_RA_Subject.csv
#         (produced by 03_VE800_colonization_summary.R - run that first)
#
# Outputs (Results/07_colonization_predictors_ordinalNet/):
#   phy_glom_list.rds, baseline_data_list.rds, coeff_list.rds, perf_list.rds,
#   ordinalNet_parallel.log, Heatmap_Summary_Ordinal_net.pdf,
#   <strain>/Barplot_<strain>_<dataset>.pdf
# =============================================================================

source(here::here("R/utils.R"))
library(phyloseq)
library(tidyverse)
library(ordinalNet)
library(doParallel)      # also attaches foreach and parallel
library(ComplexHeatmap)
library(circlize)
library(gtools)
library(Cairo)
library(ggprism)
library(wesanderson)

# ---- Parameters -------------------------------------------------------------
results_folder <- create_results_dir("Results/07_colonization_predictors_ordinalNet")

prev_frac   <- 0.1    # keep taxa present in >= 10% of samples ...
abun_cutoff <- 1e-3   # ... AND contributing >= 0.1% of total abundance
min_class_n <- 3      # colonization classes with < 3 subjects are merged into a neighbour
tax_levels  <- c("Phylum", "Class", "Order", "Genus")


# =============================================================================
# 1. Microbiome: baseline melanoma samples, prevalence/abundance filtering
# =============================================================================
phy_tax <- readRDS(here("Processed_data/phy_mic.rds"))
phy_tax <- subset_samples(phy_tax, Cohort %in% c("Melanoma"))
phy_tax <- subset_samples(phy_tax, Visit.Name.Norm %in% c("Baseline"))
phy_tax <- prune_taxa(taxa_sums(phy_tax) > 0, phy_tax)
taxa_names(phy_tax) <- make.names(taxa_names(phy_tax))

# Prevalence (number of samples with non-zero counts) and relative total abundance
prevalence <- apply(otu_table(phy_tax),
                    MARGIN = ifelse(taxa_are_rows(phy_tax), 1, 2),
                    FUN = function(x) sum(x > 0))
prevdf <- data.frame(OTU            = names(prevalence),
                     Prevalence     = prevalence,
                     TotalAbundance = taxa_sums(phy_tax))
prevdf$Prop <- prevdf$TotalAbundance / sum(prevdf$TotalAbundance)

keep_taxa <- prevdf$OTU[prevdf$Prevalence >= prev_frac * nsamples(phy_tax) &
                          prevdf$Prop >= abun_cutoff]
phy_fil <- prune_taxa(keep_taxa, phy_tax)
phy_fil <- prune_taxa(taxa_sums(phy_fil) > 0, phy_fil)

# Agglomerate to higher taxonomic levels; taxa renamed "<Level>_<Name>"
phy_glom <- lapply(setNames(tax_levels, tax_levels), function(lev) {
  phy_g <- tax_glom(phy_fil, lev)
  taxa_names(phy_g) <- paste0(lev, "_", as.character(tax_table(phy_g)[, lev]))
  prune_taxa(taxa_sums(phy_g) > 0, phy_g)
})
phy_glom[["Species"]] <- phy_fil
saveRDS(phy_glom, file.path(results_folder, "phy_glom_list.rds"))


# =============================================================================
# 2. Import functions for each baseline data set (wide: one row per subject)
# =============================================================================

# Microbiome feature table (subjects x taxa) at a given taxonomic level
import_mic_data <- function(phy_glom, tax_lev) {
  message("Assembling input data @ ", tax_lev)
  psmelt(phy_glom[[tax_lev]]) %>%
    pivot_wider(id_cols = Subject_ID, names_from = OTU, values_from = Abundance) %>%
    data.frame()
}

# Serum metabolites (SCFA or BA): baseline only, log10-transformed
import_serum_metabolites <- function(file, exclude) {
  read.csv(here("Processed_data", file)) %>%
    filter(Time %in% "Baseline", !Analyte %in% exclude) %>%
    select(Subject_ID, Analyte, Abun) %>%
    mutate(Abun = log10(Abun)) %>%
    pivot_wider(id_cols = Subject_ID, names_from = Analyte, values_from = Abun)
}

data_list <- list(
  SCFA         = import_serum_metabolites("SCFA_serum_melanoma.csv",
                                          exclude = c("Hexanoic acid")),
  BA           = import_serum_metabolites("BA_serum_melanoma.csv",
                                          exclude = c("Dehydrolithocholic Acid",
                                                      "Isodeoxycholic Acid",
                                                      "Taurolithocholic Acid",
                                                      "Alloiso_Isolithocholic Acid")),
  Mic_Phylum   = import_mic_data(phy_glom, "Phylum"),
  Mic_Class    = import_mic_data(phy_glom, "Class"),
  Mic_Order    = import_mic_data(phy_glom, "Order"),
  Mic_Genus    = import_mic_data(phy_glom, "Genus"),
  Mic_Species  = import_mic_data(phy_glom, "Species")
)
saveRDS(data_list, file.path(results_folder, "baseline_data_list.rds"))


# =============================================================================
# 3. VE800 colonization outcome (mean prevalence per subject and strain)
# =============================================================================
colon_file <- here("Results/03_VE800_colonization_summary/VE800_Mean_RA_Subject.csv")
if (!file.exists(colon_file)) stop("Run scripts/03_VE800_colonization_summary.R first: ", colon_file, " not found")
prev_sp_dt <- read.csv(colon_file) %>%
  filter(!organism %in% c("VE800-06")) %>%
  transmute(Subject_ID, mean_abun, mean_prev, Species = organism) %>%
  # Map thirds onto the quarter scale: 2/3 -> 0.5, 1/3 -> 0.25
  # (tolerance-based comparison instead of exact floating-point equality)
  mutate(mean_prev = case_when(near(mean_prev, 2/3, tol = 1e-6) ~ 0.5,
                               near(mean_prev, 1/3, tol = 1e-6) ~ 0.25,
                               TRUE ~ mean_prev))

ve800_sp <- unique(prev_sp_dt$Species)

# Sanity check: colonization categories per strain
print(table(prev_sp_dt$Species, prev_sp_dt$mean_prev))


# =============================================================================
# 4. Helper functions (used by both model fitting and plotting)
# =============================================================================

# Merge classes with fewer than `min_n` observations into a neighbouring class.
# `lev_order` gives the ordinal order of the classes. Repeatedly takes the LAST
# under-sized class and merges it into the class below it (or the class above
# it if it is the first class) until all classes have >= min_n observations.
merge_rare_classes <- function(y, lev_order, min_n = 3) {
  y <- as.character(y)
  repeat {
    tab <- table(factor(y, levels = lev_order))
    tab <- tab[tab > 0]
    if (all(tab >= min_n) || length(tab) < 2) break   # guard against an endless loop
    i      <- max(which(tab < min_n))
    target <- if (i == 1) names(tab)[2] else names(tab)[i - 1]
    y[y == names(tab)[i]] <- target
  }
  y
}

# Join a feature table with one strain's colonization, keep complete cases.
# Returns a data.frame with Subject_ID as rownames and `mean_prev` as last column.
build_model_data <- function(feature_dt, ve_sp, prev_sp_dt) {
  colo_dt <- prev_sp_dt %>%
    filter(Species %in% ve_sp) %>%
    select(Subject_ID, mean_prev)
  dt <- feature_dt %>% left_join(colo_dt, by = "Subject_ID")
  dt <- data.frame(dt[complete.cases(dt), ])
  rownames(dt) <- dt$Subject_ID
  dt$Subject_ID <- NULL
  dt
}

# Fit ordinalNet models for ONE strain across all data sets.
# Returns list(coeff = <non-zero coefficients>, perf = <CV summaries>, errors = <messages>)
fit_species <- function(ve_sp, data_list, prev_sp_dt, min_class_n = 3) {
  out <- list(coeff = list(), perf = list(), errors = list())
  
  for (dt_set in names(data_list)) {
    key <- paste0(ve_sp, "_", dt_set)
    message("Fitting: ", key)
    
    res <- tryCatch({
      t_data <- build_model_data(data_list[[dt_set]], ve_sp, prev_sp_dt)
      
      # Ordered response with rare classes merged
      y_lev   <- as.character(sort(unique(t_data$mean_prev)))
      y_class <- merge_rare_classes(t_data$mean_prev, y_lev, min_class_n)
      y_var   <- factor(y_class, levels = intersect(y_lev, y_class), ordered = TRUE)
      print(table(y_var))
      
      x_mat <- as.matrix(t_data[, setdiff(names(t_data), "mean_prev"), drop = FALSE])
      
      # 3-fold CV over the lambda path (parallel/proportional-odds terms only).
      # printProgress is off because worker output only goes to the log file.
      set.seed(4000)
      fit_cv <- ordinalNetCV(x_mat, y_var, family = "cumulative", link = "logit",
                             nFolds = 3, maxiterOut = 1000, reverse = FALSE,
                             nLambda = 50, lambdaMinRatio = 1e-4,
                             parallelTerms = TRUE, nonparallelTerms = FALSE,
                             printProgress = FALSE)
      
      sum_cv_dt <- data.frame(summary(fit_cv))   # one row per fold
      sum_cv_dt$dataset <- dt_set
      sum_cv_dt$Species <- ve_sp
      
      # Final lambda = mean of the best lambda across the 3 folds
      lambda_opt <- mean(sum_cv_dt$lambda)
      
      # Refit on all subjects at lambda_opt
      set.seed(4096)
      model_fit <- ordinalNet(x_mat, y_var, family = "cumulative", link = "logit",
                              reverse = FALSE, lambdaVals = lambda_opt,
                              parallelTerms = TRUE, nonparallelTerms = FALSE)
      
      # Long table of non-zero coefficients
      coeff_dt <- data.frame(coef(model_fit, matrix = TRUE))
      coeff_dt$Predictors <- rownames(coeff_dt)
      coeff_dt <- coeff_dt %>%
        pivot_longer(cols = !Predictors, names_to = "logit", values_to = "coeff") %>%
        mutate(coeff = as.numeric(coeff)) %>%
        filter(abs(coeff) > 0) %>%
        mutate(dataset = dt_set, Species = ve_sp)
      
      list(coeff = coeff_dt, perf = sum_cv_dt)
    }, error = function(e) e)
    
    if (inherits(res, "error")) {
      out$errors[[key]] <- conditionMessage(res)
    } else {
      out$coeff[[key]] <- res$coeff
      out$perf[[key]]  <- res$perf
    }
  }
  out
}


# =============================================================================
# 5. Fit models in parallel across VE800 strains
# =============================================================================
# Seeds are set inside each fit, so results are identical to a sequential run.
# To run sequentially instead, replace the cluster setup with registerDoSEQ().
n_cores <- max(1, min(length(ve800_sp), parallel::detectCores() - 1, na.rm = TRUE))
log_file <- file.path(results_folder, "ordinalNet_parallel.log")
cl <- makeCluster(n_cores, outfile = log_file)   # worker messages go to the log
registerDoParallel(cl)
clusterExport(cl, c("merge_rare_classes", "build_model_data"))

fit_results <- foreach(ve_sp = ve800_sp,
                       .packages = c("ordinalNet", "dplyr", "tidyr"),
                       .inorder = TRUE) %dopar% {
                         fit_species(ve_sp, data_list, prev_sp_dt, min_class_n)
                       }

stopCluster(cl)
registerDoSEQ()

# Combine per-strain results (names: "<strain>_<dataset>")
coeff_list <- do.call(c, lapply(fit_results, `[[`, "coeff"))
perf_list  <- do.call(c, lapply(fit_results, `[[`, "perf"))
fit_errors <- do.call(c, lapply(fit_results, `[[`, "errors"))

if (length(fit_errors) > 0) {
  warning(length(fit_errors), " model(s) failed; see `fit_errors` and ", log_file)
  print(fit_errors)
}

saveRDS(coeff_list, file.path(results_folder, "coeff_list.rds"))
saveRDS(perf_list,  file.path(results_folder, "perf_list.rds"))


# =============================================================================
# 6. Heatmap of selected predictors (direction of association)
# =============================================================================
# With parallel terms only, each predictor has the same coefficient for every
# cumulative logit, so dropping the `logit` column and de-duplicating leaves
# one row per predictor / data set / strain.
sum_coeff_dt <- bind_rows(coeff_list) %>%
  filter(!Predictors %in% "(Intercept)", abs(coeff) > 1e-6) %>%
  distinct(Predictors, coeff, dataset, Species)

# reverse = FALSE models P(Y <= j): a positive coefficient shifts subjects towards
# LOWER colonization, so the sign is flipped (+1 = positive association with
# colonization, -1 = negative, 0 = not selected).
imp_dt_w <- sum_coeff_dt %>%
  mutate(coeff = ifelse(coeff > 0, -1, 1)) %>%
  pivot_wider(id_cols = c(Predictors, dataset), names_from = Species,
              values_from = coeff, values_fill = 0)

mat_logic <- imp_dt_w %>% select(-Predictors, -dataset) %>% as.matrix()
rownames(mat_logic) <- imp_dt_w$Predictors
mat_logic <- mat_logic[, mixedsort(colnames(mat_logic)), drop = FALSE]

split_rows <- factor(imp_dt_w$dataset,
                     levels = c("Mic_Phylum", "Mic_Class", "Mic_Order", "Mic_Genus",
                                "Mic_Species", "SCFA", "BA"),
                     labels = c("Phylum", "Class", "Order", "Genus",
                                "Species", "SCFA", "BA"))

# Readable row names: drop level prefixes and make.names() artefacts
clean_feature_names <- function(x) {
  x <- gsub("^(Genus|Species|Phylum|Class|Order)_", "", x)
  x <- gsub("\\.\\.", " ", x)
  x <- gsub("\\.", " ", x)
  gsub("^X ", "", x)
}
rownames(mat_logic) <- clean_feature_names(rownames(mat_logic))

# Append the Class of each VE800 strain to its column name
# Copied into Processed_data/ by 01_prepare_processed_data.R (falls back to the
# original location in data/ if you have the raw data locally)
tax_file <- here("Processed_data/VE800_taxonomy.csv")
if (!file.exists(tax_file)) tax_file <- here("data/VE800_taxonomy.csv")
ve800_tax_dt <- read.csv(tax_file)
lin_tax_dt <- str_split_fixed(ve800_tax_dt$Lineage, pattern = ";", n = 7) %>% data.frame()
colnames(lin_tax_dt) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
lin_tax_dt[] <- lapply(lin_tax_dt, function(x) gsub("^.{0,4}", "", x))   # strip rank prefix
ve800_tax_dt <- cbind(ve800_tax_dt, lin_tax_dt)

strain_class <- ve800_tax_dt$Class[match(colnames(mat_logic), ve800_tax_dt$StrainID)]
colnames(mat_logic) <- paste0(colnames(mat_logic), " (", strain_class, ")")

# Highlight VE800-01 to -04 in red
fontcolors <- ifelse(grepl("VE800-01|VE800-02|VE800-03|VE800-04", colnames(mat_logic)),
                     "red", "black")

col_mat <- colorRamp2(c(-1, 0, 1), c("#ff4c4c", "white", "#046C9A"))

ht1 <- Heatmap(mat_logic, name = "Coeff", column_title = NA,
               row_split = split_rows,
               row_gap = unit(2, "mm"),
               border = TRUE,
               row_title_rot = 0,
               column_title_rot = 90,
               column_names_rot = 90,
               column_names_side = "top",
               row_title_gp = gpar(fontface = "bold"),
               column_names_gp = gpar(fontsize = 12, col = fontcolors),
               cluster_rows = FALSE,
               cluster_row_slices = TRUE,
               show_row_dend = FALSE,
               row_names_side = "left", color_space = "LAB",
               col = col_mat,
               cluster_columns = TRUE,
               width = 2,
               cell_fun = function(j, i, x, y, w, h, col) {   # cell borders
                 grid.rect(x = x, y = y, width = w, height = h,
                           gp = gpar(col = "black", fill = NA))
               },
               row_names_max_width = max_text_width(rownames(mat_logic),
                                                    gp = gpar(fontsize = 12)),
               column_names_max_height = max_text_width(colnames(mat_logic),
                                                        gp = gpar(fontsize = 12)),
               show_heatmap_legend = FALSE,
               na_col = "white")

lgd <- Legend(at = c(-1, 0, 1), title = "Relation",
              labels = c("Negative", "NS", "Positive"),
              legend_gp = gpar(fill = c("#ff4c4c", "white", "#046C9A")))

cairo_pdf(file.path(results_folder, "Heatmap_Summary_Ordinal_net.pdf"),
          height = 13, width = 7, fallback_resolution = 600)
draw(ht1, annotation_legend_list = list(lgd))
dev.off()


# =============================================================================
# 7. Box plots of selected predictors by colonization class (per strain)
# =============================================================================
prev_labels <- c("0" = "No", "0.25" = "Low", "0.5" = "Partial",
                 "0.75" = "Moderate", "1" = "Complete")
lev_class <- unname(prev_labels)

for (sp in unique(sum_coeff_dt$Species)) {
  
  sp_folder <- file.path(results_folder, sp)
  dir.create(sp_folder, showWarnings = FALSE, recursive = TRUE)
  sp_dt <- sum_coeff_dt %>% filter(Species %in% sp)
  
  for (dt_set in unique(sp_dt$dataset)) {
    message("Plotting: ", sp, " / ", dt_set)
    
    t_data   <- build_model_data(data_list[[dt_set]], sp, prev_sp_dt)
    pred_vec <- sp_dt$Predictors[sp_dt$dataset %in% dt_set]
    sel_dt   <- t_data[, names(t_data) %in% c(pred_vec, "mean_prev"), drop = FALSE]
    
    # Label colonization classes and merge rare ones (same rule as the models)
    sel_dt$mean_prev <- merge_rare_classes(unname(prev_labels[as.character(sel_dt$mean_prev)]),
                                           lev_class, min_class_n)
    
    bar_dt <- sel_dt %>%
      pivot_longer(cols = -mean_prev) %>%
      mutate(mean_prev = factor(mean_prev, levels = lev_class))
    
    dt_set_lab <- gsub("_", " ", gsub("^Mic_", "", dt_set))
    
    pdf(file.path(sp_folder, paste0("Barplot_", sp, "_", dt_set, ".pdf")),
        height = 4, width = 5)
    for (var_name in unique(bar_dt$name)) {
      bar_p <- ggplot(bar_dt %>% filter(name %in% var_name),
                      aes(x = mean_prev, y = value, fill = mean_prev)) +
        geom_boxplot(outlier.colour = NA, alpha = 0.6) +
        geom_point(colour = "black", alpha = 0.5, size = 2,
                   position = position_jitter(w = 0.1, h = 0),
                   shape = 21, stroke = 1) +
        facet_wrap(~name, scales = "free", ncol = 4) +
        scale_fill_manual(name = "Colonization", values = wes_palette("IsleofDogs2")) +
        ylab("Baseline Abundance/Proportion") +
        xlab(paste0(sp, "-Colonization")) +
        theme_prism() +
        theme(axis.text.x  = element_text(size = 10, face = "bold", angle = 60, hjust = 1, vjust = 1),
              axis.text.y  = element_text(size = 12, face = "bold"),
              axis.title.y = element_text(size = 12, face = "bold"),
              axis.title.x = element_text(size = 12, face = "bold"),
              strip.text.x = element_text(size = 10, face = "bold"),
              legend.title = element_text(size = 10),
              legend.text  = element_text(size = 10)) +
        ggtitle(dt_set_lab)
      print(bar_p)
    }
    dev.off()
  }
}


