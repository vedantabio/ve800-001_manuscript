# =============================================================================
# run_all.R - run the full pipeline in order
# -----------------------------------------------------------------------------
# Usage (from the project root, or open VE800_analysis.Rproj in RStudio):
#     Rscript run_all.R
#
# Step 1 (01_prepare_processed_data.R) needs the raw data in data/, which is
# NOT in the repository. It is skipped automatically when data/ is missing,
# and the analyses then run from the files already in Processed_data/.
# =============================================================================

library(here)

prepare_script <- "01_prepare_processed_data.R"

analysis_scripts <- c(
  "02_microbiome_diversity_and_taxa.R",
  "03_VE800_colonization_summary.R",
  "04_serum_SCFA_analysis.R",
  "05_serum_bile_acid_analysis.R",
  "06_serum_cytokine_analysis.R",
  "07_colonization_predictors_ordinalNet.R"   # needs output of 03
)

run_script <- function(f) {
  message("\n========== ", f, " ==========")
  t0 <- Sys.time()
  # Each script runs in its own environment so objects don't leak between them
  source(here("scripts", f), local = new.env(), echo = FALSE)
  message("Done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}

if (dir.exists(here("data"))) {
  run_script(prepare_script)
} else {
  message("data/ not found - skipping ", prepare_script, " and using Processed_data/")
}

invisible(lapply(analysis_scripts, run_script))
message("\nAll done. Figures and tables are in Results/")
