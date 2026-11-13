# =============================================================================
# install_packages.R - one-off installation of every package the pipeline uses
# =============================================================================

cran_pkgs <- c(
  "here", "readxl", "tidyverse", "dplyr", "tidyr", "ggplot2", "purrr", "stringr", "zoo",
  "gtools", "vegan", "ggprism", "lemon", "nlme", "rstatix", "emmeans",
  "circlize", "scales", "ggrepel", "viridis", "ordinalNet", "doParallel",
  "Cairo", "wesanderson", "remotes", "BiocManager"
)
missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing)

# Bioconductor
bioc_pkgs <- c("phyloseq", "ComplexHeatmap")
missing <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(missing)) BiocManager::install(missing, update = FALSE)

# GitHub (get.samp(), used by 02_microbiome_diversity_and_taxa.R)
if (!requireNamespace("yingtools2", quietly = TRUE)) remotes::install_github("ying14/yingtools2")
