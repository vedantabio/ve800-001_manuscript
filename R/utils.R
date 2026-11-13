# =============================================================================
# R/utils.R
# -----------------------------------------------------------------------------
# Shared helper functions sourced by every script:
#
#     source(here::here("R/utils.R"))
#
# Paths are resolved with here::here(), so scripts can be run from any working
# directory as long as they live inside this project (the .Rproj file marks
# the project root).
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# =============================================================================
# Folders
# =============================================================================

#' Create (if needed) a folder under the project root and return its absolute
#' path, e.g. create_results_dir("Results", "Timeline").
create_results_dir <- function(...) {
  path <- here::here(...)
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}


# =============================================================================
# Plotting helpers
# =============================================================================

# Colour palettes used throughout
col_resp <- c("PD" = "#ff5e1e", "SD" = "#0099b0")   # Week-8 response
col_eng  <- c("HE" = "#009031", "LE" = "#7f1f98")   # Engraftment (HE/LE)

#' Bold axis text/titles used on top of ggprism::theme_prism() in most plots.
#' @param size font size
#' @param angle x-axis label angle
theme_bold_axes <- function(size = 14, angle = 45) {
  theme(axis.text.x  = element_text(size = size, angle = angle, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y  = element_text(size = size, face = "bold"),
        axis.title.y = element_text(size = size, face = "bold"),
        axis.title.x = element_text(size = size, face = "bold"))
}

#' Save a ggplot (or anything with a print method) to PDF.
save_pdf <- function(plot, file, height, width, ...) {
  pdf(file, height = height, width = width, ...)
  print(plot)
  invisible(dev.off())
}
