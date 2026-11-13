# =============================================================================
# 01_prepare_processed_data.R
# -----------------------------------------------------------------------------
# Builds EVERY file in Processed_data/ from the raw data in data/.
# This is the only script that needs the raw data; all analysis scripts
# (02-07) read from Processed_data/ only.
#
#   Part | Raw input                              | Output (Processed_data/)
#   -----+----------------------------------------+------------------------------------
#    A   | One Codex species profiles             | phy_mic.rds, sample_metadata.csv
#    B   | VE800 targeted marker panel            | VE303_abun.csv
#    C   | Serum SCFA + bile acids ( xlsx)  | SCFA_serum_melanoma.csv,
#        |                                        | BA_serum_melanoma.csv
#    D   | Serum cytokine panel ( xlsx)     | serum_cytokines_melanoma.csv
#    E   | VE800 strain taxonomy (reference)      | VE800_taxonomy.csv (copied)
#
# All parts share the same clinical metadata: Week-8 response (SD / PD) and
# VE800 engraftment (HE / LE) per subject.
#
# =============================================================================

source(here::here("R/utils.R"))
suppressPackageStartupMessages({
  library(readxl)
  library(phyloseq)
  library(stringr)
  library(zoo)      # na.locf0(): carry the last non-NA header forward
  library(gtools)   # mixedsort()
})

steps <- c("microbiome", "marker_panel", "serum_metabolites", "cytokines", "taxonomy")

out_dir <- create_results_dir("Processed_data")

# ---- Raw input files  -----------------------------------
raw_files <- list(
  engraftment     = "data/Engraftment_analysis_metadata/2024-02-26 5 Engraftment categories.xlsx",
  clinical        = "data/Microbiome data/VE800_clinical_metadata_complete_clean.xlsx",
  mgs_metadata    = "data/Microbiome data/2024-04-26 Metadata_joined_old_new.csv",
  one_codex       = "data/Microbiome data/One_Codex_Results_VE800_Ph1_20220124.csv",
  marker_panel    = "data/Microbiome data/VE800_Ph1_Extended_Marker_Data_20220202.csv",
  serum_scfa      = "data/Input_Metabolomics_data/SCFA/VEDA-0102-21TASA_Serum_CDT.xlsx",
  serum_ba        = "data/Input_Metabolomics_data/Bile_acids/VEDA-0101-21TASA_Results_Serum.xlsx",
  serum_cytokines = "data/VE800_cytokines_2021_02_26_Response_int.xlsx",
  ve800_taxonomy  = "data/VE800_taxonomy.csv"
)


# =============================================================================
# Shared helpers: clinical metadata and visit names
# =============================================================================

#' Engraftment categories per subject. Subject_ID numeric; "Colorectal" -> "CRC".
load_engraftment <- function() {
  eng_dt <- read_excel(here(raw_files$engraftment))
  names(eng_dt)[1] <- "Subject_ID"
  eng_dt$Subject_ID <- as.numeric(eng_dt$Subject_ID)
  eng_dt$Cohort[eng_dt$Cohort == "Colorectal"] <- "CRC"
  eng_dt
}

#' Clinical metadata (all cohorts) joined with engraftment
#' (joined on all shared columns: Subject_ID and Cohort).
load_clinical_metadata <- function() {
  read_excel(here(raw_files$clinical)) %>%
    left_join(load_engraftment())
}

#' Subject-level outcomes for the serum data:
#' Subject_ID (character), Response_Wk8 (PD / SD), Engraftment (HE / LE), Cohort.
load_outcome_metadata <- function() {
  meta_dt <- load_clinical_metadata() %>%
    select(Subject_ID, Response_Wk8, Engraftment, Cohort)
  meta_dt$Response_Wk8 <- gsub("PD/Off TX", "PD", meta_dt$Response_Wk8)   # treated as PD
  # Alphabetical relabelling: "High Engraft" -> "HE", "Low Engraft" -> "LE"
  meta_dt$Engraftment <- factor(meta_dt$Engraftment)
  levels(meta_dt$Engraftment) <- c("HE", "LE")
  meta_dt$Subject_ID <- as.character(meta_dt$Subject_ID)
  meta_dt
}

#' Metagenomic (MGS) sample metadata joined with clinical + engraftment data.
#' One row per MGS sample; rownames are the MGS sample IDs.
load_mgs_sample_metadata <- function() {
  meta_dt <- read.csv(here(raw_files$mgs_metadata))
  names(meta_dt)[1] <- "Subject_ID"
  meta_dt <- meta_dt %>% inner_join(load_clinical_metadata())
  rownames(meta_dt) <- meta_dt$MGS.Sample.ID
  meta_dt
}

#' Clinical visit codes -> study days, so serum time points line up with the
#' microbiome: C1D1 -> D1, C1D14 -> D14, C2D1 -> D30, C3D1 -> D58, C4D1 -> D86.
#' (Order matters: "C1" is stripped first.)
recode_visit_to_day <- function(visit) {
  visit <- gsub("C1", "", visit)
  visit <- gsub("C2D1", "D30", visit)
  visit <- gsub("C3D1", "D58", visit)
  gsub("C4D1", "D86", visit)
}


# =============================================================================
# A. Microbiome: One Codex species profiles -> phyloseq object
# =============================================================================

#' Replace NA at each rank with "uncl_<Rank>_of_<name of the rank above>".
#' Works left-to-right so the parent is already filled when a child is set.
fill_unclassified_taxa <- function(tax_matrix) {
  tax_matrix <- as.matrix(tax_matrix)
  for (i in 2:ncol(tax_matrix)) {
    for (j in which(is.na(tax_matrix[, i]))) {
      for (k in i:ncol(tax_matrix)) {
        if (is.na(tax_matrix[j, k])) {
          tax_matrix[j, k] <- paste0("uncl_", colnames(tax_matrix)[k], "_of_",
                                     tax_matrix[j, i - 1])
        }
      }
    }
  }
  as.data.frame(tax_matrix)
}

prepare_microbiome <- function() {
  meta_dt <- load_mgs_sample_metadata()
  
  # Species-level rows with non-zero abundance, for samples that have metadata
  one_codex_data <- read.csv(here(raw_files$one_codex)) %>%
    filter(sample_name %in% meta_dt$MGS.Sample.ID,
           tax_rank %in% "species",
           !is.na(abundance),
           abundance != 0) %>%
    distinct()
  
  # Abundance table: species x samples
  rel_dt <- one_codex_data %>%
    select(sample_name, species_name, abundance) %>%
    pivot_wider(names_from = sample_name, values_from = abundance, values_fill = 0) %>%
    as.data.frame()
  rel_dt$species_name <- make.names(rel_dt$species_name)
  rownames(rel_dt) <- rel_dt$species_name
  rel_dt$species_name <- NULL
  message("One Codex: ", ncol(rel_dt), " samples, ", nrow(rel_dt), " species")
  
  # Taxonomy table
  tax_dt <- one_codex_data[, c("superkingdom_name", "phylum_name", "class_name",
                               "order_name", "family_name", "genus_name", "species_name")]
  names(tax_dt) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_dt <- unique(tax_dt)
  # Absiella sp. MP454B12 is reported under two orders; keep a single entry
  tax_dt <- tax_dt[-which(tax_dt$Species == "Absiella sp. MP454B12" &
                            tax_dt$Order == "Clostridiales"), ]
  rownames(tax_dt) <- make.names(as.character(tax_dt$Species))
  
  # phyloseq object without empty taxa / samples
  phy <- phyloseq(otu_table(rel_dt, taxa_are_rows = TRUE),
                  tax_table(as.matrix(tax_dt)),
                  sample_data(meta_dt))
  phy <- prune_taxa(taxa_sums(phy) > 0, phy)
  phy <- prune_samples(sample_sums(phy) > 0, phy)
  
  # Clean taxonomy (Family dropped: mostly empty in the One Codex export)
  tax_mat <- tax_table(phy)@.Data
  tax_mat <- tax_mat[, c("Kingdom", "Phylum", "Class", "Order", "Genus", "Species")]
  tax_mat[tax_mat == ""] <- NA
  tax_clean <- fill_unclassified_taxa(tax_mat)
  tax_clean$Order[grep("Candidatus Gastranaerophilales", tax_clean$Order)] <-
    "uncl_Order_of_Candidatus Melainabacteria"
  tax_table(phy) <- tax_table(as.matrix(tax_clean))
  
  saveRDS(phy, file.path(out_dir, "phy_mic.rds"))
  write.csv(data.frame(sample_data(phy)),
            file = file.path(out_dir, "sample_metadata.csv"), row.names = FALSE)
}


# =============================================================================
# B. VE800 marker panel (strain detection / relative abundance per sample)
# =============================================================================
prepare_marker_panel <- function() {
  meta_dt <- load_mgs_sample_metadata()
  
  m_panel_dt <- read.csv(here(raw_files$marker_panel)) %>%
    mutate(sample_name = gsub("AB.fastq.gz", "", filename)) %>%   # "<sample>AB.fastq.gz"
    filter(sample_name %in% meta_dt$MGS.Sample.ID) %>%
    inner_join(meta_dt, by = c("sample_name" = "MGS.Sample.ID")) %>%
    filter(!is.na(Subject_ID))
  
  # File name is historical; the file contains VE800 data
  write.csv(m_panel_dt, file.path(out_dir, "VE303_abun.csv"), row.names = FALSE)
}


# =============================================================================
# C. Serum metabolomics: short-chain fatty acids and bile acids
# =============================================================================

#' Read one vendor metabolomics export into a tidy long table:
#' Subject_ID, Analyte, Abun, Time, Unit (+ outcome metadata).
#' "BLOQ" (below the lower limit of quantification) -> that analyte's LLOQ.
#' Unscheduled visits are removed.
read_serum_metabolite <- function(path, outcome_meta) {
  dt <- read_xlsx(here(path))
  dt$Abun <- suppressWarnings(as.numeric(as.character(dt$Result)))
  idx_bloq <- grep("BLOQ", dt$Result)
  dt$Abun[idx_bloq] <- dt$LLOQ[idx_bloq]
  
  dt %>%
    mutate(Time       = `Visit Name`,
           Subject_ID = `Subject or Animal ID`) %>%
    filter(!Time %in% c("UNSCHED")) %>%
    mutate(Time = recode_visit_to_day(Time)) %>%
    select(Subject_ID, Analyte, Abun, Time, Unit) %>%
    left_join(outcome_meta)
}

prepare_serum_metabolites <- function() {
  meta_dt <- load_outcome_metadata()
  # write.csv keeps row names, so these files have a leading "X" column
  write.csv(read_serum_metabolite(raw_files$serum_scfa, meta_dt),
            file.path(out_dir, "SCFA_serum_melanoma.csv"))
  write.csv(read_serum_metabolite(raw_files$serum_ba, meta_dt),
            file.path(out_dir, "BA_serum_melanoma.csv"))
}


# =============================================================================
# D. Serum cytokines (multiplex panel, 3 plates)
# -----------------------------------------------------------------------------
# The  file has a two-row header (cytokine / measure) and per-plate LLOQ
# rows mixed in with the sample rows.
#   1. Build a single header "<Cytokine>_<Measure>"
#   2. Split per-plate LLOQ rows from sample rows
#   3. Per cytokine keep concentration (Conc), within-range flag (LOQ), MFI flag
#   4. Long format; below detection (NA) -> plate LLOQ
#   5. Clinical metadata, visit recoding
#   6. log2 fold change vs baseline; drop measurements with MFI >= 20
# =============================================================================
prepare_cytokines <- function() {
  path <- here(raw_files$serum_cytokines)
  
  # ---- 1. Two-row header ----
  # Row 1: cytokine names (merged cells -> "...N" -> NA -> filled forward)
  head1 <- read_excel(path, col_names = TRUE) %>%
    names() %>%
    str_replace("\\.\\.\\.", NA_character_)
  head1 <- zoo::na.locf0(head1)
  # Row 2: measure names (strip "...N" suffix)
  head2 <- read_excel(path, skip = 1, col_names = TRUE) %>%
    names() %>%
    str_replace("\\.\\.\\..*", "")
  h1      <- head1[seq_along(head2)]
  headers <- ifelse(!is.na(h1) & !is.na(head2), paste(h1, head2, sep = "_"), head2)
  
  raw_data <- read_excel(path, skip = 1, col_names = TRUE)
  colnames(raw_data) <- headers
  
  # ---- 2. LLOQ rows (Plate column contains "Plate ...") vs sample rows ----
  idx_lloq <- grep("Plate", raw_data$Plate)
  lloq_dt  <- raw_data[idx_lloq, ]
  cyto_dt  <- raw_data[-idx_lloq, ]
  
  lloq_dt <- lloq_dt %>% select(where(~ !any(is.na(.x))))
  colnames(lloq_dt) <- gsub("_.*", "", colnames(lloq_dt))
  lloq_dt   <- t(lloq_dt)
  lloq_dt   <- lloq_dt[2:nrow(lloq_dt), ]          # drop the "Plate" row
  cytokines <- rownames(lloq_dt)
  lloq_dt   <- data.frame(lloq_dt, stringsAsFactors = FALSE)
  lloq_dt   <- cbind(Cytokines = cytokines, data.frame(sapply(lloq_dt, as.numeric)))
  names(lloq_dt)[2:4] <- paste0("Plate ", 1:3)
  lloq_dt <- lloq_dt %>%
    pivot_longer(cols = starts_with("Plate"), names_to = "Plate", values_to = "LLOQ")
  
  # ---- 3. Columns belonging to each cytokine (cytokines start at column 10) ----
  all_cyto <- unique(head1[10:length(head1)])
  message("Cytokines in panel: ", length(all_cyto))
  all_cyto <- gsub("IL-12 \\(", "", all_cyto)    # avoid regex special characters
  all_cyto <- gsub("\\)", "", all_cyto)
  
  sel_cyto <- unique(unlist(lapply(all_cyto, grep, x = names(cyto_dt), value = TRUE)))
  sel_dt   <- cyto_dt[, c(names(cyto_dt)[1:9], sel_cyto)]
  names(sel_dt) <- gsub("_Conc.*",     "_Conc", names(sel_dt))
  names(sel_dt) <- gsub("_Within.*",   "_LOQ",  names(sel_dt))
  names(sel_dt) <- gsub(" Concerns.*", "",      names(sel_dt))
  
  sel_dt[, 10:ncol(sel_dt)] <- sapply(sel_dt[, 10:ncol(sel_dt)], as.numeric)
  sel_dt <- sel_dt[, colSums(is.na(sel_dt)) < nrow(sel_dt)]   # drop all-NA columns
  
  # Add empty MFI columns so every cytokine has Conc / LOQ / MFI
  LOQ     <- grep("LOQ", names(sel_dt), value = TRUE)
  MFI     <- gsub("_LOQ", "_MFI", LOQ)
  new_MFI <- MFI[!MFI %in% names(sel_dt)]
  mfi_dt  <- data.frame(matrix(NA, nrow = nrow(sel_dt), ncol = length(new_MFI)))
  names(mfi_dt) <- new_MFI
  sel_dt <- cbind(sel_dt, mfi_dt)
  sel_dt <- cbind(sel_dt[, 1:9], sel_dt[, mixedsort(names(sel_dt)[10:ncol(sel_dt)])])
  sel_dt[, 10:ncol(sel_dt)] <- sapply(sel_dt[, 10:ncol(sel_dt)], as.numeric)
  sel_dt$Plate <- paste0("Plate ", sel_dt$Plate)
  
  # % of samples per plate within the quantifiable range
  lod_prev_dt <- sel_dt[, c("Plate", grep("LOQ", names(sel_dt), value = TRUE))]
  plate_cnts  <- lod_prev_dt %>% count(Plate)
  lod_prev    <- lod_prev_dt %>%
    group_by(Plate) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))
  prop_mat <- 100 * sweep(as.matrix(t(lod_prev[, 2:ncol(lod_prev)])), 2, plate_cnts$n, FUN = "/")
  prop_dt  <- data.frame(Cytokines = rownames(prop_mat), prop_mat)
  prop_dt$Cytokines <- gsub("_LOQ", "", prop_dt$Cytokines)
  names(prop_dt)[2:4] <- paste0("Plate ", 1:3)
  prop_dt <- pivot_longer(prop_dt, cols = paste0("Plate ", 1:3),
                          values_to = "Prop", names_to = "Plate")
  
  # ---- 4. Long format: one row per sample x cytokine ----
  cyto_vec <- grep("Conc", names(sel_dt), value = TRUE)
  LOQ      <- grep("LOQ",  names(sel_dt), value = TRUE)
  MFI      <- grep("MFI",  names(sel_dt), value = TRUE)
  to_long <- function(keep, drop, value_name, suffix) {
    sel_dt %>%
      select(-all_of(drop)) %>%
      pivot_longer(cols = all_of(keep), names_to = "Cytokines", values_to = value_name) %>%
      mutate(Cytokines = gsub(suffix, "", Cytokines))
  }
  sel_dt_final <- to_long(cyto_vec, c(LOQ, MFI),      "Abundance", "_Conc") %>%
    inner_join(to_long(LOQ,      c(cyto_vec, MFI), "LOQ",       "_LOQ")) %>%
    inner_join(to_long(MFI,      c(cyto_vec, LOQ), "MFI",       "_MFI")) %>%
    inner_join(lloq_dt) %>%
    inner_join(prop_dt)
  
  sel_dt_final$Abundance <- as.numeric(sel_dt_final$Abundance)
  sel_dt_final$MFI <- 100 * sel_dt_final$MFI
  # NOTE: positional renaming - depends on the column order of the vendor file
  names(sel_dt_final)[4:9] <- c("Time", "PatientID", "Patient_sel",
                                "Wk8_resp", "Wk16_resp", "Time_wks")
  sel_dt_final$MFI[is.na(sel_dt_final$MFI)] <- 0                  # missing MFI -> 0
  idx_na <- is.na(sel_dt_final$Abundance)                         # below detection -> LLOQ
  sel_dt_final$Abundance[idx_na] <- sel_dt_final$LLOQ[idx_na]
  names(sel_dt_final)[5] <- "Subject_ID"
  
  # ---- 5-6. Metadata, visit recoding, fold change, QC filter ----
  cyt_serum_dt <- sel_dt_final %>%
    left_join(load_outcome_metadata()) %>%
    rename(Visit = Time) %>%
    filter(!Visit %in% c("UNSCHED")) %>%
    mutate(Time = recode_visit_to_day(Visit)) %>%
    select(Plate, Subject_ID, Time, Cytokines, Abundance, LLOQ,
           Response_Wk8, Engraftment, Cohort, MFI) %>%
    group_by(Subject_ID, Cytokines) %>%
    mutate(Abundance_frac = log2(Abundance / Abundance[Time == "Baseline"])) %>%
    filter(MFI < 20)
  
  write.csv(cyt_serum_dt, file.path(out_dir, "serum_cytokines_melanoma.csv"))
}


# =============================================================================
# E. VE800 strain taxonomy (reference file used by the ordinal-net analysis)
# =============================================================================
prepare_taxonomy <- function() {
  file.copy(here(raw_files$ve800_taxonomy), file.path(out_dir, "VE800_taxonomy.csv"),
            overwrite = TRUE)
}


# =============================================================================
# Run
# =============================================================================
step_fun <- list(microbiome        = prepare_microbiome,
                 marker_panel      = prepare_marker_panel,
                 serum_metabolites = prepare_serum_metabolites,
                 cytokines         = prepare_cytokines,
                 taxonomy          = prepare_taxonomy)

for (s in steps) {
  message("---- Preparing: ", s)
  step_fun[[s]]()
}
message("Processed data written to ", out_dir)
