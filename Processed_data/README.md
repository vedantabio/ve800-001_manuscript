# Processed_data

Analysis-ready files produced by `scripts/01_prepare_processed_data.R` from the raw data. These are the only inputs the analysis scripts (02–07) need, so the analyses can be rerun without the raw data.

| File | Part of script 01 | Contents |
|---|---|---|
| `phy_mic.rds` | A. Microbiome | phyloseq object: One Codex species relative abundance, cleaned taxonomy, sample + clinical metadata (all cohorts) |
| `sample_metadata.csv` | A. Microbiome | Sample table of `phy_mic.rds` |
| `VE303_abun.csv` | B. Marker panel | VE800 marker-panel detection and relative abundance per sample and strain (file name is historical) |
| `SCFA_serum_melanoma.csv` | C. Serum metabolites | Serum short-chain fatty acids, long format |
| `BA_serum_melanoma.csv` | C. Serum metabolites | Serum bile acids, long format |
| `serum_cytokines_melanoma.csv` | D. Cytokines | Serum cytokines with LLOQ, MFI flag and log2 fold change vs baseline |
| `VE800_taxonomy.csv` | E. Taxonomy (copied) | Lineage of each VE800 strain (reference file, used by `07_colonization_predictors_ordinalNet.R`) |

Serum files share these columns: `Subject_ID`, `Time` (Baseline, D1, D14, D30, D58, D86), `Response_Wk8` (SD / PD), `Engraftment` (HE / LE), `Cohort`.
