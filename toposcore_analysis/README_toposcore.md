# VE800-001 TOPOSCORE Analysis

Calculates the TOPOSCORE for each VE800-001 metagenomic (MGS) stool sample and
models how it changes over time by cohort and overall clinical response.

The TOPOSCORE comes from:

> Derosa L, Iebba V, Silva C ...
> Custom scoring based on ecological topology of gut microbiota associated with
> cancer immunotherapy outcome. Cell, 187, 3373-3389.e16.
> <https://www.cell.com/cell/fulltext/S0092-8674(24)00538-5>

## Contents

- [Quick start](#quick-start)
- [Repository layout](#repository-layout)
- [Inputs](#inputs)
- [Workflow](#workflow)
- [TOPOSCORE calculation](#toposcore-calculation)
- [Timepoints](#timepoints)
- [Statistical analysis](#statistical-analysis)
- [Outputs](#outputs)
- [Notes](#notes)

## Quick start

Open `ve800-001_toposcore_ms_only.rmd` (project root) in RStudio and click
**Knit**, or run from the project root:

```sh
Rscript -e 'rmarkdown::render("ve800-001_toposcore_ms_only.rmd")'
```

**Required packages:** `tidyverse`, `readxl`, `knitr`, `here`, `rmarkdown`,
`lme4`, `lmerTest`, `emmeans`, `downloadthis`

Paths are resolved with `here()` from the project root.

## Repository layout

```
ve800-001_manuscript/
├── ve800-001_toposcore_ms_only.rmd    # analysis report (source)
├── ve800-001_toposcore_ms_only.html   # knitted report
└── toposcore_analysis/
    ├── README_toposcore.md
    ├── input/                         # input data (not modified)
    └── output/                        # generated tables and figures
```

## Inputs

Files in `toposcore_analysis/input/` read by the report:

| File | Description |
|------|-------------|
| `VE800_metaphlan_metaphlan-db_combined_reports.xlsx` | MetaPhlAn relative abundances, one column per sample. Sheets used: `species matches`, `possible matches`, and `unmatched species`. |
| `2026-10-07_VE800_metadata_for_paper.csv` | One row per MGS sample: subject, visit, timepoint, cohort, and clinical response. |

The toposcore species (paper supplemental table `mmc2_toposcore_species.xlsx`)
were matched to MetaPhlAn species names in Excel before import. Possible
matches are species whose names likely changed between MetaPhlAn database
versions.

### MetaPhlAn profiling

Metagenomes were profiled with the nf-core/taxprofiler workflow on Seqera:

- Seqera workflow ID: `5XETXkGHwD4tHu`
- Taxprofiler [v2.0.1](https://github.com/nf-core/taxprofiler/commit/70ecc15e49b4f1fcf79d876643b5d14b65c66178)
- MetaPhlAn 4.1.1
- Database `s3://vedantabio/databases/metaphlan/mpa_vJan26_CHOCOPhlAnSGB/`
  ([SegataLab FTP](https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan26_CHOCOPhlAnSGB_202605.tar))

## Workflow

1. **Import species abundances.** Read the match sheets, pivot sample columns
   to long format, and derive `MGS.Sample.ID` from each report filename.
2. **Join clinical metadata** by `MGS.Sample.ID` and apply exclusions:
   - MetaPhlAn files with no matching sample in the metadata
   - samples with a missing visit name
   - end-of-treatment (EOT) samples
3. **Summarize inclusion.** Subjects and metagenome files by cohort and
   response; toposcore species by SIG and match type.
4. **Compute the TOPOSCORE** for each sample (`calc_toposcores`).
5. **Summarize** TOPOSCORE over time by cohort and response (plot and table).
6. **Model** TOPOSCORE with linear mixed effects models (see
   [Statistical analysis](#statistical-analysis)).
7. **Export** tables and figures to `output/`. The knitted report has a
   download button above each figure.

Handled in the metadata file before import:

- Off-treatment (Off TX) subjects are recoded to `PD`.
- Subjects with `Overall_Response` of `NE` or `NE Restart` are not included.

## TOPOSCORE calculation

```
TOPOSCORE = (#SIG2/46 − #SIG1/38 + 1) / 2
```

- **SIG1:** 37 species associated with poor immune checkpoint therapy (ICT)
  response
- **SIG2:** 45 species associated with good ICT response

Species count only when present (relative abundance > 0). Akkermansia
(`Akkermansia SGB 9226`) also counts in every sample, including when absent:

| Akkermansia relative abundance | Counted as |
|---|---|
| > 0% and ≤ 4.799% (normal) | SIG2 |
| > 4.799% (high) | SIG1 |
| 0% (absent) | SIG1 |

Because Akkermansia adds one species to either signature, the paper's
denominators (45 and 37) are increased to 46 and 38, keeping the score
between 0 and 1. A score of 0 means all SIG1 species and no SIG2 species are
counted; 1 means the reverse.

> **Note:** the paper applies the Akkermansia step only to "gray zone"
> patients. Here it is applied to every sample.

## Timepoints

All summaries, plots, and models use `Visit.Name.Norm` ("Timepoint" in
figures and tables), which converts each cycle-day visit to study day:

| Visit.Name | Timepoint |
|---|---|
| Screen | Screen |
| Baseline | Baseline |
| C1D1 | D1 |
| C1D3 | D3 |
| C1D7 | D7 |
| C1D14 | D14 |
| C1D21 | D21 |
| C2D1 | D30 |
| C3D1 | D58 |
| C4D1 | D86 |
| C5D1 | D114 |
| C6D1 | D142 |
| C7D1 | D170 |
| C8D1 | D198 |
| C10D1 | D254 |

## Statistical analysis

Both analyses are linear mixed effects models (`lme4`/`lmerTest`, REML) with a
random intercept per subject. Cohorts are too small for random slopes.
Contrasts use `emmeans` with Satterthwaite df. Samples with a missing cohort
are excluded.

### 1. Change from Baseline, by cohort

| | |
|---|---|
| **Model** | `TOPOSCORE ~ Visit.Name.Norm + (1 \| Subject_ID)`, fit per cohort, Baseline as reference |
| **Data** | Timepoints with ≥ 3 subjects in the cohort (Screen included) |
| **Tests** | Type III F-test of timepoint; each timepoint vs Baseline |
| **Adjustment** | `dunnettx` within each cohort (not across cohorts) |

### 2. Melanoma: timepoint × overall response

| | |
|---|---|
| **Model** | `TOPOSCORE ~ Visit.Name.Norm * Overall_Response + (1 \| Subject_ID)`, most common response as reference |
| **Data** | Melanoma subjects with known response; timepoints where every response group has ≥ 3 subjects |
| **Tests** | Type III F-tests; change from Baseline within each response group (`dunnettx` within group); response difference at each timepoint (unadjusted); difference in change from Baseline between groups (`dunnettx` across timepoints) |

### Limitations

- Melanoma response groups are small. Check the subject counts per cell in
  the report and treat estimates from small cells as descriptive.
- The report flags singular fits (subject variance estimated as zero), where
  the model behaves like an ordinary linear model.

## Outputs

Files written to `toposcore_analysis/output/`:

| File | Description |
|------|-------------|
| `toposcores.csv` | One row per sample: subject, timepoint, cohort, response, SIG1/SIG2 counts, TOPOSCORE |
| `toposcore_summary.csv` | n, mean, SD, median, and IQR by cohort and timepoint |
| `toposcore_over_time_by_cohort_and_response.png` | Subject trajectories and mean ± SE by response, one panel per cohort |
| `lme_toposcore_contrasts_vs_baseline.csv` | Analysis 1: change from Baseline by cohort and timepoint (Dunnett-adjusted) |
| `lme_toposcore_change_from_baseline.png` | Analysis 1: change from Baseline with 95% CI, one panel per cohort |
| `lme_melanoma_response_contrasts.csv` | Analysis 2: all three contrast sets; the `comparison` column identifies each |
| `lme_melanoma_toposcore_by_response.png` | Analysis 2: model-estimated TOPOSCORE by timepoint and response, with observed points |

## Notes

- **Visit levels** are set in `visit_name_levels` and
  `visit_name_norm_levels` (`import-clinical` chunk). If the metadata contains
  a visit not in these lists, knitting stops with an error naming it; add it
  to the matching list.
