```markdown
# Malaria burden and geographic distance jointly predict HbS differentiation in Uganda

## Table of Contents
- [Overview](#overview)
- [Repository structure](#repository-structure)
- [Data description](#data-description)
- [Installation and environment setup](#installation-and-environment-setup)
- [Workflow 1: Data ingestion and subregion summaries](#workflow-1-data-ingestion-and-subregion-summaries)
- [Workflow 2: Differentiation and distance matrices](#workflow-2-differentiation-and-distance-matrices)
- [Workflow 3: Primary joint inference (MMRR)](#workflow-3-primary-joint-inference-mmrr)
- [Workflow 4: Secondary analyses (PCoA and dbRDA)](#workflow-4-secondary-analyses-pcoa-and-dbrda)
- [Workflow 5: Sensitivity analyses (Mantel and Rousset)](#workflow-5-sensitivity-analyses-mantel-and-rousset)
- [Workflow 6: Figures, tables, and report rendering](#workflow-6-figures-tables-and-report-rendering)
- [Running the full pipeline](#running-the-full-pipeline)
- [How to cite](#how-to-cite)
- [License](#license)
- [Contact information](#contact-information)

---

## Overview
This repository contains code, configuration templates, and documentation to reproduce the analysis for:

**Malaria burden and geographic distance jointly predict HbS differentiation in Uganda**

We quantify locus-specific spatial structure at the sickle-cell variant (HbS) in **HBB** using haemoglobin genotype data aggregated to 15 Ugandan subregions. We then test whether pairwise differentiation aligns with malaria-burden dissimilarity (PfPR\(_{2–10}\)) after accounting for geographic separation.

Objectives:
- Summarize subregional genotype composition (HbAA, HbAS, HbSS) and HbS allele frequency.
- Quantify locus-specific differentiation using global and pairwise Weir–Cockerham \(F_{ST}\).
- Construct genetic, geographic, and malaria-burden dissimilarity matrices.
- Estimate independent associations of geographic distance and PfPR dissimilarity with HbS differentiation using multiple matrix regression with randomization (MMRR) (primary inference).
- Provide secondary multivariate summaries (PCoA, dbRDA) and sensitivity analyses (Mantel-family tests, Rousset IBD regression).
- Produce manuscript-ready tables and figures.

---

## Repository structure
```

uganda-hbs-malaria-distance/
├── README.md                    ← this file
├── CITATION.cff                 ← repository citation metadata
├── LICENSE                      ← dual license notice (MIT for code, CC BY 4.0 for non-code content)
├── .gitignore                   ← ignore patterns (raw data, outputs, caches)
│
├── .github/
│   └── workflows/               ← optional CI workflows
│
├── code/
│   ├── config/
│   │   ├── config.yml           ← analysis parameters (permutations, bootstrap, seeds)
│   │   └── paths.example.yml    ← example input paths (copy to paths.yml locally)
│   │
│   ├── pipeline/
│   │   ├── _targets.R           ← targets pipeline definition
│   │   └── globals.R            ← global options, seeds, package loads
│   │
│   ├── functions/               ← reusable functions (analysis library)
│   │   ├── io.R
│   │   ├── validate.R
│   │   ├── hwe_fis.R
│   │   ├── genetics_fst.R
│   │   ├── distances.R
│   │   ├── mantel.R
│   │   ├── mmrr.R
│   │   ├── ordination.R
│   │   ├── rousset.R
│   │   ├── plotting.R
│   │   ├── tables.R
│   │   └── utils.R
│   │
│   └── scripts/                 ← readable “chapter” scripts mirroring Methods
│       ├── 00_setup.R
│       ├── 01_ingest_clean.R
│       ├── 02_subregion_summaries.R
│       ├── 03_hwe_fis.R
│       ├── 04_fst.R
│       ├── 05_distance_matrices.R
│       ├── 06_mantel_partial.R
│       ├── 07_mmrr.R
│       ├── 08_pcoa_dbrda.R
│       ├── 09_rousset_ibd.R
│       ├── 10_figures.R
│       ├── 11_tables.R
│       └── 12_report_render.R
│
├── data/
│   ├── README.md                ← what goes where in data/
│   ├── raw/                     ← input genotype data (not tracked by git)
│   ├── external/                ← PfPR and polygons (not tracked by git by default)
│   └── metadata/
│       ├── data_dictionary.md
│       ├── subregions_lookup.csv
│       └── provenance.md
│
├── outputs/
│   ├── README.md                ← description of outputs
│   ├── figures/
│   ├── tables/
│   ├── models/
│   ├── diagnostics/
│   └── logs/
│
├── paper/                       ← optional manuscript assets
│   ├── README.md
│   ├── references/
│   │   └── references.bib
│   └── supplementary/
│
├── renv/                        ← reproducible R environment
└── renv.lock

````

---

## Data description

### 1) Genotype data (required)
Place genotype data in `data/raw/` (not tracked by git).

Supported formats:
- Individual-level CSV/TSV (preferred), or
- Aggregated subregion counts

**Individual-level minimum required columns**
- `subregion` (raw or canonical label)
- `genotype` (`HbAA`, `HbAS`, `HbSS`)

**Aggregated minimum required columns**
- `subregion`
- `n_HbAA`, `n_HbAS`, `n_HbSS`

Optional but useful:
- `year` (collection year, enables time-aligned PfPR sensitivity checks)
- `district` or `facility_id` (supports diagnostics for within-subregion mixture)

### 2) Subregion polygons (required for mapping, optional for core inference)
Place polygons in `data/external/` (not tracked by git by default).
- Shapefile or GeoJSON
- Must include a subregion name field that matches the 15 canonical subregion names (or can be mapped via `data/metadata/subregions_lookup.csv`)

### 3) Malaria burden (required)
Place malaria burden data in `data/external/` (not tracked by git by default).
- CSV with:
  - `subregion`
  - `pfpr_2_10` (mean PfPR\(_{2–10}\))

If year-specific malaria summaries are available, include:
- `year`
- `pfpr_2_10`

---

## Installation and environment setup

### Option A: Reproducible environment with `renv` (recommended)
From the repository root:
```r
install.packages("renv")
renv::restore()
````

### Option B: Manual install (minimum)

```r
install.packages(c(
  "targets","tarchetypes","yaml",
  "dplyr","readr","tidyr","tibble",
  "sf","geosphere",
  "vegan","ape",
  "ggplot2","patchwork"
))
```

---

## Workflow 1: Data ingestion and subregion summaries

Purpose: Load genotype data, harmonize subregion labels, compute genotype proportions and HbS allele frequency.

Run:

```r
source("code/scripts/01_ingest_clean.R")
source("code/scripts/02_subregion_summaries.R")
```

Main outputs:

* `outputs/tables/subregion_genotypes.csv`
* `outputs/tables/subregion_allele_freq.csv`

---

## Workflow 2: Differentiation and distance matrices

Purpose: Estimate global and pairwise Weir–Cockerham (F_{ST}), then construct (D_{gen}), (D_{geo}), (D_{mal}).

Run:

```r
source("code/scripts/04_fst.R")
source("code/scripts/05_distance_matrices.R")
```

Main outputs:

* `outputs/tables/fst_global.csv`
* `outputs/tables/fst_pairwise_matrix.csv`
* `outputs/models/dist_matrices.rds`

---

## Workflow 3: Primary joint inference (MMRR)

Purpose: Estimate independent associations of geographic distance and malaria-burden dissimilarity with genetic differentiation.

Run:

```r
source("code/scripts/07_mmrr.R")
```

Main outputs:

* `outputs/tables/mmrr_coefficients.csv`
* `outputs/models/mmrr_fit.rds`

---

## Workflow 4: Secondary analyses (PCoA and dbRDA)

Purpose: Summarize genetic distance structure (PCoA) and test constraints by PfPR and spatial covariates (dbRDA).

Run:

```r
source("code/scripts/08_pcoa_dbrda.R")
```

Main outputs:

* `outputs/tables/dbrda_terms.csv`
* `outputs/models/pcoa.rds`
* `outputs/models/dbrda.rds`

---

## Workflow 5: Sensitivity analyses (Mantel and Rousset)

Purpose: Provide theory- and diagnostic-based checks.

Notes:

* Mantel and partial Mantel tests are treated as descriptive sensitivity analyses.
* Rousset regression provides a theory-linked isolation-by-distance sensitivity check with bootstrap confidence intervals.

Run:

```r
source("code/scripts/06_mantel_partial.R")
source("code/scripts/09_rousset_ibd.R")
```

Main outputs:

* `outputs/tables/mantel_results.csv`
* `outputs/tables/rousset_ibd.csv`
* `outputs/figures/rousset_ibd.png`

---

## Workflow 6: Figures, tables, and report rendering

Purpose: Produce manuscript-ready artifacts.

Run:

```r
source("code/scripts/10_figures.R")
source("code/scripts/11_tables.R")
```

Main outputs:

* `outputs/figures/` (PNG/PDF figures)
* `outputs/tables/` (CSV tables)

If you render a report (Quarto/markdown), run:

```r
source("code/scripts/12_report_render.R")
```

---

## Running the full pipeline

The recommended way to reproduce everything is to run the `targets` pipeline:

```r
library(targets)
tar_make()
```

To rerun from scratch:

```r
tar_destroy(destroy = "all")
tar_make()
```

---

## How to cite

See `CITATION.cff` for repository citation metadata.

If you archive a release on Zenodo, add the DOI to `CITATION.cff` and cite the Zenodo record in the manuscript.

---

## License

This repository uses a dual-license model (see `LICENSE`):

* Code (analysis scripts, pipeline, and functions) is licensed under the MIT License.
* Documentation, manuscript text, and figures are licensed under CC BY 4.0 unless otherwise noted.

---

## Contact information

George Paasi
Makerere University
Email: georgepaasi8@gmail.com

```
::contentReference[oaicite:0]{index=0}
```
