# Malaria burden and geographic distance jointly predict HbS differentiation in Uganda

This repository contains a reproducible analysis workflow for a cross-sectional ecological spatial study of HbS (HBB) differentiation across 15 Ugandan subregions, relating locus-specific genetic differentiation to geographic separation and malaria-burden dissimilarity (PfPR2–10).

## What’s included
- A `targets` pipeline scaffold (`code/_targets.R`) that mirrors the Methods sections:
  - subregion summaries, heterozygosity, HWE, F_IS
  - global and pairwise F_ST (placeholder implementation in this scaffold)
  - distance matrices: D_gen, D_geo, D_mal
  - Mantel and partial Mantel (sensitivity)
  - MMRR (primary inference)
  - PCoA and optional dbRDA (secondary)
  - Rousset IBD regression with bootstrap CI (sensitivity)
- A small toy dataset under `data/toy/` so the pipeline runs end-to-end.

## Quick start (toy run)
1. Open R (>= 4.4.1).
2. Install required packages:
```r
install.packages(c("targets", "yaml", "purrr"))
```
3. Run the pipeline:
```r
library(targets)
tar_make()
```
Outputs appear in `outputs/`.

## Using your real data
Place your real files in `data/raw/` (kept out of Git). Update paths and parameters in `code/config.yml`.

Expected inputs (recommended):
- `genotype_counts.csv` with subregion-level counts for HbAA, HbAS, HbSS (or adapt the ingest step for individual-level genotypes).
- subregion centroids (WGS84) or polygons.
- PfPR2–10 mean per subregion (for the time window you report).

See `data/metadata/data_dictionary.md` for required fields and examples.

## Reproducibility notes
- This scaffold avoids heavyweight dependencies by using base R implementations for distance calculations and permutation tests.
- The F_ST implementation included here is a simple placeholder. For the manuscript, replace it with a Weir–Cockerham implementation (e.g., via `hierfstat` or `pegas`) and document the exact method used.

## Suggested citation
See `CITATION.cff`.
