
# scripts/01_preprocess.R
# Step 1 — Data preprocessing: clean inputs, build genind/genpop, export QC tables

message(">>> Step 1: Preprocessing starting...")

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(janitor)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(adegenet)
  library(poppr)
})

# ---------- helpers ----------
stop_if_missing <- function(df, cols, label="columns") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop("Missing ", label, ": ", paste(miss, collapse=", "))
}

tolower_nms <- function(x) {
  # clean_names() already lowercases, but keep a safe helper
  x %>% janitor::make_clean_names()
}

# ---------- load config ----------
cfg_path <- here::here("config","config.yml")
if (!file.exists(cfg_path)) stop("Missing config/config.yml. Run scripts/00_setup.R first.")
cfg <- yaml::read_yaml(cfg_path)

# normalize expected column names in config to clean_names() style
strata_cols_raw <- cfg$genotype$strata_columns %||% c("Region","Region1","Area","District")
id_col_raw      <- cfg$genotype$id_column %||% "SampleID"
allele_sep      <- cfg$genotype$allele_sep %||% "/"
missing_char    <- cfg$genotype$missing_char %||% "NA"
ploidy          <- cfg$genotype$ploidy %||% 2

strata_cols <- janitor::make_clean_names(strata_cols_raw)
id_col      <- janitor::make_clean_names(id_col_raw)

# ---------- read genotype file ----------
geno_path_default <- here::here("data","raw","final_individual_genotypes.csv")
geno_path <- geno_path_default
if (!file.exists(geno_path)) stop("Genotype file not found at: ", geno_path, "\nPut your CSV in data/raw/.")

message("Reading genotypes: ", geno_path)
geno_raw <- readr::read_csv(geno_path, show_col_types = FALSE) %>% janitor::clean_names()

# ---------- sanity checks ----------
# Make sure required columns exist (any of them present in file)
stop_if_missing(geno_raw, c(id_col), "ID column")
present_strata <- intersect(strata_cols, names(geno_raw))
if (!length(present_strata)) stop("None of the expected strata columns found. Expected any of: ",
                                  paste(strata_cols, collapse=", "))

# ---------- pick loci columns ----------
non_locus <- unique(c(id_col, present_strata))
candidate_loci <- setdiff(names(geno_raw), non_locus)

# keep columns that look like genotype strings: contain allele separator or equal to missing char or are not constant
looks_like_geno <- function(v) {
  if (is.numeric(v)) return(FALSE)
  v <- as.character(v)
  any(grepl(allele_sep, v, fixed = TRUE), na.rm = TRUE) ||
    any(v %in% missing_char) ||
    (length(na.omit(unique(v))) > 1)
}

loci_cols <- candidate_loci[sapply(geno_raw[candidate_loci], looks_like_geno)]
if (!length(loci_cols)) stop("No genotype loci detected. Check allele_sep in config (currently '", allele_sep, "').")

message("Detected ", length(loci_cols), " loci; ID: ", id_col, "; strata: ", paste(present_strata, collapse=", "))

# ---------- build working frame ----------
df <- geno_raw %>%
  dplyr::select(all_of(c(id_col, present_strata, loci_cols)))

# replace explicit missing_char strings with NA
for (lc in loci_cols) {
  df[[lc]][df[[lc]] %in% missing_char] <- NA
  df[[lc]] <- as.character(df[[lc]])
}

# drop duplicated rows (exact duplicates across loci)
dups <- duplicated(df[c(loci_cols)])
if (any(dups)) {
  message("Dropping ", sum(dups), " duplicated genotype rows (exact duplicate across loci).")
  df <- df[!dups, , drop = FALSE]
}

# ---------- create genind ----------
# pop grouping: prefer 'region' if present, else first present strata
pop_col <- if ("region" %in% names(df)) "region" else present_strata[[1]]
pop_factor <- as.factor(df[[pop_col]])

message("Building genind (ploidy=", ploidy, ", sep='", allele_sep, "') grouped by ", pop_col, "...")
genind_obj <- adegenet::df2genind(df[loci_cols],
                                  ploidy = ploidy,
                                  sep = allele_sep,
                                  NA.char = NA,
                                  pop = pop_factor)

# attach strata so we can switch groupings later
strata_df <- tibble::as_tibble(df[, present_strata, drop=FALSE])
adegenet::strata(genind_obj) <- strata_df

# ---------- save genind ----------
dir.create(here::here("outputs","models"), recursive = TRUE, showWarnings = FALSE)
saveRDS(genind_obj, here::here("outputs","models","genind_obj.rds"))
message("Saved genind_obj to outputs/models/genind_obj.rds")

# ---------- make genpop by Region (if available) ----------
make_and_save_genpop <- function(gi, strat_name, out_name) {
  if (!(strat_name %in% names(adegenet::strata(gi)))) {
    message("Strata '", strat_name, "' not found; skipping genpop for ", strat_name)
    return(invisible(NULL))
  }
  gi2 <- gi
  adegenet::setPop(gi2) <- as.formula(paste0("~", strat_name))
  gp <- adegenet::genind2genpop(gi2)
  saveRDS(gp, here::here("outputs","models", out_name))
  message("Saved ", out_name, " with ", nPop(gp), " populations.")
  invisible(gp)
}

gp_region  <- make_and_save_genpop(genind_obj, "region",  "genpop_by_region.rds")
gp_region1 <- make_and_save_genpop(genind_obj, "region1", "genpop_by_region1.rds")

# ---------- QC: missingness ----------
# per-locus missingness (share of NAs)
locus_miss <- sapply(df[loci_cols], function(x) mean(is.na(x)))
locus_miss_df <- tibble::tibble(locus = names(locus_miss), missing_rate = as.numeric(locus_miss)) %>%
  arrange(desc(missing_rate))

# per-individual missingness
ind_miss <- apply(df[loci_cols], 1, function(r) mean(is.na(r)))
ind_miss_df <- tibble::tibble(
  !!id_col := df[[id_col]],
  pop = df[[pop_col]],
  missing_rate = as.numeric(ind_miss)
) %>% arrange(desc(missing_rate))

# export QC tables
readr::write_csv(locus_miss_df, here::here("outputs","tables","qc_missingness_by_locus.csv"))
readr::write_csv(ind_miss_df,  here::here("outputs","tables","qc_missingness_by_individual.csv"))

message("Wrote QC tables to outputs/tables/.")

# ---------- optional: coordinate & malaria checks ----------
safe_read <- function(p) { if (file.exists(p)) suppressMessages(readr::read_csv(p, show_col_types = FALSE)) else NULL }

coords_path  <- here::here(cfg$coords$file %||% "data/raw/coords.csv")
malaria_path <- here::here(cfg$malaria$file %||% "data/raw/malaria_burden.csv")

coords <- safe_read(coords_path)
malaria <- safe_read(malaria_path)

if (!is.null(coords)) {
  coords <- janitor::clean_names(coords)
  rcol <- janitor::make_clean_names(cfg$coords$region_column %||% "Region")
  lonc <- janitor::make_clean_names(cfg$coords$lon_column %||% "Longitude")
  latc <- janitor::make_clean_names(cfg$coords$lat_column %||% "Latitude")
  stop_if_missing(coords, c(rcol, lonc, latc), "coords columns")
  # coverage report
  cov <- df %>%
    distinct(!!rlang::sym(pop_col)) %>%
    mutate(has_coord = !!rlang::sym(pop_col) %in% coords[[rcol]])
  readr::write_csv(cov, here::here("outputs","tables","coverage_regions_in_coords.csv"))
  message("Wrote coverage report: outputs/tables/coverage_regions_in_coords.csv")
}

if (!is.null(malaria)) {
  malaria <- janitor::clean_names(malaria)
  rcol <- janitor::make_clean_names(cfg$malaria$region_column %||% "Region")
  vcol <- janitor::make_clean_names(cfg$malaria$value_column %||% "Malaria_Parasites_Percent")
  stop_if_missing(malaria, c(rcol, vcol), "malaria columns")
  covm <- df %>%
    distinct(!!rlang::sym(pop_col)) %>%
    mutate(has_malaria = !!rlang::sym(pop_col) %in% malaria[[rcol]])
  readr::write_csv(covm, here::here("outputs","tables","coverage_regions_in_malaria.csv"))
  message("Wrote coverage report: outputs/tables/coverage_regions_in_malaria.csv")
}

# ---------- save processed sample sheet ----------
sample_sheet <- df %>%
  dplyr::select(all_of(c(id_col, present_strata))) %>%
  dplyr::rename(sample_id = all_of(id_col))
readr::write_csv(sample_sheet, here::here("data","processed","sample_sheet.csv"))
message("Saved data/processed/sample_sheet.csv")

message(">>> Step 1 complete. You can now proceed to structure/AMOVA/IBD/RDA steps.")
