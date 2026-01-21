
# scripts/04_ibd_ibe.R
# Step 4 — Isolation by Distance (IBD) and Isolation by Environment (IBE) using Mantel & partial Mantel tests
# - Builds pop-level genetic distance (Nei)
# - Builds geographic distance from coords (km)
# - Builds environmental distance from malaria burden (% parasites)
# - Runs Mantel(IBD), Mantel(IBE), and optional partial Mantel tests
# - Saves tidy tables and publication-ready figures

message(">>> Step 4: IBD & IBE analysis starting...")

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(janitor)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(vegan)
  library(adegenet)
  library(poppr)
  library(geosphere)
})

.of <- function(...) here::here("outputs", ...)

# ---------- Load config and models ----------
cfg <- yaml::read_yaml(here::here("config","config.yml"))

genind_rds <- here::here("outputs","models","genind_obj.rds")
if (!file.exists(genind_rds)) stop("Missing genind object. Run scripts/01_preprocess.R first.")
genind_obj <- readRDS(genind_rds)

# Ensure grouping at Region level if present
st <- adegenet::strata(genind_obj)
st_names <- if (!is.null(st)) names(st) %>% janitor::make_clean_names() else character(0)
if ("region" %in% st_names) {
  adegenet::setPop(genind_obj) <- ~ region
  message("Grouping set to Region.")
} else if (nPop(genind_obj) == 1L) {
  stop("No Region strata found and only one population is set. Provide Region in Step 1.")
}

# ---------- Load coords ----------
coords_path <- here::here(cfg$coords$file %||% "data/raw/coords.csv")
if (!file.exists(coords_path)) stop("Coordinates file not found at: ", coords_path)
coords_raw <- suppressMessages(readr::read_csv(coords_path, show_col_types = FALSE)) %>%
  janitor::clean_names()

# Try to find region/lon/lat columns from config or best guess
region_col <- janitor::make_clean_names(cfg$coords$region_column %||% "Region")
lon_col    <- janitor::make_clean_names(cfg$coords$lon_column %||% "Longitude")
lat_col    <- janitor::make_clean_names(cfg$coords$lat_column %||% "Latitude")

guess_col <- function(df, candidates) {
  hit <- intersect(janitor::make_clean_names(candidates), names(df))
  if (length(hit)) hit[[1]] else NA_character_
}

if (!(region_col %in% names(coords_raw))) {
  region_col <- guess_col(coords_raw, c("Region","Pop","Population","Location","Area"))
}
if (!(lon_col %in% names(coords_raw))) {
  lon_col <- guess_col(coords_raw, c("Longitude","Long","Lon","X"))
}
if (!(lat_col %in% names(coords_raw))) {
  lat_col <- guess_col(coords_raw, c("Latitude","Lat","Y"))
}

stopifnot(region_col %in% names(coords_raw),
          lon_col %in% names(coords_raw),
          lat_col %in% names(coords_raw))

coords <- coords_raw %>%
  select(region = all_of(region_col),
         longitude = all_of(lon_col),
         latitude  = all_of(lat_col)) %>%
  mutate(
    longitude = suppressWarnings(as.numeric(longitude)),
    latitude  = suppressWarnings(as.numeric(latitude))
  ) %>%
  filter(!is.na(longitude) & !is.na(latitude)) %>%
  distinct(region, .keep_all = TRUE)

if (nrow(coords) < 2) stop("Not enough coordinate rows after cleaning.")

# ---------- Load malaria burden ----------
malaria_path <- here::here(cfg$malaria$file %||% "data/raw/malaria_burden.csv")
if (!file.exists(malaria_path)) stop("Malaria burden file not found at: ", malaria_path)
malaria <- suppressMessages(readr::read_csv(malaria_path, show_col_types = FALSE)) %>%
  janitor::clean_names()

m_region <- janitor::make_clean_names(cfg$malaria$region_column %||% "Region")
m_value  <- janitor::make_clean_names(cfg$malaria$value_column  %||% "Malaria_Parasites_Percent")

stopifnot(m_region %in% names(malaria),
          m_value  %in% names(malaria))

malaria <- malaria %>%
  select(region = all_of(m_region), malaria = all_of(m_value)) %>%
  mutate(malaria = suppressWarnings(as.numeric(malaria))) %>%
  filter(!is.na(malaria)) %>%
  distinct(region, .keep_all = TRUE)

# ---------- Harmonize Regions across sources ----------
# Use Region levels present in genind populations
pops <- levels(pop(genind_obj))
if (is.null(pops)) pops <- unique(as.character(pop(genind_obj)))
pops <- as.character(pops)

coords <- coords %>% filter(region %in% pops)
malaria <- malaria %>% filter(region %in% pops)

if (length(intersect(coords$region, malaria$region)) < 2) {
  stop("Insufficient overlap between coords and malaria regions after filtering by genind populations.")
}

# Reorder all by common region list
common_regions <- sort(intersect(intersect(pops, coords$region), malaria$region))

coords <- coords %>% filter(region %in% common_regions) %>% arrange(factor(region, levels = common_regions))
malaria <- malaria %>% filter(region %in% common_regions) %>% arrange(factor(region, levels = common_regions))

# Subset genind to just those populations (if needed)
gi <- genind_obj
gi <- gi[pop(gi) %in% common_regions, ]
adegenet::setPop(gi) <- droplevels(pop(gi))

# Build genpop for those regions
gp <- adegenet::genind2genpop(gi)

# ---------- Genetic distance (pop-level) ----------
# Nei's distance between populations
genetic_dist <- adegenet::dist.genpop(gp, method = 1)  # Nei 1972
# Ensure matrix with dimnames == regions
genetic_mat <- as.matrix(genetic_dist)
# Reindex/ensure ordering
genetic_mat <- genetic_mat[common_regions, common_regions]

# ---------- Geographic distance (km) ----------
geo_mat <- geosphere::distm(coords[, c("longitude","latitude")], fun = geosphere::distHaversine) / 1000
rownames(geo_mat) <- coords$region
colnames(geo_mat) <- coords$region
geo_mat <- geo_mat[common_regions, common_regions]

# ---------- Environmental (malaria) distance ----------
env_mat <- as.matrix(dist(malaria$malaria, method = "euclidean"))
rownames(env_mat) <- malaria$region
colnames(env_mat) <- malaria$region
env_mat <- env_mat[common_regions, common_regions]

# ---------- Save distance matrices ----------
readr::write_csv(as.data.frame(as.table(genetic_mat)), .of("tables","dist_genetic_nei_region_long.csv"))
readr::write_csv(as.data.frame(as.table(geo_mat)),     .of("tables","dist_geographic_km_region_long.csv"))
readr::write_csv(as.data.frame(as.table(env_mat)),     .of("tables","dist_environment_malaria_region_long.csv"))

# ---------- Mantel tests ----------
set.seed(cfg$seed %||% 1234)
m_ibd <- vegan::mantel(as.dist(genetic_mat), as.dist(geo_mat), method = "pearson", permutations = 999)
m_ibe <- vegan::mantel(as.dist(genetic_mat), as.dist(env_mat), method = "pearson", permutations = 999)

# Optional partial Mantel controlling for geography/environment
m_partial_env <- try(vegan::mantel.partial(as.dist(genetic_mat), as.dist(env_mat), as.dist(geo_mat), permutations = 999), silent = TRUE)
m_partial_geo <- try(vegan::mantel.partial(as.dist(genetic_mat), as.dist(geo_mat), as.dist(env_mat), permutations = 999), silent = TRUE)

# ---------- Save results table ----------
res <- dplyr::bind_rows(
  tibble::tibble(Test = "IBD (genetic ~ geographic)", r = unname(m_ibd$statistic), p = unname(m_ibd$signif)),
  tibble::tibble(Test = "IBE (genetic ~ malaria)", r = unname(m_ibe$statistic), p = unname(m_ibe$signif)),
  tibble::tibble(Test = "Partial (gen~malaria | geo)", r = if (inherits(m_partial_env,"try-error")) NA_real_ else unname(m_partial_env$statistic),
                 p = if (inherits(m_partial_env,"try-error")) NA_real_ else unname(m_partial_env$signif)),
  tibble::tibble(Test = "Partial (gen~geo | malaria)", r = if (inherits(m_partial_geo,"try-error")) NA_real_ else unname(m_partial_geo$statistic),
                 p = if (inherits(m_partial_geo,"try-error")) NA_real_ else unname(m_partial_geo$signif))
)
readr::write_csv(res, .of("tables","ibd_ibe_mantel_results.csv"))

# ---------- Plots: scatter with regression & annotated Mantel r,p ----------
# Make tidy long frames for scatterplots (upper triangle only to avoid duplicates)
upper_tri_idx <- function(m) {
  ut <- upper.tri(m, diag = FALSE)
  which(ut, arr.ind = TRUE)
}

make_scatter_df <- function(Mx, My, xlab, ylab) {
  idx <- upper_tri_idx(Mx)
  tibble::tibble(
    X = Mx[idx],
    Y = My[idx]
  ) %>% dplyr::filter(is.finite(X) & is.finite(Y)) %>%
    dplyr::mutate(xlab = xlab, ylab = ylab)
}

df_ibd <- make_scatter_df(geo_mat, genetic_mat, "Geographic distance (km)", "Genetic distance (Nei)")
df_ibe <- make_scatter_df(env_mat, genetic_mat, "Environmental distance (|Δ malaria|)", "Genetic distance (Nei)")

annot <- function(lbl, r, p) ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
                                               label = sprintf("%s\nr = %.3f, p = %.3f", lbl, r, p))

p_ibd <- ggplot(df_ibd, aes(x = X, y = Y)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  annot("Mantel (IBD)", unname(m_ibd$statistic), unname(m_ibd$signif)) +
  labs(title = "Isolation by Distance (IBD)", x = unique(df_ibd$xlab), y = unique(df_ibd$ylab)) +
  theme_minimal()

p_ibe <- ggplot(df_ibe, aes(x = X, y = Y)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  annot("Mantel (IBE)", unname(m_ibe$statistic), unname(m_ibe$signif)) +
  labs(title = "Isolation by Environment (IBE)", x = unique(df_ibe$xlab), y = unique(df_ibe$ylab)) +
  theme_minimal()

ggsave(.of("figures","ibd_scatter.png"), p_ibd, width = 7, height = 5, dpi = 300)
ggsave(.of("figures","ibe_scatter.png"), p_ibe, width = 7, height = 5, dpi = 300)

message(">>> Step 4 complete. See outputs/tables and outputs/figures for results.")
