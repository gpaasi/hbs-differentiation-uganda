
# scripts/06_gdm.R
# Step 6 — Generalized Dissimilarity Modeling (GDM) with safe fallback to GAM
# - Builds pop-level genetic, geographic, and environmental (malaria) distances
# - Tries GDM (monotone I-splines); if not available/fails, fits a GAM on pairwise distances
# - Exports tidy tables and clear figures

message(">>> Step 6: GDM starting...")

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(janitor)
  library(ggplot2)
  library(vegan)
  library(adegenet)
  library(poppr)
  library(geosphere)
})

.of <- function(...) here::here("outputs", ...)

# ---------- Load config and objects ----------
cfg <- yaml::read_yaml(here::here("config","config.yml"))

genind_rds <- here::here("outputs","models","genind_obj.rds")
if (!file.exists(genind_rds)) stop("Missing genind object. Run scripts/01_preprocess.R first.")
genind_obj <- readRDS(genind_rds)

# Ensure Region grouping if present
st <- adegenet::strata(genind_obj)
st_names <- if (!is.null(st)) names(st) %>% janitor::make_clean_names() else character(0)
if ("region" %in% st_names) {
  adegenet::setPop(genind_obj) <- ~ region
  message("Pop set to Region")
}

# ---------- Load coords & malaria ----------
coords_path <- here::here(cfg$coords$file %||% "data/raw/coords.csv")
malaria_path <- here::here(cfg$malaria$file %||% "data/raw/malaria_burden.csv")
stopifnot(file.exists(coords_path), file.exists(malaria_path))

coords <- suppressMessages(readr::read_csv(coords_path, show_col_types = FALSE)) %>% janitor::clean_names()
malaria <- suppressMessages(readr::read_csv(malaria_path, show_col_types = FALSE)) %>% janitor::clean_names()

rcol <- janitor::make_clean_names(cfg$coords$region_column %||% "Region")
lonc <- janitor::make_clean_names(cfg$coords$lon_column %||% "Longitude")
latc <- janitor::make_clean_names(cfg$coords$lat_column %||% "Latitude")

m_region <- janitor::make_clean_names(cfg$malaria$region_column %||% "Region")
m_value  <- janitor::make_clean_names(cfg$malaria$value_column  %||% "Malaria_Parasites_Percent")

stopifnot(all(c(rcol,lonc,latc) %in% names(coords)), all(c(m_region,m_value) %in% names(malaria)))

coords <- coords %>%
  transmute(region = .data[[rcol]],
            longitude = suppressWarnings(as.numeric(.data[[lonc]])),
            latitude = suppressWarnings(as.numeric(.data[[latc]]))) %>%
  filter(!is.na(longitude) & !is.na(latitude)) %>%
  distinct(region, .keep_all = TRUE)

malaria <- malaria %>%
  transmute(region = .data[[m_region]], malaria = suppressWarnings(as.numeric(.data[[m_value]]))) %>%
  filter(!is.na(malaria)) %>%
  distinct(region, .keep_all = TRUE)

# Harmonize region lists across sources
pops <- as.character(levels(pop(genind_obj)))
if (is.null(pops)) pops <- unique(as.character(pop(genind_obj)))
common_regions <- Reduce(intersect, list(pops, coords$region, malaria$region))
if (length(common_regions) < 3) stop("Need at least 3 overlapping regions across genetics/coords/malaria.")

coords <- coords %>% filter(region %in% common_regions) %>% arrange(factor(region, levels = common_regions))
malaria <- malaria %>% filter(region %in% common_regions) %>% arrange(factor(region, levels = common_regions))

# Subset genind
gi <- genind_obj[pop(genind_obj) %in% common_regions, ]
adegenet::setPop(gi) <- droplevels(pop(gi))
gp <- adegenet::genind2genpop(gi)

# ---------- Distances ----------
# Genetic (Nei 1972)
genetic_mat <- as.matrix(adegenet::dist.genpop(gp, method = 1))
genetic_mat <- genetic_mat[common_regions, common_regions]

# Geographic (km, Haversine)
geo_mat <- geosphere::distm(coords[, c("longitude","latitude")], fun = geosphere::distHaversine) / 1000
rownames(geo_mat) <- coords$region; colnames(geo_mat) <- coords$region
geo_mat <- geo_mat[common_regions, common_regions]

# Environmental (malaria diff)
env_mat <- as.matrix(dist(malaria$malaria, method = "euclidean"))
rownames(env_mat) <- malaria$region; colnames(env_mat) <- malaria$region
env_mat <- env_mat[common_regions, common_regions]

# ---------- Pairwise long table ----------
upper_idx <- which(upper.tri(genetic_mat), arr.ind = TRUE)
pairs_df <- tibble::tibble(
  region_i = rownames(genetic_mat)[upper_idx[,1]],
  region_j = colnames(genetic_mat)[upper_idx[,2]],
  gen_dist = genetic_mat[upper_idx],
  geo_km   = geo_mat[upper_idx],
  mal_diff = env_mat[upper_idx]
) %>%
  filter(is.finite(gen_dist) & is.finite(geo_km) & is.finite(mal_diff))

readr::write_csv(pairs_df, .of("tables","gdm_pairs_input.csv"))

# ---------- Try GDM ----------
gdm_ok <- FALSE
gdm_err <- NULL

try({
  suppressPackageStartupMessages(library(gdm))
  # Build site-level environmental data
  site_env <- coords %>% left_join(malaria, by = "region") %>%
    transmute(site = region, Longitude = longitude, Latitude = latitude, malaria = malaria)

  # Attempt to format using distance response (bioFormat = 3)
  sp_table <- gdm::formatsitepair(
    response = genetic_mat,
    bioFormat = 3,
    siteColumn = "site",
    XColumn = "Longitude",
    YColumn = "Latitude",
    predData = site_env
  )
  gdm_model <- gdm::gdm(sp_table, geo = TRUE)
  gdm_ok <- TRUE

  # Save summary
  capture.output(summary(gdm_model), file = .of("tables","gdm_summary.txt"))
  # Variable importance (if available)
  vi <- try(gdm::gdm.varImp(sp_table, geo = TRUE, nPerm = 50, parallel = FALSE), silent = TRUE)
  if (!inherits(vi, "try-error")) {
    readr::write_csv(as.data.frame(vi$explained), .of("tables","gdm_varimp_explained.csv"))
    readr::write_csv(as.data.frame(vi$permutations), .of("tables","gdm_varimp_permutations.csv"))
  }

  # Plot I-spline partials if supported
  png(.of("figures","gdm_partial_curves.png"), width = 900, height = 600, res = 120)
  plot(gdm_model, plot.layout = c(1,2))
  dev.off()
}, silent = TRUE)

if (!gdm_ok) {
  message("GDM failed/not available; error suppressed. Falling back to GAM on pairwise distances.")
  # ---------- Fallback: GAM on pairwise distances ----------
  suppressPackageStartupMessages(library(mgcv))

  # Transform response to stabilize variance
  y <- pairs_df$gen_dist
  # Fit GAM with thin-plate regression splines
  gam_fit <- mgcv::gam(y ~ s(geo_km, k = 5) + s(mal_diff, k = 5), data = pairs_df, method = "REML")

  # Save summary
  sink(.of("tables","gam_gdm_fallback_summary.txt")); print(summary(gam_fit)); sink()

  # Approx variable importance via drop-in-deviance (R^2 reduction)
  r2_full <- summary(gam_fit)$r.sq
  # Remove each smoother
  gam_geo_only <- mgcv::gam(y ~ s(geo_km, k = 5), data = pairs_df, method = "REML")
  gam_env_only <- mgcv::gam(y ~ s(mal_diff, k = 5), data = pairs_df, method = "REML")
  r2_geo <- summary(gam_geo_only)$r.sq
  r2_env <- summary(gam_env_only)$r.sq
  varimp <- tibble::tibble(
    model = c("full","geo_only","env_only"),
    r2 = c(r2_full, r2_geo, r2_env)
  )
  readr::write_csv(varimp, .of("tables","gam_var_importance_r2.csv"))

  # Partial dependence curves
  new_geo <- data.frame(geo_km = seq(min(pairs_df$geo_km), max(pairs_df$geo_km), length.out = 200),
                        mal_diff = median(pairs_df$mal_diff, na.rm = TRUE))
  new_env <- data.frame(mal_diff = seq(min(pairs_df$mal_diff), max(pairs_df$mal_diff), length.out = 200),
                        geo_km = median(pairs_df$geo_km, na.rm = TRUE))
  pd_geo <- tibble::tibble(geo_km = new_geo$geo_km, fit = predict(gam_fit, newdata = new_geo, type = "response"))
  pd_env <- tibble::tibble(mal_diff = new_env$mal_diff, fit = predict(gam_fit, newdata = new_env, type = "response"))
  readr::write_csv(pd_geo, .of("tables","gam_partial_geo.csv"))
  readr::write_csv(pd_env, .of("tables","gam_partial_env.csv"))

  # Plots
  p1 <- ggplot(pairs_df, aes(x = geo_km, y = gen_dist)) +
    geom_point(alpha = 0.6) +
    geom_line(data = pd_geo, aes(x = geo_km, y = fit), linewidth = 1) +
    theme_minimal() +
    labs(title = "GAM fallback: effect of geographic distance", x = "Geographic distance (km)", y = "Genetic distance (Nei)")
  ggsave(.of("figures","gam_partial_geo.png"), p1, width = 7, height = 5, dpi = 300)

  p2 <- ggplot(pairs_df, aes(x = mal_diff, y = gen_dist)) +
    geom_point(alpha = 0.6) +
    geom_line(data = pd_env, aes(x = mal_diff, y = fit), linewidth = 1) +
    theme_minimal() +
    labs(title = "GAM fallback: effect of malaria difference", x = "|Δ malaria burden|", y = "Genetic distance (Nei)")
  ggsave(.of("figures","gam_partial_env.png"), p2, width = 7, height = 5, dpi = 300)

  # Observed vs fitted
  pairs_df$fit <- predict(gam_fit, type = "response")
  p3 <- ggplot(pairs_df, aes(x = fit, y = gen_dist)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = "GAM fallback: observed vs fitted", x = "Fitted", y = "Observed")
  ggsave(.of("figures","gam_obs_vs_fit.png"), p3, width = 6, height = 5, dpi = 300)
}

message(">>> Step 6 complete. See outputs/tables and outputs/figures for results.")
