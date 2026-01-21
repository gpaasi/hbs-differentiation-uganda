
# scripts/05_rda.R
# Step 5 — Redundancy Analysis (RDA) for malaria burden
# - Aligns individuals to Region-level malaria
# - Builds allele matrix with NA imputation
# - Runs RDA, permutation tests, variance explained
# - Saves triplot and tables
# - Computes Moran's I on RDA axes with km-threshold neighbors

message(">>> Step 5: RDA (malaria burden) starting...")

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(readr)
  library(dplyr)
  library(janitor)
  library(vegan)
  library(adegenet)
  library(ggplot2)
  library(spdep)
  library(geosphere)
  library(tibble)
  library(tidyr)
})

.of <- function(...) here::here("outputs", ...)

# ---------- Load config and genind ----------
cfg <- yaml::read_yaml(here::here("config","config.yml"))
genind_rds <- here::here("outputs","models","genind_obj.rds")
if (!file.exists(genind_rds)) stop("Missing genind object. Run scripts/01_preprocess.R first.")
gi0 <- readRDS(genind_rds)

# Ensure pop = Region if available
st <- adegenet::strata(gi0)
st_names <- if (!is.null(st)) names(st) %>% janitor::make_clean_names() else character(0)
if ("region" %in% st_names) {
  adegenet::setPop(gi0) <- ~ region
  message("Pop grouping set to Region.")
} else {
  message("Region strata not found; using current pop grouping.")
}

# ---------- Allele matrix (individual x loci) with NA imputation ----------
X <- adegenet::tab(gi0, NA.method = "mean")  # imputes NAs per locus mean
ind_ids <- rownames(X)
pop_vec <- as.character(pop(gi0))
stopifnot(length(pop_vec) == nrow(X))

# ---------- Load malaria burden (Region-level) ----------
mal_path <- here::here(cfg$malaria$file %||% "data/raw/malaria_burden.csv")
if (!file.exists(mal_path)) stop("Malaria burden file not found at: ", mal_path)
mal_raw <- suppressMessages(readr::read_csv(mal_path, show_col_types = FALSE)) %>% janitor::clean_names()

m_region <- janitor::make_clean_names(cfg$malaria$region_column %||% "Region")
m_value  <- janitor::make_clean_names(cfg$malaria$value_column  %||% "Malaria_Parasites_Percent")

stopifnot(m_region %in% names(mal_raw), m_value %in% names(mal_raw))

malaria <- mal_raw %>%
  transmute(region = .data[[m_region]], malaria = suppressWarnings(as.numeric(.data[[m_value]]))) %>%
  filter(!is.na(malaria)) %>%
  distinct(region, .keep_all = TRUE)

# Map malaria to individuals via their Region/pop
env_df <- tibble(sample_id = ind_ids, region = pop_vec) %>%
  left_join(malaria, by = "region")

# Drop individuals with missing malaria
keep_idx <- which(!is.na(env_df$malaria))
if (length(keep_idx) < 3) stop("Too few individuals with malaria data after join.")
X <- X[keep_idx, , drop = FALSE]
env_df <- env_df[keep_idx, , drop = FALSE]

# ---------- Optional: add geographic covariates for partial RDA ----------
coords_path <- here::here(cfg$coords$file %||% "data/raw/coords.csv")
coords_ok <- file.exists(coords_path)
if (coords_ok) {
  cr <- suppressMessages(readr::read_csv(coords_path, show_col_types = FALSE)) %>% janitor::clean_names()
  rcol <- janitor::make_clean_names(cfg$coords$region_column %||% "Region")
  lonc <- janitor::make_clean_names(cfg$coords$lon_column %||% "Longitude")
  latc <- janitor::make_clean_names(cfg$coords$lat_column %||% "Latitude")
  if (all(c(rcol,lonc,latc) %in% names(cr))) {
    cr <- cr %>% transmute(region = .data[[rcol]],
                           longitude = suppressWarnings(as.numeric(.data[[lonc]])),
                           latitude  = suppressWarnings(as.numeric(.data[[latc]]))) %>%
                 filter(!is.na(longitude) & !is.na(latitude)) %>%
                 distinct(region, .keep_all = TRUE)
    env_df <- env_df %>% left_join(cr, by = "region")
  } else {
    coords_ok <- FALSE
    message("Coords columns not found; skipping geographic covariates.")
  }
} else {
  message("Coords file not found; skipping geographic covariates.")
}

# ---------- RDA model ----------
# Basic RDA with malaria; if lon/lat available, also run partial RDA conditioning on them
rda_fit <- vegan::rda(X ~ malaria, data = env_df)
summary_rda <- summary(rda_fit)

# Significance tests
set.seed(cfg$seed %||% 1234)
a_overall <- anova.cca(rda_fit, permutations = 999)
a_terms   <- anova.cca(rda_fit, by = "terms", permutations = 999)
a_axes    <- anova.cca(rda_fit, by = "axis",  permutations = 999)
rsq_adj   <- vegan::RsquareAdj(rda_fit)

# Save tables
readr::write_csv(
  tibble(term = rownames(a_terms), F = a_terms$F, p = a_terms$`Pr(>F)`),
  .of("tables","rda_terms.csv")
)
readr::write_csv(
  tibble(axis = paste0("RDA", seq_along(a_axes$F)), F = a_axes$F, p = a_axes$`Pr(>F)`),
  .of("tables","rda_axes.csv")
)
readr::write_csv(
  tibble(R2 = unname(rsq_adj$r.squared), R2_adj = unname(rsq_adj$adj.r.squared)),
  .of("tables","rda_rsquare_adjusted.csv")
)

# Site scores (individuals) and biplot scores (constraints)
site_scores <- as.data.frame(scores(rda_fit, display = "sites", choices = 1:2, scaling = 2)) %>%
  tibble::rownames_to_column("sample_id") %>%
  rename(RDA1 = RDA1, RDA2 = RDA2) %>%
  left_join(env_df, by = "sample_id")

bp <- as.data.frame(scores(rda_fit, display = "bp", scaling = 2))
bp$label <- rownames(bp)

readr::write_csv(site_scores, .of("tables","rda_site_scores.csv"))
readr::write_csv(bp,          .of("tables","rda_biplot_scores.csv"))

# ---------- Plot: Triplot (sites + constraint arrow) ----------
p_triplot <- ggplot(site_scores, aes(x = RDA1, y = RDA2, color = region)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(alpha = 0.8, size = 2) +
  # arrow for malaria
  geom_segment(data = subset(bp, label == "malaria"),
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.15, "cm")), inherit.aes = FALSE) +
  geom_text(data = subset(bp, label == "malaria"),
            aes(x = RDA1, y = RDA2, label = "Malaria"),
            nudge_x = 0.05, nudge_y = 0.05, inherit.aes = FALSE) +
  labs(title = "RDA triplot: genetic ~ malaria burden",
       subtitle = paste0("Overall p=",
                         formatC(a_overall$`Pr(>F)`[1], format = "f", digits = 3),
                         " | Adj R²=",
                         formatC(unname(rsq_adj$adj.r.squared), format = "f", digits = 3)),
       x = "RDA1 (scaling 2)", y = "RDA2 (scaling 2)", color = "Region") +
  theme_minimal()
ggsave(.of("figures","rda_triplot.png"), p_triplot, width = 7, height = 5, dpi = 300)

# ---------- Optional partial RDA conditioning on lon/lat ----------
if (coords_ok && all(c("longitude","latitude") %in% names(env_df))) {
  if (sum(!is.na(env_df$longitude) & !is.na(env_df$latitude)) >= 3) {
    rda_partial <- vegan::rda(X ~ malaria + Condition(longitude + latitude), data = env_df)
    a_overall_p <- anova.cca(rda_partial, permutations = 999)
    rsq_adj_p   <- vegan::RsquareAdj(rda_partial)
    readr::write_csv(
      tibble(R2 = unname(rsq_adj_p$r.squared), R2_adj = unname(rsq_adj_p$adj.r.squared),
             p_overall = a_overall_p$`Pr(>F)`[1]),
      .of("tables","rda_partial_summary.csv")
    )
  }
}

# ---------- Spatial autocorrelation: Moran's I on RDA1 ----------
# Build km-based neighbors using geodesic distances and a threshold from config
if (coords_ok && all(c("longitude","latitude") %in% names(env_df))) {
  # replicate coords per individual (already joined by region)
  xy <- env_df %>% select(longitude, latitude) %>% as.data.frame()
  # km distance matrix
  dk <- geosphere::distm(as.matrix(xy), fun = geosphere::distHaversine) / 1000
  thr <- cfg$spatial$distance_threshold_km %||% 100
  # adjacency matrix within threshold (excluding self)
  W <- (dk > 0 & dk <= thr) * 1
  diag(W) <- 0
  # create listw from weights matrix
  lw <- spdep::mat2listw(W, style = "W")
  # Moran's I on RDA1 scores
  rda1 <- site_scores$RDA1
  if (length(rda1) == nrow(W) && !anyNA(rda1)) {
    mi <- spdep::moran.test(rda1, lw, alternative = "greater")
    moran_tab <- tibble(statistic = unname(mi$estimate[["Moran I statistic"]]),
                        expectation = unname(mi$estimate[["Expectation"]]),
                        variance = unname(mi$estimate[["Variance"]]),
                        p_value = mi$p.value)
    readr::write_csv(moran_tab, .of("tables","rda1_moransI.csv"))
  } else {
    message("Skipping Moran's I: NA in RDA1 or mismatch lengths.")
  }
}

message(">>> Step 5 complete. See outputs/tables and outputs/figures for results.")
