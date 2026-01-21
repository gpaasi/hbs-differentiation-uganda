
# scripts/03_amova.R
# Step 3 — AMOVA (robust, permutation-tested) with automatic hierarchy detection and clear outputs

message(">>> Step 3: AMOVA starting...")

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(adegenet)
  library(poppr)
  library(ade4)
})

.of <- function(...) here::here("outputs", ...)

# ---------- Load config and genind ----------
cfg <- yaml::read_yaml(here::here("config","config.yml"))
genind_rds <- here::here("outputs","models","genind_obj.rds")
if (!file.exists(genind_rds)) stop("Missing genind object. Run scripts/01_preprocess.R first.")
genind_obj <- readRDS(genind_rds)

# ---------- Detect available hierarchy ----------
st <- adegenet::strata(genind_obj)
if (is.null(st)) stop("No strata found on genind_obj. Ensure Step 1 attached Region/Region1/District strata.")

st_names <- names(st) %>% tolower()
has_region  <- "region"  %in% st_names
has_region1 <- "region1" %in% st_names
has_district <- "district" %in% st_names

message("Strata available: ", paste(names(st), collapse=", "))

# Helper to build formula by presence
get_hier_formula <- function() {
  f <- NULL
  if (has_region1 && has_region && has_district) {
    f <- ~ region1/region/district
  } else if (has_region1 && has_region) {
    f <- ~ region1/region
  } else if (has_region) {
    f <- ~ region
  } else {
    # fallback to current pop (if set) or the first strata column
    if (!is.null(pop(genind_obj))) return(~ Pop)
    # set temporary pop to first strata
    adegenet::setPop(genind_obj) <<- as.formula(paste0("~", names(st)[1]))
    f <- as.formula(paste0("~", names(st)[1]))
  }
  return(f)
}

hier_formula <- get_hier_formula()
message("Hierarchy formula selected: ", deparse(hier_formula))

# ---------- Distance choice (robust to memory) ----------
# Try a light-weight distance (bitwise) first; if that fails, fallback to Euclidean on allele table
dist_obj <- NULL
dist_note <- NULL
tryCatch({
  dist_obj <- poppr::bitwise.dist(genind_obj, percent = FALSE, missing_match = TRUE, scale_missing = TRUE)
  dist_note <- "bitwise.dist (scale_missing=TRUE)"
}, error = function(e) {
  message("bitwise.dist failed: ", conditionMessage(e), " — falling back to Euclidean distance on tab(genind).")
  X <- adegenet::tab(genind_obj, NA.method = "mean")
  dist_obj <<- dist(X)
  dist_note <<- "Euclidean dist on allele table (NA.method='mean')"
})

message("Using distance: ", dist_note)

# ---------- Run AMOVA ----------
set.seed(cfg$seed %||% 1234)

amv <- NULL
amv_test <- NULL

# poppr.amova can take a genind + hier formula and optional dist
amv <- poppr::poppr.amova(genind_obj, hier = hier_formula, dist = dist_obj, nperm = 0)  # compute once
# permutation test (ade4::randtest method for amova objects)
amv_test <- ade4::randtest(amv, nrepet = 999)

# ---------- Extract results ----------
# Variance components
vc <- as.data.frame(amv$componentsofcovariance)
vc$Source <- rownames(vc)
vc_long <- vc %>%
  dplyr::mutate(Percent = 100 * Sigma / sum(Sigma)) %>%
  dplyr::select(Source, Sigma, Percent, Df)

readr::write_csv(vc_long, .of("tables","amova_variance_components.csv"))

# Phi-statistics (analogous to F-statistics)
phi <- as.data.frame(amv$statphi)
phi$Phi <- rownames(phi)
readr::write_csv(phi, .of("tables","amova_phi_stats.csv"))

# Permutation p-values
perm_df <- data.frame(
  Observed = amv_test$obs,
  Pvalue = amv_test$pvalue
)
readr::write_csv(perm_df, .of("tables","amova_permutation_results.csv"))

# ---------- Plots ----------
# 1) Variance partition bar
p_var <- ggplot(vc_long, aes(x = reorder(Source, -Percent), y = Percent)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = -0.3, size = 3) +
  labs(title = "AMOVA: Variance Partition", x = "Source", y = "Variance Explained (%)",
       subtitle = paste0("Distance: ", dist_note)) +
  theme_minimal()
ggsave(.of("figures","amova_variance_partition.png"), p_var, width = 7, height = 5, dpi = 300)

# 2) Permutation distribution
perm_plot_df <- data.frame(Simulated = amv_test$sim)
p_perm <- ggplot(perm_plot_df, aes(x = Simulated)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = amv_test$obs, linetype = "dashed") +
  labs(title = "AMOVA Permutation Test",
       subtitle = paste("Observed =", round(amv_test$obs, 4), "| p =", signif(amv_test$pvalue, 3)),
       x = "Simulated statistic", y = "Count") +
  theme_minimal()
ggsave(.of("figures","amova_permutation_hist.png"), p_perm, width = 7, height = 5, dpi = 300)

# ---------- Also run single-level AMOVA by Region (if Region exists) ----------
if (has_region) {
  message("Running additional single-level AMOVA by Region...")
  amv_region <- try(poppr::poppr.amova(genind_obj, hier = ~ region, dist = dist_obj, nperm = 0), silent = TRUE)
  if (!inherits(amv_region, "try-error")) {
    vc_r <- as.data.frame(amv_region$componentsofcovariance)
    vc_r$Source <- rownames(vc_r)
    vc_r <- vc_r %>% mutate(Percent = 100*Sigma/sum(Sigma)) %>% select(Source, Sigma, Percent, Df)
    readr::write_csv(vc_r, .of("tables","amova_variance_components_region.csv"))
  }
}

message(">>> Step 3 complete. See outputs/tables and outputs/figures for results.")
