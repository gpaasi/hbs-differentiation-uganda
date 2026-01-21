
# scripts/07_dbmem.R
# Step 7 — Spatial eigenvector mapping (dbMEM) to model spatial structure explicitly
# - Builds Region-level allele matrix (Hellinger transformed)
# - Generates Moran's Eigenvector Maps with adespatial::dbmem
# - Forward selection of MEMs (999 perms) using adjusted R² threshold
# - Fits final RDA(Y ~ selected MEMs);
# - Optional variation partitioning with malaria burden
# - Saves tidy tables and clear figures

message(">>> Step 7: dbMEM (spatial eigenvectors) starting...")

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(janitor)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(adegenet)
  library(geosphere)
})

.of <- function(...) here::here("outputs", ...)

# ---------- Load config and genind ----------
cfg <- yaml::read_yaml(here::here("config","config.yml"))
genind_rds <- here::here("outputs","models","genind_obj.rds")
if (!file.exists(genind_rds)) stop("Missing genind object. Run scripts/01_preprocess.R first.")
gi0 <- readRDS(genind_rds)

# Ensure Region grouping if available
st <- adegenet::strata(gi0)
st_names <- if (!is.null(st)) names(st) %>% janitor::make_clean_names() else character(0)
if ("region" %in% st_names) {
  adegenet::setPop(gi0) <- ~ region
  message("Pop set to Region.")
} else {
  message("Region strata not found; using current pop grouping.")
}

# Build genpop and allele matrix (pop x alleles), then Hellinger-transform
gp <- adegenet::genind2genpop(gi0)
A <- gp@tab
A[is.na(A)] <- 0
# Hellinger transform is recommended for community-like multivariate data
Y <- vegan::decostand(A, method = "hellinger")

regions <- rownames(Y)
if (is.null(regions) || length(regions) < 3) stop("Need at least 3 regions/populations.")

# ---------- Load coordinates (Region-level) ----------
coords_path <- here::here(cfg$coords$file %||% "data/raw/coords.csv")
if (!file.exists(coords_path)) stop("Coordinates file not found at: ", coords_path)
coords <- suppressMessages(readr::read_csv(coords_path, show_col_types = FALSE)) %>% janitor::clean_names()

rcol <- janitor::make_clean_names(cfg$coords$region_column %||% "Region")
lonc <- janitor::make_clean_names(cfg$coords$lon_column %||% "Longitude")
latc <- janitor::make_clean_names(cfg$coords$lat_column %||% "Latitude")
stopifnot(all(c(rcol,lonc,latc) %in% names(coords)))

coords <- coords %>%
  transmute(region = .data[[rcol]],
            longitude = suppressWarnings(as.numeric(.data[[lonc]])),
            latitude  = suppressWarnings(as.numeric(.data[[latc]]))) %>%
  filter(!is.na(longitude) & !is.na(latitude)) %>%
  distinct(region, .keep_all = TRUE)

# Harmonize order with Y
coords <- coords %>% filter(region %in% regions) %>% arrange(factor(region, levels = regions))
Y <- Y[coords$region, , drop = FALSE]
stopifnot(identical(coords$region, rownames(Y)))

# ---------- Create dbMEM eigenvectors ----------
# Ensure adespatial is available
if (!requireNamespace("adespatial", quietly = TRUE)) {
  stop("Package 'adespatial' is required for dbMEM. Please run scripts/00_setup.R and add 'adespatial' to the package list, or install manually: install.packages('adespatial')")
}
library(adespatial)

# dbMEM: uses a threshold equal to the longest edge of the Minimum Spanning Tree if thresh = NULL
mem <- adespatial::dbmem(as.matrix(coords[, c("longitude","latitude")]), silent = TRUE, mem.autocor = "positive")
MEM <- as.data.frame(mem)
rownames(MEM) <- coords$region

# Save eigenvalues and vectors
eig_vals <- attributes(mem)$values
ev_df <- tibble::tibble(component = paste0("MEM", seq_along(eig_vals)), eigenvalue = eig_vals)
readr::write_csv(ev_df, .of("tables","dbmem_eigenvalues.csv"))
MEM_out <- MEM %>% tibble::rownames_to_column("region")
readr::write_csv(MEM_out, .of("tables","dbmem_vectors_by_region.csv"))

# Scree plot of eigenvalues
p_scree <- ggplot(ev_df, aes(x = seq_along(eigenvalue), y = eigenvalue)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title = "dbMEM eigenvalues (positive spatial autocorrelation)",
       x = "Component", y = "Eigenvalue")
ggsave(.of("figures","dbmem_scree.png"), p_scree, width = 7, height = 5, dpi = 300)

# ---------- Forward selection of MEMs ----------
# Global RDA with all MEMs to get adjusted R² threshold
rda_all <- vegan::rda(Y ~ ., data = MEM)
R2adj_all <- vegan::RsquareAdj(rda_all)$adj.r.squared

# Forward selection (999 perms) to avoid overfitting
fs <- adespatial::forward.sel(Y, MEM, adjR2thresh = R2adj_all, nperm = 999)
# fs$order returns selected column indices in selection order
sel_idx <- sort(fs$order)
if (length(sel_idx) == 0) {
  message("No MEMs selected under adjusted R² threshold; will keep the first MEM as a minimal spatial term.")
  sel_idx <- 1
}
sel_names <- colnames(MEM)[sel_idx]

sel_table <- tibble::tibble(order = seq_along(sel_idx),
                            MEM = sel_names,
                            R2Cum = fs$R2Cum[order(fs$order)][seq_along(sel_idx)],
                            F = fs$F[order(fs$order)][seq_along(sel_idx)],
                            p = fs$pvalue[order(fs$order)][seq_along(sel_idx)])
readr::write_csv(sel_table, .of("tables","dbmem_forward_selection.csv"))

# ---------- Final RDA with selected MEMs ----------
MEMsel <- MEM[, sel_names, drop = FALSE]
rda_sel <- vegan::rda(Y ~ ., data = MEMsel)
a_overall <- anova.cca(rda_sel, permutations = 999)
axes_test <- anova.cca(rda_sel, by = "axis", permutations = 999)
rsq_adj <- vegan::RsquareAdj(rda_sel)

# Save summaries
readr::write_csv(
  tibble::tibble(R2 = rsq_adj$r.squared, R2_adj = rsq_adj$adj.r.squared,
                 p_overall = a_overall$`Pr(>F)`[1]),
  .of("tables","dbmem_rda_summary.csv")
)
axes_df <- tibble::tibble(axis = paste0("RDA", seq_along(axes_test$F)),
                          F = axes_test$F, p = axes_test$`Pr(>F)`)
readr::write_csv(axes_df, .of("tables","dbmem_rda_axes.csv"))

# Biplot of first two canonical axes (sites only; MEM arrows can be cluttered)
site_sc <- as.data.frame(scores(rda_sel, display = "sites", choices = 1:2, scaling = 2)) %>%
  tibble::rownames_to_column("region") %>%
  dplyr::rename(RDA1 = RDA1, RDA2 = RDA2)
p_sites <- ggplot(site_sc, aes(RDA1, RDA2, label = region)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 2) + ggrepel::geom_text_repel(size = 3) +
  theme_minimal() + labs(title = "dbMEM-constrained RDA: sites (Region)",
                         subtitle = paste0("Adj R²=", round(rsq_adj$adj.r.squared,3),
                                           ", p=", signif(a_overall$`Pr(>F)`[1],3)))
ggsave(.of("figures","dbmem_rda_sites.png"), p_sites, width = 7, height = 5, dpi = 300)

# ---------- Optional: Variation partitioning with malaria burden ----------
malaria_path <- here::here(cfg$malaria$file %||% "data/raw/malaria_burden.csv")
if (file.exists(malaria_path)) {
  mal <- suppressMessages(readr::read_csv(malaria_path, show_col_types = FALSE)) %>% janitor::clean_names()
  m_region <- janitor::make_clean_names(cfg$malaria$region_column %||% "Region")
  m_value  <- janitor::make_clean_names(cfg$malaria$value_column  %||% "Malaria_Parasites_Percent")
  if (all(c(m_region,m_value) %in% names(mal))) {
    mal <- mal %>%
      transmute(region = .data[[m_region]], malaria = suppressWarnings(as.numeric(.data[[m_value]]))) %>%
      filter(!is.na(malaria)) %>%
      distinct(region, .keep_all = TRUE)
    # align to regions
    mal <- mal %>% filter(region %in% regions) %>% arrange(factor(region, levels = regions))
    mal <- mal %>% filter(region %in% coords$region)
    # reorder Y to matched region set again
    keep <- intersect(rownames(Y), mal$region)
    if (length(keep) >= 3) {
      Yp <- Y[keep, , drop = FALSE]
      Xspace <- MEMsel[keep, , drop = FALSE]
      Xenv <- data.frame(malaria = mal$malaria[match(keep, mal$region)])
      vp <- vegan::varpart(Yp, ~ malaria, Xspace, data = Xenv)
      # Save fractions
      png(.of("figures","varpart_space_malaria.png"), width = 800, height = 600, res = 120)
      plot(vp, bg = c("skyblue","salmon"), Xnames = c("Malaria","Space"),
           id.size = 0.8, cutoff = 0.001, cex = 1.0)
      dev.off()
      # Numeric fractions (adjusted R2)
      # Extract adjusted R2 components via redundancy analysis
      # [a] unique malaria, [b] shared, [c] unique space; [d] residual
      # Compute using vegan::RsquareAdj with partial RDA fits
      rda_env <- vegan::rda(Yp ~ malaria, data = Xenv)
      rda_space <- vegan::rda(Yp ~ ., data = Xspace)
      rda_env_space <- vegan::rda(Yp ~ malaria + ., data = cbind(Xenv, Xspace))
      R2_env <- vegan::RsquareAdj(rda_env)$adj.r.squared
      R2_space <- vegan::RsquareAdj(rda_space)$adj.r.squared
      R2_both <- vegan::RsquareAdj(rda_env_space)$adj.r.squared
      a <- max(0, R2_env - (R2_both - R2_space))  # unique env
      c <- max(0, R2_space - (R2_both - R2_env))  # unique space
      b <- max(0, R2_both - a - c)                # shared
      d <- max(0, 1 - R2_both)                    # residual
      frac <- tibble::tibble(component = c("env_unique","shared","space_unique","residual"),
                             fraction = c(a,b,c,d))
      readr::write_csv(frac, .of("tables","varpart_space_malaria_fractions.csv"))
    }
  }
}

# ---------- Quick maps for first two MEMs ----------
if (ncol(MEM) >= 2) {
  mem12 <- coords %>%
    dplyr::select(region, longitude, latitude) %>%
    bind_cols(MEM[,1:2])
  p1 <- ggplot(mem12, aes(longitude, latitude)) +
    geom_point(aes(size = abs(MEM1), color = MEM1)) +
    scale_size_continuous(range = c(2,6)) +
    theme_minimal() + labs(title = "Spatial pattern: MEM1")
  p2 <- ggplot(mem12, aes(longitude, latitude)) +
    geom_point(aes(size = abs(MEM2), color = MEM2)) +
    scale_size_continuous(range = c(2,6)) +
    theme_minimal() + labs(title = "Spatial pattern: MEM2")
  ggsave(.of("figures","mem1_map.png"), p1, width = 6, height = 5, dpi = 300)
  ggsave(.of("figures","mem2_map.png"), p2, width = 6, height = 5, dpi = 300)
}

message(">>> Step 7 complete. See outputs/tables and outputs/figures for results.")
