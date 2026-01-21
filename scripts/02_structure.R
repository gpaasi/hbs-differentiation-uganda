
# scripts/02_structure.R
# Step 2 — Population structure: PCA, DAPC (with safe fallbacks), pairwise Fst, NJ tree

message(">>> Step 2: Structure analysis starting...")

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(adegenet)
  library(hierfstat)
  library(ape)
  library(poppr)
})

# ---------- Load config and objects ----------
cfg <- yaml::read_yaml(here::here("config","config.yml"))
genind_rds <- here::here("outputs","models","genind_obj.rds")
if (!file.exists(genind_rds)) stop("Missing genind object. Run scripts/01_preprocess.R first.")
genind_obj <- readRDS(genind_rds)

# Try to load genpop by Region; if absent, build on the fly
genpop_rds <- here::here("outputs","models","genpop_by_region.rds")
if (file.exists(genpop_rds)) {
  genpop_obj <- readRDS(genpop_rds)
} else {
  message("genpop_by_region.rds not found; creating from genind using Region (if available).")
  if ("region" %in% names(adegenet::strata(genind_obj))) {
    gi <- genind_obj
    adegenet::setPop(gi) <- ~region
    genpop_obj <- adegenet::genind2genpop(gi)
  } else {
    genpop_obj <- adegenet::genind2genpop(genind_obj)
  }
  saveRDS(genpop_obj, genpop_rds)
}

# Helper to make safe file paths
.of <- function(...) here::here("outputs", ...)

# ---------- PCA (memory-aware) ----------
message("Running PCA (genlight + glPca if possible)...")
pca_ok <- TRUE
pca_err <- NULL
pca_scores <- NULL

tryCatch({
  # Convert to genlight to save memory
  gl <- try(as.genlight(genind_obj), silent = TRUE)
  if (inherits(gl, "try-error")) stop("as.genlight conversion failed; falling back to dudi.pca().")
  # choose nf reasonably
  nf <- min(20, nInd(genind_obj) - 1L, 50L)
  gpca <- glPca(gl, nf = nf, parallel = FALSE)
  ind_scores <- as.data.frame(gpca$scores)
  ind_scores$sample_id <- indNames(genind_obj)
  ind_scores$pop <- as.character(pop(genind_obj))
  # save
  readr::write_csv(ind_scores, .of("tables","pca_individual_scores.csv"))
  # plot PC1 vs PC2
  p_pca <- ggplot(ind_scores, aes(x = PC1, y = PC2, color = pop)) +
    geom_point(alpha = 0.8) +
    guides(color = guide_legend(title = "Population")) +
    labs(title = "PCA (PC1 vs PC2)", x = "PC1", y = "PC2") +
    theme_minimal()
  ggsave(.of("figures","pca_pc1_pc2.png"), p_pca, width = 7, height = 5, dpi = 300)
}, error = function(e) {
  pca_ok <<- FALSE
  pca_err <<- conditionMessage(e)
})

if (!pca_ok) {
  message("glPca path failed: ", pca_err, " — trying dudi.pca(tab(...)).")
  # Fallback: base PCA on allele count table with NAs imputed by mean
  X <- adegenet::tab(genind_obj, NA.method = "mean")
  # choose nf safely
  nf <- min(10, ncol(X)-1L, 50L)
  dp <- ade4::dudi.pca(X, scannf = FALSE, nf = nf)
  ind_scores <- as.data.frame(dp$li) %>% tibble::rownames_to_column("sample_id")
  ind_scores$pop <- as.character(pop(genind_obj))
  readr::write_csv(ind_scores, .of("tables","pca_individual_scores.csv"))
  p_pca <- ggplot(ind_scores, aes(x = Axis1, y = Axis2, color = pop)) +
    geom_point(alpha = 0.8) +
    guides(color = guide_legend(title = "Population")) +
    labs(title = "PCA (Axis1 vs Axis2)", x = "Axis1", y = "Axis2") +
    theme_minimal()
  ggsave(.of("figures","pca_axis1_axis2.png"), p_pca, width = 7, height = 5, dpi = 300)
}

# ---------- DAPC (robust, non-interactive) ----------
message("Running DAPC...")

dapc_ok <- TRUE
dapc_err <- NULL

# Safe n.pca choice: limit to min(50, nInd-1, 5% of loci)
n_loci <- nLoc(genind_obj)
n_ind  <- nInd(genind_obj)
n_pca  <- max(2L, min( min(50L, n_ind - 1L), max(2L, floor(0.05 * n_loci)) ))

# Try a data-driven clustering; if it fails, fall back to Region grouping
dapc_res <- NULL
grp_used <- NULL

tryCatch({
  # cap max clusters to avoid "more centers than distinct points"
  # derive an upper bound from unique rows of allele table
  X_small <- adegenet::tab(genind_obj, NA.method = "mean")
  uniq_n <- nrow(unique(round(X_small, 6)))
  kmax <- max(2L, min(10L, uniq_n - 1L))
  if (kmax < 2L) stop("Insufficient unique points for clustering; will use population grouping instead.")
  fc <- adegenet::find.clusters(genind_obj, n.pca = n_pca, max.n.clust = kmax, n.clust = NULL, choose.n.clust = FALSE, criterion = "bic", quiet = TRUE)
  grp_used <- fc$grp
  dapc_res <- adegenet::dapc(genind_obj, grp_used, n.pca = n_pca, n.da = min(length(levels(grp_used)) - 1L, 2L))
}, error = function(e) {
  dapc_ok <<- FALSE
  dapc_err <<- conditionMessage(e)
})

if (!dapc_ok || is.null(dapc_res)) {
  message("Clustering DAPC failed or not possible (", dapc_err, "). Falling back to Region grouping if available.")
  if (!is.null(pop(genind_obj))) {
    grp_used <- pop(genind_obj)
    dapc_res <- adegenet::dapc(genind_obj, grp_used, n.pca = n_pca, n.da = min(nPop(genind_obj) - 1L, 2L))
  } else {
    stop("No grouping available for DAPC.")
  }
}

# Save DAPC object and plots
saveRDS(dapc_res, .of("models","dapc_result.rds"))
members <- as.data.frame(dapc_res$posterior) %>% tibble::rownames_to_column("sample_id")
readr::write_csv(members, .of("tables","dapc_membership_posterior.csv"))

p_dapc <- scatter(dapc_res, scree.da = TRUE, posi.da = "bottomright", ellipse = 0.5, bg = "white")
# Save as PDF and PNG
ggsave(.of("figures","dapc_scatter.png"), plot = p_dapc, width = 7, height = 5, dpi = 300)

# ---------- Pairwise and global Fst (hierfstat) ----------
message("Computing pairwise and global Fst (Weir & Cockerham)...")

# hierfstat format
hdat <- adegenet::genind2hierfstat(genind_obj)
# pairwise WC Fst
pwfst <- hierfstat::pairwise.WCfst(hdat)
readr::write_csv(as.data.frame(pwfst), .of("tables","pairwise_fst_wc84.csv"))

# global F-statistics
wc <- hierfstat::wc(hdat)
glob <- data.frame(Fst_overall = wc$Fst["overall"])
readr::write_csv(glob, .of("tables","global_fst_wc84.csv"))

# Heatmap of pairwise Fst
pw_long <- reshape2::melt(as.matrix(pwfst), varnames = c("Pop1","Pop2"), value.name = "Fst")
p_fst <- ggplot(pw_long, aes(Pop1, Pop2, fill = Fst)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Pairwise Fst (WC84)", x = "Population", y = "Population")
ggsave(.of("figures","pairwise_fst_heatmap.png"), p_fst, width = 7, height = 6, dpi = 300)

# ---------- NJ tree on population distances ----------
message("Building NJ tree on genpop Nei's distance...")

# Nei's distance between populations
D <- adegenet::dist.genpop(genpop_obj, method = 1)  # 1 = Nei's 1972
tree <- ape::nj(as.dist(D))

# save tree plot
png(.of("figures","nj_tree_populations.png"), width = 900, height = 700, res = 120)
plot(tree, main = "Neighbor-Joining Tree (Nei's distance, populations)", cex = 0.8)
ape::add.scale.bar()
dev.off()

# save Newick
ape::write.tree(tree, file = .of("tables","nj_tree_populations.newick"))

message(">>> Step 2 complete. Outputs saved in outputs/figures and outputs/tables.")
