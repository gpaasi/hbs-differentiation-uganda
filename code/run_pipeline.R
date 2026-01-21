library(yaml)
source_files <- list.files('code/functions', full.names = TRUE, pattern='\\.R$')
for (f in source_files) source(f)

cfg <- yaml::read_yaml('code/config.yml')
set_seed(cfg$project$seed)

counts <- read_counts(cfg$inputs$genotype_counts_csv)
pfpr <- read_pfpr(cfg$inputs$pfpr_csv)
cent <- read_centroids(cfg$inputs$centroids_csv)

summ <- summarise_subregions(counts)
summ <- add_hwe(summ)
summ$HWE_q <- bh_fdr(summ$HWE_p)

fst <- pairwise_fst_placeholder(summ)
D_geo <- geo_distance_matrix(cent)
D_mal <- mal_distance_matrix(pfpr)

m_geo <- mantel_test(fst, D_geo, perms=cfg$parameters$permutations_mantel)
m_mal <- mantel_test(fst, D_mal, perms=cfg$parameters$permutations_mantel)
mm <- mmrr(fst, D_geo, D_mal, perms=cfg$parameters$permutations_mmrr)
rs <- rousset_regression(fst, D_geo, bootstrap=cfg$parameters$bootstrap_rousset)

outdir <- cfg$outputs$dir
dir.create(file.path(outdir, 'tables'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir, 'models'), recursive=TRUE, showWarnings=FALSE)

save_csv(summ, file.path(outdir, 'tables', 'subregion_summary_hwe.csv'))
saveRDS(fst, file.path(outdir, 'models', 'pairwise_fst.rds'))
saveRDS(mm, file.path(outdir, 'models', 'mmrr_fit.rds'))
saveRDS(rs, file.path(outdir, 'models', 'rousset_fit.rds'))

message('Done. Outputs in: ', outdir)
