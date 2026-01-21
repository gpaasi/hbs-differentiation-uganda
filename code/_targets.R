library(targets)
library(yaml)
library(purrr)

purrr::walk(list.files('code/functions', full.names = TRUE, pattern='\\.R$'), source)

cfg <- yaml::read_yaml('code/config.yml')
set_seed(cfg$project$seed)

tar_option_set(packages = c('yaml','purrr'))

list(
  tar_target(config, cfg),
  tar_target(genotype_counts, read_counts(config$inputs$genotype_counts_csv), format='rds'),
  tar_target(pfpr, read_pfpr(config$inputs$pfpr_csv), format='rds'),
  tar_target(centroids, read_centroids(config$inputs$centroids_csv), format='rds'),

  tar_target(subregion_summary, summarise_subregions(genotype_counts), format='rds'),
  tar_target(subregion_hwe, { d <- add_hwe(subregion_summary); d$HWE_q <- bh_fdr(d$HWE_p); d }, format='rds'),

  tar_target(pairwise_fst, pairwise_fst_placeholder(subregion_summary), format='rds'),
  tar_target(global_fst, global_fst_placeholder(pairwise_fst)),

  tar_target(D_geo, geo_distance_matrix(centroids), format='rds'),
  tar_target(D_mal, mal_distance_matrix(pfpr), format='rds'),
  tar_target(D_gen, pairwise_fst, format='rds'),

  tar_target(mantel_geo, mantel_test(D_gen, D_geo, perms=config$parameters$permutations_mantel), format='rds'),
  tar_target(mantel_mal, mantel_test(D_gen, D_mal, perms=config$parameters$permutations_mantel), format='rds'),
  tar_target(partial_geo, partial_mantel(D_gen, D_geo, D_mal, perms=config$parameters$permutations_mantel), format='rds'),
  tar_target(partial_mal, partial_mantel(D_gen, D_mal, D_geo, perms=config$parameters$permutations_mantel), format='rds'),

  tar_target(mmrr_fit, mmrr(D_gen, D_geo, D_mal, perms=config$parameters$permutations_mmrr), format='rds'),

  tar_target(pcoa_fit, pcoa(D_gen), format='rds'),
  tar_target(dbrda_fit, { pred <- merge(pfpr, centroids, by='subregion'); dbrda_optional(D_gen, pred) }, format='rds'),

  tar_target(rousset_fit, rousset_regression(D_gen, D_geo,
    neg_to_zero=config$parameters$fst_negative_to_zero,
    log_base=config$parameters$rousset_log_base,
    bootstrap=config$parameters$bootstrap_rousset), format='rds'),

  tar_target(write_outputs, {
    outdir <- config$outputs$dir
    dir.create(file.path(outdir, 'tables'), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(outdir, 'models'), recursive=TRUE, showWarnings=FALSE)

    save_csv(subregion_hwe, file.path(outdir, 'tables', 'subregion_summary_hwe.csv'))
    save_rds(pairwise_fst, file.path(outdir, 'models', 'pairwise_fst.rds'))
    save_rds(mmrr_fit, file.path(outdir, 'models', 'mmrr_fit.rds'))
    save_rds(rousset_fit, file.path(outdir, 'models', 'rousset_fit.rds'))

    list(ok=TRUE, outdir=outdir)
  }, format='rds')
)
