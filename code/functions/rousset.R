rousset_regression <- function(Fst_mat, Dgeo, neg_to_zero=TRUE, log_base='e', bootstrap=2000) {
  F <- Fst_mat
  if (neg_to_zero) F[F < 0] <- 0
  tr <- F / (1 - F + 1e-12)
  y <- upper_tri_vec(tr)
  xdist <- upper_tri_vec(Dgeo)
  x <- if (log_base == '10') log10(xdist + 1e-12) else log(xdist + 1e-12)

  fit <- stats::lm(y ~ x)
  slope <- stats::coef(fit)[['x']]
  r2 <- summary(fit)$r.squared

  n_pairs <- length(y)
  boot_slopes <- numeric(bootstrap)
  for (b in seq_len(bootstrap)) {
    idx <- sample.int(n_pairs, replace=TRUE)
    fb <- stats::lm(y[idx] ~ x[idx])
    boot_slopes[b] <- stats::coef(fb)[2]
  }
  ci <- stats::quantile(boot_slopes, probs=c(0.025, 0.975), na.rm=TRUE)

  list(slope=slope, r2=r2, ci=as.numeric(ci))
}
