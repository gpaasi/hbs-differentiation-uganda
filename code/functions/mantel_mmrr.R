mantel_test <- function(A, B, perms=9999, method='pearson') {
  stopifnot(all(dim(A) == dim(B)))
  vA <- upper_tri_vec(A)
  vB <- upper_tri_vec(B)
  obs <- stats::cor(vA, vB, method=method, use='pairwise.complete.obs')
  n <- nrow(A)
  perm_stats <- numeric(perms)
  for (k in seq_len(perms)) {
    idx <- sample.int(n)
    Bp <- B[idx, idx]
    perm_stats[k] <- stats::cor(vA, upper_tri_vec(Bp), method=method, use='pairwise.complete.obs')
  }
  p <- (sum(abs(perm_stats) >= abs(obs)) + 1) / (perms + 1)
  list(r=obs, p=p)
}

partial_mantel <- function(A, B, C, perms=9999) {
  stopifnot(all(dim(A) == dim(B)), all(dim(A) == dim(C)))
  vA <- upper_tri_vec(A); vB <- upper_tri_vec(B); vC <- upper_tri_vec(C)
  rA <- stats::residuals(stats::lm(vA ~ vC))
  rB <- stats::residuals(stats::lm(vB ~ vC))
  obs <- stats::cor(rA, rB, use='pairwise.complete.obs')
  n <- nrow(A)
  perm_stats <- numeric(perms)
  for (k in seq_len(perms)) {
    idx <- sample.int(n)
    Bp <- B[idx, idx]
    vBp <- upper_tri_vec(Bp)
    rBp <- stats::residuals(stats::lm(vBp ~ vC))
    perm_stats[k] <- stats::cor(rA, rBp, use='pairwise.complete.obs')
  }
  p <- (sum(abs(perm_stats) >= abs(obs)) + 1) / (perms + 1)
  list(r=obs, p=p)
}

mmrr <- function(Dgen, Dgeo, Dmal, perms=10000) {
  y <- upper_tri_vec(Dgen)
  x1 <- upper_tri_vec(zscore_mat(Dgeo))
  x2 <- upper_tri_vec(zscore_mat(Dmal))
  fit <- stats::lm(y ~ x1 + x2)
  coefs <- stats::coef(fit)
  r2 <- summary(fit)$r.squared

  n <- nrow(Dgen)
  perm_b <- matrix(NA_real_, nrow=perms, ncol=2, dimnames=list(NULL, c('x1','x2')))
  perm_r2 <- numeric(perms)
  for (k in seq_len(perms)) {
    idx <- sample.int(n)
    Dp <- Dgen[idx, idx]
    yp <- upper_tri_vec(Dp)
    fp <- stats::lm(yp ~ x1 + x2)
    perm_b[k,] <- stats::coef(fp)[c('x1','x2')]
    perm_r2[k] <- summary(fp)$r.squared
  }
  p_x1 <- (sum(abs(perm_b[,'x1']) >= abs(coefs['x1'])) + 1)/(perms+1)
  p_x2 <- (sum(abs(perm_b[,'x2']) >= abs(coefs['x2'])) + 1)/(perms+1)
  p_r2 <- (sum(perm_r2 >= r2) + 1)/(perms+1)

  list(coef=coefs, r2=r2, p=list(x1=p_x1, x2=p_x2, r2=p_r2))
}
