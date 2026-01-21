assert_has_cols <- function(df, cols) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

upper_tri_vec <- function(mat) mat[upper.tri(mat, diag = FALSE)]

zscore_mat <- function(mat) {
  v <- upper_tri_vec(mat)
  mu <- mean(v, na.rm = TRUE)
  sdv <- stats::sd(v, na.rm = TRUE)
  if (sdv == 0) return(mat * 0)
  (mat - mu) / sdv
}

set_seed <- function(seed) {
  set.seed(seed)
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion")
}
