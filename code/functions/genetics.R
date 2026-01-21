summarise_subregions <- function(counts) {
  counts$N <- counts$HbAA + counts$HbAS + counts$HbSS
  counts$pS <- (counts$HbAS + 2*counts$HbSS) / (2*counts$N)
  counts$pA <- 1 - counts$pS
  counts$Ho <- counts$HbAS / counts$N
  counts$He <- 2*counts$pS*counts$pA
  counts$FIS <- ifelse(counts$He > 0, 1 - counts$Ho / counts$He, NA_real_)
  counts$prop_AA <- counts$HbAA / counts$N
  counts$prop_AS <- counts$HbAS / counts$N
  counts$prop_SS <- counts$HbSS / counts$N
  counts
}

hwe_chisq <- function(row) {
  N <- row[["N"]]
  pS <- row[["pS"]]
  pA <- 1 - pS
  exp <- c(N*(pA^2), N*(2*pA*pS), N*(pS^2))
  obs <- c(row[["HbAA"]], row[["HbAS"]], row[["HbSS"]])
  if (any(exp < 5)) return(list(stat=NA_real_, p=NA_real_))
  stat <- sum((obs-exp)^2/exp)
  p <- stats::pchisq(stat, df=1, lower.tail=FALSE)
  list(stat=stat, p=p)
}

add_hwe <- function(df) {
  res <- lapply(seq_len(nrow(df)), function(i) hwe_chisq(df[i,]))
  df$HWE_chisq <- vapply(res, `[[`, numeric(1), "stat")
  df$HWE_p <- vapply(res, `[[`, numeric(1), "p")
  df
}

bh_fdr <- function(p) stats::p.adjust(p, method = "BH")

pairwise_fst_placeholder <- function(df) {
  subs <- df$subregion
  p <- df$pS
  n <- length(subs)
  mat <- matrix(0, n, n, dimnames=list(subs, subs))
  for (i in 1:(n-1)) for (j in (i+1):n) {
    pbar <- (p[i] + p[j]) / 2
    Ht <- 2*pbar*(1-pbar)
    Hs <- (2*p[i]*(1-p[i]) + 2*p[j]*(1-p[j])) / 2
    fst <- ifelse(Ht > 0, (Ht - Hs)/Ht, 0)
    mat[i,j] <- fst; mat[j,i] <- fst
  }
  mat
}

global_fst_placeholder <- function(pairwise_mat) mean(upper_tri_vec(pairwise_mat), na.rm=TRUE)
