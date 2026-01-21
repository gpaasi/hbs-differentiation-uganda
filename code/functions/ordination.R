pcoa <- function(D) {
  fit <- stats::cmdscale(D, k=2, eig=TRUE, add=TRUE)
  list(points=as.data.frame(fit$points), eig=fit$eig)
}

dbrda_optional <- function(Dgen, predictors_df) {
  if (!requireNamespace('vegan', quietly = TRUE)) {
    return(list(available=FALSE, message='Package vegan not installed; skipping dbRDA.'))
  }
  D <- stats::as.dist(Dgen)
  predictors <- predictors_df
  rownames(predictors) <- predictors$subregion
  predictors$subregion <- NULL
  mod <- vegan::capscale(D ~ ., data=predictors, add=TRUE)
  anova_terms <- vegan::anova.cca(mod, by='term', permutations=9999)
  list(available=TRUE, model=mod, anova_terms=anova_terms)
}
