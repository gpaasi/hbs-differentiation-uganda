\
    # scripts/utils_genetics.R
    # Genetics-specific helpers

    make_genind_from_wide <- function(df, loci_cols, strata = NULL, ploidy = 2, sep = "/", NA.char = "NA") {
      stopifnot(all(loci_cols %in% names(df)))
      X <- df[, loci_cols, drop = FALSE]
      adegenet::df2genind(
        X, ploidy = ploidy, sep = sep, NA.char = NA.char,
        pop = if (!is.null(strata)) as.factor(df[[strata]]) else NULL
      )
    }
