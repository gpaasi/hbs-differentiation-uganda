\
    # scripts/utils_io.R
    # Small helpers for safe reading/writing

    quiet_read_csv <- function(path, ...) {
      stopifnot(file.exists(path))
      readr::read_csv(path, show_col_types = FALSE, progress = FALSE, ...)
    }

    write_csv_safely <- function(df, path, ...) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(df, path, na = "", ...)
      message("Wrote: ", path)
      invisible(path)
    }
