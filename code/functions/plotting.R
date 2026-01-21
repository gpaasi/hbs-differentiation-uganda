save_csv <- function(df, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  utils::write.csv(df, path, row.names=FALSE)
}

save_rds <- function(obj, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  saveRDS(obj, path)
}
