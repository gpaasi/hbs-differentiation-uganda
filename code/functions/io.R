read_counts <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  assert_has_cols(df, c("subregion","HbAA","HbAS","HbSS"))
  df$HbAA <- as.integer(df$HbAA)
  df$HbAS <- as.integer(df$HbAS)
  df$HbSS <- as.integer(df$HbSS)
  df
}

read_pfpr <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  assert_has_cols(df, c("subregion","pfpr_mean"))
  df$pfpr_mean <- as.numeric(df$pfpr_mean)
  df
}

read_centroids <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  assert_has_cols(df, c("subregion","lat","lon"))
  df$lat <- as.numeric(df$lat)
  df$lon <- as.numeric(df$lon)
  df
}
