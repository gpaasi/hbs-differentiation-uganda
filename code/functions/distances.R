haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371.0088
  to_rad <- function(x) x * pi / 180
  lat1 <- to_rad(lat1); lon1 <- to_rad(lon1)
  lat2 <- to_rad(lat2); lon2 <- to_rad(lon2)
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  a <- sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  R * c
}

geo_distance_matrix <- function(centroids) {
  assert_has_cols(centroids, c('subregion','lat','lon'))
  subs <- centroids$subregion
  n <- length(subs)
  mat <- matrix(0, n, n, dimnames=list(subs, subs))
  for (i in 1:(n-1)) for (j in (i+1):n) {
    d <- haversine_km(centroids$lat[i], centroids$lon[i], centroids$lat[j], centroids$lon[j])
    mat[i,j] <- d; mat[j,i] <- d
  }
  mat
}

mal_distance_matrix <- function(pfpr_df) {
  assert_has_cols(pfpr_df, c('subregion','pfpr_mean'))
  subs <- pfpr_df$subregion
  n <- length(subs)
  mat <- matrix(0, n, n, dimnames=list(subs, subs))
  for (i in 1:(n-1)) for (j in (i+1):n) {
    d <- abs(pfpr_df$pfpr_mean[i] - pfpr_df$pfpr_mean[j])
    mat[i,j] <- d; mat[j,i] <- d
  }
  mat
}
