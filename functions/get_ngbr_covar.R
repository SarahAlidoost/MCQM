get_ngbr_covar <- function(n_ngbr,
                           var_name,
                           grd_obs,
                           grd_prd,
                           prediction) {
  if (prediction) {
    ngbr <-
      getNeighbours(
        grd_obs,
        grd_prd,
        size = c(n_ngbr + 1),
        var = var_name,
        prediction = prediction,
        min.dist = 1
      )
  } else{
    ngbr <-
      getNeighbours(
        grd_obs,
        size = c(n_ngbr + 1),
        var = var_name,
        prediction = prediction,
        min.dist = 1
      )
  }
  return (as.matrix(ngbr@data[, -c(1)]))
}