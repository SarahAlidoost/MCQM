estimate_inv_marginal_cdf <- function(x) {
  #' Get inverse of marginal distributions or invCDF.
  #' A monotone cubic spline is used to make it continous.
  Rnx <- rank(x) / (length(x) + 1)
  return(splinefun(Rnx, x, method = "monoH.FC"))
}