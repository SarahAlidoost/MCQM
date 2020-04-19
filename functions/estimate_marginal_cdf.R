estimate_marginal_cdf <- function(x) {
  #' Get marginal distributions or CDF.
  #' A monotone cubic spline is used to make it continous.
  Rnx <- rank(x) / (length(x) + 1)
  return(splinefun(x, Rnx, method = "monoH.FC"))
}