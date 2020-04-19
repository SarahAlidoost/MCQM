estimate_Rvm <- function (vindata) {
  #' Estimate multivariate copula.
  vineDim <- ncol(vindata)
  vineFit <- vineCopula(as.integer(vineDim))
  rvm <-
    RVineCopSelect(vindata, familyset = c(1:5), vineFit@RVM$Matrix)
  return(vineCopula(rvm))
}