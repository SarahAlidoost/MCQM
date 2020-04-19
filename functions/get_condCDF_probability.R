get_condCDF_probability <- function(Var, multiVar, cop) {
  #' Get Conditional probabilities.
  #' probability is the value of conditional CDF.
  p_cond <- rep(NA, nrow(multiVar))
  for (i in 1:nrow(multiVar)) {
    condVar <- multiVar[i, ]
    condSecVine <- get_condPDF_density(condVar, cop)
    ePred <-
      integrate(function(x)
        (condSecVine(x)),0+.Machine$double.eps,Var[i],subdivisions=10000L)
    p_cond[i] <- ePred$value
  }
  return(p_cond)
}