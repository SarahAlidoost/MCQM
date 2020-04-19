get_condQuantile <- function(multiVar, P, cop) {
  #' Get conditional quantiles.
  condQ <- rep(NA, nrow(multiVar))
  
  for (i in 1:nrow(multiVar)) {
    condVar <- multiVar[i, ]
    condSecVine <- get_condPDF_density(condVar, cop)
    p <- P[i]
    
    xVals <- attr(condSecVine, "xVals")
    density <- condSecVine(xVals)
    nx <- length(xVals)
    int <- cumsum(c(0, diff(xVals) * (0.5 * diff(density) + density[-nx])))
    lower <- max(which(int <= p))
    m <-
      (density[lower + 1] - density[lower]) / (xVals[lower + 1] - xVals[lower])
    b <- density[lower]
    xRes <- -b / m + sign(m) * sqrt(b ^ 2 / m ^ 2 + 2 * (p - int[lower]) /
                                      m)
    condQ[i] <- xVals[lower] + xRes
  }
  return(condQ)
}