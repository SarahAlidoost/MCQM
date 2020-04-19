get_condPDF_density <- function(condVar, cop) {
  #' Get conditional densities.
  n <- 1000
  rat <- 50:1 %x% c(1e-06, 1e-05, 1e-04, 0.001)
  xVals <- unique(sort(c(rat, 1 - rat, 1:(n - 1) / n)))
  nx <- length(xVals)
  repCondVar <- matrix(condVar,
                       ncol = length(condVar),
                       nrow = nx,
                       byrow = T)
  density <- dCopula(cbind(xVals, repCondVar), cop)
  density <- c(max(0, 2 * density[1] - density[2]),
               density,
               max(0, 2 * density[nx] - density[nx - 1]))
  linAppr <- approxfun(c(0, xVals, 1), density)
  int <- sum(diff(c(0, xVals, 1)) * (0.5 * diff(density) +
                                       density[-(nx + 2)]))
  condVineFun <- function(u)
    linAppr(u) / int
  attr(condVineFun, "xVals") <- c(0, xVals, 1)
  
  return(condVineFun)
}