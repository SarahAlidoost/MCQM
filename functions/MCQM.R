MCQM <- function (u_covar1,
                  u_covar2,
                  u_var,
                  Cop_Fun1,
                  Cop_Fun2) {
  if (Cop_Fun1@dimension == 2) {
    #' Bivariate copula quantile mapping.
    p <- dduCopula(cbind(u_covar1, u_var), Cop_Fun1)
    id <- which(p > 1 | p < 0)
    if (length(id) > 0)
      (p[id] <- max(p[-c(id)]))
    v_var = invdduCopula(u = u_covar2, Cop_Fun2, y = p)
  } else {
    #' Multivariate copula quantile mapping.
    p <- get_condCDF_probability(u_var, u_covar1, Cop_Fun1)
    id <- which(p > 1)
    if (length(id) > 0)
      (p[id] <- c(0.999))
    id <- which(p < 0)
    if (length(id) > 0)
      (p[id] <- c(0.001))
    v_var = get_condQuantile(u_covar2, p, Cop_Fun2)
  }
  return(v_var)
}