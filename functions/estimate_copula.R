estimate_copula <- function (u1, u2) {
  #' estimate bivariate copula.
  #' BiCopSelect: etimate copula.
  #' copulaFromFamilyIndex: get the family.
  #' five families are selected.
  flist <- c(1:5)
  test_cor <- cor(u1, u2, method = c("kendall"))
  if ((test_cor < -0.90) || (test_cor > 0.90))
    (flist <- c(1:3))
  if (test_cor < 0)
    (flist <- flist[-c(3, 4)])
  
  Bicop <- BiCopSelect(u1, u2, familyset = flist, rotations = FALSE)
  cop <- copulaFromFamilyIndex(Bicop$family, Bicop$par, Bicop$par2)
  return(cop)
}