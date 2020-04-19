#' This script is used to implement
#' multivariate copula quantile mapping (MCQM) for bias correction
#' Variables are:
#' obs: data and locations of observations
#' prd: data and locations of predictions
#' ws:  weather station (or x)
#' era: era_interim (or y)
#' dem: digital elevation model (or z)


# Set directories
library_path = "PATH_TO/R/win-library/3.5"
functions_path = "PATH_TO/MCQM/functions"
data_path = "PATH_TO/MCQM/example_data"
output_path = "PATH_TO/MCQM"


# Load libraries
.libPaths(library_path)
library(sp)
library(gstat)
library(VineCopula)
library(copula)
library(spcopula)
# Remove local variables
rm(library_path)


# Source functions
for (file_name in list.files(functions_path, pattern = "[.][R]")) {
  source(file.path(functions_path, file_name))
}
# Remove local variables
rm(functions_path, file_name)


# Load input
file_name <- list.files(data_path, pattern = "[.][RData]")
load(file.path(data_path, file_name))
# Remove local variables
rm(data_path, file_name)


# Estimate marginal distributions
#' marginal_cdf: u = F(x)
#' inv_marginal_cdf: x = F^-1(u)
inv_marginal_cdf_x <-
  estimate_inv_marginal_cdf(grd_obs@data$ws_obs)
marginal_cdf_x <-  estimate_marginal_cdf(grd_obs@data$ws_obs)
grd_obs@data$u_obs <- marginal_cdf_x(grd_obs@data$ws_obs)

#' marginal_cdf: v = F(y)
marginal_cdf_y <-
  estimate_marginal_cdf(c(grd_obs@data$era_obs, grd_prd@data$era_prd))
grd_obs@data$v_obs <- marginal_cdf_y(grd_obs@data$era_obs)
grd_prd@data$v_prd <- marginal_cdf_y(grd_prd@data$era_prd)

#' marginal_cdf: w = F(z)
marginal_cdf_z <-
  estimate_marginal_cdf(c(grd_obs@data$dem_obs, grd_prd@data$dem_prd))
grd_obs@data$w_obs <- marginal_cdf_z(grd_obs@data$dem_obs)
grd_prd@data$w_prd <- marginal_cdf_z(grd_prd@data$dem_prd)
# Remove local variables
rm(marginal_cdf_x, marginal_cdf_y, marginal_cdf_z)


# Run QM
#' One dimentional Quantile Mapping.
#' u = F(x)
#' v = F(y)
#' x = F-1(v)
grd_prd@data$QM_x <- inv_marginal_cdf_x(grd_prd@data$v_prd)


# Run MCQM-I
#' Covariate = height from DEM (or z).
#' u = F(x)
#' v = F(y)
#' w = F(z)
#' u = C^-1( C(v|w) |w)
cop_wv <-
  estimate_copula(grd_obs@data$w_obs, grd_obs@data$v_obs)
cop_wu <-
  estimate_copula(grd_obs@data$w_obs, grd_obs@data$u_obs)
MCQM1_u <- MCQM(grd_prd@data$w_prd,
                grd_prd@data$w_prd,
                grd_prd@data$v_prd,
                cop_wv,
                cop_wu)
grd_prd@data$MCQM1_x <- inv_marginal_cdf_x(MCQM1_u)
# Remove local variables
rm(MCQM1_u, cop_wv, cop_wu)


# Get covariates to be used in MCQM-II and MCQM-III
#' Covariate = one nearest neighbours (or ui).
grd_obs@data$ui_obs <-
  get_ngbr_covar(1, "u_obs", grd_obs, grd_prd, FALSE)
grd_obs@data$vi_obs <-
  get_ngbr_covar(1, "v_obs", grd_obs, grd_prd, FALSE)

grd_prd@data$ui_prd <-
  get_ngbr_covar(1, "u_obs", grd_obs, grd_prd, TRUE)
grd_prd@data$vi_prd <-
  get_ngbr_covar(1, "v_obs", grd_obs, grd_prd, TRUE)


# Run MCQM-II
#' Covariate = nearest neighbours (or -i).
#' u = F(x)
#' v = F(y)
#' u-i = F(x-i)
#' u = C^-1( C(v|v-i) |u-i)
cop_viv <-
  estimate_copula(grd_obs@data$vi_obs, grd_obs@data$v_obs)
cop_uiu <-
  estimate_copula(grd_obs@data$ui_obs, grd_obs@data$u_obs)
MCQM2_u <- MCQM(grd_prd@data$vi_prd,
                grd_prd@data$ui_prd,
                grd_prd@data$v_prd,
                cop_viv,
                cop_uiu)
grd_prd@data$MCQM2_x <- inv_marginal_cdf_x(MCQM2_u)
# Remove local variables
rm(MCQM2_u, cop_viv, cop_uiu)


# Run MCQM-III.
#' Covariates = nearest neighbours (or -i) and height (or z).
#' u = F(x)
#' v = F(y)
#' w = F(z)
#' u-i = F(x-i)
#' u = C^-1( C(v|v-i,w) |u-i,w)
cop_viwv <- estimate_Rvm(cbind(grd_obs@data$v_obs,
                               grd_obs@data$w_obs,
                               grd_obs@data$vi_obs))
cop_uiwu <- estimate_Rvm(cbind(grd_obs@data$u_obs,
                               grd_obs@data$w_obs,
                               grd_obs@data$ui_obs))
MCQM3_u <- MCQM(
  cbind(grd_prd@data$w_prd, grd_prd@data$vi_prd),
  cbind(grd_prd@data$w_prd, grd_prd@data$ui_prd),
  grd_prd@data$v_prd,
  cop_viwv,
  cop_uiwu
)
grd_prd@data$MCQM3_x <- inv_marginal_cdf_x(MCQM3_u)
# Remove local variables
rm(inv_marginal_cdf_x, MCQM3_u, cop_viwv, cop_uiwu)


# Save output
save.image(paste(output_path, "/MCQM_output.RData", sep = ""))
