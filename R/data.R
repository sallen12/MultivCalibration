#' Gridded wind speed data
#'
#' @description
#' A subset of gridded 10m wind speed forecasts and observations from the EUMETNET's
#' post-processing benchmark dataset.
#'
#' @format
#' An object of type list containing 6 elements:
#' \describe{
#'   \item{obs}{Array of dimension (1045, 33, 32) containing the observed
#'   wind speed for each date and grid point.}
#'   \item{fc_ifs}{Array of dimension (1045, 33, 32, 11) containing the 11 IFS
#'   ensemble members at each date and grid point.}
#'   \item{fc_ecc}{Array of dimension (1045, 33, 32, 11) containing the 11
#'   ensemble members obtained using ensemble copula coupling at each date
#'   and grid point.}
#'   \item{fc_ss}{Array of dimension (1045, 33, 32, 11) containing the 11
#'   ensemble members obtained using the Schaake shuffle at each date
#'   and grid point.}
#'   \item{coord}{Longitude and Latitude coordinates of the 1056 grid points.}
#'   \item{valtime}{Vector of dates for which the forecasts and observations
#'   are available.}
#' }
#'
#' @details
#' Forecasts and observations are available on a grid over central Europe comprised
#' of 33 longitudes and 32 latitudes. These coordinates are available in \code{coord},
#' with the first column containing the longitudes and the second column containing
#' the latitudes.
#'
#' The forecasts are reforecasts from ECMWF's Integrated Forecasting System (IFS)
#' with 11 ensemble members. Reforecasts are generated for the five years prior to
#' each Monday and Thursday in 2017 and 2018. There are 209 Mondays and Thursdays
#' in 2017 and 2018, resulting in 5*209=1045 forecast cases.
#'
#' The observations are ERA5 reanalysis fields on the same gridded domain.
#'
#' Ensemble forecasts are also provided that correspond to the output of two
#' statistical post-processing methods applied to the IFS reforecasts. Both
#' methods use an ensemble model output statistics framework to post-process the
#' wind speed forecasts at each grid point - it is assumed that the wind speed
#' follows a truncated logistic distribution. Quantiles from the truncated logistic
#' distribution are then reordered using a relevant dependence template. Forecasts
#' are provided for two dependence templates, obtained using ensemble copula coupling
#' (ECC) and the Schaake shuffle (SS).
#'
#' @usage data("wind_dat")
#'
#' @source https://github.com/EUPP-benchmark/climetlab-eumetnet-postprocessing-benchmark
#'
"wind_dat"
