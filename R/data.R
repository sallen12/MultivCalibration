#' Gridded wind speed data
#'
#' A subset of gridded 10m wind speed forecasts and observations from the EUMETNET's
#' post-processing benchmark dataset. Forecasts and observations are available on a
#' grid over central Europe comprised of 33 longitudes and 32 latitudes.
#'
#' The forecasts are reforecasts from ECMWF's Integrated Forecasting System (IFS)
#' with 11 ensemble members. Reforecasts are generated for the five years prior to
#' each Monday and Thursday in 2017 and 2018. There are 209 Mondays and Thursdays
#' in 2017 and 2018, resulting in 5*209=1045 forecast cases.
#'
#' @usage data("wind_dat")
#'
#' @format ## `wind_dat`
#'
#' A list with 5 elements:
#' \describe{
#'   \item{fc}{Array of dimensions (33, 32, 5, 11, 209) containing the (re)forecast
#'   wind speed for each longitude, latitude, year, ensemble member and time.}
#'   \item{obs}{Array of dimensions (33, 32, 5, 209) containing the observed
#'   wind speed for each longitude, latitude, year, and time.}
#'   \item{lats}{Vector of latitudes of the grid points.}
#'   \item{lons}{Vector of longitudes of the grid points.}
#'   \item{times}{Days in the two year period for which the reforecasts are created.}
#' }
#' @source https://github.com/EUPP-benchmark/climetlab-eumetnet-postprocessing-benchmark
#'
#' @examples
#' data("wind_dat")
"wind_dat"