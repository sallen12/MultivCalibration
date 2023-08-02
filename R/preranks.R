#' Pre-rank functions to assess multivariate calibration
#'
#' Calculate pre-ranks of multivariate forecasts and observations
#'
#' @param y multivariate observation (numeric vector of length d).
#' @param x samples from multivariate forecast distribution (numeric matrix with d rows).
#' @param prerank the pre-rank function to be used. This is either a string from
#'  a list of possible options (see details below), or a function.
#' @param return_rank logical specifying whether the rank should be returned
#'  (rather than the vector of pre-ranks).
#' @param ... additional arguments to the pre-rank function.
#'
#' @return
#' Rank of the pre-rank transformed observation among the forecast sample (if
#' \code{return_rank = T}), or a vector of pre-ranks corresponding to the observation
#' and sample members (if \code{return_rank = F}).
#'
#' @references
#'
#' Gneiting, T., Stanberry, L. I., Grimit, E. P., Held, L. and N. A. Johnson (2008):
#' `Assessing probabilistic forecasts of multivariate quantities, with an application to ensemble predictions of surface winds'.
#' \emph{Test} 17, 211-235.
#' \doi{10.1007/s11749-008-0114-x}
#'
#' Thorarinsdottir, T. L., Scheuerer, M. and C. Heinz (2016):
#' `Assessing the calibration of high-dimensional ensemble forecasts using rank histograms'.
#' \emph{Journal of computational and graphical statistics} 25, 105-122.
#' \doi{10.1080/10618600.2014.977447}
#'
#' Knüppel, M., Krüger, F., & Pohle, M. O. (2022):
#' `Score-based calibration testing for multivariate forecast distributions'.
#' \emph{arXiv preprint}. arXiv:2211.16362.
#' \doi{10.48550/arXiv.2211.16362}
#'
#' Allen, S., Ziegel, J. and D. Ginsbourger (2023):
#' `Assessing the calibration of multivariate probabilistic forecasts'.
#' \emph{arXiv preprint}.
#'
#' @author Sam Allen
#'
#' @details
#' When assessing multivariate calibration, it is common to convert the multivariate
#' forecasts and observations into univariate objects, and then apply univariate
#' methods. In the context of multivariate calibration, the function used to
#' perform this transformation is often called a pre-rank function. The function
#' \code{get_prerank()} can be used to apply this transformation, and
#' to extract the rank of the transformed observation among the transformed
#' samples from the forecast (i.e. ensemble members).
#'
#' The argument \code{prerank} specifies which pre-rank function to use. This can
#' either be a string corresponding to one of several in-built options, or it can
#' be a user-specified function.
#'
#' The in-built pre-rank functions currently available are the multivariate rank
#' (\code{prerank = "multivariate_rank"}), the average rank (\code{"average_rank"}),
#' the band-depth rank (\code{"band_depth"}), the mean (\code{"mean"}), the variance
#' (\code{"variance"}), the energy score (\code{"energy_score"}), and the
#' fraction of threshold exceedances (\code{"fte_rank"}).
#' Pre-rank functions will later be added for the variogram, isotropy, and
#' minimum spanning tree. See references for details.
#'
#' If \code{prerank} is a function, it should convert a vector of dimension d, to
#' a single numeric value. Checks are in place to ensure this is satisfied. The
#' \code{prerank} function could also take additional inputs, in which case these
#' inputs should be included as additional arguments in \code{get_prerank}.
#' See examples below.
#'
#' @examples
#' d <- 5
#' M <- 10
#'
#' y <- as.vector(mvtnorm::rmvnorm(1, rep(0, d)))
#' x <- t(mvtnorm::rmvnorm(M, rep(0, d)))
#'
#' get_prerank(y, x, prerank = "average_rank", return_rank = F)
#' get_prerank(y, x, prerank = "average_rank")
#'
#' get_prerank(y, x, prerank = "variance", return_rank = F)
#' get_prerank(y, x, prerank = "variance")
#'
#' get_prerank(y, x, prerank = "FTE", t = 0.5, return_rank = F)
#' get_prerank(y, x, prerank = "FTE", t = 0.5) # ties are resolved at random
#' get_prerank(y, x, prerank = "FTE", t = 0.5)
#'
#'
#' ### use sapply to apply to several forecast cases
#' n <- 1000
#' y <- t(mvtnorm::rmvnorm(n, rep(0, d)))
#'
#' # calibrated forecasts
#' x <- array(t(mvtnorm::rmvnorm(n*M, rep(0, d))), c(d, M, n))
#' mvranks <- sapply(1:n, function(i) get_prerank(y[, i], x[, , i], prerank = "mean"))
#' barplot(table(mvranks)) # observation is equally likely to take each rank
#'
#' # miscalibrated mean
#' x <- array(t(mvtnorm::rmvnorm(n*M, rep(-0.5, d))), c(d, M, n))
#' mvranks <- sapply(1:n, function(i) get_prerank(y[, i], x[, , i], prerank = "mean"))
#' barplot(table(mvranks))
#' # forecast's under-estimate the mean, so the observation often has a higher rank
#' # when evaluated using the mean pre-rank function
#'
#' # miscalibrated variance
#' x <- array(t(mvtnorm::rmvnorm(n*M, sigma = 0.7*diag(d))), c(d, M, n))
#' mvranks <- sapply(1:n, function(i) get_prerank(y[, i], x[, , i], prerank = "variance"))
#' barplot(table(mvranks)) # observation is more likely to take a higher rank
#' # forecast's under-estimate the variance, so the observation often has a higher variance rank
#' # when evaluated using the variance pre-rank function
#'
#' # custom pre-rank function
#' x <- array(t(mvtnorm::rmvnorm(n*M, rep(0, d))), c(d, M, n))
#' prerank <- function(x) mean((x - mean(x))^3) # function to quantify skew
#' mvranks <- sapply(1:n, function(i) get_prerank(y[, i], x[, , i], prerank = prerank))
#' barplot(table(mvranks))
#'
#' prerank <- function(x, q) mean((x - mean(x))^q) # function to quantify q-th centered moment
#' mvranks <- sapply(1:n, function(i) get_prerank(y[, i], dat[, , i], prerank = prerank, q = 3))
#' barplot(table(mvranks))
#'
#' @name preranks
NULL


#' @rdname preranks
#' @export
get_prerank <- function(y, x, prerank, return_rank = TRUE, ...) {
  check_inputs(y = y, x = x, prerank = prerank, ...)
  if (is.function(prerank)) {
    custom_rank(y, x, prerank, return_rank, ...)
  } else if (prerank == "multivariate_rank") {
    mv_rank(y, x, return_rank)
  } else if (prerank == "average_rank") {
    av_rank(y, x, return_rank)
  } else if (prerank == "band_depth") {
    bd_rank(y, x, return_rank)
  } else if (prerank == "mean") {
    mean_rank(y, x, return_rank)
  } else if (prerank == "variance") {
    var_rank(y, x, return_rank)
  } else if (prerank == "energy_score") {
    es_rank(y, x, return_rank)
  } else if (prerank == "FTE") {
    fte_rank(y, x, t, return_rank)
  } else if (prerank == "variogram") {
    vg_rank(y, x, w, p, std, return_rank)
  }
}

# multivariate rank
mv_rank <- function(y, x, return_rank = TRUE) {
  S <- cbind(y, x)
  M <- ncol(S)
  rho <- sapply(1:M, function(i) sum(sapply(1:M, function(j) all(S[, j] <= S[, i]))))
  names(rho) <- c("obs", sprintf("ens%d", 1:ncol(x)))
  if (return_rank) {
    rank_y <- rank(rho, ties.method = "random")[1]
    return(unname(rank_y))
  } else {
    return(rho)
  }
}

# average rank
av_rank <- function(y, x, return_rank = TRUE) {
  S <- cbind(y, x)
  d <- length(y)
  C <- sapply(1:d, function(l) rank(S[l, ]))
  rho <- rowMeans(C)
  names(rho) <- c("obs", sprintf("ens%d", 1:ncol(x)))
  if (return_rank) {
    rank_y <- rank(rho, ties.method = "random")[1]
    return(unname(rank_y))
  } else {
    return(rho)
  }
}

# band-depth
bd_rank <- function(y, x, return_rank = TRUE) {
  S <- cbind(y, x)
  d <- length(y)
  M <- ncol(S)
  c <- sapply(1:d, function(l) rank(S[l, ]))
  rho <- rowMeans((M - c)*(c - 1))
  names(rho) <- c("obs", sprintf("ens%d", 1:ncol(x)))
  if (return_rank) {
    rank_y <- rank(rho, ties.method = "random")[1]
    return(unname(rank_y))
  } else {
    return(rho)
  }
}

# mean
mean_rank <- function(y, x, return_rank = TRUE) {
  custom_rank(y, x, prerank = mean, return_rank)
}

# variance
var_rank <- function(y, x, return_rank = TRUE) {
  custom_rank(y, x, prerank = var, return_rank)
}

# energy score
es_rank <- function(y, x, return_rank = TRUE) {
  custom_rank(y, x, prerank = scoringRules::es_sample, return_rank, dat = x)
}

# fraction of threshold exceedances
fte_rank <- function(y, x, t, return_rank = TRUE) {
  custom_rank(y, x, prerank = sum, return_rank, t = t)
}

# custom pre-rank function
custom_rank <- function(y, x, prerank, return_rank = TRUE, ...) {
  g_y <- prerank(y, ...)
  g_x <- apply(x, 2, prerank, ...)
  rho <- c(g_y, g_x)
  names(rho) <- c("obs", sprintf("ens%d", 1:ncol(x)))
  if (return_rank) {
    rank_y <- rank(rho, ties.method = "random")[1]
    return(unname(rank_y))
  } else {
    return(rho)
  }
}

################################################################################
# helper functions

# input checks (adapted from scoringRules)
check_inputs <- function (y, x, prerank, ...) {
  if (!is.numeric(y)) stop("'y' is not numeric")
  if (!is.numeric(x)) stop("'x' is not numeric")
  if (!is.vector(y)) stop("'y' is not a vector")
  if (!is.matrix(x)) stop("'x' is not a matrix ")
  if (length(y) != dim(x)[1]) stop("Dimensions of 'y' and 'x' do not match")

  if (is.function(prerank)) {
    g_y <- prerank(y, ...)
    g_x <- apply(x, 2, prerank, ...)
    if (!is.numeric(g_y) || !is.numeric(g_x)) stop("The pre-rank function returns non-numeric values")
    if (length(g_y) > 1) stop("The pre-rank function does not return a single value")
  } else{
    admissible_preranks <- c("multivariate_rank", "average_rank", "band_depth",
                             "mean", "variance", "energy_score", "FTE")
    if (!(prerank %in% admissible_preranks)) {
      stop(paste("'prerank' must be one of:", paste(admissible_preranks, collapse = ", ")))
    }
    if (prerank == "FTE") {
      if (is.null(t)) stop("The FTE pre-rank function requires an additional argument 't'")
      if (!is.numeric(t)) stop("'t' is not numeric")
      if (length(t) > 1) stop("'t' must be a single numeric value")
    }
  }
}
