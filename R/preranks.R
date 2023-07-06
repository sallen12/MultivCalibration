#' Pre-rank functions for multivariate calibration
#'
#' Calulate pre-ranks of multivariate forecasts and observations
#'
#' @param y multivariate observation (numeric vector of length d).
#' @param dat samples from multivariate forecast distribution (numeric matrix with d rows).
#' @param prerank the pre-rank function to be used. See details for a list of options.
#' @param return_rank logical specifying whether the rank should be returned
#'  (rather than the vector of pre-ranks). Default is TRUE.
#' @param t threshold for the FTE pre-rank function (single numeric value).
#'
#' @return
#' Rank of the pre-rank transformed observation among the forecast sample (if
#' \code{prerank = F}), or a vector of pre-ranks corresponding to the observation
#' and sample members (if \code{prerank = T}).
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
#' Allen, S., Ziegel, J. and D. Ginsbourger (2023):
#' `Assessing the calibration of multivariate probabilistic forecasts'.
#' \emph{arXiv preprint}.
#'
#' @author Sam Allen
#'
#' @details
#' The pre-rank functions currently included are the multivariate rank
#' (\code{prerank = "multivariate_rank}), the average rank (\code{"average_rank"}),
#' the band-depth rank (\code{"band_depth"}), the mean (\code{"mean"}), the variance
#' (\code{"variance"}), the energy score (\code{"energy_score"}), and the
#' fraction of threshold exceedances (\code{fte_rank}). See references for details.
#'
#' Pre-rank functions will also be added for the variogram, isotropy, and
#' minimum spanning tree.
#'
#' @examples
#' d <- 5
#' M <- 10
#'
#' y <- as.vector(mvtnorm::rmvnorm(1, rep(0, d)))
#' dat <- t(mvtnorm::rmvnorm(M, rep(0, d)))
#'
#' get_prerank(y, dat, prerank = "average_rank", return_rank = F)
#' get_prerank(y, dat, prerank = "average_rank")
#'
#' get_prerank(y, dat, prerank = "variance", return_rank = F)
#' get_prerank(y, dat, prerank = "variance")
#'
#' get_prerank(y, dat, prerank = "FTE", t = 0.5, return_rank = F)
#' get_prerank(y, dat, prerank = "FTE", t = 0.5) # ties are resolved at random
#' get_prerank(y, dat, prerank = "FTE", t = 0.5)
#'
#'
#' ### use sapply to apply to several forecast cases
#' n <- 1000
#' y <- t(mvtnorm::rmvnorm(n, rep(0, d)))
#'
#' # calibrated forecasts
#' dat <- array(t(mvtnorm::rmvnorm(n*M, rep(0, d))), c(d, M, n))
#' mvranks <- sapply(1:n, function(i) get_prerank(y[, i], dat[, , i], prerank = "mean"))
#' barplot(table(mvranks)) # observation is equally likely to take each rank
#'
#' # miscalibrated mean
#' dat <- array(t(mvtnorm::rmvnorm(n*M, rep(-0.5, d))), c(d, M, n))
#' mvranks <- sapply(1:n, function(i) get_prerank(y[, i], dat[, , i], prerank = "mean"))
#' barplot(table(mvranks))
#' # forecast's under-estimate the mean, so the observation often has a higher rank
#' # when evaluated using the mean pre-rank function
#'
#' # miscalibrated variance
#' dat <- array(t(mvtnorm::rmvnorm(n*M, sigma = 1.5*diag(d))), c(d, M, n))
#' mvranks <- sapply(1:n, function(i) get_prerank(y[, i], dat[, , i], prerank = "variance"))
#' barplot(table(mvranks)) # observation is more likely to take a higher rank
#' # forecast's under-estimate the variance, so the observation often has a higher variance rank
#' # when evaluated using the variance pre-rank function
#'
#' @name preranks
NULL


#' @rdname preranks
#' @export
get_prerank <- function(y, dat, prerank, return_rank = TRUE, t = NULL) {
  check_inputs(list(y = y, dat = dat, prerank = prerank, t = t))
  if (prerank == "multivariate_rank") {
    mv_rank(y, dat, return_rank)
  } else if (prerank == "average_rank") {
    av_rank(y, dat, return_rank)
  } else if (prerank == "band_depth") {
    bd_rank(y, dat, return_rank)
  } else if (prerank == "mean") {
    mean_rank(y, dat, return_rank)
  } else if (prerank == "variance") {
    var_rank(y, dat, return_rank)
  } else if (prerank == "energy_score") {
    es_rank(y, dat, return_rank)
  } else if (prerank == "FTE") {
    fte_rank(y, dat, t, return_rank)
  }
}

# multivariate rank
mv_rank <- function(y, dat, return_rank = FALSE) {
  S <- cbind(y, dat)
  M <- ncol(S)
  rho <- sapply(1:M, function(i) sum(sapply(1:M, function(j) all(S[, j] <= S[, i]))))
  names(rho) <- c("obs", sprintf("ens%d", 1:ncol(dat)))
  if (return_rank) {
    rank_y <- rank(rho, ties.method = "random")[1]
    return(unname(rank_y))
  } else {
    return(rho)
  }
}

# average rank
av_rank <- function(y, dat, return_rank = FALSE) {
  S <- cbind(y, dat)
  d <- length(y)
  C <- sapply(1:d, function(l) rank(S[l, ]))
  rho <- rowMeans(C)
  names(rho) <- c("obs", sprintf("ens%d", 1:ncol(dat)))
  if (return_rank) {
    rank_y <- rank(rho, ties.method = "random")[1]
    return(unname(rank_y))
  } else {
    return(rho)
  }
}

# band-depth
bd_rank <- function(y, dat, return_rank = FALSE) {
  S <- cbind(y, dat)
  d <- length(y)
  M <- ncol(S)
  c <- sapply(1:d, function(l) rank(S[l, ]))
  rho <- rowMeans((M - c)*(c - 1))
  names(rho) <- c("obs", sprintf("ens%d", 1:ncol(dat)))
  if (return_rank) {
    rank_y <- rank(rho, ties.method = "random")[1]
    return(unname(rank_y))
  } else {
    return(rho)
  }
}

# mean
mean_rank <- function(y, dat, return_rank = FALSE) {
  g_y <- mean(y)
  g_dat <- colMeans(dat)
  rho <- c(g_y, g_dat)
  names(rho) <- c("obs", sprintf("ens%d", 1:ncol(dat)))
  if (return_rank) {
    rank_y <- rank(rho, ties.method = "random")[1]
    return(unname(rank_y))
  } else {
    return(rho)
  }
}

# variance
var_rank <- function(y, dat, return_rank = FALSE) {
  g_y <- var(y)
  g_dat <- apply(dat, 2, var)
  rho <- c(g_y, g_dat)
  names(rho) <- c("obs", sprintf("ens%d", 1:ncol(dat)))
  if (return_rank) {
    rank_y <- rank(rho, ties.method = "random")[1]
    return(unname(rank_y))
  } else {
    return(rho)
  }
}

# energy score
es_rank <- function(y, dat, return_rank = FALSE) {
  g_y <- scoringRules::es_sample(y, dat)
  g_dat <- apply(dat, 2, scoringRules::es_sample, dat = dat)
  rho <- c(g_y, g_dat)
  names(rho) <- c("obs", sprintf("ens%d", 1:ncol(dat)))
  if (return_rank) {
    rank_y <- rank(rho, ties.method = "random")[1]
    return(unname(rank_y))
  } else {
    return(rho)
  }
}

# fraction of threshold exceedances
fte_rank <- function(y, dat, t, return_rank = FALSE) {
  g_y <- sum(y > t)
  g_dat <- colSums(dat > t)
  rho <- c(g_y, g_dat)
  names(rho) <- c("obs", sprintf("ens%d", 1:ncol(dat)))
  if (return_rank) {
    rank_y <- rank(rho, ties.method = "random")[1]
    return(unname(rank_y))
  } else {
    return(rho)
  }
}


################################################################################
# helper functions

# input checks (from scoringRules)
check_inputs <- function(input) {
  admissible_preranks <- c("multivariate_rank", "average_rank", "band_depth",
                           "mean", "variance", "energy_score", "FTE")
  if (!(input$prerank %in% admissible_preranks)) {
    stop(paste("'prerank' must be one of:", paste(admissible_preranks, collapse = ", ")))
  }
  if (input$prerank == "FTE") {
    if (!is.numeric(input$t)) stop("'t' is not numeric")
    if (length(t) > 1) stop("'t' must be a single numeric value")
  }
  if (!is.numeric(input$y)) stop("'y' is not numeric")
  if (!is.numeric(input$dat)) stop("'dat' is not numeric")
  if (!is.vector(input$y)) stop("'y' is not a vector")
  if (!is.matrix(input$dat)) stop("'dat' is not a matrix ")
  if (length(input$y) != dim(input$dat)[1]) stop("Dimensions of 'y' and 'dat' do not match")
}
