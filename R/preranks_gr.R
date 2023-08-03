#' Pre-rank functions to assess multivariate calibration of gridded forecasts
#'
#' Calculate multivariate ranks of gridded forecasts and observations
#'
#' @param y multivariate observation (numeric matrix with d1 columns and d2 rows).
#' @param x samples from multivariate forecast distribution (numeric array with dimension (d1, d2, M)).
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
#' Allen, S., Ziegel, J. and D. Ginsbourger (2023):
#' `Assessing the calibration of multivariate probabilistic forecasts'.
#' \emph{arXiv preprint}. arXiv:2307.05846.
#' \doi{10.48550/arXiv.2307.05846}
#'
#' @author Sam Allen
#'
#' @details
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
#' d1 <- 5
#' d2 <- 6
#' M <- 10
#'
#' y <- matrix(rnorm(d1*d2), nrow = d1, ncol = d2)
#' x <- array(rnorm(d1*d2*M), c(d1, d2, M))
#'
#' get_prerank_gr(y, x, prerank = "average_rank", return_rank = F)
#' get_prerank(as.vector(y), matrix(x, nrow = d1*d2, ncol = M), prerank = "average_rank", return_rank = F)
#'
#' get_prerank_gr(y, x, prerank = "variance", return_rank = F)
#' get_prerank_gr(y, x, prerank = "variance")
#'
#' get_prerank_gr(y, x, prerank = "variogram", h = c(0, 1), return_rank = F)
#' get_prerank_gr(y, x, prerank = "variogram", h = c(0, -1), return_rank = F)
#' get_prerank_gr(y, x, prerank = "variogram", h = c(0, 1), std = F, return_rank = F)
#'
#' get_prerank_gr(y, x, prerank = "isotropy", return_rank = F)
#' get_prerank_gr(y, x, prerank = "isotropy")
#' get_prerank_gr(y, x, prerank = "isotropy", h = 2)
#'
#' @name preranks_gr
NULL

#' @rdname preranks_gr
#' @export
get_prerank_gr <- function(y, x, prerank, return_rank = TRUE, ...) {
  check_inputs_gr(y = y, x = x, prerank = prerank, return_rank, ...)
  varargs <- list(...)
  if (!is.function(prerank)) {
    if (prerank %in% c("multivariate_rank", "average_rank", "band_depth",
                       "mean", "variance", "energy_score", "FTE")) {
      y <- as.vector(y)
      dim(x) <- c(prod(dim(x)[1:2]), dim(x)[3])
    } else if (prerank == "variogram") {
      if (is.null(varargs$p)) varargs$p <- 2
      if (is.null(varargs$std)) varargs$std <- TRUE
    } else if (prerank == "isotropy") {
      if (is.null(varargs$h)) varargs$h <- 1
      if (is.null(varargs$p)) varargs$p <- 2
    }
  }
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
    fte_rank(y, x, return_rank, t = varargs$t)
  } else if (prerank == "variogram") {
    vg_rank_mat(y, x, return_rank, h = varargs$h, p = varargs$p, std = varargs$std)
  } else if (prerank == "isotropy") {
    iso_rank(y, x, return_rank, h = varargs$h)
  }
}

# variogram
vg_rank_mat <- function(y, x, h, p, std, return_rank = TRUE) {
  custom_rank(y, x, prerank = vario_func_mat, return_rank, h = h, p = p, std = std)
}

# variogram helper function
vario_func_mat <- function(x, h, p, std) {
  if (is.matrix(h)) {
    g_x <- -sum(apply(h, 1, vario_mat, y = x, p = p))
  } else {
    g_x <- -vario_mat(x, h, p)
  }
  if (std) g_x <- g_x/var(as.vector(x))
  return(g_x)
}

# isotropy
iso_rank <- function(y, x, h, return_rank = TRUE) {
  custom_rank(y, x, prerank = iso_func, return_rank, h = h, p = p)
}

# isotropy helper function
iso_func <- function(x, h, p) {
  h <- h*rbind(c(0, 1), c(1, 0), c(1, 1), c(-1, 1))
  g_x_vec <- apply(h, 1, vario_mat, y = x, p = p)
  g_x <- ((g_x_vec[1] - g_x_vec[2])/(g_x_vec[1] + g_x_vec[2]))^2 +
    ((g_x_vec[3] - g_x_vec[4])/(g_x_vec[3] + g_x_vec[4]))^2
  return(-g_x)
}

# custom pre-rank function
custom_rank_gr <- function(y, x, prerank, return_rank = TRUE, ...) {
  g_y <- prerank(y, ...)
  g_x <- apply(x, 3, prerank, ...)
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
check_inputs_gr <- function (y, x, prerank, return_rank, ...) {
  if (!is.logical(return_rank)) stop("'return_rank' is not a logical")
  if (!is.numeric(y)) stop("'y' is not numeric")
  if (!is.numeric(x)) stop("'x' is not numeric")
  if (!is.matrix(y)) stop("'y' is not a matrix")
  if (!is.array(x)) stop("'x' is not an array ")
  if (!all.equal(dim(y), dim(x)[1:2])) stop("Dimensions of 'y' and 'x' do not match")

  if (is.function(prerank)) {
    g_y <- prerank(y, ...)
    g_x <- apply(x, 3, prerank, ...)
    if (!is.numeric(g_y) || !is.numeric(g_x)) stop("The pre-rank function returns non-numeric values")
    if (length(g_y) > 1) stop("The pre-rank function does not return a single value")
  } else{
    varargs <- list(...)
    admissible_preranks <- c("multivariate_rank", "average_rank", "band_depth", "mean",
                             "variance", "energy_score", "FTE", "variogram", "isotropy")
    if (!(prerank %in% admissible_preranks)) {
      stop(paste("'prerank' must be one of:", paste(admissible_preranks, collapse = ", ")))
    }
    if (prerank == "FTE") {
      if (is.null(varargs$t)) stop("The FTE pre-rank function requires an additional argument 't'")
      if (!is.numeric(varargs$t)) stop("'t' is not numeric")
      if (length(varargs$t) > 1) stop("'t' must be a single numeric value")
    } else if (prerank == "variogram") {
      if (is.null(varargs$h)) stop("The variogram pre-rank function requires an additional argument 'h'")
      if (!is.numeric(varargs$h)) stop("'h' is not numeric")
      if (!is.vector(varargs$h) || !is.matrix(varargs$h)) stop("'h' must be either a vector or a matrix")
      if (is.vector(varargs$h)) {
        if (length(varargs$h) != 2) stop("'h' must be a vector of length 2")
      } else if (is.matrix(varargs$h)) {
        if (dim(varargs$h)[2] != 2) stop("'h' must be a matrix with 2 columns")
      }

      if (!is.null(varargs$p)) {
        if (!is.numeric(varargs$p)) stop("'p' is not numeric")
        if (length(varargs$p) > 1) stop("'p' must be a single numeric value")
        if (varargs$p <= 0) stop(paste("'p' must be a positive number, got ", varargs$p))
      }

      if (!is.null(varargs$std)) {
        if (!is.logical(varargs$std)) stop("'std' must be a logical")
      }
    } else if (prerank == "isotropy") {
      if (!is.null(varargs$h)) {
        if (!is.numeric(varargs$h)) stop("'h' is not numeric")
        if (length(varargs$h) > 1) stop("'h' must be a single numeric value")
        if (varargs$h - as.integer(varargs$h) != 0)
        message("'h' is not an integer. Using the integer part of 'h' when calculating the isotropy")
      }

      if (!is.null(varargs$p)) {
        if (!is.numeric(varargs$p)) stop("'p' is not numeric")
        if (length(varargs$p) > 1) stop("'p' must be a single numeric value")
        if (varargs$p <= 0) stop(paste("'p' must be a positive number, got ", varargs$p))
      }
    }
  }
}
