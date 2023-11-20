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
#' @seealso \link{preranks} for pre-rank functions for multivariate forecasts
#'
#' @details
#' The argument \code{prerank} specifies which pre-rank function to use. This can
#' either be a string corresponding to one of several in-built options, or it can
#' be a user-specified function.
#'
#' The in-built pre-rank functions currently available are the multivariate rank
#' (\code{prerank = "multivariate_rank"}), the average rank (\code{"average_rank"}),
#' the band-depth rank (\code{"band_depth"}), the mean (\code{"mean"}), the variance
#' (\code{"variance"}), the energy score (\code{"energy_score"}), the
#' fraction of threshold exceedances (\code{"FTE"}), the variogram (\code{"variogram"}),
#' and the isotropy of the variogram (\code{"isotropy"}). See references for details.
#'
#' If \code{prerank} is a function, it should convert a matrix to
#' a single numeric value. Checks are in place to ensure this is satisfied. The
#' \code{prerank} function could also take additional inputs, in which case these
#' inputs should be included as additional arguments in \code{get_prerank_gr()}.
#' See examples below.
#'
#' The FTE pre-rank requires a threshold parameter \code{t}, which must be a single real value.
#' The isotropy pre-rank function requires that a single integer \code{h} is provided,
#' denoting the lag at which the variograms are compared.
#'
#' The variogram pre-rank function requires either an additional argument \code{h}
#' that specifies the multivariate lag at which to calculate the variogram, or an argument
#' \code{w} which is an array of non-negative weights assigned to the different pairs of
#' dimensions. The lag \code{h} should be a vector containing two integers, or a matrix with two
#' columns, in which case the variogram is calculated for each row of the given matrix
#' and the sum of the variograms is returned.
#'
#' Assuming \code{y} is a matrix of dimension (d1, d2), the weight matrix \code{w} must be
#' an array with dimensions (d1, d2, d1, d2). The element in the (i, j, k, l) entry contains
#' the weight between the grid point (i, j) and the grid point (k, l). If both \code{h}
#' and \code{w} are provided, then only \code{h} will be used.
#'
#' The variogram and isotropy pre-rank functions additionally require an argument \code{p} specifying
#' the exponent to be used in the variogram. The variogram also takes an input \code{std}
#' specifying whether the variogram should be standardised by the variance across the dimensions.
#' The exponent \code{p} is a positive real number, whereas \code{std} is a logical.
#' The default is \code{p = 2} and \code{std = TRUE}.
#'
#' @examples
#' d1 <- 5
#' d2 <- 6
#' M <- 10
#'
#' y <- matrix(rnorm(d1*d2), nrow = d1, ncol = d2)
#' x <- array(rnorm(d1*d2*M), c(d1, d2, M))
#'
#' # create array with weights that are inversely proportional to the distance
#' w <- array(NA, c(d1, d2, d1, d2))
#' for (i in 1:d1) {
#'   for (j in 1:d2) {
#'     w[i, j, , ] <- outer(1:d1, 1:d2, FUN = function(k, l) 1 / sqrt((i - k)^2 + (j - l)^2))
#'   }
#' }
#' w[is.infinite(w)] <- 0
#'
#' get_prerank_gr(y, x, prerank = "average_rank", return_rank = FALSE)
#' get_prerank(as.vector(y), matrix(x, nrow = d1*d2, ncol = M),
#'             prerank = "average_rank", return_rank = FALSE)
#'
#' get_prerank_gr(y, x, prerank = "variance", return_rank = FALSE)
#' get_prerank_gr(y, x, prerank = "variance")
#'
#' get_prerank_gr(y, x, prerank = "variogram", h = c(0, 1), return_rank = FALSE)
#' get_prerank_gr(y, x, prerank = "variogram", h = c(0, -1), return_rank = FALSE)
#' get_prerank_gr(y, x, prerank = "variogram", h = c(0, 1), std = FALSE, return_rank = FALSE)
#' get_prerank_gr(y, x, prerank = "variogram", h = rbind(c(0, 1), c(1, 0)))
#' get_prerank_gr(y, x, prerank = "variogram", w = w)
#'
#' get_prerank_gr(y, x, prerank = "isotropy", return_rank = FALSE)
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
    vg_rank_mat(y, x, return_rank, h = varargs$h, w = varargs$w, p = varargs$p, std = varargs$std)
  } else if (prerank == "isotropy") {
    iso_rank(y, x, return_rank, h = varargs$h, p = varargs$p)
  }
}

# variogram
vg_rank_mat <- function(y, x, h = NULL, w = NULL, p = 2, std = TRUE, return_rank = TRUE) {
  if (is.null(h)) {
    custom_rank_gr(y, x, prerank = vario_func_mat_w, return_rank, w = w, p = p, std = std)
  } else {
    custom_rank_gr(y, x, prerank = vario_func_mat_h, return_rank, h = h, p = p, std = std)
  }
}

# variogram helper function - lag
vario_func_mat_h <- function(x, h, p, std) {
  if (is.matrix(h)) {
    g_x <- -sum(apply(h, 1, vario_mat, y = x, p = p))
  } else {
    g_x <- -vario_mat(x, h, p)
  }
  if (std) g_x <- g_x/var(as.vector(x))
  return(g_x)
}

# variogram helper function - weight matrix
vario_func_mat_w <- function(x, w, p, std) {
  x <- as.vector(x)
  w <- matrix(w, nrow = length(x), ncol = length(x))
  g_x <- -vario(x, w, p)
  if (std) g_x <- g_x/var(x)
  return(g_x)
}

# isotropy
iso_rank <- function(y, x, h, p, return_rank = TRUE) {
  custom_rank_gr(y, x, prerank = iso_func, return_rank, h = h, p = p)
}

# isotropy helper function
iso_func <- function(x, h, p) {
  if (length(h) > 1) {
    g_x <- sapply(h, iso_func, x = x, p = p)
    g_x <- sum(g_x)
  } else {
    h <- h*rbind(c(0, 1), c(1, 0), c(1, 1), c(-1, 1))
    g_x_vec <- apply(h, 1, vario_mat, y = x, p = p)
    w_hv <- 1/(2*(g_x_vec[1]^2)/(nrow(x)*(ncol(x) - h)) + 2*(g_x_vec[2]^2)/(ncol(x)*(nrow(x) - h)))
    w_di <- 1/(2*(g_x_vec[3]^2 + g_x_vec[4]^2)/((ncol(x) - h)*(nrow(x) - h)))
    g_x_hv <- w_hv*(g_x_vec[1] - g_x_vec[2])^2
    g_x_di <- w_di*(g_x_vec[3] - g_x_vec[4])^2
    g_x <- -max(g_x_hv, g_x_di)
  }
  return(g_x)
}

# custom pre-rank function
custom_rank_gr <- function(y, x, prerank, return_rank = TRUE, ...) {
  g_y <- prerank(y, ...)
  g_x <- apply(x, 3, prerank, ...)
  rho <- c(g_y, g_x)
  names(rho) <- c("obs", sprintf("ens%d", 1:dim(x)[3]))
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
      if (is.null(varargs$h) && is.null(varargs$w)) stop("The variogram pre-rank function requires an additional argument 'h' or 'w'")

      if (!is.null(varargs$h)) {
        if (!is.numeric(varargs$h)) stop("'h' is not numeric")
        if (!is.vector(varargs$h) && !is.matrix(varargs$h)) stop("'h' must be either a vector or a matrix")
        if (is.vector(varargs$h)) {
          if (length(varargs$h) != 2) stop("'h' must be a vector of length 2")
        } else if (is.matrix(varargs$h)) {
          if (dim(varargs$h)[2] != 2) stop("If 'h' is a matrix, it must have 2 columns")
        }
      }

      if (!is.null(varargs$w)) {
        if (!is.array(varargs$w)) stop("'w' is not an array")
        if (!is.numeric(varargs$w)) stop("'w' must be a numeric")
        if (any(varargs$w < 0)) stop("The weight matrix 'w' contans negative entries")
        if (length(dim(varargs$w)) != 4) stop("'w' must be an array with four dimensions")
        if ((dim(varargs$w)[1] != nrow(y)) || dim(varargs$w)[3] != nrow(y))
          stop("The 1st and 3rd dimensions of 'w' should match the number of rows in 'y'")
        if ((dim(varargs$w)[2] != ncol(y)) || dim(varargs$w)[4] != ncol(y))
          stop("The 2nd and 4th dimensions of 'w' should match the number of columns in 'y'")
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
        if (!is.vector(varargs$h)) stop("'h' must be a numeric vector")
        if (any(varargs$h != as.integer(varargs$h)))
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
