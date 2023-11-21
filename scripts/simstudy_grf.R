##### Gaussian random field simulation study

set.seed(1000)


##### set up

# devtools::install_github("sallen12/WeightedForecastVerification")
library(geoR)
library(ggplot2)
library(WeightedForecastVerification)
library(MultivCalibration)

p <- q <- 30
d <- p*q
sig2 <- 1
phi <- 1
n <- 1000
M <- 20

# create array with weights that decrease exponentially with the distance between grid points
w_mat <- array(NA, c(p, q, p, q))
for (i in 1:p) {
  for (j in 1:q) {
    w_mat[i, j, , ] <- outer(1:p, 1:q, FUN = function(k, l) exp(-sqrt((i - k)^2 + (j - l)^2)))
  }
}


##### utility functions

# function to plot random fields
plot_grf <- function(coords, z) {
  grf_df <- data.frame(x = coords[, 1], y = coords[, 2], z)
  ggplot(grf_df) + geom_raster(aes(x = x, y = y, fill = z)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_distiller(limits = c(-3, 3), palette = "RdBu") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
}

# wrapper to get data frame of ranks corresponding to each pre-rank function
get_ranks <- function(y, x, w_mat, t = 1, h = 1, parallel = TRUE) {
  n <- nrow(y)

  get_ranks_i <- function(y, x, w_mat, t, h) {
    ranks_i <- c(get_prerank_gr(y, x, prerank = "average_rank"),
                 get_prerank_gr(y, x, prerank = "band_depth"),
                 get_prerank_gr(y, x, prerank = "mean"),
                 get_prerank_gr(y, x, prerank = "variance"),
                 get_prerank_gr(y, x, prerank = "variogram", w = w_mat),
                 get_prerank_gr(y, x, prerank = "FTE", t = t),
                 get_prerank_gr(y, x, prerank = "isotropy", h = h))
    names(ranks_i) <- c("avr", "bdr", "mea", "var", "dep", "fte", "iso")
    return(ranks_i)
  }

  if (parallel) {
    library(foreach)
    library(doParallel)
    library(doSNOW)
    cl <- makeCluster(6)
    registerDoSNOW(cl)
    pb <- txtProgressBar(min = 1, max = n, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    output <- foreach(i = 1:n, .combine=rbind, .packages = "MultivCalibration", .options.snow = opts) %dopar% {
      get_ranks_i(y[i, , ], x[i, , , ], w_mat, t, h)
    }
    rank_df <- data.frame(output)
    close(pb)
    stopCluster(cl)
  } else {
    output <- sapply(1:n, function(i) {
      print(i)
      get_ranks_i(y[i, , ], x[i, , , ], w_mat, t, h)})
    rank_df <- data.frame(t(output))
  }

  return(rank_df)
}

# wrapper to plot and save multivariate rank histograms
plot_ranks <- function(rank_df, filename = NULL, ylabs = NULL, ymax = 0.2) {
  if (is.null(ylabs)) ylabs <- rep(" ", ncol(rank_df))

  plot_list <- lapply(1:ncol(rank_df), function(i)
    pit_hist(rank_df[[i]], M + 1, ylab = ylabs[i], ymax = ymax, xlab = NULL, yticks = F, xticks = F))

  plot_list$ncol <- 1
  allplot <- do.call(gridExtra::grid.arrange, plot_list)

  if (!is.null(filename)) {
    width <- 1.6
    height <- 1.4
    ggsave(filename, allplot, width = width, height = height*ncol(rank_df))
  }
}


filedir <- "scripts/fig_"
ylabs <- c("Average rank", "Band-depth", "Mean", "Variance", "Dependence", "FTE", "Isotropy")


##### simulate observations

y <- grf(d, grid = "reg", cov.pars = c(sig2, phi), nsim = n)
coords <- y$coords*(d - 1)
y <- aperm(array(y$data, c(p, q, n)), c(3, 1, 2))
plot_grf(coords, as.vector(y[10, , ]))

# Note: y[1, i, j] returns the (i, j) coordinate of the 1st field


#####  simulate forecasts

### Type 1: errors in the mean

# a
mu_x <- rep(-0.25, d)
x <- grf(d, grid = "reg", cov.pars = c(sig2, phi), mean = mu_x, nsim = n*M)
x <- aperm(array(x$data, c(p, q, M, n)), c(4, 1, 2, 3))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, ylabs = ylabs)

# b
mu_x <- rep(0.25, d)
x <- grf(d, grid = "reg", cov.pars = c(sig2, phi), mean = mu_x, nsim = n*M)
x <- aperm(array(x$data, c(p, q, M, n)), c(4, 1, 2, 3))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, ylabs = ylabs)


### Type 2: errors in the variance

# a
sig2_x <- 0.85
x <- grf(d, grid = "reg", cov.pars = c(sig2_x, phi), nsim = n*M)
x <- aperm(array(x$data, c(p, q, M, n)), c(4, 1, 2, 3))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "2a", ".pdf"), ylabs = ylabs)

# b
sig2_x <- 1.25
x <- grf(d, grid = "reg", cov.pars = c(sig2_x, phi), nsim = n*M)
x <- aperm(array(x$data, c(p, q, M, n)), c(4, 1, 2, 3))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "2b", ".pdf"))


### Type 3: errors in the correlation

# a
phi_x <- 0.5
x <- grf(d, grid = "reg", cov.pars = c(sig2, phi_x), nsim = n*M)
x <- aperm(array(x$data, c(p, q, M, n)), c(4, 1, 2, 3))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "2c", ".pdf"))

# b
phi_x <- 2
x <- grf(d, grid = "reg", cov.pars = c(sig2, phi_x), nsim = n*M)
x <- aperm(array(x$data, c(p, q, M, n)), c(4, 1, 2, 3))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "2d", ".pdf"))


### Type 4: errors in the isotropy

# a
x <- grf(d, grid = "reg", cov.pars = c(sig2, phi), aniso.pars = c(0, 1.1), nsim = n*M)
x <- aperm(array(x$data, c(p, q, M, n)), c(4, 1, 2, 3))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "2e", ".pdf"))

# b
y <- t(grf(d, grid = "reg", cov.pars = c(sig2, phi), aniso.pars = c(0, 1.1), nsim = n)$data)
x <- grf(d, grid = "reg", cov.pars = c(sig2, phi), nsim = n*M)
x <- aperm(array(x$data, c(p, q, M, n)), c(4, 1, 2, 3))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "2f", ".pdf"))


