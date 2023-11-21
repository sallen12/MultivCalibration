##### Multivariate normal simulation study

set.seed(1000)


##### set up

# devtools::install_github("sallen12/WeightedForecastVerification")
library(MASS)
library(WeightedForecastVerification)
library(MultivCalibration)

d <- 10         # number of dimensions
n <- 10000      # number of repetitions
M <- 20         # number of samples in the forecast distributions
sig2 <- 1       # scale parameter in data generating process
phi <- 1        # correlation/range parameter in data generating process

# create matrix with weights that decrease exponentially with the distance between dimensions
w_mat <- exp(-1/abs(outer(1:d, 1:d, FUN = "-")))


##### utility functions

# wrapper to get data frame of ranks corresponding to each pre-rank function
get_ranks <- function(y, x, w_mat) {
  n <- nrow(y)
  rank_df <- data.frame(mvr = sapply(1:n, function(i) get_prerank(y[i, ], x[i, , ], prerank = "multivariate_rank")),
                        avr = sapply(1:n, function(i) get_prerank(y[i, ], x[i, , ], prerank = "average_rank")),
                        bdr = sapply(1:n, function(i) get_prerank(y[i, ], x[i, , ], prerank = "band_depth")),
                        esr = sapply(1:n, function(i) get_prerank(y[i, ], x[i, , ], prerank = "energy_score")),
                        mea = sapply(1:n, function(i) get_prerank(y[i, ], x[i, , ], prerank = "mean")),
                        var = sapply(1:n, function(i) get_prerank(y[i, ], x[i, , ], prerank = "variance")),
                        dep = sapply(1:n, function(i) get_prerank(y[i, ], x[i, , ], w = w_mat, prerank = "variogram")))
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
ylabs <- c("Multivariate", "Average rank", "Band-depth", "Energy score", "Mean", "Variance", "Dependence")


##### simulate observations

mu_y <- rep(0, d)
Sig_y <- sapply(1:d, function(i) sapply(1:d, function(j) sig2*exp(-abs(i-j)/phi)))
y <- mvrnorm(n, mu = mu_y, Sigma = Sig_y)


##### simulate forecasts

### Type 1: errors in the mean

# a
mu_x <- rep(-0.25, d)
x <- replicate(M, mvrnorm(n, mu = mu_x, Sigma = Sig_y))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "1a", ".pdf"), ylabs = ylabs)

# b
mu_x <- rep(0.25, d)
x <- replicate(M, mvrnorm(n, mu = mu_x, Sigma = Sig_y))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "1b", ".pdf"))


### Type 2: errors in the variance

# a
sig2_x <- 0.85
Sig_x <- Sig_y*sig2_x
x <- replicate(M, mvrnorm(n, mu = mu_y, Sigma = Sig_x))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "1c", ".pdf"))

# b
sig2_x <- 1.25
Sig_x <- Sig_y*sig2_x
x <- replicate(M, mvrnorm(n, mu = mu_y, Sigma = Sig_x))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "1d", ".pdf"))


### Type 3: errors in the correlation

# a
phi_x <- 0.5
Sig_x <- sapply(1:d, function(i) sapply(1:d, function(j) sig2*exp(-abs(i-j)/phi_x)))
x <- replicate(M, mvrnorm(n, mu = mu_y, Sigma = Sig_x))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "1e", ".pdf"))

# b
phi_x <- 2
Sig_x <- sapply(1:d, function(i) sapply(1:d, function(j) sig2*exp(-abs(i-j)/phi_x)))
x <- replicate(M, mvrnorm(n, mu = mu_y, Sigma = Sig_x))
rank_df <- get_ranks(y, x, w_mat)
plot_ranks(rank_df, filename = paste0(filedir, "1f", ".pdf"))

