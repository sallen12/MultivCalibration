##### Gaussian random field simulation study



### set up

library(geoR)
library(ggplot2)
library(Rcpp)

source("utility_funcs.R")
sourceCpp("variogram_func.cpp") # not used in vg1_rank

d <- 2500
sig2 <- 1
phi <- 3
n <- 10000
M <- 20

nbhd_size <- 1 # neighbourhood size for variogram prerank
coords <- expand.grid(1:sqrt(d), 1:sqrt(d))
w_mat <- sapply(1:d, function(i) sapply(1:d, function(j) as.numeric(max(abs(coords[i, ] - coords[j, ])) <= nbhd_size)))

### observations

y <- grf(d, grid = "reg", cov.pars = c(sig2, phi), nsim = n)
coords <- y[[1]]
y <- t(y$data)

plot_grf(coords, y[100, ])


### forecasts

## Type 1: errors in the mean

# a
mu_x <- rep(-0.5, d)
x <- grf(d, grid = "reg", cov.pars = c(sig2, phi), mean = mu_x, nsim = n*M)
x <- aperm(array(x$data, c(d, M, n)), c(3, 1, 2))
plot_rank(y, x, w_mat, ylab = T, fignum = "2a")

# b
mu_x <- rep(0.5, d)
x <- grf(d, grid = "reg", cov.pars = c(sig2, phi), mean = mu_x, nsim = n*M)
x <- aperm(array(x$data, c(d, M, n)), c(3, 1, 2))
plot_rank(y, x, w_mat, fignum = "2b")


## Type 2: errors in the variance

# a
sig2_x <- 0.65
x <- grf(d, grid = "reg", cov.pars = c(sig2_x, phi), nsim = n*M)
x <- aperm(array(x$data, c(d, M, n)), c(3, 1, 2))
plot_rank(y, x, w_mat, fignum = "2c")

# b
sig2_x <- 1.35
x <- grf(d, grid = "reg", cov.pars = c(sig2_x, phi), nsim = n*M)
x <- aperm(array(x$data, c(d, M, n)), c(3, 1, 2))
plot_rank(y, x, w_mat, fignum = "2d")


## Type 3: errors in the correlation

# a
phi_x <- 1.5
x <- grf(d, grid = "reg", cov.pars = c(sig2, phi_x), nsim = n*M)
x <- aperm(array(x$data, c(d, M, n)), c(3, 1, 2))
plot_rank(y, x, w_mat, fignum = "2e")

# b
phi_x <- 5
x <- grf(d, grid = "reg", cov.pars = c(sig2, phi_x), nsim = n*M)
x <- aperm(array(x$data, c(d, M, n)), c(3, 1, 2))
plot_rank(y, x, w_mat, fignum = "2f")

