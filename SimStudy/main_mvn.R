##### Multivariate normal simulation study



### set up

library(MASS)
library(ggplot2)

source("utility_funcs.R")
Rcpp::sourceCpp("variogram_func.cpp")

d <- 10
sig2 <- 1
phi <- 3
n <- 10000
M <- 20

lag <- 1 # lag for variogram prerank
w_mat <- matrix(as.numeric(abs(outer(1:d, 1:d, FUN = "-")) <= lag), nrow = d)

### observations

mu_y <- rep(0, d)
Sig_y <- matrix(NA, nrow = d, ncol = d)
for (i in 1:d) {
  for (j in 1:d) {
    Sig_y[i, j] <- sig2*exp(-abs(i-j)/phi)
  }
}  

y <- mvrnorm(n, mu = mu_y, Sigma = Sig_y)


### forecasts

## Type 1: errors in the mean

# a
mu_x <- rep(-0.5, d)
x <- replicate(M, mvrnorm(n, mu = mu_x, Sigma = Sig_y))
plot_rank(y, x, w_mat, ylab = T, fignum = "1a")

# b
mu_x <- rep(0.5, d)
x <- replicate(M, mvrnorm(n, mu = mu_x, Sigma = Sig_y))
plot_rank(y, x, w_mat, fignum = "1b")


## Type 2: errors in the variance

# a
sig2_x <- 0.65
Sig_x <- Sig_y*sig2_x
x <- replicate(M, mvrnorm(n, mu = mu_y, Sigma = Sig_x))
plot_rank(y, x, w_mat, fignum = "1c")

# b
sig2_x <- 1.35
Sig_x <- Sig_y*sig2_x
x <- replicate(M, mvrnorm(n, mu = mu_y, Sigma = Sig_x))
plot_rank(y, x, w_mat, fignum = "1d")


## Type 3: errors in the correlation

# a
phi_x <- 1.5
Sig_x <- matrix(NA, nrow = d, ncol = d)
for (i in 1:d) {
  for (j in 1:d) {
    Sig_x[i, j] <- sig2*exp(-abs(i-j)/phi_x)
  }
}  
x <- replicate(M, mvrnorm(n, mu = mu_y, Sigma = Sig_x))
plot_rank(y, x, w_mat, fignum = "1e")

# b
phi_x <- 5
Sig_x <- matrix(NA, nrow = d, ncol = d)
for (i in 1:d) {
  for (j in 1:d) {
    Sig_x[i, j] <- sig2*exp(-abs(i-j)/phi_x)
  }
}  
x <- replicate(M, mvrnorm(n, mu = mu_y, Sigma = Sig_x))
plot_rank(y, x, w_mat, fignum = "1f")

