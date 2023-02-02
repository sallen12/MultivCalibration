##### Multivariate normal simulation study

### set up

library(MASS)
library(ggplot2)
library(Rcpp)

source("~/R/MultivCal/SimStudy/prerank_funcs.R")
sourceCpp("~/R/MultivCal/SimStudy/variogram_func.cpp") # not used in vg1_rank

# function to plot and save results
plot_rank <- function(y, x, ylab = F, fignum = NULL) {
  
  rank_df <- data.frame(mvr = sapply(1:n, function(i) mv_rank(y[i, ], x[i, , ])),
                        avr = sapply(1:n, function(i) av_rank(y[i, ], x[i, , ])),
                        bdr = sapply(1:n, function(i) bd_rank(y[i, ], x[i, , ])),
                        esr = sapply(1:n, function(i) es_rank(y[i, ], x[i, , ])),
                        loc = sapply(1:n, function(i) mean_rank(y[i, ], x[i, , ])),
                        sc = sapply(1:n, function(i) var_rank(y[i, ], x[i, , ])),
                        dep = sapply(1:n, function(i) vg1_rank(y[i, ], x[i, , ])))
  
  if (ylab) {
    lab_vec <- c("Multivariate", "Average", "Band-depth", "Energy score", 
                 "Location", "Scale", "Dependence")
  } else {
    lab_vec <- rep(" ", ncol(rank_df))
  }
  
  plot_mvr <- plot_hist(rank_df$mvr, M + 1, ylab = lab_vec[1], axis = F)
  plot_avr <- plot_hist(rank_df$avr, M + 1, ylab = lab_vec[2], axis = F)
  plot_bdr <- plot_hist(rank_df$bdr, M + 1, ylab = lab_vec[3], axis = F)
  plot_esr <- plot_hist(rank_df$esr, M + 1, ylab = lab_vec[4], axis = F)
  plot_loc <- plot_hist(rank_df$loc, M + 1, ylab = lab_vec[5], axis = F)
  plot_sc <- plot_hist(rank_df$sc, M + 1, ylab = lab_vec[6], axis = F)
  plot_dep <- plot_hist(rank_df$dep, M + 1, ylab = lab_vec[7], axis = F)
  
  if (!is.null(fignum)) {
    filedir <- "~/R/MultivCal/SimStudy/Figures/fig_1"
    ggsave(paste(filedir, fignum, "i.pdf", sep=""), plot_mvr, width=2.15, height=1.8, device="pdf")
    ggsave(paste(filedir, fignum, "ii.pdf", sep=""), plot_avr, width=2.15, height=1.8, device="pdf")
    ggsave(paste(filedir, fignum, "iii.pdf", sep=""), plot_bdr, width=2.15, height=1.8, device="pdf")
    ggsave(paste(filedir, fignum, "iv.pdf", sep=""), plot_esr, width=2.15, height=1.8, device="pdf")
    ggsave(paste(filedir, fignum, "v.pdf", sep=""), plot_loc, width=2.15, height=1.8, device="pdf")
    ggsave(paste(filedir, fignum, "vi.pdf", sep=""), plot_sc, width=2.15, height=1.8, device="pdf")
    ggsave(paste(filedir, fignum, "vii.pdf", sep=""), plot_dep, width=2.15, height=1.8, device="pdf")
  }
}

d <- 10
sig2 <- 1
phi <- 3
n <- 10000
M <- 50


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
plot_rank(y, x, ylab = T, fignum = "a")


# b
mu_x <- rep(0.5, d)
x <- replicate(M, mvrnorm(n, mu = mu_x, Sigma = Sig_y))
plot_rank(y, x, fignum = "b")


## Type 2: errors in the variance

# a
sig2_x <- 0.65
Sig_x <- Sig_y*sig2_x
x <- replicate(M, mvrnorm(n, mu = mu_y, Sigma = Sig_x))
plot_rank(y, x, fignum = "c")

# b
sig2_x <- 1.35
Sig_x <- Sig_y*sig2_x
x <- replicate(M, mvrnorm(n, mu = mu_y, Sigma = Sig_x))
plot_rank(y, x, fignum = "d")


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
plot_rank(y, x, fignum = "e")

# b
phi_x <- 5
Sig_x <- matrix(NA, nrow = d, ncol = d)
for (i in 1:d) {
  for (j in 1:d) {
    Sig_x[i, j] <- sig2*exp(-abs(i-j)/phi_x)
  }
}  
x <- replicate(M, mvrnorm(n, mu = mu_y, Sigma = Sig_x))
plot_rank(y, x, fignum = "f")


