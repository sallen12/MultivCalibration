##### Application to wind speed fields



### set up

library(maps)
library(ggplot2)

source("utility_funcs.R")
Rcpp::sourceCpp("variogram_func.cpp")


##### load data (lead time = 5 days)

x <- R.matlab::readMat("CaseStudy/WS.mat")[[1]]

n_train <- 2149
n_test <- 1432
n_ens <- 10
n_locs <- 1353
lons <- -21:19
lats <- 37:69

tr_data <- x[, 1:n_train, ]
ts_data <- x[, n_train + (1:n_test), ]
rm(x)

tr_obs <- tr_data[1, , ]
ts_obs <- ts_data[1, , ]

tr_ens <- tr_data[2:11, , ]
ts_ens <- ts_data[2:11, , ]
rm(tr_data, ts_data)

tr_ens_mean <- apply(tr_ens, c(2, 3), mean)
ts_ens_mean <- apply(ts_ens, c(2, 3), mean)
tr_ens_sd <- apply(tr_ens, c(2, 3), sd)
ts_ens_sd <- apply(ts_ens, c(2, 3), sd)


## plot on map

plot_map(lons, lats, z = colMeans(tr_obs))

################################################################################
## perform post-processing

source("postprocessing.R")

## univariate post-processing

get_pars(tr_obs, tr_ens_mean, tr_ens_sd)


### multivariate post-processing

rank_mat <- apply(ts_ens, c(2, 3), rank) # ensemble dependence template

## sample from post-processed distribution

# wrapper to sample from post-processed forecast distributions
get_ens <- function(ens_mean, ens_sd, par_mat) {
  n_test <- nrow(ens_mean)
  n_locs <- ncol(ens_mean)
  pp_mu <- t(sapply(1:n_test, function(i) rowSums(cbind(1, ens_mean[i, ]) * t(par_mat[1:2, ]))))
  pp_sig <- t(sapply(1:n_test, function(i) rowSums(cbind(1, ens_sd[i, ]) * t(par_mat[3:4, ]^2))))
  ens <- lapply(1:n_test, function(i) sapply(1:n_locs, function(j) get_sample(pp_mu[i, j], pp_sig[i, j])))
  ens <- aperm(simplify2array(ens), c(1, 3, 2))
  return(ens)
}

ts_com_ens <- get_ens(ts_ens_mean, ts_ens_sd, par_mat, n_test) # 
ts_ecc_ens <- apply_ecc(ts_com_ens, rank_mat) # apply ECC
ts_ind_ens <- apply_ecc(ts_com_ens) # apply randomisation


# plot example forecasts
day <- sample(1:n_test, 1) # 1254
plot_map(lons, lats, z = ts_obs[day, ])
plot_map(lons, lats, z = ts_ecc_ens[10, day, ])
plot_map(lons, lats, z = ts_com_ens[1, day, ])
plot_map(lons, lats, z = ts_ind_ens[3, day, ])



################################################################################
## univariate calibration

# function to get univariate ranks
get_univ_ranks <- function(y, ens) {
  univ_ranks <- sapply(1:nrow(y), function(i) {
    sapply(1:ncol(y), function(j) {
      rank(c(y[i, j], ens[, i, j]), ties.method = "random")[1]
    })
  })
}

univ_ranks <- get_univ_ranks(ts_obs, ts_ens)
plot_hist(univ_ranks, n_bins = 11, title = "Raw")
univ_ranks <- get_univ_ranks(ts_obs, ts_ecc_ens)
plot_hist(univ_ranks, n_bins = 11, title = "ECC") 
univ_ranks <- get_univ_ranks(ts_obs, ts_com_ens)
plot_hist(univ_ranks, n_bins = 11, title = "Com.")
univ_ranks <- get_univ_ranks(ts_obs, ts_ind_ens)
plot_hist(univ_ranks, n_bins = 11, title = "Ind.")



################################################################################
## multivariate calibration

## standardise prior to multivariate evaluation

# function to standardise fields
standardise_fields <- function(x, mean_vec, sd_vec) {
  mean_mat <- matrix(mean_vec, nrow = nrow(x), ncol = ncol(x), byrow = T)
  sd_mat <- matrix(sd_vec, nrow = nrow(x), ncol = ncol(x), byrow = T)
  x_std <- (x - mean_vec)/sd_vec
  return(x_std)
}

tr_obs_mn <- colMeans(tr_obs)
tr_obs_sd <- apply(tr_obs, 2, sd)
ts_obs <- standardise_fields(ts_obs, tr_obs_mn, tr_obs_sd)
for (i in 1:n_ens) {
  ts_ens[i, , ] <- standardise_fields(ts_ens[i, , ], tr_obs_mn, tr_obs_sd)
  ts_ecc_ens[i, , ] <- standardise_fields(ts_ecc_ens[i, , ], tr_obs_mn, tr_obs_sd)
  ts_com_ens[i, , ] <- standardise_fields(ts_com_ens[i, , ], tr_obs_mn, tr_obs_sd)
  ts_ind_ens[i, , ] <- standardise_fields(ts_ind_ens[i, , ], tr_obs_mn, tr_obs_sd)
}


## get weight matrix for variogram prerank

nbhd_size <- 1 # neighbourhood size
coords <- expand.grid(lons, lats)
w_mat <- sapply(1:n_locs, function(i) sapply(1:n_locs, function(j) as.numeric(max(abs(coords[i, ] - coords[j, ])) <= nbhd_size)))


## plot and save multivariate rank histograms

plot_rank(ts_obs, aperm(ts_ens, c(2, 3, 1)), w_mat, ylab = T) # raw ensemble
plot_rank(ts_obs, aperm(ts_ecc_ens, c(2, 3, 1)), w_mat, ylab = T) # ECC
plot_rank(ts_obs, aperm(ts_com_ens, c(2, 3, 1)), w_mat, ylab = T) # comonotonic copula
plot_rank(ts_obs, aperm(ts_ind_ens, c(2, 3, 1)), w_mat, ylab = T) # independence copula
  

