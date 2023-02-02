################################################################################
###### utility functions to perform post-processing of wind speed fields #######
################################################################################

# function to calculate the CRPS of a truncated normal distribution from parameters
objective_crps <- function(par, y, ens_mn, ens_sd, obj = T) {
  m <- cbind(1, ens_mn) %*% par[1:2]
  s <- cbind(1, ens_sd) %*% (par[3:4]^2)
  objective <- scoringRules::crps_tnorm(y, m, s)
  if (obj) {
    objective <- mean(objective)
  }
  return(objective)
}

# function to optimise the CRPS
optim_crps <- function(y, ens_mn, ens_sd) {
  opt_par <- optim(fn = objective_crps, par = c(0, 1, 1, 1), y = y, ens_mn = ens_mn, ens_sd = ens_sd)
  return(opt_par$par)
}

# function to get post-processing model parameters
get_pars <- function(y, tr_mn, tr_sd, load_pars = T, save_pars = F) {
  if (!load_pars) {
    par_mat <- sapply(1:ncol(y), function(i) {
      print(i)
      optim_crps(y[, i], tr_mn[, i], tr_sd[, i])
    })
    if (save_pars) {
      save("pp_pars")
    }
  } else {
    load("pp_pars.RData")
  }
  return(par_mat)
}

# function to sample from a truncated normal distribution
get_sample <- function(mean, sd, M = n_ens){
  quants <- (1:M)/(M + 1)
  ens <- crch::qtnorm(quants, mean = mean, sd = sd, left = 0)
  return(ens)
}

# function to apply ensemble copula coupling
apply_ecc <- function(ref_ens, rank_template = NULL) {
  ens <- array(NA, dim(ref_ens))
  n_ens <- dim(ens)[1]
  n_test <- dim(ens)[2]
  n_locs <- dim(ens)[3]
  if (!is.null(rank_template)) {
    for(i in 1:n_test){
      for(j in 1:n_locs){
        rank_ord <- rank_mat[, i, j]
        ens[ , i, j] <- ref_ens[rank_ord, i, j]
      }
    }
  } else {
    for(i in 1:n_test){
      for(j in 1:n_locs){
        ens[ , i, j] <- ref_ens[sample(1:n_ens), i, j]
      }
    }
  }
  return(ens)
}

