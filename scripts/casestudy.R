##### Application to EU MetNet wind speed reforecast fields

set.seed(1000)


##### set up

# devtools::install_github("sallen12/WeightedForecastVerification")
#devtools::install_github("AlexanderHenzi/epit")
library(ggplot2)
library(WeightedForecastVerification)
library(MultivCalibration)


################################################################################
## load data

data("wind_dat", package = "MultivCalibration")
list2env(wind_dat, globalenv())
rm(wind_dat)

n <- 1045
p <- 32
q <- 33
d <- p*q
M <- 11

# transform fields to be in the correct format
obs <- array(obs, c(n, p, q))
fc_ifs <- array(fc_ifs, c(n, p, q, M))
fc_ecc <- array(fc_ecc, c(n, p, q, M))
fc_ss <- array(fc_ss, c(n, p, q, M))

# create array with weights that decrease exponentially with the distance between grid points
w_mat <- array(NA, c(p, q, p, q))
for (i in 1:p) {
  for (j in 1:q) {
    w_mat[i, j, , ] <- outer(1:p, 1:q, FUN = function(k, l) exp(-sqrt((i - k)^2 + (j - l)^2)))
  }
}


################################################################################
## plot example forecasts

# function to plot and save fields over Europe
plot_map <- function(lons, lats, z, filename = NULL, ymin = 0, ymax = 15, title = NULL){
  if (is.matrix(z)) z <- as.vector(z)

  world <- map_data("world")
  ind <- (world$long >= min(lons) & world$long < max(lons)) & (world$lat >= min(lats) & world$lat <= max(lats))
  world <- world[ind, ]

  df <- data.frame(lat = lats, lon = lons, z = z)
  plot_obj <- ggplot() + geom_tile(data = df, aes(lon, lat, fill = z), color = "white") +
    borders("world") +
    coord_fixed(ylim = range(lats), xlim = range(lons)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradient(low = "white", high = "red", limits = c(ymin, ymax),
                        guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    theme_void() + theme(legend.title = element_blank(), legend.position = "bottom",
                         panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                         legend.key.width = unit(0.3, "in")) +
    ggtitle(title)


  if (!is.null(filename)) {
    ggsave(filename, plot_obj, width = 2, height = 2.5)
  }

  return(plot_obj)

}

day <- sample(1:n, 1) # 1008
plot_map(coord[, 2], coord[, 1], z = obs[day, , ], filename = "scripts/fig_3a.pdf")
plot_map(coord[, 2], coord[, 1], z = fc_ifs[day, , , 1], filename = "scripts/fig_3b.pdf")
plot_map(coord[, 2], coord[, 1], z = fc_ecc[day, , , 1], filename = "scripts/fig_3c.pdf")
plot_map(coord[, 2], coord[, 1], z = fc_ss[day, , , 1], filename = "scripts/fig_3d.pdf")


################################################################################
## evaluation

##### univariate calibration

# wrapper to get univariate ranks
get_univ_ranks <- function(y, x) {
  ranks <- sapply(1:nrow(y), function(i) sapply(1:ncol(x), function(j) get_prerank_gr(y, x, prerank = function(z) z[i, j])))
  return(ranks)
}

univ_ranks_ifs <- sapply(1:n, function(i) get_univ_ranks(obs[i, , ], fc_ifs[i, , , ]))
univ_ranks_ecc <- sapply(1:n, function(i) get_univ_ranks(obs[i, , ], fc_ecc[i, , , ]))

pit_hist(as.vector(univ_ranks_ifs), ymax = 0.31, xlab = NULL, ylab = NULL, xticks = F, yticks = F)
ggsave("scripts/fig_4a.pdf", width = 3, height = 2.5)
pit_hist(as.vector(univ_ranks_ecc), ymax = 0.31, xlab = NULL, ylab = NULL, xticks = F, yticks = F)
ggsave("scripts/fig_4b.pdf", width = 3, height = 2.5)


##### multivariate calibration


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
plot_ranks <- function(rank_df, filename = NULL, ylabs = NULL, ymax = 0.31) {
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


## IFS
rank_df_ifs <- get_ranks(obs, fc_ifs, w_mat, t = 6, parallel = TRUE)
plot_ranks(rank_df_ifs, filename = paste0(filedir, "5a", ".pdf"), ylabs = ylabs)

## ECC
rank_df_ecc <- get_ranks(obs, fc_ecc, w_mat, t = 6, parallel = TRUE)
plot_ranks(rank_df_ecc, filename = paste0(filedir, "5b", ".pdf"))

## SS
rank_df_ss <- get_ranks(obs, fc_ss, w_mat, t = 6, parallel = TRUE)
plot_ranks(rank_df_ss, filename = paste0(filedir, "5c", ".pdf"))


##### e-values

alpha <- 0.05

get_evals <- function(r, lag, m, n0 = 1, strategy = "betabinom") {
  evals <- epit::e_rank_histogram(r = r, h = lag, m = m, options = list(n0 = n0), strategy = strategy)$evalues_h
  evals <- epit:::evalue_combine_h(lapply(evals, function(x) x$e))
  return(evals)
}

plot_evals <- function(mv_ranks_ifs, mv_ranks_ecc, mv_ranks_ss, filename = NULL) {
  mv_ranks <- list()

  mv_evals_ifs <- get_evals(r = mv_ranks_ifs, lag = 5, m = M, n0 = 20)
  mv_evals_ecc <- get_evals(r = mv_ranks_ecc, lag = 5, m = M, n0 = 20)
  mv_evals_ss <- get_evals(r = mv_ranks_ss, lag = 5, m = M, n0 = 20)

  df <- data.frame(x = 1:n, e = c(mv_evals_ifs, mv_evals_ecc, mv_evals_ss),
                   mth = rep(c(" IFS", "ECC", "SS"), each = n))
  plt_prerank <- ggplot() + geom_line(data = df, aes(x = x, y = e, col = mth)) +
    geom_hline(aes(yintercept = 3*(exp(1)*log(5))/alpha), lty = "dotted") +
    geom_vline(aes(xintercept = (1:5)*(n/5)), col = "grey") +
    scale_y_continuous(name = "Cumulative e-value", trans = "log",
                       limits = c(1e-2, 1e12),
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_continuous(name = "Day",
                       breaks = seq(200, 1400, 200),
                       expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.title = element_blank(),
          legend.justification = c(0, 1), legend.position = c(0.01, 0.99)) +
    guides(col = guide_legend(nrow = 1))

  if (!is.null(filename)) {
    ggsave(filename, plt_prerank, width = 3.5, height = 3.3)
  }

  return(plt_prerank)
}

plot_evals(rank_df_ifs$bdr, rank_df_ecc$bdr, rank_df_ss$bdr, filename = "scripts/fig_6a.pdf")
plot_evals(rank_df_ifs$var, rank_df_ecc$var, rank_df_ss$var, filename = "scripts/fig_6b.pdf")
plot_evals(rank_df_ifs$dep, rank_df_ecc$dep, rank_df_ss$dep, filename = "scripts/fig_6c.pdf")


##### spearman's rank correlation coefficient

get_srcc <- function(rank_df) {
  k <- ncol(rank_df)
  srcc <- sapply(1:k, function(i) sapply(1:k, function(j) cor(rank_df[[i]], rank_df[[j]])))
  rownames(srcc) <- colnames(srcc) <- colnames(rank_df)
  return(srcc)
}

round(get_srcc(rank_df_ifs), 2) ## Table 2
round(get_srcc(rank_df_ecc), 2)
round(get_srcc(rank_df_ss), 2)


