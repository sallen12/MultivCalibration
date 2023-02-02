################################################################################
############### utility functions to perform simulation studies ################
################################################################################

################################################################################
## prerank functions

# function to get rank from preranks
get_rank <- function(preranks, obs_pos = 1) {
  sum(preranks[obs_pos] > preranks[-obs_pos]) + 1 + 
    sample(0:sum(preranks[obs_pos] == preranks[-obs_pos]), 1)
}

# multivariate rank
mv_rank <- function(y, x) {
  dat <- cbind(y, x)
  M <- ncol(dat)
  rank <- sapply(1:M, function(i) sum(sapply(1:M, function(j) all(dat[, j] <= dat[, i]))))
  rank_y <- get_rank(rank)
  return(rank_y)
}

# average rank
av_rank <- function(y, x) {
  dat <- cbind(y, x)
  d <- length(y)
  M <- ncol(dat)
  c <- sapply(1:M, function(i) sapply(1:d, function(l) sum(dat[l, ] <= dat[l, i])))
  rank <- colMeans(c)
  rank_y <- get_rank(rank)
  return(rank_y)
}

# band-depth rank
bd_rank <- function(y, x) {
  dat <- cbind(y, x)
  M <- ncol(dat)
  c <- sapply(1:M, function(i) sapply(1:d, function(l) sum(dat[l, ] <= dat[l, i])))
  rank <- rowMeans(sapply(1:d, function(l) 
    c[l, ]*(M - c[l, ]) + (c[l, ] - 1)*sapply(1:M, function(i) sum(c[l, ] == c[l, i]))))
  rank_y <- get_rank(rank)
  return(rank_y)
}

# mean rank
mean_rank <- function(y, x) {
  g_y <- mean(y)
  g_x <- colMeans(x)
  rank_y <- get_rank(c(g_y, g_x))
  return(rank_y)
}

# variance rank
var_rank <- function(y, x) {
  g_y <- var(y)
  g_x <- apply(x, 2, var)
  rank_y <- get_rank(c(g_y, g_x))
  return(rank_y)
}

# variogram rank
vg_rank <- function(y, x, w_mat = NULL, p = 0.5) {
  if (is.null(w_mat)) {
    d <- length(y)
    w_mat <- matrix(1, nrow = d, ncol = d)
  }
  g_y <- vario(y, w_mat, p)
  g_x <- apply(x, 2, vario, w_mat, p)
  rank_y <- get_rank(c(g_y, g_x))
  return(rank_y)
}

# energy score rank
es_rank <- function(y, x) {
  g_y <- scoringRules::es_sample(y, x)
  g_x <- apply(x, 2, function(z) scoringRules::es_sample(z, x))
  rank_y <- get_rank(c(g_y, g_x))
  return(rank_y)
}

################################################################################
## plot functions

# function to plot rank histograms
plot_hist <- function(rank, n_bins, title = NULL, ylab = "Rel. Freq.", ymax = 0.2, axis = T){
  
  rank_freq <- sapply(1:n_bins, function(i) mean(rank == i, na.rm = T))
  df <- data.frame(freq = rank_freq, rank = as.factor(1:n_bins))
  df$freq[df$freq > ymax] <- ymax
  
  hist_plot <- ggplot(df, aes(x = rank, y = freq)) + geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 1/n_bins), col="red", lty = "dashed") + 
    scale_x_discrete(name = "Rank") +
    scale_y_continuous(name = ylab, limits = c(0, ymax), expand = c(0, 0)) +
    theme_bw() +
    theme(legend.title = element_blank(), panel.grid = element_blank()) +
    ggtitle(title)
  if (!axis) {
    hist_plot <- hist_plot + theme(axis.text = element_blank(), axis.title.x = element_blank())
  }
  return(hist_plot)
}

# function to plot and save results
plot_rank <- function(y, x, w_mat = NULL, ylab = F, fignum = NULL) {
  
  # calculate weight matrix for variogram prerank
  if (is.null(w_mat)) {
    d <- ncol(y)
    w_mat <- matrix(1, nrow = d, ncol = d)
  }
  
  rank_df <- data.frame(mvr = sapply(1:n, function(i) mv_rank(y[i, ], x[i, , ])),
                        avr = sapply(1:n, function(i) av_rank(y[i, ], x[i, , ])),
                        bdr = sapply(1:n, function(i) bd_rank(y[i, ], x[i, , ])),
                        esr = sapply(1:n, function(i) es_rank(y[i, ], x[i, , ])),
                        loc = sapply(1:n, function(i) mean_rank(y[i, ], x[i, , ])),
                        sc = sapply(1:n, function(i) var_rank(y[i, ], x[i, , ])),
                        dep = sapply(1:n, function(i) vg_rank(y[i, ], x[i, , ], w_mat)))
  
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
    filedir <- "Figures/fig_"
    width <- 1.6
    height <- 1.4
    ggsave(paste(filedir, fignum, "i.pdf", sep=""), plot_mvr, width=width, height=height, device="pdf")
    ggsave(paste(filedir, fignum, "ii.pdf", sep=""), plot_avr, width=width, height=height, device="pdf")
    ggsave(paste(filedir, fignum, "iii.pdf", sep=""), plot_bdr, width=width, height=height, device="pdf")
    ggsave(paste(filedir, fignum, "iv.pdf", sep=""), plot_esr, width=width, height=height, device="pdf")
    ggsave(paste(filedir, fignum, "v.pdf", sep=""), plot_loc, width=width, height=height, device="pdf")
    ggsave(paste(filedir, fignum, "vi.pdf", sep=""), plot_sc, width=width, height=height, device="pdf")
    ggsave(paste(filedir, fignum, "vii.pdf", sep=""), plot_dep, width=width, height=height, device="pdf")
  } else {
    gridExtra::grid.arrange(plot_mvr, plot_avr, plot_bdr, plot_esr, plot_loc, plot_sc, plot_dep, ncol = 1)
  }
}

# function to plot random fields
plot_grf <- function(coords, z) {
  grf_df <- data.frame(coords, z)
  ggplot(grf_df) + geom_raster(aes(x = x, y = y, fill = z)) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_distiller(limits = c(-3, 3), palette = "RdBu") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          legend.position = "none")
}