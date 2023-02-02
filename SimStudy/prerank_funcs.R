## pre-rank functions

mv_rank <- function(y, x) {
  dat <- cbind(y, x)
  M <- ncol(dat)
  rank <- sapply(1:M, function(i) sum(sapply(1:M, function(j) all(dat[, j] <= dat[, i]))))
  rank_y <- sum(rank[1] > rank[2:M]) + 1 + sample(0:sum(rank[1] == rank[2:M]), 1)
  return(rank_y)
}

av_rank <- function(y, x) {
  dat <- cbind(y, x)
  d <- length(y)
  M <- ncol(dat)
  c <- sapply(1:M, function(i) sapply(1:d, function(l) sum(dat[l, ] <= dat[l, i])))
  rank <- colMeans(c)
  rank_y <- sum(rank[1] > rank[2:M]) + 1 + sample(0:sum(rank[1] == rank[2:M]), 1)
  return(rank_y)
}

bd_rank <- function(y, x) {
  dat <- cbind(y, x)
  M <- ncol(dat)
  c <- sapply(1:M, function(i) sapply(1:d, function(l) sum(dat[l, ] <= dat[l, i])))
  rank <- rowMeans(sapply(1:d, function(l) 
    c[l, ]*(M - c[l, ]) + (c[l, ] - 1)*sapply(1:M, function(i) sum(c[l, ] == c[l, i]))))
  rank_y <- sum(rank[1] > rank[2:M]) + 1 + sample(0:sum(rank[1] == rank[2:M]), 1)
  return(rank_y)
}

mean_rank <- function(y, x) {
  g_y <- mean(y)
  g_x <- colMeans(x)
  rank_y <- sum(g_y > g_x) + 1 + sample(0:sum(g_y == g_x), 1)
  return(rank_y)
}

var_rank <- function(y, x) {
  g_y <- var(y)
  g_x <- apply(x, 2, var)
  rank_y <- sum(g_y > g_x) + 1 + sample(0:sum(g_y == g_x), 1)
  return(rank_y)
}

vg_rank <- function(y, x, p = 0.5) {
  g_y <- vsC(y, p)
  g_x <- apply(x, 2, vsC, p)
  rank_y <- sum(g_y > g_x) + 1 + sample(0:sum(g_y == g_x), 1)
  return(rank_y)
}

vg1_rank <- function(y, x, p = 0.5) {
  g_y <- sum(abs(diff(y))^p)
  g_x <- apply(x, 2, function(z) sum(abs(diff(z))^p))
  rank_y <- sum(g_y > g_x) + 1 + sample(0:sum(g_y == g_x), 1)
  return(rank_y)
}

es_rank <- function(y, x) {
  g_y <- scoringRules::es_sample(y, x)
  g_x <- apply(x, 2, function(z) scoringRules::es_sample(z, x))
  rank_y <- sum(g_y > g_x) + 1 + sample(0:sum(g_y == g_x), 1)
  return(rank_y)
}

## plot functions

plot_hist <- function(rank, n_bins, title = NULL, ylab = "Rel. Freq.", ymax = 0.2, axis = T){
  
  rank_freq <- sapply(1:n_bins, function(i) mean(rank == i, na.rm = T))
  df <- data.frame(freq = rank_freq, rank = as.factor(1:n_bins))
  
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