light="#DCBCBC"
light_highlight="#C79999"
mid="#B97C7C"
mid_highlight="#A25050"
dark="#8F2727"
dark_highlight="#7C0000"

plot_group_params <- function(params, key, max_group_idx) {
  x_plots <- floor(sqrt(max_group_idx))
  y_plots <- ceiling(sqrt(max_group_idx))
  
  par(mfrow=c(x_plots, y_plots))

  for (group_idx in 1:max_group_idx) {
    min <- min(params[[key]][,group_idx])
    max <- max(params[[key]][,group_idx])
    breaks <- seq(min, max, (max - min) / 25)
    
    hist(params[[key]][,group_idx], breaks=breaks,
         col=c_dark, border=c_dark_highlight, main=paste(key, "=", group_idx), 
         xlab=key, xlim=c(min, max), yaxt='n', ylab="")
  }
}

plot_group_params_comp <- function(grouped_params, hier_params, key, max_group_idx) {
  x_plots <- floor(sqrt(max_group_idx))
  y_plots <- ceiling(sqrt(max_group_idx))
  
  par(mfrow=c(x_plots, y_plots))

  for (group_idx in 1:max_group_idx) {
    min <- min(c(grouped_params[[key]][,group_idx], hier_params[[key]][,group_idx]))
    max <- max(c(grouped_params[[key]][,group_idx], hier_params[[key]][,group_idx]))
    breaks <- seq(min, max, (max - min) / 25)
    
    hist(grouped_params[[key]][,group_idx], breaks=breaks,
         col=c_dark, border=c_dark_highlight, main=paste(key, "=", group_idx), 
         xlab=key, xlim=c(min, max), yaxt='n', ylab="")
    hist(hier_params[[key]][,group_idx], breaks=breaks,
         col=c_light, border=c_light_highlight, add=T)
  }
}

plot_ppc <- function(params, data) {
  breaks_delta <- 1.0 / length(data$y)
  breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

  ppc_hist <- hist(params$p_hat_ppc, breaks=breaks, plot=F)
  y_max <- max(ppc_hist$density)

  ave <- sum(data$y) / length(data$y)
  
  x <- hist(ave, breaks=breaks, freq=F,
       col="black", border="black", main="", 
       xlab="p_hat_ppc", xlim=c(-breaks_delta, 1 + breaks_delta),
       yaxt='n', ylab="", ylim=c(0, y_max + 1))
  plot(ppc_hist, freq=F, col=c_dark, border=c_dark_highlight, add=T)
}

plot_ppc_comp <- function(grouped_params, hier_params, data) {
  breaks_delta <- 1.0 / length(data$y)
  breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

  g_ppc_hist <- hist(grouped_params$p_hat_ppc, breaks=breaks, plot=F)
  h_ppc_hist <- hist(hier_params$p_hat_ppc, breaks=breaks, plot=F)
  y_max <- max(c(g_ppc_hist$density, h_ppc_hist$density))

  ave <- sum(data$y) / length(data$y)
  
  x <- hist(ave, breaks=breaks, freq=F,
       col="black", border="black", main="", 
       xlab="p_hat_ppc", xlim=c(-breaks_delta, 1 + breaks_delta),
       yaxt='n', ylab="", ylim=c(0, y_max + 1))
  plot(g_ppc_hist, freq=F, col=c_dark, border=c_dark_highlight, add=T)
  plot(h_ppc_hist, freq=F, col=c_light, border=c_light_highlight, add=T)
}

plot_group_ppc <- function(params, data, key, max_group_idx) {
  x_plots <- floor(sqrt(max_group_idx))
  y_plots <- ceiling(sqrt(max_group_idx))
  
  par(mfrow=c(x_plots, y_plots))

  for (group_idx in 1:max_group_idx) {
    group_counts <- data$y[data[[key]] == group_idx]
    breaks_delta <- 1.0 / length(group_counts)
    breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

    ppc_hist <- hist(params[[paste("p_hat_", key, "_ppc", sep="")]][,group_idx], 
                     breaks=breaks, plot=F)
    y_max <- max(ppc_hist$density)

    ave <- sum(group_counts) / length(group_counts)
    
    x <- hist(ave, breaks=breaks, freq=F,
         col="black", border="black", main=paste(key, "=", group_idx), 
         xlab="p_hat_ppc", xlim=c(-breaks_delta, 1 + breaks_delta),
         yaxt='n', ylab="", ylim=c(0, y_max + 1))
    plot(ppc_hist, freq=F, col=c_dark, border=c_dark_highlight, add=T)
  }
}

plot_group_ppc_comp <- function(grouped_params, hier_params, data, key, max_group_idx) {
  x_plots <- floor(sqrt(max_group_idx))
  y_plots <- ceiling(sqrt(max_group_idx))
  
  par(mfrow=c(x_plots, y_plots))

  for (group_idx in 1:max_group_idx) {
    group_counts <- data$y[data[[key]] == group_idx]
    breaks_delta <- 1.0 / length(group_counts)
    breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

    name <- paste("p_hat_", key, "_ppc", sep="")
    g_ppc_hist <- hist(grouped_params[[name]][,group_idx], breaks=breaks, plot=F)
    h_ppc_hist <- hist(hier_params[[name]][,group_idx], breaks=breaks, plot=F)
    y_max <- max(c(g_ppc_hist$density, h_ppc_hist$density))

    ave <- sum(group_counts) / length(group_counts)
    
    x <- hist(ave, breaks=breaks, freq=F,
         col="black", border="black", main=paste(key, "=", group_idx), 
         xlab="p_hat_ppc", xlim=c(-breaks_delta, 1 + breaks_delta),
         yaxt='n', ylab="", ylim=c(0, y_max + 1))
    plot(g_ppc_hist, freq=F, col=c_dark, border=c_dark_highlight, add=T)
    plot(h_ppc_hist, freq=F, col=c_light, border=c_light_highlight, add=T)
  }
}
