xacf <- function(xdat, thresh, lag.max = 10, plot = TRUE) {
  qlev <- mean(thresh >= xdat)
  xr <- rank(xdat, ties.method = "random", na.last = TRUE) /
    (length(na.omit(xdat)))

  lag.max <- as.integer(lag.max)
  xacf_res <- matrix(nrow = lag.max, ncol = 4L)
  nr <- length(xr)
  for (i in 1:lag.max) {
    chi <- mev::xdep.chi(
      xdat = cbind(xr[-(1:i)], xr[-seq(nr, length.out = i, by = -1)]),
      qlev = qlev,
      margtrans = "none",
      plot = FALSE
    )
    xacf_res[i, ] <- c(i, chi$coef, chi$lower, chi$upper)
  }
  colnames(xacf_res) <- c("lag", "coef", "lower", "upper")
  xacf_res <- as.data.frame(xacf_res)
  class(xacf_res) <- c("mev_xacf", "data.frame")
  if (isTRUE(plot)) {
    plot(xacf_res)
  }
  return(invisible(xacf_res))
}

plot.mev_xacf <- function(x, ...) {
  plot(
    x$lag,
    x$coef,
    type = "h",
    ylim = c(0, 1),
    xlab = "lag",
    ylab = "tail correlation",
    bty = "l"
  )
}
