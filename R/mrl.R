#' Automated mean residual life plots
#'
#' This function implements the automated proposal from Section 2.2 of Langousis et al. (2016) for mean residual life plots. It returns the threshold that minimize the weighted mean square error and moment estimators for the scale and shape parameter based on weighted least squares.
#'
#' The procedure consists in estimating the usual mean residual life as a function of the threshold, and looking for an order statistic or threshold value above which the fit is more or less linear.
#'
#' @references Langousis, A., A. Mamalakis, M. Puliga and R. Deidda (2016).
#' \emph{Threshold detection for the generalized Pareto distribution:
#' Review of representative methods and application to the
#' NOAA NCDC daily rainfall database}, Water Resources Research, \strong{52}, 2659--2681.
#'
#' @return a list containing
#' \itemize{
#' \item \code{thresh}: candidate threshold vector
#' \item \code{thresh0}: selected threshold
#' \item \code{scale}: scale parameter estimate
#' \item \code{shape}: shape parameter estimate
#' \item \code{mrl}: empirical mean excess values
#' \item \code{xdat}: ordered observations
#' \item \code{intercept}: intercept for mean excess value at chosen threshold
#' \item \code{slope}: slope for mean excess value at chosen threshold
#' \item \code{tmanual}: logical; \code{TRUE} if the user passed a vector of thresholds
#' }
#'
#' @param xdat [numeric] vector of observations
#' @param thresh [numeric] vector of thresholds; if missing, uses all order statistics from the 20th largest until \code{kmax} as candidates
#' @param kmax [integer] maximum number of order statistics
#' @param plot [logical] if \code{TRUE} (default), return a plot of the mean residual life plot with the fitted slope
#' and the chosen threshold
#' @param ... additional arguments, currently ignored
#' @export
#' @examples
#' thselect.mrl(rgp(n = 100))
thselect.mrl <- function(
  xdat,
  thresh,
  kmax,
  plot = TRUE,
  ...
) {
  k <- 9
  # expected conditional exceedances based
  # on at least k+1 obs (default to 10, hardcoded)
  xdat <- sort(xdat[is.finite(xdat)], decreasing = TRUE)
  if (!missing(thresh)) {
    stopifnot(isTRUE(all(is.finite(thresh))))
    thresh <- sort(thresh)
    nt <- length(thresh)
    xdat <- xdat[xdat > thresh[1]]
    if (isTRUE(findInterval(x = xdat[10], vec = thresh) < nt)) {
      stop(
        "Not enough observations for reliable estimation:\n users must provide at least 10 exceedances over largest threshold."
      )
    }
    user_thresh <- TRUE
  } else {
    if (!missing(kmax)) {
      kmax <- as.integer(kmax)
      stopifnot(kmax <= length(xdat))
      xdat <- xdat[seq_len(kmax)]
    }
    # Use order statistics
    thresh <- xdat[-(1:20)]
    nt <- length(thresh)
    user_thresh <- FALSE
  }
  n <- length(xdat)
  # At least 20 observations
  stopifnot(n >= 20)
  # Compute running mean and variance efficiently
  meanX <- cumsum(xdat) / seq_along(xdat)
  # Running variance
  varX <- cumsum(
    (xdat[-1] - meanX[-1]) *
      (xdat[-1] - meanX[-n])
  ) /
    seq.int(1L, to = n - 1, by = 1L)
  # Exclude 10 largest observations to ensure that
  # we have not too much variability
  excu <- (meanX[-n] - xdat[-1])[-(1:k)]
  weights <- (seq_len(n - 1) / varX)[-(1:k)]
  xk <- xdat[-(1:(k + 1))]
  # Containers
  mse <- numeric(nt)
  coefs <- matrix(0, nrow = nt, ncol = 2)
  # Could use rank-one update, but is improvement worth it?
  for (i in seq_along(thresh)) {
    fit <- lm(excu ~ xk, weights = weights, subset = xk > thresh[i])
    mse[i] <- weighted.mean(
      x = fit$residuals^2,
      w = weights[seq_len(nobs(fit))]
    )
    coefs[i, ] <- coef(fit)
  }
  # plot(x = thresh,
  #      y = sqrt(mse),
  #      type = "l",
  #      bty = "l")
  index <- which.min(mse)
  thresh0 <- thresh[index]
  shape <- coefs[index, 2] / (1 + coefs[index, 2])
  scale <- coefs[index, 1] * (1 - shape) + shape * xdat[k + index]
  ret <- list(
    thresh = thresh,
    thresh0 = thresh0,
    scale = scale,
    shape = shape,
    mrl = excu,
    xdat = xdat,
    intercept = coefs[index, 1],
    slope = coefs[index, 2],
    tmanual = user_thresh
  )
  class(ret) <- "mev_thselect_automrl"
  if (plot) {
    plot(ret)
  }
  return(invisible(ret))
}
#' @export
plot.mev_thselect_automrl <- function(x, ...) {
  k <- 9L # hardcoded
  plot(
    x = x$xdat[-(1:(k + 1))],
    y = x$mrl,
    ylab = "mean excess value",
    xlab = "observation",
    pch = 20,
    bty = "l"
  )
  if (isTRUE(x$tmanual)) {
    rug(x$thresh)
  }
  abline(v = x$thresh0, lty = 2)
  abline(a = x$intercept, b = x$slope)
}

#' @export
print.mev_thselect_automrl <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  cat("Threshold selection method: mean residual life plot\n")
  cat(
    "Langousis' weighted regression line: intercept",
    round(x$intercept, digits),
    "and slope",
    round(x$slope, digits),
    "\n"
  )
  cat("Selected threshold:", round(x$thresh0, digits), "\n")
  return(invisible(NULL))
}

#' Automated mean residual life plots
#'
#' This function implements the automated proposal from
#' Section 2.2 of Langousis et al. (2016)
#' for mean residual life plots. It returns the threshold
#' that minimize the weighted mean square error and
#' moment estimators for the scale and shape parameter
#' based on weighted least squares.
#'
#' The procedure consists in estimating the usual
#'
#' @references Langousis, A., A. Mamalakis, M. Puliga and R. Deidda (2016).
#' \emph{Threshold detection for the generalized Pareto distribution:
#' Review of representative methods and application to the
#' NOAA NCDC daily rainfall database}, Water Resources Research, \strong{52}, 2659--2681.
#'
#' @return a list containing
#' \itemize{
#' \item \code{thresh}: selected threshold
#' \item \code{scale}: scale parameter estimate
#' \item \code{shape}: shape parameter estimate
#' }
#'
#' @param xdat [numeric] vector of observations
#' @param thresh [numeric] vector of thresholds; if missing, uses all order statistics from the 20th largest until \code{kmax} as candidates
#' @param kmax [integer] maximum number of order statistics
#' @param plot [logical] if \code{TRUE} (default), return a plot of the mean residual life plot with the fitted slope
#' and the chosen threshold
#' @param ... additional arguments, currently ignored
#' @export
#'
#' @keywords internal
automrl <- function(
  xdat,
  kmax,
  thresh,
  plot = TRUE,
  ...
) {
  .Deprecated(new = "thselect.mrl", package = "mev")
  thselect.mrl(
    xdat = xdat,
    thresh = thresh,
    kmax = kmax,
    plot = plot,
    ...
  )
}

#' Mean residual life plot
#'
#' Computes mean of sample exceedances over a range of thresholds or for a pre-specified number of largest order statistics, and returns a plot with 95\% Wald-based confidence intervals as a function of either the threshold or the number of exceedances. The main purpose is the plotting method, which generates the so-called mean residual life plot. The latter should be approximately linear over the threshold for a generalized Pareto distribution
#'
#' @references Davison, A.C. and R.L. Smith (1990). Models for Exceedances over High Thresholds (with discussion), \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, \bold{52}(3), 393--442.
#' @export
#' @param xdat vector of sample observations
#' @param thresh vector of thresholds
#' @param kmin integer giving the minimum number of exceedances; ignored if \code{thresh} is provided. Default to 10
#' @param kmax integer giving the maximum number of exceedances; ignored if \code{thresh} is provided. Default to sample size.
#' @param plot logical; if \code{TRUE}, call the plot method
#' @param level double giving the level of confidence intervals for the plot, default to 0.95
#' @param xlab string indicating whether to use thresholds (\code{thresh}) or number of largest order statistics (\code{nexc}) for the x-axis
#' @param type string whether to plot pointwise confidence intervals using segments (\code{"ptwise"}) or using dashed lines (\code{"band"})
#' @param ... additional arguments, currently ignored
#' @return an invisible list with mean sample exceedances and standard deviation, number of exceedances, threshold
#' @examples
#' tstab.mrl(
#'  xdat = rgp(n = 100, shape = -0.5),
#'  xlab = "thresh",
#'  kmax = 50)
#' tstab.mrl(
#'  rexp(100),
#'  thresh = qexp(seq(0, 0.9, by = 0.01)))
tstab.mrl <- function(
  xdat,
  thresh,
  kmin = 10L,
  kmax = length(xdat),
  plot = TRUE,
  level = 0.95,
  xlab = c("thresh", "nexc"),
  type = c("band", "ptwise"),
  ...
) {
  # expected conditional exceedances based
  # on at least k+1 obs (default to 10, hardcoded)
  xdat <- as.numeric(xdat)
  xdat <- sort(xdat[is.finite(xdat)], decreasing = TRUE)
  if (!missing(thresh)) {
    stopifnot(isTRUE(all(is.finite(thresh))))
    thresh <- sort(thresh)
    nt <- length(thresh)
    xdat <- xdat[xdat > thresh[1]]
    user_thresh <- TRUE
  } else {
    if (!missing(kmax)) {
      kmax <- as.integer(kmax)
      kmin <- as.integer(kmin)
      stopifnot(kmax <= length(xdat))
    } else {
      kmax <- length(xdat)
    }
    stopifnot(kmax > kmin, kmin > 1)
    xdat <- xdat[seq_len(kmax + 1L)]
    # Use order statistics
    thresh <- xdat[-(1:kmin)]
    nexc <- kmin:kmax
    nt <- length(thresh)
    user_thresh <- FALSE
  }
  n <- length(xdat)
  if (!user_thresh) {
    # Compute running mean and variance efficiently
    meanX <- cumsum(xdat) / seq_along(xdat)
    # Running variance
    varX <- cumsum(
      (xdat[-1] - meanX[-1]) *
        (xdat[-1] - meanX[-n])
    ) /
      seq.int(1L, to = n - 1, by = 1L)
    # Exclude 10 largest observations to ensure that
    # we have not too much variability
    mrl <- (meanX[-n] - xdat[-1])[-(1:(kmin - 1))]
    sdmrl <- sqrt(varX)[-(1:(kmin - 1))]
  } else {
    mrl <- nexc <- sdmrl <- numeric(nt)
    for (i in seq_along(thresh)) {
      exc <- xdat[xdat > thresh[i]] - thresh[i]
      nexc[i] <- length(exc)
      mrl[i] <- mean(exc)
      sdmrl[i] <- sd(exc)
    }
  }
  ret <- list(
    nexc = nexc,
    thresh = thresh,
    mrl = mrl,
    sd = sdmrl
  )
  class(ret) <- "mev_tstab_mrl"
  if (plot) {
    plot(ret, level = level, xlab = xlab, type = type)
  }
  return(invisible(ret))
}

#' Mean residual life parameter stability plot
#'
#' @param xlab [string]; whether to plot mean residual life plot as a function of threshold value of number of exceedances
#' @param level [numeric] level of Wald confidence intervals
#' @param type [string] whether to plot pointwise confidence intervals using segments (\code{"ptwise"}) or using dashed lines (\code{"band"})
#' @param x object resulting from a call to \code{tstab.mrl}
#' @param ... additional arguments, currently ignored
#' @return \code{NULL}; use to produce plots
#' @export
plot.mev_tstab_mrl <- function(
  x,
  xlab = c("thresh", "nexc"),
  level = 0.95,
  type = c("band", "ptwise"),
  ...
) {
  xlab <- match.arg(xlab)
  type <- match.arg(type)
  os <- isTRUE(xlab == "nexc")
  if (isTRUE(os)) {
    xv <- x$nexc
    xlab = "number of exceedances"
  } else {
    xv <- x$thresh
    xlab <- "threshold"
  }
  stopifnot(
    is.numeric(level),
    length(level) == 1L,
    level < 1,
    level > 0
  )
  alpha <- 1 - level
  lower <- x$mrl + qnorm(alpha / 2) * x$sd / sqrt(x$nexc)
  upper <- x$mrl + qnorm(1 - alpha / 2) * x$sd / sqrt(x$nexc)
  if (type == "ptwise") {
    plot(
      x = xv,
      y = x$mrl,
      ylab = "mean excess value",
      xlab = xlab,
      ylim = c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE)),
      pch = 20,
      bty = "l"
    )
    segments(
      x0 = xv,
      y0 = lower,
      y1 = upper,
      lwd = 0.1
    )
  } else {
    matplot(
      x = xv,
      y = cbind(x$mrl, lower, upper),
      ylab = "mean excess value",
      xlab = xlab,
      ylim = c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE)),
      type = "l",
      lty = c(1, 2, 2),
      bty = "l",
      col = c("black", "grey", "grey")
    )
  }
  return(invisible(NULL))
}
