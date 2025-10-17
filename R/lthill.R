#' Lower-trimmed Hill estimator for the shape parameter
#'
#' Given a sample of Pareto-tailed samples (positive tail index),
#' compute the lower-trimmed Hill estimator. If \eqn{k0=k}, the estimator reduces to Hill's estimator for the shape index
#' @param xdat [numeric] vector of positive observations
#' @param k0 [integer] vector of number of largest order statistics, no greater than \code{k}
#' @param k [integer] number of order statistics for the threshold
#' @param sorted [logical] if \code{TRUE}, data are assumed to be sorted in decreasing order.
#' @param ... additional arguments for other routines (notably \code{vectorize})
#' @export
#' @references Bladt, M., Albrecher, H. & Beirlant, J. (2020) \emph{Threshold selection and trimming in extremes}. Extremes, 23, 629-665 . \doi{10.1007/s10687-020-00385-0}
#' @return a scalar with the shape parameter estimate if \code{k0} is a scalar, otherwise a data frame with columns \code{k0} for the number of exceedances and \code{shape} for the tail index.
#' @examples
#' # Pareto sample
#' n <- 200
#' xdat <- 10/(1 - runif(n)) - 10
#' shape.lthill(xdat = xdat, k = 100, k0 = 5:100)
#' @importFrom grDevices grey.colors
shape.lthill <- function(
  xdat,
  k,
  k0 = k,
  # vectorize = FALSE,
  sorted = FALSE,
  ...
) {
  args <- list(...)
  k0min <- 5L
  if (isTRUE(args$vectorize)) {
    k0min <- 1L
    k0 <- 1:k
  }
  sorted <- isTRUE(sorted)
  if (length(k) != 1L) {
    stop(
      "Invalid argument \"k\": must be an integer of length 1."
    )
  }
  if (!sorted) {
    xdat <- sort(xdat, decreasing = TRUE, na.last = NA)
  }
  k0 <- as.integer(k0)
  if (length(k0) > 1L) {
    k0o <- sort(k0)
    k0 <- k0o[length(k0)]
  } else {
    k0o <- k0
  }
  k <- as.integer(k)
  n <- length(xdat)
  stopifnot(k0 <= k, k0 >= 1, k < n)
  k0o <- k0o[k0o >= k0min & k0o <= k]
  # Threshold
  th <- xdat[k + 1]
  # Observations should be Pareto-tailed, and positive
  stopifnot(th > 0)
  vectorize <- length(k0) == 1
  if (!isTRUE(vectorize)) {
    return(
      (mean(log(xdat[1:k0])) - log(th)) /
        (1 + ifelse(k0 < k, sum(1 / ((k0 + 1):k)), 0))
    )
  } else {
    y <- log(xdat[1:k0])
    denom <- sapply(1:k0, function(j) {
      1 + sum(1 / ((j + 1):k))
    })
    if (k0 == k) {
      denom[k0] <- 1
    }
    shape <- (cumsum(y) / (1:k0) - log(th)) / denom
    df <- data.frame(k0 = k0o, shape = shape[1:k0 %in% k0o])
    attr(df, which = "k") <- k
    return(df)
  }
}

#' Trimmed Hill estimator for the shape parameter
#'
#' Given a sample of Pareto-tailed samples (positive tail index),
#' compute the trimmed Hill estimator. If \eqn{k0=k}, the estimator
#' reduces to Hill's estimator for the shape index
#' @param xdat [numeric] vector of positive observations
#' @param k0 [integer] number of largest order statistics, strictly less than \code{k}
#' @param k [integer] number of order statistics for the threshold
#' @param sorted [logical] if \code{TRUE}, data are assumed to be sorted in decreasing order
#' @export
#' @return a scalar with the shape parameter estimate
#' @references Bhattacharya, S., Kallitsis, M. and S. Stoev, (2019) Data-adaptive trimming of the Hill estimator and detection of outliers in the extremes of heavy-tailed data. Electronic Journal of Statistics 13, 1872–1925
shape.trimhill <- function(xdat, k, k0, sorted = FALSE) {
  sorted <- isTRUE(sorted)
  if (length(k0) != 1L) {
    stop(
      "Invalid argument \"k0\": must be an integer of length 1."
    )
  }
  if (length(k) != 1L) {
    stop(
      "Invalid argument \"k\": must be an integer of length 1."
    )
  }
  if (!sorted) {
    xdat <- sort(xdat, decreasing = TRUE, na.last = NA)
  }
  n <- length(xdat)
  k0 <- as.integer(k0)
  k <- as.integer(k)
  stopifnot(k0 >= 0, k0 < k, k < n)
  # Threshold
  th <- xdat[k + 1]
  # Observations should be Pareto-tailed, and positive
  stopifnot(th > 0)
  logy <- log(xdat[(k0 + 1):k]) - log(th)
  (k0 + 1) / (k - k0) * logy[1] + sum(logy[-1]) / (k - k0)
}


#' Threshold stability plots for left-truncated Hill estimators
#'
#' Given a vector of exceedances and some potential choices of \eqn{k} for the threshold, compute the left-truncated Hill estimators for each value of \code{k} and use these to compute the variance and slope of the estimator
#' @inheritParams shape.lthill
#' @param which [string] the type of plot, showing the left-truncated Hill plot on the log, the log of the variance of the estimator, or the log slope
#' @param ... additional parameters for color, etc. to be passed to plot
#' @param log [logical] if \code{TRUE} (default), shows the Hill plot on the log-scale
#' @return an invisible list with lthill, order statistics, the log variance and the log scale.
#' @references Bladt, M., Albrecher, H. & Beirlant, J. (2020) \emph{Threshold selection and trimming in extremes}. Extremes, 23, 629-665 . \doi{10.1007/s10687-020-00385-0}
#' @export
#' @examples
#' xdat <- 10/(1 - runif(n = 1000)) - 10
#' tstab.lthill(xdat = xdat, k = c(50,100,200))
tstab.lthill <- function(
  xdat,
  k,
  which = c("lthill", "var", "slope"),
  log = TRUE,
  ...
) {
  args <- list(...)
  which <- match.arg(which, several.ok = TRUE)
  n <- length(xdat)
  nmax <- min(n, 1000L)
  k <- sort(as.integer(k))
  nk <- length(k)
  kmax <- k[nk]
  if (!is.null(args$col)) {
    cols <- args$col
    if (!length(cols) %in% c(1L, nk)) {
      stop(paste(
        "Optional argument \"col\" must be a scalar or a vector of length",
        nk
      ))
    }
  } else {
    cols <- grey.colors(n = nk, start = 0, end = 0.8)
  }
  stopifnot(kmax < length(xdat), k[1] > 10)
  # old.par <- par(no.readonly = TRUE)
  # on.exit(par(old.par))
  # par(mfrow = c(length(which), 1))
  # Hill-type plot
  xdat <- sort(xdat, decreasing = TRUE, na.last = NA)
  if ("lthill" %in% which) {
    lt <- matrix(NA, ncol = nk, nrow = k[nk] - 4L)
    for (i in seq_along(k)) {
      ks <- k[i]
      lt[1:(ks - 4), i] <- shape.lthill(
        xdat = xdat,
        k = ks,
        k0 = 5:ks,
        sorted = TRUE
      )$shape
    }
    matplot(
      x = 5:kmax,
      y = lt,
      log = ifelse(isTRUE(log), 'x', ''),
      col = cols,
      type = "l",
      ylab = "shape",
      xlab = "number of exceedances (log scale)",
      bty = "l"
    )
    lt <- t(lt)
  } else {
    lt <- NULL
  }
  var_plot <- isTRUE("var" %in% which)
  slope_plot <- isTRUE("slope" %in% which)
  if (var_plot | slope_plot) {
    os <- 10:kmax
    logslope <- logvar <- numeric(kmax - 9)
    for (ki in 1:(kmax - 9)) {
      kv <- 10 + ki - 1
      lth_s <- shape.lthill(
        xdat = xdat,
        k = kv,
        k0 = kv,
        vectorize = TRUE,
        sorted = TRUE
      )$shape
      if (var_plot) {
        logvar[ki] <- log(var(lth_s))
      }
      if (slope_plot) {
        logslope[ki] <- log(abs(coef(lm(lth_s ~ seq_len(kv)))[2]))
      }
    }
  } else {
    os <- NULL
  }
  if (var_plot) {
    plot(
      os,
      logvar,
      type = "l",
      log = ifelse(isTRUE(log), 'x', ''),
      ylab = "log variance",
      xlab = "number of order statistics",
      bty = "l"
    )
  } else {
    logvar <- NULL
  }
  if (slope_plot) {
    plot(
      os,
      logslope,
      log = ifelse(isTRUE(log), 'x', ''),
      type = "l",
      ylab = "log of slope",
      xlab = "number of order statistics",
      bty = "l"
    )
  } else {
    logslope <- NULL
  }
  invisible(
    list(
      lthill = lt,
      os,
      logvar,
      logslope
    )
  )
}

#' Lower truncated Hill threshold selection
#'
#' Given a sample of positive data with Pareto tail, the algorithm computes the optimal number of order statistics that minimizes the variance of the average left truncated tail index estimator, and uses the relationship to the Hill estimator for the Hall class of distributions to derive the optimal number (minimizing the asymptotic mean squared error) of the Hill estimator. The default value for the second order regular variation index is taken to be \eqn{\rho=-1}.
#' @param xdat [vector] positive vector of exceedances
#' @param kmin [int] minimum number of exceedances
#' @param kmax [int] maximum number of exceedances for the estimation of the shape parameter.
#' @param rho [double] scalar for the second order regular variation index, a negative number.
#' @param test [logical] if \code{TRUE}, computes the goodness-of-fit statistic for the model using Monte Carlo
#' @param nsim [int] number of replications for Monte Carlo test, used only if \code{test=TRUE}.
#' @param level [double] confidence level for test
#'
#' @return a list with the number of order statistics for the Hill estimator, \code{k0} and the corresponding shape estimate \code{shape}, the average lower-trimmed Hill estimator \code{shape._lth} and the number of order statistics upon which the latter is based, \code{k0_lth}.
#' @references Bladt, M., Albrecher, H. & Beirlant, J. (2020) \emph{Threshold selection and trimming in extremes}. Extremes, 23, 629-665 . \doi{10.1007/s10687-020-00385-0}
#' @export
thselect.bab <- function(
  xdat,
  kmin = floor(0.2 * length(xdat)),
  kmax = length(xdat) - 1L,
  rho = -1,
  test = FALSE,
  nsim = 999L,
  level = 0.95
) {
  range <- as.integer(c(kmin, kmax))
  stopifnot(
    length(rho) == 1L,
    rho < 0,
    range[1] >= 10,
    range[2] < length(xdat)
  )
  xdat <- sort(xdat, decreasing = TRUE, na.last = NA)
  Cst <- 0.502727
  fcst <- function(p) {
    #TODO add dependence to expint
    if (length(p) == 1L & isTRUE(all.equal(p, -1, check.attributes = FALSE))) {
      e1 <- 0.2193839344
      e1mp <- 0.04890051071
      e1m2p <- 0.01304838109
    } else {
      if (!requireNamespace("expint", quietly = TRUE)) {
        stop(
          "Package \"expint\" must be installed to use this function.",
          call. = FALSE
        )
      }
      e1 <- expint::expint(1)
      e1mp <- expint::expint(1 - p)
      e1m2p <- expint::expint(1 - 2 * p)
    }
    (1 - exp(1 - 2 * p) * (1 - 2 * p) * e1m2p - exp(2 - 2 * p) * e1mp^2) /
      (p^2 * (1 - p)^2) +
      2 /
        (p^2 * (1 - p)) *
        (exp(2 - p) * e1mp * e1 - 1 + exp(1 - p) * (1 - p) * e1mp) +
      (1 - exp(1) * e1 - exp(2) * e1^2) / p^2
  }
  # Check function plot matches the figure in the article
  # ps <- seq(-5,0, length.out = 101)
  # fcst(ps)
  varshape <- shape.lth <- numeric(length = diff(range))
  iter <- 0L
  kseq <- seq(from = range[1], to = range[2], by = 1L)
  for (i in kseq) {
    iter <- iter + 1L
    shape <- shape.lthill(
      xdat = xdat,
      k = i,
      vectorize = TRUE,
      sorted = TRUE
    )$shape
    shape.lth[iter] <- mean(shape)
    varshape[iter] <- var(shape)
  }
  k0_lth <- kseq[which.min(varshape)]
  k0 <- floor(k0_lth * (Cst / ((1 - rho)^2 * fcst(rho)))^(-1 / (1 - 2 * rho)))
  if (isTRUE(test)) {
    level <- as.numeric(level)[1]
    stopifnot(is.finite(level), level < 1, level > 0)
    if (level < 0.5) {
      warning("Level of confidence interval is low.")
    }

    shape_lthill <- shape.lthill(
      xdat,
      k = k0_lth,
      k0 = k0_lth,
      vectorize = TRUE,
      sorted = TRUE
    )$shape
    Tstat <- shape_lthill[-1] / shape_lthill[-k0_lth]
    nsim <- as.integer(nsim)
    stopifnot(nsim >= 19L)
    # Constants
    omega <- function(k) {
      k <- as.integer(k)
      1:(k - 1) * (1 + rev(cumsum(1 / (k:2))))
    }
    om <- omega(k0_lth)
    omr <- c(om[-length(om)] / om[-1], 1)
    R <- matrix(nrow = nsim, ncol = k0_lth - 1)
    for (i in seq_len(nsim)) {
      pp <- cumsum(rexp(k0_lth + 1L))
      lpr <- log(pp[-(k0_lth + 1L)] / pp[k0_lth + 1])
      R[i, ] <- omr * (1 + lpr[-1] / cumsum(lpr)[-length(lpr)])
    }
    envel <- boot::envelope(mat = R, level = level)
    CI <- envel$overall
    # CI <- apply(R, 2, quantile, prob = c(alpha_st/2, 1-alpha_st/2))
    size <- mean(Tstat < CI[1, ] | Tstat > CI[2, ])
    Rst <- (Tstat - CI[1, ]) / (CI[2, ] - CI[1, ])
  } else {
    Rst <- NULL
  }

  res <- list(
    k0 = k0,
    k0_lth = k0_lth,
    thresh0 = xdat[k0],
    shape = mev::shape.hill(xdat, k = k0),
    shape_lth = shape[k0_lth],
    stat = Rst
  )
  class(res) <- "mev_thselect_bab"
  return(invisible(res))
}

#' @export
plot.mev_thselect_bab <- function(x, ...) {
  if (is.null(x$stat)) {
    stop("Rerun \"thselect.bab\" with argument \"test=TRUE\".")
  }
  plot(
    x = 2:(length(x$stat) + 1),
    y = x$stat,
    xlab = "b",
    ylab = "standardized ratio statistic",
    ylim = c(min(min(x$stat), 0), max(max(x$stat), 1)),
    bty = "l",
    type = "l",
    col = ifelse(x$stat < 0 | x$stat > 1, "grey", "black"),
    panel.first = {
      abline(h = c(0, 1), lty = 2)
    }
  )
}

#' @export
print.mev_thselect_bab <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  cat(
    "Threshold selection method: Bladt, Albrecher and Beirlant (2020)\n"
  )
  cat("Selected threshold:", round(x$thresh0, digits), "\n")
  cat("Number of exceedances:", round(x$k0, digits), "\n")
  cat("Shape estimate (Hill):", round(x$shape, digits), "\n\n")
  cat(
    "Number of exceedances of intermediate sequence (lower-trimmed Hill):",
    round(x$k0_lth, digits),
    "\n"
  )
  cat("Shape estimate (lower-trimmed Hill):", round(x$shape_lth, digits), "\n")

  return(invisible(NULL))
}
