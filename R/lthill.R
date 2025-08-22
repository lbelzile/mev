#' Left-trimmed Hill estimator for the shape parameter
#'
#' Given a sample of Pareto-tailed samples (positive tail index),
#' compute the left-trimmed Hill estimator. If \eqn{k0=k}, the estimator
#' reduces to Hill's estimator for the shape index
#' @param xdat [numeric] vector of positive observations
#' @param k0 [integer] scalar number of largest order statistics, no greater than \code{k}
#' @param vectorize [logical] if \code{TRUE}, computes for each value of
#' @param k [integer] number of order statistics for the threshold
#' @param sorted [logical] if \code{TRUE}, data are assumed to be sorted in decreasing order
#' @return a scalar with the shape parameter estimate
#' @examples
#' # Pareto sample
#' n <- 200
#' xdat <- 10/(1 - runif(n)) - 10
#' shape.lthill(xdat = xdat, k = 100, k0 = 1:100)
shape.lthill <- function(xdat, k, k0 = k, vectorize = FALSE, sorted = FALSE) {
  sorted <- isTRUE(sorted)
  if (length(k) != 1L) {
    stop(
      "Invalid argument \"k\": must be an integer of length 1."
    )
  }
  if (length(k0) != 1L) {
    stop(
      "Invalid argument \"k0\": must be an integer of length 1."
    )
  }
  if (!sorted) {
    xdat <- sort(xdat, decreasing = TRUE, na.last = NA)
  }
  n <- length(xdat)
  k0 <- as.integer(k0)
  k <- as.integer(k)
  stopifnot(k0 <= k, k0 >= 1, k < n)
  # Threshold
  th <- xdat[k + 1]
  # Observations should be Pareto-tailed, and positive
  stopifnot(th > 0)
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
    return((cumsum(y) / (1:k0) - log(th)) / denom)
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
#' @return a scalar with the shape parameter estimate
#' @references Bhattacharya, S., Kallitsis, M. and S. Stoev, (2019) Data-adaptive trimming of the Hill estimator and detection of outliers in the extremes of heavy-tailed data. Electronic Journal of Statistics 13, 1872–1925
shape.trimedhill <- function(xdat, k, k0, sorted = FALSE) {
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


lthill.diag <- function(
  xdat,
  k,
  which = c("lthill", "var", "slope"),
  ...
) {
  which <- match.arg(which, several.ok = TRUE)
  n <- length(xdat)
  nmax <- min(n, 1000L)
  k <- sort(as.integer(k))
  kmax <- k[length(k)]
  stopifnot(kmax < length(xdat), k[1] > 10)
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(mfrow = c(length(which), 1))
  # Hill-type plot
  xdat <- sort(xdat, decreasing = TRUE, na.last = NA)
  if ("lthill" %in% which) {
    lt <- matrix(NA, ncol = length(k), nrow = max(k))
    for (i in seq_along(k)) {
      ks <- k[i]
      lt[1:ks, i] <- shape.lthill(
        xdat = xdat,
        k = ks,
        vectorize = TRUE,
        sorted = TRUE
      )
    }
    matplot(
      x = 1:kmax,
      y = lt,
      type = "l",
      ylab = "shape",
      xlab = "number of order statistics"
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
        vectorize = TRUE,
        sorted = TRUE
      )
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
      ylab = "log variance",
      xlab = "number of order statistics"
    )
  } else {
    logvar <- NULL
  }
  if (slope_plot) {
    plot(
      os,
      logslope,
      type = "l",
      ylab = "log variance",
      xlab = "number of order statistics"
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
#' Given a sample of positive data with Pareto tail, the algorithm computes the optimal number of order statistics that minimizes the variance of the average left truncated tail index estimator, and uses the relationship to the Hill estimator for the Hall class of distributions to derive the optimal number (minimizing the asymptotic mean squared error) of the Hill estimator. The default value for the second order regular variation index is taken to be \eqn{rho=-1}.
#' @param xdat vector of observations.
#' @param range [integer] vector of length 2 containing the minimum and maximum number of observations for the estimation of the shape parameter.
#' @param rho second order regular variation index, a negative number.
#'
#' @return a list with the number of order statistics for the Hill estimator, \code{k0} and the corresponding shape estimate \code{shape.hill}, the average left-trimmed Hill estimator \code{shape.lth} and the number of order statistics upon which the latter is based, \code{nexc}.
#' @references Bladt, M., Albrecher, H. & Beirlant, J. (2020) \emph{Threshold selection and trimming in extremes}. Extremes 23, 629-665 . \doi{10.1007/s10687-020-00385-0}
thselect.bab <- function(
  xdat,
  range = c(floor(0.2 * length(xdat)), length(xdat) - 1L),
  rho = -1,
  ...
) {
  range <- sort(as.integer(range))
  stopifnot(
    length(range) == 2L,
    length(rho) == 1L,
    rho < 0,
    range[1] > 10,
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
    )
    shape.lth[iter] <- mean(shape)
    varshape[iter] <- var(shape)
  }
  klth <- kseq[which.min(varshape)]
  k0 <- floor(klth * (Cst / ((1 - rho)^2 * fcst(rho)))^(-1 / (1 - 2 * rho)))
  list(
    k0 = k0,
    nexc = klth,
    shape.hill = mev::shape.hill(xdat, k = k0),
    shape.lth = shape[klth]
  )
}
