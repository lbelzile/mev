#' Threshold selection via SAMSEE
#'
#' Smooth asymptotic mean squared error estimator
#' of Schneider et al. (2021) for threshold selection.
#' The implementation uses a second-order regular variation index of -1
#'
#' @references Schneider, L.F., Krajina, A. & Krivobokova, T. (2021). \emph{Threshold selection in univariate extreme value analysis}, Extremes, \bold{24}, 881–913 \doi{10.1007/s10687-021-00405-7}
#' @param xdat vector of positive exceedances
#' @return a list with elements
#' \describe{
#' \item{\code{k0}}{optimal number of exceedances}
#' \item{\code{shape}}{Hill estimator of the tail index}
#' \item{\code{thresh0}}{selected threshold}
#' }
#' @export
thselect.samsee <- function(xdat) {
  # Keep log exceedances and sort
  xdat <- xdat[is.finite(xdat) & xdat > 0]
  n <- length(xdat)
  logdat <- sort(log(xdat), decreasing = TRUE)
  cumlogdat <- cumsum(logdat)
  k <- 1:(n - 1)
  # Hill estimator for k=2, ..., n
  gamma_hill <- cumlogdat[k] / k - logdat[k + 1]
  # Compute square log spacing
  Mn <- vector(mode = "numeric", length = length(gamma_hill))
  for (i in seq_along(Mn)) {
    Mn[i] <- mean((logdat[1:k[i]] - logdat[k[i] + 1])^2)
  }
  # Compute de Vries and generalized jackknife estimators
  gamma_v <- 0.5 * Mn / gamma_hill
  gamma_gj <- 2 * gamma_v - gamma_hill
  # Bias via averaging of Hill's estimator
  mean_gamma_hill <- cumsum(gamma_hill) / k

  # Minimum value for K > k
  Kmin <- 5L
  # Compute MSE of gamma_v using bias term
  AD <- vector(mode = "numeric", length = length(gamma_hill) - Kmin)
  # Compute bar_gamma(1, K) to bar_gamma(K-1,K)
  gamma_kK <- function(K, hill = gamma_hill) {
    rev(cumsum(hill[K:1])) / (K - 1:K + 1)
  }
  # Compute mean squared error approx
  for (i in seq_along(AD)) {
    K <- i + Kmin - 1
    AD[i] <- mean(
      (gamma_v[1:K] + gamma_kK(K) - mean_gamma_hill[K] - gamma_hill[1:K])^2
    )
  }
  # Compute average 'derivative' of AD around K
  # via finite differences (vectorized)
  i <- 3:(length(AD) - 2)
  D_AD <- abs(AD[i + 2] - AD[i]) /
    2 +
    abs(AD[i + 1] - AD[i]) +
    abs(AD[i - 2] - AD[i]) / 2 +
    abs(AD[i - 1] - AD[i])
  # Find optimal K
  Kstar <- which.min(D_AD) + Kmin
  # Compute SAMSEE
  SAMSEE <- gamma_gj[Kstar]^2 /
    (1:Kstar) +
    4 * (gamma_kK(Kstar) - mean_gamma_hill[Kstar])^2
  kstar <- which.min(SAMSEE[-c(1, length(SAMSEE))]) + 1L
  res <- list(
    k0 = kstar,
    shape = gamma_hill[kstar],
    thresh0 = exp(logdat[kstar + 1])
    # shape.hill = gamma_hill,
    # shape.v = gamma_v,
    # shape.gj = gamma_gj
  )
  class(res) <- "mev_thselect_samsee"
  return(invisible(res))
}

#' @export
print.mev_thselect_samsee <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  cat(
    "Threshold selection method: Schneider, Krajina, Krivobokova (2021)\n"
  )
  cat("Smooth asymptotic mean squared error estimator\n")
  cat("Selected threshold:", round(x$thresh0, digits), "\n")
  cat("Number of exceedances:", round(x$k0, digits), "\n")
  cat("Shape estimate:", round(x$shape, digits), "\n")
  return(invisible(NULL))
}


#' Minimum distance threshold selection procedure
#'
#'
#' @references Clauset, A., Shalizi, C.R. and Newman, M.E.J. (2009). \emph{Power-Law Distributions in Empirical Data}. SIAM Review. Society for Industrial and Applied Mathematics, \bold{51}, 661-703, \doi{10.1137/070710111}
#' @param xdat vector of positive exceedances
#' @return a list with components
#' \describe{
#' \item{\code{k0}}{order statistic corresponding to threshold (number of exceedances)}
#' \item{\code{shape}}{Hill's estimator of the tail index based on k0 exceedances}
#' \item{\code{thresh0}}{numerical value of the threshold, the n-k0+1 order statistic of the original sample}
#' }
#' @export
thselect.mdps <- function(xdat) {
  xdat <- xdat[is.finite(xdat) & xdat > 0]
  n <- length(xdat)
  xdat <- sort(xdat, decreasing = TRUE)
  logdat <- log(xdat)
  cumlogdat <- cumsum(logdat)
  kseq <- 2:n
  # Hill estimator for k=2, ..., n
  hill <- cumlogdat[kseq - 1] / (kseq - 1) - logdat[kseq]
  # browser()
  Dk <- numeric(length = n - 1L)
  for (k in seq_along(Dk)) {
    j <- 1:(k - 1)
    powS <- (xdat[j] / xdat[k])^(-1 / hill[k])
    Dk[k] <- max(powS - (j - 1) / (k - 1), j / (k - 1) - powS)
  }
  # plot(log(Dk))
  k0 <- which.min(Dk)
  res <- list(
    k0 = k0,
    shape = hill[k0],
    thresh0 = xdat[k0]
  )
  class(res) <- "mev_thselect_mdps"
  return(invisible(res))
}

#' @export
print.mev_thselect_mdps <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  cat(
    "Threshold selection method: Clauset, Shalizi and Newman (2009)\n"
  )
  cat("Minimum distance procedure (MDPS)\n")
  cat("Selected threshold:", round(x$thresh0, digits), "\n")
  cat("Number of exceedances:", round(x$k0, digits), "\n")
  cat("Shape estimate:", round(x$shape, digits), "\n")
  return(invisible(NULL))
}


#' Prediction error C-criterion threshold selection method
#'
#' This function computes the non-robust Pareto prediction error of Dupuis and Victoria-Feser (2003), termed C-criterion, for the Hill estimator of the shape parameter. The threshold returned is the value of the threshold, taken from order statistics, that minimizes the average prediction error.
#'
#' @references Dupuis, D.J. and M.-P. Victoria-Feser (2003). A Prediction Error Criterion for Choosing the Lower Quantile in Pareto Index Estimation, University of Geneva, technical report, \url{https://archive-ouverte.unige.ch/unige:5789}.
#'
#' @param xdat vector of observations
#' @param kmax maximum number of order statistics to consider. Default to sample size if left unspecified.
#' @return a list with the number of exceedances \code{k}, the chosen threshold \code{thresh0} and the corresponding Hill estimator shape estimate \code{shape}.
#' @export
thselect.pec <- function(
  xdat,
  kmax
) {
  xdat <- sort(xdat[is.finite(xdat)])
  if (missing(kmax)) {
    kmax <- length(xdat)
  } else {
    kmax <- as.integer(kmax)
    stopifnot(kmax <= length(xdat))
  }
  stopifnot(xdat[kmax] > 0)
  kseq <- 10:kmax
  hill <- mev::shape.hill(
    xdat = xdat,
    k = kseq
  )$shape
  predcrit <- numeric(length(kseq))
  for (i in seq_along(kseq)) {
    k <- kseq[i]
    x0 <- xdat[k] + 1L
    #w <- 1 / cumsum(1 / (k:1)^2)
    w <- (k + 1 - 1:k) / (1:k)
    predcrit[i] <- -1 +
      mean(
        w *
          (hill[i]^2 *
            ((log(xdat[1:k]) - log(x0)) +
              1 / hill[i] * (log((k + 1 - 1:k)) - log(k + 1)))^2 +
            2 * (log((k + 1 - 1:k)) - log(k + 1))^2)
      )
  }
  # plot(kseq, predcrit)
  i0 <- which.min(predcrit)
  res <- list(
    k0 = kseq[i0],
    thresh0 = xdat[length(xdat) - kseq[i0]],
    shape = hill[i0]
  )
  class(res) <- "mev_thselect_pec"
  return(invisible(res))
}

#' @export
print.mev_thselect_pec <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  cat(
    "Threshold selection method: Dupuis and Victoria-Feser (2003)\n"
  )
  cat("Prediction error C-criterion (non-robust version)\n")
  cat("Selected threshold:", round(x$thresh0, digits), "\n")
  cat("Number of exceedances:", round(x$k0, digits), "\n")
  cat("Shape estimate:", round(x$shape, digits), "\n")
  return(invisible(NULL))
}


#' Threshold selection by shape mean square error minimization
#'
#' Use a semiparametric bootstrap to calculate the mean squared error
#' of the shape parameter using maximum likelihood for different thresholds, and return the one that minimize the mean squared error.
#'
#' @param xdat vector of observations
#' @param thresh vector of thresholds
#' @param B number of bootstrap replications
#' @export
#' @return an object of class \code{mev_thselect_cbm} containing
#' \itemize{
#' \item{\code{thresh}: ordered vector of candidate thresholds}
#' \item{\code{thresh0}: selected threshold}
#' \item{\code{shape}: shape parameter coefficient estimates at each threshold}
#' \item{\code{nexc}: number of exceedances at each threshold}
#' \item{\code{bias}: vector of bootstrap bias estimates}
#' \item{\code{var}: vector of bootstrap variance estimates}
#' \item{\code{mse}: vector of mean squared error bootstrap estimates}
#' }
#' @references Caers, J., Beirlant, J. and Maes, M.A. (1999). Statistics for Modeling Heavy Tailed Distributions in Geology: Part I. Methodology. \emph{Mathematical Geology}, 31, 391–410. <doi:10.1023/A:1007538624271>
#' @examples
#' set.seed(2025)
#' xdat <- rnorm(1000)
#' thresh <- qnorm(c(0.8, 0.9, 0.95))
#' thselect.cbm(xdat, thresh, B = 50)
thselect.cbm <- function(xdat, thresh, B = 100) {
  thresh <- sort(thresh)
  xdat <- sort(xdat[is.finite(xdat)])
  n <- length(xdat)
  nth <- length(thresh)
  stopifnot(thresh[nth] < xdat[n])
  B <- as.integer(B)
  stopifnot(B >= 19L)
  # empdist <- ecdf(xdat)(thresh)
  nexc <- sapply(thresh, function(u) {
    sum(xdat > u)
  })
  pu <- 1 - nexc / n
  boot_shape <- matrix(nrow = B, ncol = nth)
  shapes <- bias <- boot_mean <- varia <- numeric(nth)
  for (i in seq_along(thresh)) {
    gpd_fit <- mev::fit.gpd(xdat = xdat, thresh = thresh[i])
    shapes[i] <- coef(gpd_fit)['shape']
    for (b in 1:B) {
      boot_u <- runif(n)
      above <- which(boot_u > pu[i])
      boot_samp <- mev::qgp(
        p = (boot_u[above] - pu[i]) / (1 - pu[i]),
        scale = coef(gpd_fit)['scale'],
        shape = coef(gpd_fit)['shape']
      )
      boot_shape[b, i] <- coef(mev::fit.gpd(boot_samp, thresh = 0))['shape']
    }
  }
  # Compute MSE
  bootmean <- colMeans(boot_shape)
  bias <- bootmean - shapes
  varia <- apply(boot_shape, 2, var)
  mse <- bias^2 + varia
  ret <- list(
    thresh = thresh,
    thresh0 = thresh[which.min(mse)],
    nexc = nexc,
    shape = shapes,
    bias = bias,
    var = varia,
    mse = mse
  )
  class(ret) <- "mev_thselect_cbm"
  return(invisible(ret))
}

#' @export
plot.mev_thselect_cbm <- function(x, ...) {
  plot(
    x$thresh,
    x$mse,
    type = "b",
    ylim = c(0, max(x$mse)),
    yaxs = "i",
    xlab = "threshold",
    ylab = "mean squared error",
    panel.first = {
      abline(v = x$thresh0, lty = 3, col = "gray50")
    }
  )
  lines(x$thresh, x$bias^2, lty = 2)
  points(x$thresh, x$bias^2, pch = 2)
  lines(x$thresh, x$var, lty = 3)
  points(x$thresh, x$var, pch = 3)
  return(invisible(NULL))
}

#' @export
print.mev_thselect_cbm <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  thind <- which(x$thresh %in% x$thresh0)
  cat(
    "Threshold selection method: Caers, Beirlant and Maes (1999)\n"
  )
  cat("Bootstrap mean square error minimization\n")
  cat("Selected threshold:", round(x$thresh0, digits), "\n")
  cat("Number of exceedances:", round(x$nexc[thind], digits), "\n")
  cat("Shape estimate:", round(x$shape[thind], digits), "\n")
  return(invisible(NULL))
}
#' Threshold selection via minimization of the weighted Cramér-von Mises
#'
#' For a Pareto-type sample, return the threshold that
#' minimizes a weighted Cramér-von Mises criterion for the
#' exponential sample with scale \eqn{H_{n, n_u}} and
#' the log increments.
#'

#' @references Goegebeur , Y., Beirlant , J., and de Wet , T. (2008). Linking Pareto-Tail Kernel Goodness-of-fit Statistics with Tail Index at Optimal Threshold and Second Order Estimation. REVSTAT-Statistical Journal, 6(\bold{1}), 51–69. <doi:10.57805/revstat.v6i1.57>
thselect.wcvm <- function(xdat, k) {
  xdat <- sort(xdat, decreasing = TRUE)
  xdat <- xdat[is.finite(xdat) & xdat > 0]
  if (missing(k)) {
    shape_hill <- shape.hill(xdat)
  } else {
    shape_hill <- shape.hill(xdat, k = k)
  }
  kseq <- shape_hill$k
  nk <- length(kseq)
  logxdat <- log(xdat)
  criterion <- numeric(length = nk)
  for (i in seq_len(nk)) {
    nexc <- kseq[i]
    Hnu <- shape_hill$shape[i]
    criterion[i] <- mean(
      1:nexc /
        (nexc - 1:nexc + 1) *
        (logxdat[1:nexc] -
          logxdat[nexc + 1] +
          Hnu * log(1:nexc / (nexc + 1)))^2
    ) /
      Hnu
  }
  k0ind <- which.min(criterion)
  out <- list(
    k0 = kseq[k0ind],
    shape = shape_hill$shape[k0ind],
    thresh = xdat[kseq[k0ind] + 1],
    criterion = data.frame(k = kseq, crit = criterion)
  )
  class(out) <- "mev_thselect_wcvm"
  return(invisible(out))
}


#' @export
print.mev_thselect_wcvm <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  cat(
    "Threshold selection method: Goegebeur, Beirlant and de Weit (2008) \n Weighted Cramer-von Mises distance for log exceedances\n"
  )
  cat("Selected threshold:", round(x$thresh, digits), "\n")
  cat("Shape parameter:", round(x$shape, digits), "\n")
  cat("Number of exceedances:", round(x$k0, digits), "\n")
  return(invisible(NULL))
}


#' @export
plot.mev_thselect_wcvm <- function(x, ...) {
  args <- list(...)
  args$x <- x$criterion$k
  args$y <- log(x$criterion$crit)
  if (is.null(args$ylab)) {
    args$ylab <- "log of weighted Cramér-von Mises distance"
  }
  if (is.null(args$xlab)) {
    args$xlab <- "number of exceedances"
  }
  if (is.null(args$bty)) {
    args$bty <- "l"
  }
  if (is.null(args$type)) {
    args$type <- "l"
  }
  do.call(plot, args = args)
  abline(v = x$k0, col = "grey90", lty = 2)

  return(invisible(NULL))
}

#' Pickands' order statistics threshold selection method
#'
#' Restricting to the largest fourth of the data, returns the number of exceedances that minimizes the Kolmogorov-Smirnov statistic, i.e., the maximum absolute difference between the estimated generalized Pareto and the empirical distribution of exceedances. Relative to the paper, different estimation methods are proposed.
#'
#' @references  James Pickands III (1975). \emph{Statistical inference using extreme order statistics}, Annals of Statistics, 3(\bold{1}) 119-131, \doi{10.1214/aos/1176343003}
#' @param xdat [numeric] vector of observations
#' @param method [string] estimation method, either the quartiles of Pickands (1975), maximum likelihood, probability weighted moments or L-moments
#' @param thresh [numeric] vector of candidate thresholds. If missing, defaults to order statistics from the 10th to a quarter of the sample size.
#' @return a list with components
#' \itemize{
#' \item \code{k0}: number of exceedances
#' \item \code{thresh0}: selected threshold returned by the procedure
#' \item \code{thresh}: vector of candidate thresholds
#' \item \code{dist}; vector of Kolmogorov-Smirnoff distance
#' \item \code{method}; string for the estimation method
#' \item \code{scale}: estimated scale parameter at the chosen threshold
#' \item \code{shape}: estimated shape parameter at the chosen threshold
#' }
#' @note The quartiles estimator of Pickands is robust, but very inefficient. It is provided for historical reasons.
#' @export
thselect.pickands <- function(
  xdat,
  thresh,
  method = c("mle", "lmom", "quartiles")
) {
  method <- match.arg(method)
  xdat <- as.vector(xdat)
  xdat <- sort(xdat[is.finite(xdat)], decreasing = TRUE)
  n <- length(xdat)
  if (missing(thresh)) {
    mmax <- floor(n / 4)
    mmin <- 10L
    m_candidate <- mmax:mmin
    thresh <- xdat[m_candidate]
  } else {
    thresh <- sort(thresh)
    m_candidate <- sapply(thresh, function(th) {
      sum(xdat > th)
    })
  }
  shape <- scale <- dist <- numeric(length(m_candidate))
  for (i in seq_along(m_candidate)) {
    m <- m_candidate[i]
    samp <- xdat[seq_len(m - 1)] - thresh[i]
    if (method == "quartiles") {
      quants <- as.numeric(quantile(samp, probs = c(0.5, 0.75)))
      shape[i] <- (log(diff(quants)) - log(quants[1])) / log(2)
      scale[i] <- quants[1] * shape[i] / (2^shape[i] - 1)
    } else if (method == "mle") {
      coefs <- coef(mev::fit.gpd(xdat = samp, threshold = 0))
      shape[i] <- coefs['shape']
      scale[i] <- coefs['scale']
    } else if (method == "lmom") {
      pars <- gpd.lmom(rev(samp), sorted = TRUE, Lskew = FALSE)
      scale[i] <- pars['scale']
      shape[i] <- pars['shape']
    } #else if (method == "lmom") {
    #pars <- gpd.lmom(samp, sorted = TRUE, Lskew = TRUE)
    #scale[i] <- pars['scale']
    #shape[i] <- pars['shape']
    #}
    dist[i] <- max(abs(
      rank(samp) /
        length(samp) -
        mev::pgp(samp, scale = scale[i], shape = shape[i])
    ))
  }
  ind <- which.min(dist)
  k0 <- as.integer(m_candidate[ind])
  ret <- list(
    k0 = k0,
    thresh0 = thresh[ind],
    dist = dist,
    thresh = thresh,
    method = method,
    scale = scale[ind],
    shape = shape[ind]
  )
  class(ret) <- "mev_thselect_pickands"
  return(invisible(ret))
}

#' @export
plot.mev_thselect_pickands <- function(x, ...) {
  plot(
    x = x$thresh,
    y = x$dist,
    pch = 20,
    xlab = "threshold",
    ylab = "Kolmogorov-Smirnoff distance",
    bty = "l",
    panel.first = {
      abline(v = x$thresh0, col = "gray", lty = 2)
    }
  )
}

#' @export
print.mev_thselect_pickands <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  cat(
    "Threshold selection method: Pickands (1975)\n Kolmogorov-Smirnoff goodness-of-fit statistic\n"
  )
  cat(
    "Estimation method:",
    switch(
      x$method,
      mle = "maximum likelihood",
      lmom = "L-moments",
      quartile = "quartile"
    ),
    "\n"
  )
  cat("Selected threshold:", round(x$thresh0, digits), "\n")
  cat("Number of exceedances:", round(x$k0, digits), "\n")
  return(invisible(NULL))
}
