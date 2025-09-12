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
    shape = 1 / hill[k0],
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
