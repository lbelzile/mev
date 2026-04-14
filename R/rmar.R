#' Simulation from first-order max-autoregressive processes
#'
#' Generate data from stationary sequences for extremes for
#' non-negative shapes, following Tavares (1977) for the Gumbel case.
#'
#' The models are parametrized in terms of extremal index \eqn{\theta \in (0,1]}.
#'
#' When \code{shape = 0}, the stationary process has unit Gumbel margins.
#' When \code{shape > 0}, the margins have Frechet margins with distribution \eqn{F(x) = \exp(-x^{-1/\xi}).}
#' @param n sample size
#' @param theta extremal index, a value in (0,1]
#' @param shape non-negative shape parameter of the GEV
#' @return a vector of length \code{n} drawn from the stationary distribution.
#' @examples
#' X1 <- rmar1(n = 1000, theta = 0.5)
#' X2 <- rmar1(n = 1000, theta = 0.2, shape = 0.2)
#' par(mfrow = c(1, 2))
#' plot(X1)
#' plot(X2)
#' xacf(X1, qlev = 0.9)
#' xacf(X2, qlev = 0.9)
#' @references Valadares Tavares, L. (1977). The Exact Distribution of Extremes of a Non-Gaussian Process. \emph{Stochastic Processes and Their Applications} (2): 151-56.\doi{10.1016/0304-4149(77)90026-6.12}
#' @references Davis, Richard A., and Sidney I. Resnick (1989). Basic Properties and Prediction of Max-ARMA Processes, \emph{Advances in Applied Probability}, 21 (\bold{4}): 781–803. \doi{10.2307/1427767}.
#' @export
rmar1 <- function(n, theta, shape = 0) {
  n <- as.integer(n[1])
  stopifnot(n >= 1)
  stopifnot(length(par) == 1L, theta > 0, theta <= 1)
  Y <- numeric(n + 1)
  if (abs(shape) < 1e-6) {
    if (theta == 1) {
      return(mev::rgev(n = n, shape = 0, loc = 0, scale = 1))
    }
    par <- log(theta) - log(1 - theta)
    cst <- log(1 - theta)
    # Gumbel innovations
    Z <- mev::rgev(n = n + 1, shape = 0, loc = par, scale = 1)
    Y[1] <- Z[1] - par # std. gumbel
    for (i in seq_len(n)) {
      Y[i + 1] <- max(Y[i], Z[i + 1]) + cst
    }
    # Maximum of m variables is Weibull or GEV(loc = log1p((m-1)*exp(-(cst-par))), scale = 1, shape = 0)
    # Marginal is Gumbel
    # } else if (type == "tavares80") {
    #   stopifnot(length(par) == 2L, isTRUE(all(par > 0)))
    #   K <- (par[2] + par[1]) / par[2]
    #   Z <- rexp(n = n, rate = par[1])
    #   Y[1] <- rexp(n = 1, rate = par[2])
    #   for (i in seq_len(n)) {
    #     Y[i + 1] <- K * min(Y[i], Z[i])
    #   }
  } else {
    if (theta == 1) {
      return(mev::rgev(n = n, shape = shape, loc = 1, scale = shape))
    }
    par <- (1 - theta)^shape
    Y[1] <- mev::rgev(
      n = 1,
      loc = 1,
      scale = shape,
      shape = shape
    )
    Z <- mev::rgev(
      n,
      shape = shape,
      scale = shape * theta^(shape),
      loc = theta^shape
    )
    for (i in seq_len(n)) {
      Y[i + 1] <- max(par * Y[i], Z[i])
    }
  }
  return(Y[-1])
}
