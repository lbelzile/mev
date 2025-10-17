#' Exponential regression estimator
#'
#' This function implements the exponential regression estimator of the shape parameter for the case of Pareto tails with positive shape index.
#'
#' @param xdat vector of observations
#' @param k vector of integer, the number of largest observations to consider
#' @param method string; one of \code{bdgm} for the approach of Beirlant, Dierckx, Goegebeur and Matthys (1999)  or \code{fh} for Feuerverger and Hall (1999)
#' @param bounds vector of length 2 giving the bounds for \code{rho}, the second order parameter. Default to \eqn{\rho \in [-5, -0.5]}
#' @references Feuerverger, A. and P. Hall (1999), Estimating a tail exponent by modelling departure from a Pareto distribution, \emph{The Annals of Statistics} 27(\bold{2}), 760-781. <doi:10.1214/aos/1018031215>
#'
#' @details The second-order parameter is difficult to pin down, and while values within \eqn{[-1,0)} are most logical under Hall model, the model parameters become unidentifiable when \eqn{\rho \to 0}. The default constraint restrict \eqn{-5 <\rho < -0.5}, with the upper bound changed to \eqn{-0.25} for sample of sizes larger than 5000 observations. Users can set the value of the bounds for \eqn{\rho} via argument \code{bounds}. The optimization is initialized at the Hill estimator.
#' @export
#' @references Beirlant, J., Dierckx, G., Goegebeur, Y. G. Matthys (1999). Tail Index Estimation and an Exponential Regression Model. \emph{Extremes}, 2, 177–200 (1999). <doi:10.1023/A:1009975020370>
#'
#' @return a data frame with columns
#' \itemize{
#' \item \code{k} number of exceedances
#' \item \code{shape} estimate of the shape parameter
#' \item \code{rho} estimate of the second-order regular variation index
#' \item \code{scale} estimate of the scale parameter
#'
#' }
shape.erm <- function(
  xdat,
  k,
  method = c("bdgm", "fh"),
  bounds = NULL
) {
  method <- match.arg(method)
  k <- sort(as.integer(k))
  xdat <- sort(
    x = xdat[is.finite(xdat)],
    decreasing = TRUE
  )
  n <- length(xdat)
  kmax <- k[length(k)] + 1
  stopifnot(kmax <= n)
  xdat <- as.numeric(xdat[1:kmax])
  stopifnot(xdat[kmax] > 0)
  logdata <- log(xdat)
  Z <- 1:(kmax - 1) * (logdata[1:(kmax - 1)] - logdata[2:kmax])
  shape <- rho <- b <- numeric(length = length(k))
  if (is.null(bounds)) {
    rhobounds <- c(-5, ifelse(kmax < 5000, -0.5, -0.25))
  } else {
    bounds <- as.numeric(bounds[is.finite(bounds)])
    stopifnot(length(bounds) == 2L, isTRUE(all(bounds < 0)))
    rhobounds <- sort(bounds)
  }
  if (method == "bdgm") {
    expreg <- function(par, z) {
      k <- length(z)
      shape <- par[1]
      # This parameter tends to be very small
      bn <- exp(par[2])
      rho <- par[3]
      scale <- shape + bn * ((1:k) / (k + 1))^(-rho)
      if (isTRUE(any(scale < 0))) {
        return(1e10)
      }
      -sum(dexp(z, rate = 1 / scale, log = TRUE))
    }
    ineqfn <- function(par, z) {
      k <- length(z)
      par[1] + exp(par[2]) * (c(1, k) / (k + 1))^(-par[3])
    }
    for (i in seq_along(k)) {
      # The mean of Z provides Hill's estimator
      opt <- Rsolnp::solnp(
        pars = c(mean(Z[1:k[i]]), -2.5, mean(rhobounds)),
        fun = expreg,
        LB = c(1e-8, -50, rhobounds[1]),
        UB = c(
          2,
          3,
          ifelse(is.null(bounds) & k[i] > 5000, -0.25, rhobounds[2])
        ),
        ineqfun = ineqfn,
        ineqLB = rep(0, 2),
        ineqUB = rep(1e10, 2),
        z = Z[1:k[i]],
        control = list(trace = 0)
      )
      b[i] <- exp(opt$par[2])
      shape[i] <- opt$par[1]
      rho[i] <- opt$par[3]
    }
  } else if (method == "fh") {
    expreg <- function(par, z) {
      k <- length(z)
      shape <- par[1]
      bn <- exp(par[2])
      rho <- par[3]
      scale <- shape * exp(bn * ((1:k) / (k + 1))^(-rho))
      -sum(dexp(z, rate = 1 / scale, log = TRUE))
    }
    for (i in seq_along(k)) {
      # if(i > 1){
      #   start <- opt$pars
      #   if(start[1] < 1e-4){
      #     start[1] <- 1e-1
      #   }
      #   if(start[3] < -0.98 | start[3] > -1e-2){
      #     start[3] <- -0.1
      #   }
      # }
      start <- c(mean(Z[1:k[i]]), -2.5, mean(rhobounds))
      opt <- Rsolnp::solnp(
        pars = start,
        fun = expreg,
        LB = c(0, -25, rhobounds[1]),
        UB = c(
          2,
          3,
          ifelse(is.null(bounds) & k[i] > 5000, -0.25, rhobounds[2])
        ),
        z = Z[1:k[i]],
        control = list(trace = 0)
      )
      b[i] <- exp(opt$par[2]) / opt$par[1]
      shape[i] <- opt$par[1]
      rho[i] <- opt$par[3]
    }
  }
  data.frame(
    k = k,
    shape = shape,
    scale = b,
    rho = pmin(0, rho)
  )
}
