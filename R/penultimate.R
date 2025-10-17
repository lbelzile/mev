#' Smith's penultimate approximations
#'
#' The function takes as arguments the distribution and density functions. There are two options:
#' \code{method='bm'} yields block maxima and \code{method='pot'} threshold exceedances.
#' For \code{method='bm'}, the user should provide in such case the block sizes via the
#' argument \code{m}, whereas if \code{method='pot'}, a vector of threshold values should be
#' provided. The other argument (\code{u} or \code{m} depending on the method) is ignored.
#'
#' Alternatively, the user can provide functions \code{densF}, \code{quantF} and \code{distF} for the density,
#' quantile function and distribution functions, respectively. The user can also supply the derivative
#' of the density function, \code{ddensF}. If the latter is missing, it will be approximated using finite-differences.
#'
#'
#' @details For \code{method = "pot"}, the function computes the reciprocal hazard and its derivative on the log scale to avoid numerical overflow. Thus, the density function should have argument \code{log} and the distribution function arguments \code{log.p} and \code{lower.tail}, respectively.
#' @param family the name of the parametric family. Will be used to obtain \code{dfamily}, \code{pfamily}, \code{qfamily}
#' @param method either block maxima (\code{'bm'}) or peaks-over-threshold (\code{'pot'}) are supported
#' @param u vector of thresholds for method \code{'pot'}
#' @param qu vector of quantiles for method \code{'pot'}. Ignored if argument \code{u} is provided.
#' @param m vector of block sizes for method \code{'bm'}
#' @param returnList logical; should the arguments be returned as a list or as a matrix of parameter
#' @param ... additional arguments passed to \code{densF} and \code{distF}
#' @author Leo Belzile
#' @importFrom methods formalArgs
#' @import stats
#' @return either a vector, a matrix if either \code{length(m)>1} or \code{length(u)>1} or a list (if \code{returnList}) containing
#' \itemize{
#' \item \code{loc}: location parameters (\code{method='bm'})
#' \item \code{scale}: scale parameters
#' \item \code{shape}: shape parameters
#' \item \code{u}: thresholds (if \code{method='pot'}), percentile corresponding to threshold (if \code{method='pot'})
#' \item \code{m}: block sizes (if \code{method='bm'})
#' }
#' @references Smith, R.L. (1987). Approximations in extreme value theory. \emph{Technical report 205}, Center for Stochastic Process, University of North Carolina, 1--34.
#' @examples
#' #Threshold exceedance for Normal variables
#' qu <- seq(1,5,by=0.02)
#' penult <- smith.penult(family = "norm", ddensF=function(x){-x*dnorm(x)},
#'    method = 'pot', u = qu)
#' plot(qu, penult$shape, type='l', xlab='Quantile',
#'    ylab='Penultimate shape', ylim=c(-0.5,0))
#' #Block maxima for Gamma variables -
#' #User must provide arguments for shape (or rate)
#' m <- seq(30, 3650, by=30)
#' penult <- smith.penult(family = 'gamma', method = 'bm', m=m, shape=0.1)
#' plot(m, penult$shape, type='l', xlab='Quantile', ylab='Penultimate shape')
#'
#' #Comparing density of GEV approximation with true density of maxima
#' m <- 100 #block of size 100
#' p <- smith.penult(family='norm',
#'    ddensF=function(x){-x*dnorm(x)}, method='bm', m=m, returnList=FALSE)
#' x <- seq(1, 5, by = 0.01)
#' plot(x, m*dnorm(x)*exp((m-1)*pnorm(x,log.p=TRUE)),type='l', ylab='Density',
#' main='Distribution of the maxima of\n 100 standard normal variates')
#' lines(x, mev::dgev(x,loc=p[1], scale=p[2], shape=0),col=2)
#' lines(x, mev::dgev(x,loc=p[1], scale=p[2], shape=p[3]),col=3)
#' legend(x = 'topright',lty = c(1,1,1,1), col = c(1,2,3,4),
#'    legend = c('exact', 'ultimate', 'penultimate'), bty = 'n')
#' @export
#' @keywords internal
smith.penult <- function(
  family,
  method = c("bm", "pot"),
  u = NULL,
  qu = NULL,
  m = NULL,
  returnList = TRUE,
  ...
) {
  ellips <- list(...)
  #Compatibility condition with version 1.12 and before, via ellipsis argument
  if (!is.null(ellips$model) && length(method) == 2) {
    method <- ellips$model
  }
  # Redefine density, quantile and distribution functions from family
  if (!missing(family)) {
    densF <- paste0("d", family)
    distF <- paste0("p", family)
    quantF <- paste0("q", family)
    computeQuant <- TRUE
  } else {
    #compatibility - copy from previous
    if (any(c(is.null(ellips$densF), is.null(ellips$distF)))) {
      stop("Argument \"family\" missing.")
    } else {
      densF <- ellips$densF
      distF <- ellips$distF
      if (!is.null(ellips$quantF)) {
        quantF <- ellips$quantF
        computeQuant <- TRUE
      } else {
        computeQuant <- FALSE
      }
    }
  }
  # Matching extra arguments with additional ones passed via ellipsis
  # Which are formals of the function
  ellips$log <- ellips$log.p <- ellips$lower.tail <- NULL
  indf <- names(ellips) %in% formalArgs(densF)
  indF <- names(ellips) %in% formalArgs(distF)
  if (!is.null(quantF)) {
    indQ <- names(ellips) %in% formalArgs(quantF)
  }
  fn.arg <- ellips[which(indf * (indf == indF) == 1)]
  method <- match.arg(method)
  logRecipHaz <- FALSE
  if (
    isTRUE(all(
      "log" %in% formalArgs(densF),
      c("log.p", "lower.tail") %in% formalArgs(distF)
    ))
  ) {
    logRecipHaz <- TRUE
  }

  # Distribution function, density and density derivative
  densFn <- function(x) {
    do.call(densF, c(x = x, fn.arg))
  }
  distFn <- function(x) {
    do.call(distF, c(q = x, fn.arg))
  }
  if (is.null(ellips$ddensF)) {
    ddensFn <- function(x) {
      tol <- 6e-6
      #numDeriv::grad(densFn, x = x, method = "Richardson")
      (densFn(x + tol) - densFn(x - tol)) / (2 * tol)
    }
  } else {
    ddensF <- ellips$ddensF
    if (!inherits(ddensF, "function")) {
      stop("Invalid arguments. Please provide valid functions.")
    }
    ddensFn <- function(x) {
      do.call(ddensF, c(x, fn.arg))
    }
  }
  # not checking for concordance via numerical derivatives, but could be added Block maxima
  if (method == "bm") {
    if (is.null(m)) {
      stop("Sequence of block size must be provided.")
    }
    # Normalizing sequence
    bm <- sapply(m, function(n) {
      if (!is.null(quantF)) {
        bmroot <- do.call(quantF, c(p = exp(-1 / n), fn.arg))
      } else {
        bmroot <- uniroot(
          f = function(bn) {
            do.call(distF, c(bn, fn.arg)) - exp(-1 / n)
          },
          lower = -1e+15,
          upper = 1e+15,
          f.lower = -exp(-1 / n),
          f.upper = 1 - exp(-1 / n),
          tol = 1e-08,
          check.conv = TRUE
        )
        if (abs(bmroot$f.root) < 1e-05) {
          return(bmroot$root)
        } else {
          warning("Could not find \"bm\" using numerical root finder.")
          return(NA)
        }
      }
    })
    # Penultimate scale and shape functions
    phi <- function(x) {
      -sapply(x, function(xval) {
        distFn(xval) * log(distFn(xval)) / densFn(xval)
      })
    }
    dphi <- function(x) {
      sapply(x, function(xval) {
        -(1 + log(distFn(xval))) +
          distFn(xval) * log(distFn(xval)) * ddensFn(xval) / (densFn(xval)^2)
      })
    }
    if (returnList) {
      params <- list(loc = bm, scale = phi(bm), shape = dphi(bm), m = m)
    } else {
      params <- cbind(loc = bm, scale = phi(bm), shape = dphi(bm), m = m)
      if (nrow(params) == 1L) {
        params <- params[1, ]
      }
    }
    return(params)
  } else if (method == "pot") {
    if (is.null(u) & is.null(qu)) {
      stop("Sequence of thresholds must be provided.")
    } else if (is.null(u) & !is.null(qu)) {
      if (computeQuant) {
        u <- sapply(qu, function(p) {
          do.call(quantF, c(p = p, fn.arg))
        })
      } else {
        u <- rep(NA, length(qu))
      }
    } else if (!is.null(u) & is.null(qu)) {
      qu <- sapply(u, function(q) {
        distFn(x = q)
      })
    }
    if (logRecipHaz) {
      # Compute reciprocal hazard on log scale to avoid overflow
      phi <- function(x) {
        exp(
          do.call(
            distF,
            c(list(q = x), fn.arg, lower.tail = FALSE, log.p = TRUE)
          ) -
            do.call(densF, c(list(x = x), fn.arg, log = TRUE))
        )
      }
      dphi <- function(x) {
        sapply(x, function(xval) {
          -1 -
            sign(ddensFn(xval)) *
              exp(
                do.call(
                  distF,
                  c(q = xval, fn.arg, lower.tail = FALSE, log.p = TRUE)
                ) +
                  log(abs(ddensFn(xval))) -
                  2 * do.call(densF, c(x = xval, fn.arg, log = TRUE))
              )
        })
      }
    } else {
      phi <- function(x) {
        sapply(x, function(xval) {
          (1 - distFn(xval)) / densFn(xval)
        })
      }
      dphi <- function(x) {
        sapply(x, function(xval) {
          -1 -
            (1 - distFn(xval)) *
              ddensFn(xval) /
              (densFn(xval)^2)
        })
      }
    }
    if (returnList) {
      params <- list(u = u, scale = phi(u), shape = dphi(u), qu = qu)
    } else {
      params <- cbind(u = u, scale = phi(u), shape = dphi(u), qu = qu)
      if (nrow(params) == 1L) {
        params <- params[1, ]
      }
    }
    return(params)
  }
  ## The approximations are for \eqn{F^m(x)} for the GEV distribution and for \eqn{1-F(u+x)}{1-F(u)} for the GP distribution.
}


#' Smith's penultimate approximations
#'
#' The function takes as arguments the distribution and density functions. There are two options:
#' \code{method='bm'} yields block maxima and \code{method='pot'} threshold exceedances.
#' For \code{method='bm'}, the user should provide in such case the block sizes via the
#' argument \code{m}, whereas if \code{method='pot'}, a vector of threshold values should be
#' provided. The other argument (\code{thresh} or \code{m} depending on the method) is ignored.
#'
#' Alternatively, the user can provide functions \code{densF}, \code{quantF} and \code{distF} for the density,
#' quantile function and distribution functions, respectively. The user can also supply the derivative
#' of the density function, \code{ddensF}. If the latter is missing, it will be approximated using finite-differences.
#'
#'
#' @details For \code{method = "pot"}, the function computes the reciprocal hazard and its derivative on the log scale to avoid numerical overflow. Thus, the density function should have argument \code{log} and the distribution function arguments \code{log.p} and \code{lower.tail}, respectively.
#' @param family the name of the parametric family. Will be used to obtain \code{dfamily}, \code{pfamily}, \code{qfamily}
#' @param method either block maxima (\code{'bm'}) or peaks-over-threshold (\code{'pot'}) are supported
#' @param thresh vector of thresholds for method \code{'pot'}
#' @param qlev vector of quantile levels for method \code{'pot'}, e.g., 0.9, 0.95, ... Ignored if argument \code{thresh} is provided.
#' @param m vector of block sizes for method \code{'bm'}
#' @param ... additional arguments passed to \code{densF} and \code{distF}
#' @author Leo Belzile
#' @importFrom methods formalArgs
#' @import stats
#' @return a data frame containing
#' \itemize{
#' \item \code{loc}: location parameters (\code{method='bm'})
#' \item \code{scale}: scale parameters
#' \item \code{shape}: shape parameters
#' \item \code{thresh}: thresholds (if \code{method='pot'}), percentile corresponding to threshold (if \code{method='pot'})
#' \item \code{m}: block sizes (if \code{method='bm'})
#' }
#' @references Smith, R.L. (1987). Approximations in extreme value theory. \emph{Technical report 205}, Center for Stochastic Process, University of North Carolina, 1--34.
#' @examples
#' # Threshold exceedance for Normal variables
#' quants <- seq(1, 5, by = 0.02)
#' penult <- penultimate(
#'    family = "norm",
#'    method = 'pot',
#'    thresh = quants,
#'    ddensF = function(x){-x*dnorm(x)}, # optional argument
#'    )
#' plot(x = quants,
#'      y = penult$shape,
#'      type = 'l',
#'      xlab = 'quantile',
#'     ylab = 'Penultimate shape',
#'     ylim = c(-0.5, 0))
#' # Block maxima for Gamma variables
#' # User must provide arguments for shape (or rate), for which there is no default
#' m <- seq(30, 3650, by = 30)
#' penult <- penultimate(family = 'gamma', method = 'bm', m = m, shape = 0.1)
#' plot(x = m,
#'      y = penult$shape,
#'      type = 'l',
#'      xlab = 'quantile',
#'      ylab = 'penultimate shape')
#'
#' # Comparing density of GEV approximation with true density of maxima
#' m <- 100 # block of size 100
#' p <- penultimate(
#'   family = 'norm',
#'   ddensF = function(x){-x*dnorm(x)},
#'   method = 'bm',
#'   m = m)
#' x <- seq(1, 5, by = 0.01)
#' plot(
#'   x = x,
#'   y = m * dnorm(x) * exp((m-1) * pnorm(x, log.p = TRUE)),
#'   type = 'l',
#'   ylab = 'density',
#'   main = 'Distribution of the maxima of\n 100 standard normal variates')
#' lines(x, mev::dgev(x, loc = p$loc, scale = p$scale, shape = 0), col = 2)
#' lines(x, mev::dgev(x, loc = p$loc, scale = p$scale, shape = p$shape), col = 4)
#' legend(
#'  x = 'topright',
#'  lty = c(1, 1, 1),
#'  col = c(1, 2, 4),
#'  legend = c('exact', 'ultimate', 'penultimate'),
#'  bty = 'n')
#' @export
penultimate <- function(
  family,
  method = c("bm", "pot"),
  thresh,
  qlev,
  m,
  ...
) {
  #  Sort argument by default
  if (!missing(m)) {
    m <- sort(m)
  } else {
    m <- NULL
  }
  if (!missing(thresh)) {
    thresh <- sort(thresh)
  } else {
    thresh <- NULL
  }
  if (!missing(qlev)) {
    qlev <- sort(qlev)
  } else {
    qlev <- NULL
  }
  args <- list(...)
  args$returnList <- TRUE
  args$qu <- qlev
  args$u <- thresh
  args$m <- m
  res <- do.call(
    what = smith.penult,
    args = c(
      family = family,
      method = method,
      args
    )
  )
  res$thresh <- res$u
  res$u <- NULL
  res <- as.data.frame(res)
  class(res) <- c("data.frame", "mev_penultimate")
  return(res)
}

#' Smith's third penultimate approximation
#'
#' This function returns the density and distribution functions
#' of the 3rd penultimate approximation for extremes of Smith (1987). It requires
#' knowledge of the exact constants \eqn{\epsilon} and \eqn{\rho} described in the paper.
#'
#' Let \eqn{F}, \eqn{f} denote respectively the distribution and density functions and define the function \eqn{\phi(x)}  as
#' \deqn{\phi(x)=-\frac{F(x)\log F(x)}{f(x)}}{\phi(x)=-F(x)log F(x)/f(x)}
#' for block maxima.
#' The sequence \code{loc} corresponds to \eqn{b_n} otherwise, defined as the solution of \eqn{F(b_n)=\exp(-1/n)}{F(b_n)=exp(-1/n)}.
#'
#' The \code{scale} is given by \eqn{a_n=\phi(b_n)}, the \code{shape} as \eqn{\gamma_n=\phi'(b_n)}. These are returned by a call to \link{smith.penult}.
#'
#' For threshold exceedances, \eqn{b_n} is replaced by the sequence of thresholds \eqn{u} and we
#' take instead \eqn{\phi(x)} to be the reciprocal hazard function \eqn{\phi(x)=(1-F(x))/f(x)}{\phi(x)=(1-F(x))/f(x)}.
#'
#' In cases where the distribution function is in the maximum domain of
#' attraction of the Gumbel distribution, \eqn{\rho} is possibly undetermined and
#' \eqn{\epsilon} can be equal to \eqn{\phi(b_n)\phi''(b_n)}.
#'
#' For distributions in the maximum domain of
#' attraction of the Gumbel distribution and that are class N, it is also possible to abstract from the \eqn{\rho} parameter by substituting the function \eqn{H_{\rho}}{H(x;\rho)} by \eqn{x^3/6} without affecting the rate of convergence. This can be done by setting \code{mdaGumbel=TRUE} in the function call.
#'
#' @section Warning:
#' The third penultimate approximation does not yield a valid distribution function over the whole range of the original distribution, but is rather valid in a neighborhood of the true support of the distribution of maxima/threshold exceedance.
#' The function handles the most standard failure (decreasing distribution function and negative densities), but any oscillatory behaviour will not necessarily be captured.
#' This is inherent to the method and can be resolved by `not' evaluating the functions \eqn{F} and \eqn{f} at the faulty points.
#' @param loc location parameter returned by \code{\link{smith.penult}} or threshold vector
#' @param scale scale parameter returned by \code{\link{smith.penult}}
#' @param shape shape parameter returned by \code{\link{smith.penult}}
#' @param eps parameter vector, see \strong{Details}.
#' @param rho second-order parameter, model dependent
#' @param method one of \code{pot} for the generalized Pareto or \code{bm} for the generalized extreme value distribution
#' @param mdaGumbel logical indicating whether the function \eqn{H_{\rho}}{H(x;\rho)} should be replaced by \eqn{x^3/6}; see \strong{Details}.
#' @param ... additional parameters, currently ignored. These are used for backward compatibility due to a change in the names of the arguments.
#' @references Smith, R.L. (1987). Approximations in extreme value theory. \emph{Technical report 205}, Center for Stochastic Process, University of North Carolina, 1--34.
#' @keywords internal
#' @examples
#' #Normal maxima example from Smith (1987)
#' m <- 100 #block of size 100
#' p <- smith.penult(family='norm',
#'    ddensF=function(x){-x*dnorm(x)}, method='bm', m=m, returnList=FALSE)
#' approx <- smith.penult.fn(loc=p[1], scale=p[2], shape=p[3],
#'    eps=p[3]^2+p[3]+p[2]^2, mdaGumbel=TRUE, method='bm')
#' x <- seq(0.5,6,by=0.001)
#' #First penultimate approximation
#' plot(x, exp(m*pnorm(x, log.p=TRUE)),type='l', ylab='CDF',
#' main='Distribution of the maxima of\n 100 standard normal variates')
#' lines(x, mev::pgev(x,loc=p[1], scale=p[2], shape=0),col=2)
#' lines(x, mev::pgev(x,loc=p[1], scale=p[2], shape=p[3]),col=3)
#' lines(x, approx$F(x),col=4)
#' legend(x='bottomright',lty=c(1,1,1,1),col=c(1,2,3,4),
#'    legend=c('Exact','1st approx.','2nd approx.','3rd approx'),bty='n')
#' plot(x, m*dnorm(x)*exp((m-1)*pnorm(x,log.p=TRUE)),type='l', ylab='Density',
#' main='Distribution of the maxima of\n 100 standard normal variates')
#' lines(x, mev::dgev(x,loc=p[1], scale=p[2], shape=0),col=2)
#' lines(x, mev::dgev(x,loc=p[1], scale=p[2], shape=p[3]),col=3)
#' lines(x, approx$f(x),col=4)
#' legend(x='topright',lty=c(1,1,1,1),col=c(1,2,3,4),
#'  legend=c('Exact','1st approx.','2nd approx.','3rd approx'),bty='n')
#'
#' #Threshold exceedances
#' par <- smith.penult(family = "norm", ddensF=function(x){-x*dnorm(x)},
#' method='pot', u=4, returnList=FALSE)
#' approx <- smith.penult.fn(loc=par[1], scale=par[2], shape=par[3],
#'  eps=par[3]^2+par[3]+par[2]^2, mdaGumbel=TRUE, method='pot')
#' x <- seq(4.01,7,by=0.01)
#' #Distribution function
#' plot(x, 1-(1-pnorm(x))/(1-pnorm(par[1])),type='l', ylab='Conditional CDF',
#' main='Exceedances over 4\n for standard normal variates')
#' lines(x, mev::pgp(x, loc=par[1], scale=par[2], shape=0),col=2)
#' lines(x, mev::pgp(x, loc=par[1], scale=par[2], shape=par[3]),col=3)
#' lines(x, approx$F(x),col=4)
#' #Density
#' plot(x, dnorm(x)/(1-pnorm(par[1])),type='l', ylab='Conditional density',
#' main='Exceedances over 4\n for standard normal variates')
#' lines(x, dgp(x, loc=par[1], scale=par[2], shape=0),col=2)
#' lines(x, dgp(x, loc=par[1], scale=par[2], shape=par[3]),col=3)
#' lines(x, approx$f(x),col=4)
#' @export
smith.penult.fn <- function(
  loc,
  scale,
  shape,
  eps,
  rho = NULL,
  method = c("bm", "pot"),
  mdaGumbel = FALSE,
  ...
) {
  ellips <- list(...)
  #Compatibility condition with version 1.12 and before, via ellipsis argument
  if (!is.null(ellips$model) && length(method) == 2) {
    method <- ellips$model
  }
  bn = loc
  an = scale
  gamma = shape
  # Functions appearing in the theory of slow variation
  hrho <- function(x, rho) {
    if (rho == 0) {
      log(x)
    } else {
      (x^rho - 1) / rho
    }
  }
  H <- function(x, eta, rho) {
    H0 <- function(x, eta)
      (0.5 *
        (log(pmax(1 + x * eta, 0))^2) -
        log(pmax(1 + x * eta, 0)) +
        1 -
        1 / (1 + x * eta)) /
        eta^3
    Hm1 <- function(x, eta)
      (log(pmax(1 + x * eta, 0)) /
        (1 + x * eta) +
        log(pmax(1 + x * eta, 0)) -
        2 * (1 - 1 / (1 + x * eta))) /
        eta^3
    Hrho <- function(x, eta, rho)
      ifelse(
        (1 + x * eta) > 0,
        hrho(1 + x * eta, rho) +
          rho * hrho(1 + x * eta, -1) -
          (rho + 1) *
            log(1 + x * eta),
        0
      )
    if (rho == 0) {
      H0(x, eta)
    } else if (rho == -1) {
      Hm1(x, eta)
    } else {
      Hrho(x, eta, rho)
    }
  }
  # derivative of dH
  dH <- function(x, eta, rho) {
    ifelse(
      (1 + x * eta) > 0,
      ((1 + x * eta)^(rho - 1) +
        rho * (1 + x * eta)^(-2) -
        (rho + 1) / (1 + x * eta)) /
        (rho * (rho + 1) * eta^2),
      0
    )
  }

  ## Block maxima - GEV-like distribution functions and densities Distribution function of third penultimate approximation
  if (method == "bm") {
    if (!mdaGumbel) {
      if (is.null(rho)) {
        stop("Invalid \"rho\" parameter")
      }
      GEV3rd <- function(x) {
        # , an, bn, gamma, eps, rho
        q <- (x - bn) / an
        if (gamma == 0) p <- exp(-exp(-q) * (1 + eps * H(q, gamma, rho))) else
          p <- exp(
            -pmax(1 + gamma * q, 0)^(-1 / gamma) * (1 + eps * H(q, gamma, rho))
          )
        invalid <- which(diff(p) < 0)
        p[invalid] <- 0
        p <- pmin(1, pmax(0, p))
        return(p)
      }
      # Density of third penultimate approximation , an, bn, gamma, eps, rho
      dGEV3rd <- function(x) {
        q <- (x - bn) / an
        if (gamma == 0) stop("Not yet implemented")
        pmax(
          0,
          GEV3rd(x) /
            an *
            ((1 + eps * H(q, gamma, rho)) *
              exp((-1 / gamma - 1) * log(pmax(1 + q * gamma, 0))) -
              exp(
                (-1 / gamma) *
                  log(pmax(1 + q * gamma, 0))
              ) *
                eps *
                dH(q, gamma, rho))
        )
      }
      return(list(F = GEV3rd, f = dGEV3rd))
    } else {
      # mdaGumbel Distribution function of the third penultimate approximation, replacing H_rho(x, eta) by x^3/6 , an, bn, gamma, eps
      GEV3rda <- function(x) {
        q <- (x - bn) / an
        if (gamma == 0) p <- exp(-exp(-q) * (1 + eps * q^3 / 6)) else
          p <- exp(-pmax(1 + gamma * q, 0)^(-1 / gamma) * (1 + eps * q^3 / 6))
        invalid <- which(diff(p) < 0)
        p[invalid] <- 0
        p <- pmin(1, pmax(0, p))
        return(p)
      }
      # Density of third penultimate approximation, replacing H_rho(x, eta) by x^3/6 , an, bn, gamma, eps
      dGEV3rda <- function(x) {
        q <- (x - bn) / an
        if (gamma == 0) stop("Not yet implemented")
        pmax(
          0,
          GEV3rda(x) /
            an *
            ((1 + eps * q^3 / 6) *
              exp((-1 / gamma - 1) * log(pmax(1 + q * gamma, 0))) -
              exp(
                (-1 / gamma) *
                  log(pmax(1 + q * gamma, 0))
              ) *
                eps *
                q^2 /
                2)
        )
      }
      return(list(F = GEV3rda, f = dGEV3rda))
    }
  } else if (method == "pot") {
    if (!mdaGumbel) {
      if (is.null(rho)) {
        stop("Invalid \"rho\" parameter")
      }
      # Peaks-over-threshold - GP-like distribution functions and densities , an, bn, gamma, eps, rho
      GP3rd <- function(x) {
        q <- (x - bn) / an
        p <- 1 -
          pmax(0, (1 + gamma * q))^(-1 / gamma) * (1 + eps * H(q, gamma, rho))
        p[which(diff(p) < 0)] <- 0
        p <- pmin(1, pmax(0, p))
        p
      }
      dGP3rd <- function(x) {
        # , an, bn, gamma, eps, rho
        q <- (x - bn) / an
        pmax(
          0,
          pmax(0, (1 + gamma * q))^(-1 / gamma - 1) /
            an *
            ((1 + eps * H(q, gamma, rho)) -
              pmax(0, (1 + gamma * q)) * eps * dH(q, gamma, rho))
        )
      }
      return(list(F = GP3rd, f = dGP3rd))
    } else {
      GP3rda <- function(x) {
        # , an, bn, gamma, eps
        q <- (x - bn) / an
        p <- 1 - pmax(0, (1 + gamma * q))^(-1 / gamma) * (1 + eps * q^3 / 6)
        p[which(diff(p) < 0)] <- 0
        p <- pmin(1, pmax(0, p))
        p
      }

      dGP3rda <- function(x) {
        # , an, bn, gamma, eps
        q <- (x - bn) / an
        pmax(
          0,
          pmax(0, (1 + gamma * q))^(-1 / gamma - 1) /
            an *
            ((1 + eps * q^3 / 6) - pmax(0, (1 + gamma * q)) * eps * q^2 / 2)
        )
      }
      return(list(F = GP3rda, f = dGP3rda))
    }
  }
}


#' Transform arguments using max stability
#'
#' Given a vector of location, scale and shape parameters,
#' compute the corresponding parameters for block of size
#' \code{m} assuming a generalized extreme value distribution.
#' @param pars vector of location, scale and shape parameters
#' @param m [integer] block size
#' @param inverse [logical] whether to compute the parameters for the inverse relationship (defaults to \code{FALSE})
#' @export
#' @examples
#' maxstable(pars = maxstable(pars = c(1,2,0), m = 10), m = 10, inv = TRUE)
#' maxstable(pars = maxstable(pars = c(1,2,0.1), m = 5), m = 1/5)
maxstable <- function(pars, m = 1L, inverse = FALSE) {
  stopifnot(m > 0, length(pars) == 3L, is.logical(inverse))
  pars <- as.numeric(pars)
  if (isTRUE(all.equal(as.numeric(pars[3]), 0))) {
    isZero <- TRUE
  } else {
    isZero <- FALSE
  }
  if (isTRUE(inverse)) {
    m <- 1 / m
  }
  if (!isZero) {
    return(
      c(
        loc = pars[1] - pars[2] / pars[3] * (1 - m^pars[3]),
        scale = pars[2] * m^pars[3],
        shape = pars[3]
      )
    )
  } else {
    return(
      c(loc = pars[1] + pars[2] * log(m), scale = pars[2], shape = 0)
    )
  }
}
