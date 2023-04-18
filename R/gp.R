
#------------------------------------------------------------------------------#
# Fit GP (sigma, xi) distribution using 1D optimisation #
#------------------------------------------------------------------------------#
.gpd_1D_fit <- function(xdat, threshold, npy = 365, ydat = NULL, sigl = NULL, shl = NULL, siglink = identity, shlink = identity, siginit = NULL,
                        shinit = NULL, show = TRUE, method = "BFGS", maxit = 10000, xi.tol = 0.001, phi.input = NULL, calc.se = F, ...) {
  z <- list()
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  n <- length(xdat)
  z$trans <- FALSE
  if (is.function(threshold))
    stop("\"threshold\" cannot be a function")
  u <- rep(threshold, length.out = n)
  if (length(unique(u)) > 1)
    z$trans <- TRUE
  xdatu <- xdat[xdat > u]
  xind <- (1:n)[xdat > u]
  u <- u[xind]
  ex.data <- xdatu - u
  in2 <- sqrt(6 * var(xdat))/pi
  in1 <- mean(xdat, na.rm = TRUE) - 0.57722 * in2
  phi.init <- (2^0.1 - 1)/median(xdatu - u)

  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, length(xdatu)))
    if (is.null(siginit)) {
      siginit <- in2
    }
  } else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, length(xdatu)), ydat[xind, sigl])
    if (is.null(siginit))
      siginit <- c(in2, rep(0, length(sigl)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdatu)))
    if (is.null(shinit))
      shinit <- 0.1
  } else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdatu)), ydat[xind, shl])
    if (is.null(shinit))
      shinit <- c(0.1, rep(0, length(shl)))
  }
  init <- c(siginit, shinit)
  z$model <- list(sigl, shl)
  z$link <- deparse(substitute(c(siglink, shlink)))
  z$threshold <- threshold
  z$nexc <- length(xdatu)
  z$data <- xdatu

  GP.1D.negloglik <- function(phi) {
    # negated 1D loglikelihood
    k <- length(ex.data)
    if (phi == 0) {
      sigma.mle <- (1/k) * sum(ex.data)
      return(k * log(sigma.mle) + (1/sigma.mle) * sum(ex.data))
    }
    zz <- 1 + phi * ex.data
    if (min(zz) <= 0)
      return(1e+30)
    xi.of.phi <- mean(log(zz))
    sigma.of.phi <- xi.of.phi/phi
    if (sigma.of.phi <= 0)
      return(1e+30)
    k * log(sigma.of.phi) + (1 + 1/xi.of.phi) * sum(log(1 + phi * ex.data))
  }  #................................# end of GP.1D.negloglik()
  #
  gp.1D.grad <- function(a) {
    # gradient of negated log-likelihood
    phi <- a
    n <- length(ex.data)
    yy <- 1 + phi * ex.data
    xi.phi <- mean(log(yy))
    -n/phi + (1 + 1/xi.phi) * sum(ex.data/yy)
  }
  #
  if (!is.null(phi.input)) {
    phi.init <- phi.input
  }
  temp <- try(optim(phi.init, GP.1D.negloglik, gr = gp.1D.grad, hessian = FALSE, method = method, control = list(maxit = maxit,
                                                                                                                 ...)))
  if (inherits(temp, what = "try-error")) {
    z$conv <- 50
    return(z)
  }
  phi <- temp$par
  zz <- 1 + phi * (xdatu - u)
  if (min(zz) <= 0) {
    z$conv <- 50
    return(z)
  }
  xi <- mean(log(zz))
  sc <- xi/phi
  z$mle <- c(sc, xi)
  if (calc.se && xi > -0.5) {
    z$cov <- solve(gpd.infomat(par = c(sc, xi), dat = xdatu - u))
  }
  z$conv <- temp$convergence
  z$counts <- temp$counts
  z$nllh <- temp$value
  z$vals <- cbind(sc, xi, u)
  z$rate <- length(xdatu)/n
  if (calc.se) {
    z$se <- tryCatch(sqrt(diag(z$cov)), error = function(e) NULL)
  }
  z$n <- n
  z$npy <- npy
  z$xdata <- xdat
  if (show)
    print(z[c(4, 5, 9, 10, 7, 12, 13)])
  uu <- unique(u)
  z$mle.t <- c(z$mle[1] - z$mle[2] * uu, z$mle[2])  # sigma-xi*u replaces sigma

  if (calc.se) {
    d <- matrix(c(1, -uu, 0, 1), 2, 2, byrow = T)  # derivatives of sigma-xi*u and xi
    v <- d %*% z$cov %*% t(d)  # new VC matrix
    z$cov.t <- v
    z$se.t <- sqrt(diag(z$cov.t))  # new SEs
  }

  class(z) <- "gpd.fit"
  invisible(z)
}

#------------------------------------------------------------------------------#
# Fit GP(sigma, xi) distribution using 2D optimisation #
#------------------------------------------------------------------------------#

# Stolen from ismev; gradient function added.
#' Maximum likelihood method for the generalized Pareto Model
#'
#' Maximum-likelihood estimation for the generalized Pareto model, including generalized linear modelling of each parameter. This function was adapted by Paul Northrop to include the gradient in the \code{gpd.fit} routine from \code{ismev}.
#'
#' @param xdat numeric vector of data to be fitted.
#' @param threshold a scalar or a numeric   vector of the same length as \code{xdat}.
#' @param npy number of observations per year/block.
#' @param ydat  matrix of covariates for generalized linear modelling of the parameters (or \code{NULL} (the default) for stationary fitting). The number of rows should be the same as the length of \code{xdat}.
#' @param sigl numeric vector of integers, giving the columns of \code{ydat} that contain covariates for generalized linear modelling of the scale parameter (or \code{NULL} (the default) if the corresponding parameter is stationary).
#' @param shl numeric vector of integers, giving the columns of \code{ydat} that contain covariates for generalized linear modelling of the shape parameter (or \code{NULL} (the default) if the corresponding parameter is stationary).
#' @param siglink inverse link functions for generalized linear modelling of the scale parameter
#' @param shlink inverse link functions for generalized linear modelling of the shape parameter
#' @param siginit numeric giving initial value(s) for parameter estimates. If \code{NULL} the default is \code{sqrt(6 * var(xdat))/pi}
#' @param shinit numeric giving initial value(s) for the shape parameter estimate; if \code{NULL}, this is 0.1.  If using parameter covariates, then these values are used for the constant term, and zeros for all other terms.
#' @param show logical; if \code{TRUE} (default), print details of the fit.
#' @param method optimization method (see \code{\link{optim}} for details).
#' @param maxit  maximum number of iterations.
#' @param ...other control parameters for the optimization. These are passed to components of the \code{control} argument of \code{optim}.
#'
#' @details  For non-stationary fitting it is recommended that the covariates within the generalized linear models are (at least approximately) centered and scaled (i.e. the columns of \code{ydat} should be approximately centered and scaled).
#'
#' The form of the GP model used follows Coles (2001) Eq (4.7).  In particular, the shape parameter is defined so that positive values imply a heavy tail and negative values imply a bounded upper value.
#' @return a list with components
#' \describe{
#' \item{nexc}{scalar giving the number of threshold exceedances.}
#' \item{nllh}{scalar giving the negative log-likelihood value.}
#' \item{mle}{numeric vector giving the MLE's for the scale and shape parameters, resp.}
#'  \item{rate}{scalar giving the estimated probability of exceeding the threshold.}
#'  \item{se}{numeric vector giving the standard error estimates for the scale and shape parameter estimates, resp.}
#' \item{trans}{logical indicator for a non-stationary fit.}
#' \item{model}{list with components \code{sigl} and \code{shl}.}
#' \item{link}{character vector giving inverse link functions.}
#' \item{threshold}{threshold, or vector of thresholds.}
#' \item{nexc}{number of data points above the threshold.}
#' \item{data}{data that lie above the threshold. For non-stationary models, the data are standardized.}
#' \item{conv}{convergence code, taken from the list returned by \code{\link{optim}}. A zero indicates successful convergence.}
#' \item{nllh}{negative log likelihood evaluated at the maximum likelihood estimates.}
#' \item{vals}{matrix with three columns containing the maximum likelihood estimates of the scale and shape parameters, and the threshold, at each data point.}
#' \item{mle}{vector containing the maximum likelihood estimates.}
#' \item{rate}{proportion of data points that lie above the threshold.}
#' \item{cov}{covariance matrix.}
#' \item{se}{numeric vector containing the standard errors.}
#' \item{n}{number of data points (i.e., the length of \code{xdat}).}
#' \item{npy}{number of observations per year/block.}
#' \item{xdata}{data that has been fitted.}
#' }
#' @references Coles, S., 2001.  An Introduction to Statistical Modeling of Extreme Values.  Springer-Verlag, London.
#' @export
#' @keywords internal
.gpd_2D_fit <- function(xdat, threshold, npy = 365, ydat = NULL, sigl = NULL, shl = NULL, siglink = identity, shlink = identity, siginit = NULL,
                        shinit = NULL, show = TRUE, method = "Nelder-Mead", maxit = 10000, ...) {
  z <- list()
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  n <- length(xdat)
  z$trans <- FALSE
  if (is.function(threshold))
    stop("\"threshold\" cannot be a function")
  u <- rep(threshold, length.out = n)
  if (length(unique(u)) > 1)
    z$trans <- TRUE
  xdatu <- xdat[xdat > u]
  xind <- (1:n)[xdat > u]
  u <- u[xind]
  in2 <- sqrt(6 * var(xdat))/pi
  in1 <- mean(xdat, na.rm = TRUE) - 0.57722 * in2
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, length(xdatu)))
    if (is.null(siginit))
      siginit <- in2
  } else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, length(xdatu)), ydat[xind, sigl])
    if (is.null(siginit))
      siginit <- c(in2, rep(0, length(sigl)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdatu)))
    if (is.null(shinit))
      shinit <- 0.1
  } else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdatu)), ydat[xind, shl])
    if (is.null(shinit))
      shinit <- c(0.1, rep(0, length(shl)))
  }
  init <- c(siginit, shinit)
  z$model <- list(sigl, shl)
  z$link <- deparse(substitute(c(siglink, shlink)))
  z$threshold <- threshold
  z$nexc <- length(xdatu)
  z$data <- xdatu
  gpd.lik <- function(a) {
    sc <- siglink(sigmat %*% (a[seq(1, length.out = npsc)]))
    xi <- shlink(shmat %*% (a[seq(npsc + 1, length.out = npsh)]))
    y <- (xdatu - u)/sc
    y <- 1 + xi * y
    if (min(sc) <= 0) {
      l <- 10e5
    } else {
      if (min(y) <= 0) {
        l <- 10e5
      } else {
        l <- sum(log(sc)) + sum(log(y) * (1/xi + 1))
      }
    }
    l
  }
  #
  gp.grad <- function(a) {
    # gradient of negated log-likelihood
    sigma <- a[1]
    xi <- a[2]
    y <- xdatu - u
    n <- length(y)
    yy <- 1 + xi * y/sigma
    s1 <- -n * sigma^(-1) + (1 + xi) * sigma^(-2) * sum(y/yy)
    if (any(yy < 0)) {
      # TODO FIX WHETHER THIS MAKES SENSE
      -c(s1, 1e+30)
    } else {
      s2 <- xi^(-2) * sum(log(yy)) - (1 + 1/xi) * sigma^(-1) * sum(y/yy)
      -c(s1, s2)
    }
  }
  #
  x <- try(optim(init, gpd.lik, gr = gp.grad, hessian = TRUE, method = method, control = list(maxit = maxit, ...)))
  if (inherits(x, what = "try-error")) {
    z$conv <- 50
    return(z)
  }
  sc <- siglink(sigmat %*% (x$par[seq(1, length.out = npsc)]))
  xi <- shlink(shmat %*% (x$par[seq(npsc + 1, length.out = npsh)]))
  z$conv <- x$convergence
  z$counts <- x$counts
  z$nllh <- x$value
  z$vals <- cbind(sc, xi, u)
  if (z$trans) {
    z$data <- -log(as.vector((1 + (xi * (xdatu - u))/sc)^(-1/xi)))
  }
  z$mle <- x$par
  z$rate <- length(xdatu)/n
  z$cov <- tryCatch(solve(x$hessian), error = function(e) {
    matrix(NA, ncol(x$hessian), ncol(x$hessian))
  })  #TODO fix this
  suppressWarnings(z$se <- sqrt(diag(z$cov)))
  z$n <- n
  z$npy <- npy
  z$xdata <- xdat
  if (show) {
    if (z$trans)
      print(z[c(2, 3)])
    if (length(z[[4]]) == 1)
      print(z[4])
    print(z[c(5, 7)])
    if (!z$conv)
      print(z[c(8, 10, 11, 13)])
  }
  class(z) <- "gpd.fit"
  invisible(z)
}

#------------------------------------------------------------------------------#
# Fit GP(sigma, xi) distribution using 1D optimisation, using nlm # ... to force gradients to be very close to zero at MLE #
#------------------------------------------------------------------------------#

.GP_1D_fit_nlm <- function(xdat, threshold, init.val = NULL, calc.se = F, ...) {
  z.fit <- list()
  ex.data <- xdat[xdat > threshold] - threshold
  #
  GP.1D.negloglik <- function(phi) {
    # negated 1D loglikelihood
    k <- length(ex.data)
    if (phi == 0) {
      sigma.mle <- mean(ex.data)
      return(k * log(sigma.mle) + (1/sigma.mle) * sum(ex.data))
    }
    zz <- 1 + phi * ex.data
    if (min(zz) <= 0)
      return(1e+30)
    xi.of.phi <- mean(log(zz))
    sigma.of.phi <- xi.of.phi/phi
    if (sigma.of.phi <= 0)
      return(1e+30)
    neg.log.lik <- k * log(sigma.of.phi) + (1 + 1/xi.of.phi) * sum(log(1 + phi * ex.data))
    #
    attr(neg.log.lik, "gradient") <- -k/phi + (1 + 1/xi.of.phi) * sum(ex.data/zz)
    #
    attr(neg.log.lik, "hessian") <- k/phi^2 - xi.of.phi^(-2) * (sum(ex.data/zz))^2/k - (1 + 1/xi.of.phi) * sum(ex.data^2/zz^2)
    #
    neg.log.lik
  }  #................................# end of GP.1D.negloglik()
  #
  init <- ifelse(is.null(init.val), (2^0.1 - 1)/median(ex.data), init.val)
  f.scale <- GP.1D.negloglik(init)
  typ.size <- init
  suppressWarnings(res <- nlm(GP.1D.negloglik, init, fscale = f.scale, typsize = typ.size, iterlim = 1000, check.analyticals = FALSE,
                              ...))
  phi.hat <- res$estimate
  xi.hat <- mean(log(1 + res$estimate * ex.data))
  sigma.hat <- xi.hat/phi.hat
  z.fit$mle <- c(sigma.hat, xi.hat)
  z.fit$nllh <- res$minimum
  z.fit$gradient <- res$gradient
  z.fit$code <- res$code
  z.fit$counts["function"] <- res$iterations
  z.fit$convergence <- ifelse(res$code %in% 1:2, 0, 1)
  if (calc.se)
    z.fit$se <- sqrt(diag(solve(gpd.infomat(c(sigma.hat, xi.hat), dat = ex.data, method = "obs"))))
  invisible(z.fit)
}


#------------------------------------------------------------------------------#
# GP fitting function of Grimshaw (1993) #
#------------------------------------------------------------------------------#

#' GP fitting function of Grimshaw (1993)
#'
#' Function for estimating parameters \code{k} and \code{a} for a random sample from a GPD.
#'
#' @author Scott D. Grimshaw
#' @param x sample values
#' @return a list with the maximum likelihood estimates of components \code{a} and \code{k}
#' @keywords internal
.fit.gpd.grimshaw <- function (x)
  {
    n <- length(x)
    xq <- sort(x)
    xbar <- mean(x)
    sumx2 <- sum(x^2)/n
    x1 <- xq[1]
    xn <- xq[n]
    epsilon <- 1e-6/xbar # careful: too low epsilon leads to xi=0
    lobnd <- 2 * (x1 - xbar)/x1^2
    maxiter <- 400L
    if (lobnd >= 0) {
      lobnd <- -epsilon
    }
    hibnd <- (1/xn) - epsilon #modif
    if (hibnd <= 0) {
      hibnd <- epsilon
    }
    secderiv <- sumx2 - 2 * xbar^2
    if (secderiv > 0) {
      thzeros <- cbind(c(0, 0), c(0, 0))
      nzeros <- 2
      hlo <- (1 + sum(log(1 - lobnd * x))/n) * (sum(1/(1 - lobnd * x))/n) - 1
      if (hlo < 0) {
        thlo <- lobnd
        thhi <- -epsilon
      }   else {
        thlo <- -epsilon
        thhi <- lobnd
      }
      thzero <- lobnd
      dxold <- abs(thhi - thlo)
      dx <- dxold
      temp1 <- mean(log(1 - thzero * x))
      temp2 <- mean(1/(1 - thzero * x))
      temp3 <- mean(1/(1 - thzero * x)^2)
      h <- (1 + temp1) * (temp2) - 1
      hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3))/thzero
      j <- 1
      while (j <= maxiter) {
        c1 <- (((thzero - thhi) * hprime - h) * ((thzero - thlo) * hprime - h) >= 0)
        c2 <- (abs(2 * h) > abs(dxold * hprime))
        if (c1 + c2 >= 1) {
          dxold <- dx
          dx <- (thhi - thlo)/2
          thzero <- thlo + dx
          if (thlo == thzero) {
            j <- 1000
          }
        } else {
          dxold <- dx
          dx <- h/hprime
          temp <- thzero
          thzero <- thzero - dx
          if (temp == thzero) {
            j <- 1001
          }
        }
        if (abs(dx) < epsilon * abs(thlo + thhi)/2) {
          j <- 999
        }
        temp1 <- mean(log(1 - thzero * x))
        temp2 <- mean(1/(1 - thzero * x))
        temp3 <- mean(1/(1 - thzero * x)^2)
        h <- (1 + temp1) * (temp2) - 1
        hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3))/thzero
        if (h < 0) {
          thlo <- thzero
        }  else {
          thhi <- thzero
        }
        j <- j + 1
      }
      if (j > maxiter + 1) {
        thzeros[1, ] <- cbind(thzero, j)
      }
      hlo <- (1 + sum(log(1 - epsilon * x))/n) * (sum(1/(1 -epsilon * x))/n) - 1
      if (hlo < 0) {
        thlo <- epsilon
        thhi <- hibnd
      }  else {
        thlo <- hibnd
        thhi <- epsilon
      }
      thzero <- (hibnd + epsilon)/2
      dxold <- abs(thhi - thlo)
      dx <- dxold
      temp1 <- mean(log(1 - thzero * x))
      temp2 <- mean(1/(1 - thzero * x))
      temp3 <- mean(1/(1 - thzero * x)^2)
      h <- (1 + temp1) * (temp2) - 1
      hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3))/thzero
      j <- 1
      while (j <= maxiter) {
        c1 <- (((thzero - thhi) * hprime - h) * ((thzero - thlo) * hprime - h) >= 0)
        c2 <- (abs(2 * h) > abs(dxold * hprime))
        if (c1 + c2 >= 1) {
          dxold <- dx
          dx <- (thhi - thlo)/2
          thzero <- thlo + dx
          if (thlo == thzero) {
            j <- 1000
          }
        } else {
          dxold <- dx
          dx <- h/hprime
          temp <- thzero
          thzero <- thzero - dx
          if (temp == thzero) {
            j <- 1001
          }
        }
        if (abs(dx) < epsilon * abs(thlo + thhi)/2) {
          j <- 999
        }
        temp1 <- mean(log(1 - thzero * x))
        temp2 <- mean(1/(1 - thzero * x))
        temp3 <- mean(1/(1 - thzero * x)^2)
        h <- (1 + temp1) * (temp2) - 1
        hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3))/thzero
        if (h < 0) {
          thlo <- thzero
        } else {
          thhi <- thzero
        }
        j <- j + 1
      }
      if (j > maxiter + 1) {
        thzeros[2, ] <- cbind(thzero, j)
      }
    } else {
      thzeros <- matrix(rep(0, 8), ncol = 2)
      nzeros <- 4
      hlo <- (1 + sum(log(1 - lobnd * x))/n) * (sum(1/(1 -lobnd * x))/n) - 1
      if (hlo < 0) {
        thlo <- lobnd
        thhi <- -epsilon
      } else {
        thlo <- -epsilon
        thhi <- lobnd
      }
      thzero <- lobnd
      dxold <- abs(thhi - thlo)
      dx <- dxold
      temp1 <- mean(log(1 - thzero * x))
      temp2 <- mean(1/(1 - thzero * x))
      temp3 <- mean(1/(1 - thzero * x)^2)
      h <- (1 + temp1) * (temp2) - 1
      hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3))/thzero
      j <- 1
      while (j <= maxiter) {
        c1 <- (((thzero - thhi) * hprime - h) * ((thzero - thlo) * hprime - h) >= 0)
        c2 <- (abs(2 * h) > abs(dxold * hprime))
        if (c1 + c2 >= 1) {
          dxold <- dx
          dx <- (thhi - thlo)/2
          thzero <- thlo + dx
          if (thlo == thzero) {
            j <- 1000
          }
        }  else {
          dxold <- dx
          dx <- h/hprime
          temp <- thzero
          thzero <- thzero - dx
          if (temp == thzero) {
            j <- 1001
          }
        }
        if (abs(dx) < epsilon * abs(thlo + thhi)/2) {
          j <- 999
        }
        temp1 <- mean(log(1 - thzero * x))
        temp2 <- mean(1/(1 - thzero * x))
        temp3 <- mean(1/(1 - thzero * x)^2)
        h <- (1 + temp1) * (temp2) - 1
        hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3))/thzero
        if (h < 0) {
          thlo <- thzero
        } else {
          thhi <- thzero
        }
        j <- j + 1
      }
      if (j > maxiter + 1) {
        thzeros[1, ] <- cbind(thzero, j)
      }
      temp1 <- mean(log(1 - thzero * x))
      temp2 <- mean(1/(1 - thzero * x))
      temp3 <- mean(1/(1 - thzero * x)^2)
      hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3))/thzero
      if (hprime > 0) {
        thlo <- thzero
        thhi <- -epsilon
        thzero <- thhi
        dx <- thlo - thhi
        j <- 1

        while (j <= maxiter) {
          dx <- 0.5 * dx
          thmid <- thzero + dx
          hmid <- (1 + sum(log(1 - thmid * x))/n) * (sum(1/(1 - thmid * x))/n) - 1
          if (hmid < 0) {
            thzero <- thmid
          } else if (hmid == 0) {
            j <- 999
          }
          if (abs(dx) < epsilon * abs(thlo + thhi)/2) {
            j <- 999
          }
          j <- j + 1
        }
        if (j > maxiter + 1) {
          thzeros[2, ] <- cbind(thzero, j)
        }
      } else {
        thlo <- lobnd
        thhi <- thzero
        thzero <- thlo
        dx <- thhi - thlo
        j <- 1
        while (j <= maxiter) {
          dx <- 0.5 * dx
          thmid <- thzero + dx
          hmid <- (1 + sum(log(1 - thmid * x))/n) * (sum(1/(1 - thmid * x))/n) - 1
          if (hmid < 0) {
            thzero <- thmid
          }  else if (hmid == 0) {
            j <- 999
          }
          if (abs(dx) < epsilon * abs(thlo + thhi)/2) {
            j <- 999
          }
          j <- j + 1
        }
        if (j > maxiter + 1) {
          thzeros[2, ] <- cbind(thzero, j)
        }
      }
      hlo <- (1 + sum(log(1 - epsilon * x))/n) * (sum(1/(1 - epsilon * x))/n) - 1
      if (hlo < 0) {
        thlo <- epsilon
        thhi <- hibnd
      } else {
        thlo <- hibnd
        thhi <- epsilon
      }
      thzero <- (hibnd + epsilon)/2
      dxold <- abs(thhi - thlo)
      dx <- dxold
      temp1 <- mean(log(1 - thzero * x))
      temp2 <- mean(1/(1 - thzero * x))
      temp3 <- mean(1/(1 - thzero * x)^2)
      h <- (1 + temp1) * (temp2) - 1
      hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3))/thzero
      j <- 1
      while (j <= maxiter) {
        c1 <- (((thzero - thhi) * hprime - h) * ((thzero - thlo) * hprime - h) >= 0)
        c2 <- (abs(2 * h) > abs(dxold * hprime))
        if (c1 + c2 >= 1) {
          dxold <- dx
          dx <- (thhi - thlo)/2
          thzero <- thlo + dx
          if (thlo == thzero) {
            j <- 1000
          }
        } else {
          dxold <- dx
          dx <- h/hprime
          temp <- thzero
          thzero <- thzero - dx
          if (temp == thzero) {
            j <- 1001
          }
        }
        if (abs(dx) < epsilon * abs(thlo + thhi)/2) {
          j <- 999
        }
        temp1 <- mean(log(1 - thzero * x))
        temp2 <- mean(1/(1 - thzero * x))
        temp3 <- mean(1/(1 - thzero * x)^2)
        h <- (1 + temp1) * (temp2) - 1
        hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3))/thzero
        if (h < 0) {
          thlo <- thzero
        } else {
          thhi <- thzero
        }
        j <- j + 1
      }
      if (j > maxiter + 1) {
        thzeros[3, ] <- cbind(thzero, j)
      }
      temp1 <- mean(log(1 - thzero * x))
      temp2 <- mean(1/(1 - thzero * x))
      temp3 <- mean(1/(1 - thzero * x)^2)
      hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3))/thzero
      if (hprime > 0) {
        thlo <- thzero
        thhi <- hibnd
        thzero <- thhi
        dx <- thlo - thhi
        j <- 1
        while (j <= maxiter) {
          dx <- 0.5 * dx
          thmid <- thzero + dx
          hmid <- (1 + sum(log(1 - thmid * x))/n) * (sum(1/(1 - thmid * x))/n) - 1
          if (hmid < 0) {
            thzero <- thmid
          }  else if (hmid == 0) {
            j <- 999
          }
          if (abs(dx) < epsilon * abs(thlo + thhi)/2) {
            j <- 999
          }
          j <- j + 1
        }
        if (j > maxiter + 1) {
          thzeros[4, ] <- cbind(thzero, j)
        }
      } else {
        thlo <- epsilon
        thhi <- thzero
        thzero <- thlo
        dx <- thhi - thlo
        j <- 1
        while (j <= maxiter) {
          dx <- 0.5 * dx
          thmid <- thzero + dx
          hmid <- (1 + sum(log(1 - thmid * x))/n) * (sum(1/(1 - thmid * x))/n) - 1
          if (hmid < 0) {
            thzero <- thmid
          }  else if (hmid == 0) {
            j <- 999
          }
          if (abs(dx) < epsilon * abs(thlo + thhi)/2) {
            j <- 999
          }
          j <- j + 1
        }
        if (j > maxiter + 1) {
          thzeros[4, ] <- cbind(thzero, j)
        }
      }
    }
    thetas <- thzeros[thzeros[, 2] > maxiter + 1, ]
    nzeros <- nrow(thetas)
    proll <- rep(0, nzeros)
    mles <- matrix(rep(0, 4 * nzeros), ncol = 4)
    i <- 1
    while (i <= nzeros) {
      temp1 <- sum(log(1 - thetas[i, 1] * x))
      mles[i, 1] <- -temp1/n
      mles[i, 2] <- mles[i, 1]/thetas[i, 1]
      mles[i, 3] <- -n * log(mles[i, 2]) + (1/mles[i, 1] - 1) * temp1
      mles[i, 4] <- 999
      i <- i + 1
    }
    ind <- seq_along(mles[, 4])
    ind <- ind[mles[, 4] == 999]
    if (sum(ind) == 0) {
      nomle <- 0
    } else {
      nomle <- 1
    }
    if (nomle != 0) {
      mles <- mles[ind, ]
      outside <- which(mles[, 1] > 1 + 1e-10)
      if (length(outside) > 0) {
        mles[outside, 3] <- -Inf
      }
      mles <- rbind(mles,
                    c(1, xn, -n * log(xn), 999), # case xi=-1
                    c(0, xbar, sum(log(dexp(xq, rate = 1/xbar))), 0)) # exponential model

      nmles <- nrow(mles)
      maxlogl <- max(mles[, 3])
      ind <- order(mles[, 3])
      ind <- ind[nmles]
      k <- mles[ind, 1]
      a <- mles[ind, 2]
      conv <- 0
    } else {
      conv <- 50
      k <- NA
      a <- NA
    }
    list(k = k, a = a, conv = conv)
  }


".Zhang_Stephens_posterior" <- function(x) {
  # x: sample data from the GPD
  n <- length(x)
  x <- sort(x)
  lx <- function(b, x) {
    k <- -mean(log(1 - b * x))
    log(b/k) + k - 1
  }
  m <- 20 + floor(sqrt(n))
  b <- 1/x[n] + (1 - sqrt(m/(1:m - 0.5)))/3/x[floor(n/4 + 0.5)]
  L <- sapply(1:m, function(i) {
    n * lx(b[i], x)
  })
  w <- sapply(1:m, function(i) {
    1/sum(exp(L - L[i]))
  })
  b <- sum(b * w)
  xi <- mean(log(1 - b * x))
  sigma <- -xi/b
  mle <- c(sigma, xi)
  names(mle) <- c("scale", "shape")
  list(mle = mle)
}

".Zhang_posterior" <- function(x) {
  n <- length(x)
  x <- sort(x)
  lx <- function(b, x) {
    k <- -mean(log(1 - b * x))
    if (b == 0) {
      k - 1 - log(mean(x))
    } else {
      k - 1 + log(b/k)
    }
  }
  p <- (3:9)/10
  xp <- x[round(n * (1 - p) + 0.5)]
  m <- 20 + round(n^0.5)
  xq <- x[round(n * (1 - p * p) + 0.5)]
  k <- log(xq/xp - 1, p)
  a <- k * xp/(1 - p^k)
  a[k == 0] <- (-xp/log(p))[k == 0]
  k <- -1
  b <- w <- L <- (n - 1)/(n + 1)/x[n] - (1 - ((1:m - 0.5)/m)^k)/k/median(a)/2
  L <- sapply(1:m, function(i) n * lx(b[i], x))
  w <- sapply(1:m, function(i) 1/sum(exp(L - L[i])))
  b <- sum(b * w)
  k <- -mean(log(1 - b * x))
  mle <- c(k/b, -k)
  names(mle) <- c("scale", "shape")
  list(mle = mle)
}



#' Maximum likelihood estimate of generalized Pareto applied to threshold exceedances
#'
#' The function \code{\link[mev]{fit.gpd}} is a wrapper around \code{gp.fit}
#' @export
#' @inheritParams fit.gpd
#' @keywords internal
gp.fit <- function(xdat, threshold, method = c("Grimshaw", "auglag", "nlm", "optim", "ismev", "zs", "zhang"), show = FALSE, MCMC = NULL, fpar = NULL, warnSE = TRUE) {
  xi.tol = 1e-04
  xdat <- na.omit(xdat)
  xdatu <- xdat[xdat > threshold] - threshold
  # Optimization of model, depending on routine
  method <- match.arg(method)
  if(!is.null(fpar)){
    stopifnot(is.list(fpar),
              length(fpar) == 1L,
              names(fpar) %in% c("scale","shape"))
    if(method %in% c("zs", "zhang")){
      stop("Invalid method choice.")
    }
    method <- "fpar"
  }
  if (!is.null(MCMC) && !method %in% c("zs", "zhang"))
    warning("Ignoring argument \"MCMC\" for frequentist estimation")
  if (missing(method)) {
    method = "Grimshaw"
  }
  if (method == "nlm") {
    temp <- .Zhang_Stephens_posterior(xdatu)
    temp <- .GP_1D_fit_nlm(xdat, threshold, init.val = temp$mle[2]/temp$mle[1], gradtol = 1e-10, steptol = 1e-05, calc.se = FALSE)
    ifelse(temp$code < 2, temp$conv <- 0, temp$conv <- 50)
    # 1D max, use nlm to get gradients v close to zero
  } else if (method == "ismev") {
    temp <- .gpd_2D_fit(xdat, threshold, show = FALSE)  # ismev GP fitting function
    if (temp$conv != 0) {
      # algorithm failed to converge
      warning("Algorithm did not converge. Switching method to Grimshaw")
      temp <- .fit.gpd.grimshaw(xdatu)
      temp$mle <- c(temp$a, -temp$k)
    }
    temp <- .gpd_2D_fit(xdat, threshold, show = FALSE, siginit = temp$mle[1], shinit = temp$mle[2], method = "BFGS", reltol = 1e-30,
                        abstol = 1e-30)
  } else if (method %in% c("copt", "auglag")) {
    method <- "auglag"
    mdat <- xdat[xdat > threshold] - threshold
    maxdat <- max(mdat)
    temp <- try(alabama::constrOptim.nl(c(1,0.1),
                                        fn = gpd.ll,
                                        gr = gpd.score,
                                        hin = function(par, dat){c(par[1], par[2]+1, ifelse(par[2]<0,par[2]*maxdat + par[1],1e-4))},
                                        dat = mdat, control.outer = list(trace = FALSE),
                                        control.optim = list(fnscale = -1, trace = FALSE)))
    if(!inherits(temp, what = "try-error")){
      if(temp$convergence != 0){
        warning("Algorithm did not converge.")
        temp <- .fit.gpd.grimshaw(mdat)
        temp$mle <- c(temp$a, -temp$k)
      }
    } else{
      temp <- .fit.gpd.grimshaw(mdat)
      temp$mle <- c(temp$a, -temp$k)
    }
    temp$mle <- temp$par
    temp$nllh <- -temp$value
    nll_limxi <- -length(mdat) * log(max(mdat))
    if(temp$value <  nll_limxi){
      temp$mle <- c(max(mdat), -1+1e-7)
      temp$nllh <- -nll_limxi
    }
    temp$conv <- temp$convergence
  } else if (method == "fpar") {
    method <- "auglag"
    mdat <- xdat[xdat > threshold] - threshold
    maxdat <- max(mdat)
    param_names <- c("scale", "shape")
    wf <- param_names == names(fpar)
    if(!any(wf)){
      stop("List \"fpar\" must contain either \"scale\" or \"shape\"")
    }
    wfo <- order(c(which(!wf), which(wf)))
    fixed_param <- as.vector(unlist(fpar))
    stopifnot(length(fixed_param) == 1L)
    if(wf[2] & isTRUE(all.equal(fixed_param, 0, check.attributes = FALSE))){# shape is fixed}
      mean_xdat <- mean(xdat)
      temp <- list(par = mean_xdat,
                   value = -gpd.ll(c(mean_xdat, 0), mdat),
                   convergence = 0,
                   counts = c("function" = 1, "gradient" = 0))
    } else{
    start <- ifelse(wf[2], ifelse(fixed_param < 0, 1.1*maxdat/abs(fixed_param), mean(xdat)), 0.1)
    temp <- try(alabama::auglag(par = start, fpar = fixed_param,
                                wf = wf, wfo = wfo,
                                fn = function(par, fpar, wfo, wf, dat, ...){
                                          -gpd.ll(c(par,fpar)[wfo], dat = dat)},
                                        gr = function(par, fpar, wfo, wf, dat, ...){
                                          - gpd.score(c(par,fpar)[wfo], dat = dat)[!wf]},
                                        hin = function(par, fpar, wfo, wf, dat, ...){
                                          param <- c(par,fpar)[wfo]
                                          c(param[1], param[2]+1,
                                            ifelse(param[2]<0, param[2]*maxdat + param[1],1e-4))},
                                        dat = mdat, control.outer = list(trace = FALSE),
                                        control.optim = list(trace = FALSE)))
      if(!inherits(temp, what = "try-error")){
      if(!isTRUE(all(temp$kkt1, temp$kkt2)) && !all.equal(c(temp$par), -1)){
        warning("Optimization algorithm did not converge.")
      } else{
        temp$convergence == 0
      }
    } else{
      stop("Optimization algorithm did not converge.")
    }
  }
    temp$mle <- temp$par
    temp$param <- c(temp$par, fixed_param)[wfo]
    temp$nllh <- temp$value
    temp$conv <- temp$convergence
    # Initialize standard errors and covariance matrix
    temp$vcov <- matrix(NA)
    temp$se <- NA
    if(temp$param[2] > -0.5){
      infomat <- gpd.infomat(dat = mdat, par = temp$param, method = "obs")[!wf,!wf]
      if(isTRUE(c(infomat) > 0)){
      temp$vcov <- matrix(1/infomat)
      temp$se <- sqrt(c(temp$vcov))
      }
    }
    names(temp$mle) <- names(temp$se) <- param_names[!wf]
    names(temp$param) <- param_names
    output <- structure(list(estimate = temp$mle,
                             std.err = temp$se,
                             param = temp$param,
                             vcov = temp$vcov,
                             threshold = threshold,
                             method = method,
                             nllh = temp$nllh,
                             nat = sum(xdat > threshold),
                             pat = length(xdatu)/length(xdat),
                             convergence = ifelse(temp$conv == 0, "successful", temp$conv),
                             counts = temp$counts,
                             exceedances = xdatu,
                             wfixed = wf),
                        class = "mev_gpd")
    if (show) {
      print(output)
    }
    return(invisible(output))
  } else if (method == "optim") {
    temp <- .gpd_1D_fit(xdat, threshold, show = FALSE, xi.tol = xi.tol)  # 1D max, algebraic Hessian
    if (temp$conv != 0) {
      # algorithm failed to converge
      warning("Algorithm did not converge. Switching method to Grimshaw")
      temp <- .fit.gpd.grimshaw(xdat[xdat > threshold] - threshold)
      temp$mle <- c(temp$a, -temp$k)
    }
    temp <- .gpd_1D_fit(xdat, threshold, show = FALSE, xi.tol = xi.tol,
                        phi.input = temp$mle[2]/temp$mle[1], reltol = 1e-30, abstol = 1e-30)
    #1D max, algebraic Hessian
  } else if (method == "zs" || method == "zhang") {
    xdat_ab = xdat[xdat > threshold] - threshold
    temp <- switch(method, zs = .Zhang_Stephens_posterior(xdat_ab), zhang = .Zhang_posterior(xdat_ab))
    if (!is.null(MCMC) && MCMC != FALSE) {
      if (is.logical(MCMC) && !is.na(MCMC) && MCMC) {
        # if TRUE
        burn = 2000
        thin = 1
        niter = 10000
      } else if (is.list(MCMC)) {
        burn <- ifelse(is.null(MCMC$burnin), 2000, MCMC$burnin)
        thin <- ifelse(is.null(MCMC$thin), 1, MCMC$thin)
        niter <- ifelse(is.null(MCMC$niter), 10000, MCMC$niter)
      }
      bayespost <- switch(method, zs = Zhang_Stephens(xdat_ab, init = -temp$mle[2]/temp$mle[1], burnin = burn, thin = thin,
                                                      niter = niter, method = 1), zhang = Zhang_Stephens(xdat_ab, init = -temp$mle[2]/temp$mle[1], burnin = burn, thin = thin,
                                                                                                         niter = niter, method = 2))
      # Copying output for formatting and printing
      post.mean <- bayespost$summary[1, ]
      post.se <- sqrt(bayespost$summary[2, ])
      names(post.mean) <- names(post.se) <- c("scale", "shape")
      post <- structure(list(method = method,
                             estimate = post.mean,
                             param = post.mean,
                             threshold = threshold,
                             nat = sum(xdat > threshold),
                             pat = sum(xdat > threshold)/length(xdat),
                             approx.mean = temp$mle,
                             post.mean = post.mean,
                             post.se = post.se,
                             accept.rate = bayespost$rate,
                             thin = bayespost$thin,
                             burnin = bayespost$burnin,
                             niter = bayespost$niter,
                             exceedances = xdatu), class = c("mev_gpdbayes", "mev_gpd"))

    } else {
      post <- structure(list(method = method,
                             estimate = temp$mle,
                             param = temp$mle,
                             threshold = threshold,
                             nat = sum(xdat > threshold),
                             pat = sum(xdat > threshold)/length(xdat),
                             approx.mean = temp$mle,
                             exceedances = xdatu), class = c("mev_gpdbayes"))

    }
    if (show)
      print(post)
    return(invisible(post))
  }
  if (method == "Grimshaw") {
    pjn <- .fit.gpd.grimshaw(xdatu)  # Grimshaw (1993) function, note: k is -xi, a is sigma
    temp <- list()
    temp$mle <- c(pjn$a, -pjn$k)  # mle for (sigma, xi)
    sc <- rep(temp$mle[1], length(xdatu))
    xi <- temp$mle[2]
    temp$nllh <- sum(log(sc)) + sum(log(1 + xi * xdatu/sc) * (1/xi + 1))
    temp$conv <- pjn$conv
  }
  if(temp$mle[2] < -1){
    #Transform the solution (unbounded) to boundary - with maximum observation for scale and -1 for shape.
    temp$mle <- c(max(xdatu) + 1e-10, -1)
  }

  # Collecting observations from temp and formatting the output
  invobsinfomat <- tryCatch(solve(gpd.infomat(dat = xdatu,
                                              par = c(temp$mle[1], temp$mle[2]), method = "obs")),
                            error = function(e) {
                              "notinvert"
                            }, warning = function(w) w)
  if (any(c(isTRUE(invobsinfomat == "notinvert"), all(is.nan(invobsinfomat)), all(is.na(invobsinfomat))))) {
    if(warnSE){
    warning("Cannot calculate standard errors based on observed information")
    }
    if (!is.null(temp$se)) {
      std.errors <- diag(temp$se)
    } else {
      std.errors <- rep(NA, 2)
    }
  } else if (!is.null(temp$mle) && temp$mle[2] > -0.5 && temp$conv == 0) {
    # If the MLE was returned
    std.errors <- sqrt(diag(invobsinfomat))
  } else {
    if(warnSE){
      warning("Cannot calculate standard error based on observed information")
    }
    std.errors <- rep(NA, 2)
  }
  if (temp$mle[2] < -1 && temp$conv == 0) {
    warning("The MLE is not a solution to the score equation for \"xi < -1'")
  }
  names(temp$mle) <- names(std.errors) <- c("scale", "shape")
  output <- structure(list(estimate = temp$mle,
                           param = temp$mle,
                           std.err = std.errors,
                           vcov = invobsinfomat,
                           threshold = threshold,
                           method = method,
                           nllh = temp$nllh,
                           nat = sum(xdat > threshold),
                           pat = length(xdatu)/length(xdat),
                           convergence = ifelse(temp$conv == 0, "successful", temp$conv),
                           counts = temp$counts,
                           exceedances = xdatu,
                           wfixed =rep(FALSE, 2L)),
                      class = "mev_gpd")
  if (show) {
    print(output)
  }
  invisible(output)
}

#' Robust threshold selection of Dupuis
#'
#' The optimal bias-robust estimator (OBRE) for the generalized Pareto.
#' This function returns robust estimates and the associated weights.
#'
#' @references Dupuis, D.J. (1998). Exceedances over High Thresholds: A Guide to Threshold Selection,
#' \emph{Extremes}, \bold{1}(3), 251--261.
#'
#' @param dat a numeric vector of data
#' @param thresh threshold parameter
#' @param k bound on the influence function; the constant \code{k} is a robustness parameter
#' (higher bounds are more efficient, low bounds are more robust). Default to 4.
#' @param tol numerical tolerance for OBRE weights iterations.
#' @param show logical: should diagnostics and estimates be printed. Default to \code{FALSE}.
#' @seealso \code{\link{fit.gpd}}
#' @return a list with the same components as \code{\link{fit.gpd}},
#' in addition to
#' \itemize{
#' \item{\code{estimate}:}{optimal bias-robust estimates of the \code{scale} and \code{shape} parameters.}
#' \item{\code{weights}:}{vector of OBRE weights.}
#' }
#' @keywords internal
#' @importFrom utils head
#' @export
#' @examples
#' dat <- rexp(100)
#' .fit.gpd.rob(dat, 0.1)
.fit.gpd.rob <- function(dat, thresh, k = 4, tol = 1e-5, show = FALSE){
  k <- max(k, sqrt(2))
  ninit <- length(dat)
  gpd.score.i <- function(par, dat) {
    sigma = par[1]
    xi = as.vector(par[2])
    if (!isTRUE(all.equal(0, xi))) {
      rbind(dat * (xi + 1)/(sigma^2 * (dat * xi/sigma + 1)) - 1/sigma,
            (-dat * (1/xi + 1)/(sigma *(dat * xi/sigma + 1)) + log1p(pmax(-1, dat * xi/sigma))/xi^2))
    } else {
      rbind((dat - sigma)/sigma^2, (1/2 * (dat - 2 * sigma) * dat/sigma^2))
    }
  }
  #B-score weight function
  Wfun <- function(dat, par, A, a, k){
    pmin(1, k/sqrt(colSums((A %*% (gpd.score.i(dat = dat, par = par) - a))^2)))
  }
  afun <- function(par, A, a, k){
    upbound <- ifelse(par[2] < 0, -par[1]/par[2], Inf)
    c(integrate(f = function(x){
      gpd.score.i(par = par, dat = x)[1,] *
        Wfun(dat = x, par = par, A = A, a = a, k = k) *
        mev::dgp(x = x, loc = 0, scale = par[1], shape = par[2])
    }, lower = 0, upper = upbound)$value,
    integrate(f = function(x){
      gpd.score.i(par = par, dat = x)[2,] *
        Wfun(dat = x, par = par, A = A, a = a, k = k) *
        mev::dgp(x = x, loc = 0, scale = par[1], shape = par[2])
    }, lower = 0, upper = upbound)$value)/
      integrate(f = function(x){
        Wfun(dat = x, par = par, A = A, a = a, k = k) *
          mev::dgp(x = x, loc = 0, scale = par[1], shape = par[2])
      }, lower = 0, upper = upbound)$value
  }
  #Initialize algorithm
  # keep only exceedances
  dat <- as.vector(dat[dat > thresh] - thresh)
  # starting value is maximum likelihood estimates
  par_val <- gp.fit(dat, threshold = 0)$estimate
  a <- rep(0, 2)
  A <- t(solve(chol(gpd.infomat(par = par_val, dat = dat, method = "obs"))))
  dtheta <- rep(Inf, 2)
  u <- runif(1e4);
  niter <- 0L;
  niter_max <- 1e3L
  while(isTRUE(max(abs(dtheta/par_val)) > tol) && niter < niter_max){
    niter <- niter + 1L
    if(all(is.finite(dtheta))){
      par_val <- par_val + dtheta
    }
    #Monte-Carlo integration with antithetic variables
    xd <- mev::qgp(c(u, 1-u), loc = 0, scale = par_val[1], shape = par_val[2])
    score <- gpd.score.i(par = par_val, dat = xd)
    Wc <- Wfun(dat = xd, par = par_val, A = A, a = a, k = k)
    a <- rowSums(score %*% Wc)/ sum(Wc)
    Wcsq <- Wfun(dat = xd, par = par_val, A = A, a = a, k = k)^2
    M2e <- c(mean((score[1,]-a[1])^2*Wcsq),
             mean((score[1,]-a[1])*(score[2,]-a[2])*Wcsq),
             mean((score[2,]-a[2])^2*Wcsq))
    M2 <- matrix(c(M2e[1], M2e[2], M2e[2], M2e[3]), ncol = 2, nrow = 2, byrow = TRUE)
    #Compute Matrix A
    A <- chol(solve(M2))
    #Compute matrix M1 with new values of a, A and Delta theta
    Wc <- Wfun(dat = xd, par = par_val, A = A, a = a, k = k)
    M1e <- c(mean((score[1,]-a[1])^2*Wc),
             mean((score[1,]-a[1])*(score[2,]-a[2])*Wc),
             mean((score[2,]-a[2])^2*Wc))
    M1 <- matrix(c(M1e[1], M1e[2], M1e[2], M1e[3]), ncol = 2, nrow = 2, byrow = TRUE)
    Wgt_dat <- Wfun(dat = dat, par = par_val, A = A, a = a, k = k)
    dtheta <- c(solve(M1) %*% colMeans(t(gpd.score.i(dat = dat, par = par_val) - a) * Wgt_dat))
  }
  #Test for standard errors - Huber robust sandwich variance
  Wc <- Wfun(dat = xd, par = par_val, A = A, a = a, k = k)
  M3 <- c(mean((score[1,]-a[1])*score[1,]*Wc),
          mean((score[1,]-a[1])*score[2,]*Wc),
          #mean((score[2,]-a[2])*score[1,]*Wc),
          # should be symmetric, not quite numerically but error 10e-8
          mean((score[2,]-a[2])*score[2,]*Wc))
  #Compute inverse of -\int \partial/\partial theta \psi dF(\theta)
  Masy <- solve(matrix(c(M3[1], M3[2], M3[2], M3[3]), ncol = 2, nrow = 2, byrow = TRUE))
    #Hampel et al. (1986), eq. 4.2.13, p. 231
  vcov <- Masy %*% M1 %*% Masy / length(dat)
  colnames(vcov) <- rownames(vcov) <- c("scale","shape")
  stderr <- sqrt(diag(vcov))
  names(stderr) <- names(par_val)
  convergence <- "successful"
  conv <- 0L
  if(niter == niter_max){
    convergence <- "Algorithm did not converge: maximum number of iterations reached"
    conv <- 1L
  }
  if(!isTRUE(is.finite(gpd.ll(par_val, dat = dat)))){
    convergence <- "Solution not feasible; algorithm aborted."
    conv <- 2L
  }
  ret <- structure(list(estimate = par_val,
                        std.err = stderr,
                        vcov = vcov,
                        threshold = thresh,
                        method = "obre",
                        nllh = ifelse(conv == 2L, Inf, -as.vector(gpd.ll(par_val, dat = dat))),
                        convergence = convergence,
                        nat = length(dat),
                        pat = length(dat)/ninit,
                        counts = c("function" = niter),
                        param = par_val,
                        exceedances = dat,
                        weights = Wgt_dat), class = "mev_gpd")
  if(show){
    print(ret)
    if(convergence == "successful"){
    matw <- head(cbind("exceedances" =  ret$exceedances,
                       "weights" = ret$weights,
                       "p-value" = rank(ret$weights)/length(ret$weights))[order(ret$exceedances, decreasing = TRUE),])
    rownames(matw) <- 1:6
    matw <- matw[matw[,2] != 1,]
    if(length(matw)> 0){
      cat("Largest observations: OBRE weights\n")
      print(round(matw, digits = 3))
    }
    }
  }
  return(invisible(ret))
}
