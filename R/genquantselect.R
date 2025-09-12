#' Generalized quantile threshold selection
#'
#' The methodology proposed by Beirlant, Vynckier and Teugels (1996)
#' uses an asymptotic expansion of the mean squared error for Hill's estimator given
#' a random sample with Pareto tails and positive shape, using an exponential regression.
#' The value of \code{k} is selected to minimize the mean squared error given optimal weighting scheme. This
#' depends on the order of regular variation \eqn{\rho}, which is obtained based on the slope of the difference
#' in Hill estimators, suitably reweighted. The iterative procedure of Beirlant et al. alternates between parameter estimation
#' until convergence. It returns the Hill shape estimate, the number of higher order statistic, the parameter \code{rho} and
#' estimates of the standard error of the shape and the mean squared error, based on the ultimate parameter values.
#' Since the weights can become negative, there is no guarantee that the mean squared error estimate is positive, nor that
#' the estimated value of \eqn{\rho} is nonpositive.
#'
#' @references Beirlant, J., Vynckier, P., & Teugels, J. L. (1996). Excess Functions and Estimation of the Extreme-Value Index. \emph{Bernoulli}, 2(\bold{4}), 293–318. \doi{10.2307/3318416}
#' @param xdat [vector] sample of exceedances
#' @param maxiter [int] maximum number of iteration
#' @param tol [double] tolerance for difference in value of \eqn{k} for the fixed point
#' @param kmin [int] minimum number of exceedances for the estimator
#' @param kmax [int] maximum number of exceedances for the estimator
#' @param ... additional arguments, currently ignored
#' @return a list with components
#' \itemize{
#' \item \code{shape} the Hill estimator of the shape, based on the \code{k} largest order statistics
#' \item \code{k0} number of high order statistics for estimation of the shape using Hill's estimator
#' \item \code{rho} estimate of the second order regular variation parameter
#' \item \code{mse} mean squared error estimate of the shape parameter
#' \item \code{se} standard error of the shape parameter
#' \item \code{convergence} logical; if \code{TRUE}, indicates that the method converged to a fixed point within \code{tol} before reaching the maximum number of iterations \code{maxiter}
#' }
#' @export
#' @examples
#' # Simulate Pareto data - log(xdat) is exponential with rate 2
#' xdat <- rgp(n = 200, loc = 1, scale = 0.5, shape = 0.5)
#' (thselect.expgqt(xdat))
thselect.expgqt <- function(
  xdat,
  maxiter = 10L,
  tol = 2,
  kmin = max(10, floor(length(xdat) / 100)),
  kmax = floor(0.8 * length(xdat)),
  ...
  # rho = c("btv","fgh"),
  # plot = FALSE
) {
  xdat <- xdat[xdat > 0 & is.finite(xdat)]
  xdat <- sort(xdat, decreasing = TRUE)
  logx <- log(xdat)
  n <- length(logx)
  kmin <- as.integer(kmin)
  kmax <- as.integer(kmax)
  if (kmax <= kmin) {
    stop(
      "Invalid input: the series is too short for estimation of the tail index"
    )
  }
  # alphaf <- function(j, k) {
  #   stopifnot(j <= k)
  #   sum(1 / (k - 1:(k - j + 1) + 1))
  # }
  # betaf <- function(j, k) {
  #   stopifnot(j <= k)
  #   sum(1 / (k - 1:(k - j + 1) + 1) ^ 2)
  # }
  # cf <- function(k, weight = c("1", "2")) {
  #   stopifnot(k >= 1)
  #   weight <- match.arg(weight)
  #   if (weight == "1") {
  #     w <- 1
  #   } else{
  #     w <- (1:k) / (k + 1)
  #   }
  #   k <- as.integer(k)
  #   betas <- sapply(1:k, function(j) {
  #     betaf(j, k)
  #   })
  #   alphas <- sapply(1:k, function(j) {
  #     alphaf(j, k)
  #   })
  #   sum(w * (betas + (alphas - log(k + 1) + log(1:k)) ^ 2))
  # }
  # c1k <- function(k){cf(k = k, weight = "1")}
  # c2k <- function(k){cf(k = k, weight = "2")}
  af <- function(k, weight, shape) {
    k <- as.integer(k)
    stopifnot(length(k) == 1L, length(weight) %in% c(1, k))
    if (shape >= 0) {
      mean(weight * ((k + 1) / (1:k) - 1))
    } else {
      mean(
        weight *
          ((k + 1) /
            (1:k) +
            (-4 * (((1:k) / (k + 1))^(-shape)) + (3 - shape - 2 * shape^2)) /
              (1 + shape + 2 * shape^2))
      )
    }
  }
  bf <- function(k, weight, rho) {
    k <- as.integer(k)
    rho <- as.numeric(rho)[1]
    stopifnot(length(k) == 1L, length(weight) == k, rho <= 0)
    krho <- function(t, rho) {
      if (!isTRUE(all.equal(rho, 0))) {
        (t^rho - 1) / rho
      } else {
        log(t)
      }
    }
    (1 - rho)^2 * mean(weight * krho((1:k) / (k + 1), rho)^2)
  }
  # Weighting functions of Beirlant, Vynckier and Teugels (1996b), Bernoulli
  w1f <- function(k) {
    (1:k) / (k + 1)
  }
  w2f <- function(k) {
    log(k + 1) - log(1:k)
  }
  # Optimal weighting
  wopt <- function(k, shape, rho) {
    w1 <- w1f(k)
    w2 <- w2f(k)
    a1 <- af(k = k, shape = shape, weight = w1)
    a2 <- af(k = k, shape = shape, weight = w2)
    b1 <- bf(k = k, rho = rho, weight = w1)
    b2 <- bf(k = k, rho = rho, weight = w2)
    (w1 * (b2 - a2) + w2 * (a1 - b1)) / (a1 * b2 - b1 * a2)
  }
  # Step 1: compute weights c with w = "1" and w = "2" for each candidate k number of order statistics
  # Step 2: compute Hill's estimator of the shape
  fithill <- function(logdata, sorted) {
    logdata <- logdata[is.finite(logdata)]
    n <- length(logdata)
    if (!isTRUE(sorted)) {
      logdata <- sort(logdata, decreasing = TRUE)
    }
    cumlogdat <- cumsum(logdata)
    k <- 1:(n - 1)
    cumlogdat[k] / k - logdata[k + 1]
  }

  hillest <- fithill(logx, sorted = TRUE)
  logUH <- logx[-length(logx)] + log(hillest)
  # if(isTRUE(plot)){
  # plot(x = -log(1:(length(logx)-1)/length(logx)),
  #       y = logUH,
  #       type = "l",
  #       bty = "l",
  #       xlab = "-log(j/n)",
  #       ylab = "generalized Hill")
  # }
  shape_genhill <- function(k, logUH, weight = 1) {
    sum(weight * (log(k + 1) - log(1:k)) * (logUH[1:k] - logUH[k + 1])) /
      sum(weight * ((log(k + 1) - log(1:k))^2))
  }
  shapes <- sapply(1:(n - 2), function(k) {
    shape_genhill(k = k, logUH = logUH)
  })
  # if(isTRUE(plot)){
  # plot(-log(11:(n-2)) +log(n),
  #      shapes[-(1:10)], bty = "l", type = "l",
  #      ylab = "shape", xlab = "-log(j/n)")
  # }
  # Distance with distance to unit exponential Q-Q plot
  mse_a_fun <- function(logdata, k, shape, weights = 1) {
    # if(!isTRUE(sorted)){
    #   logdata <- sort(logdata, decreasing = TRUE)
    # } else {
    # # Poor man sanity check
    # # Data must be sorted in decreasing order
    #   stopifnot(logdata[1] >= logdata[2])
    # }
    mean(
      weights *
        (logdata[1:k] - logdata[k + 1] - shape * (log(k + 1) - log(1:k)))^2
    )
  }
  # Pick the index that minimizes the distance among candidates
  kmin <- max(ceiling(n / 100), 10)
  # Value of k0 must not be too large to estimate m0
  kcand <- kmin:floor(n^0.9)
  mse_a <- sapply(kcand, function(k) {
    mse_a_fun(logUH, k = k, shape = shapes[k])
  })
  k0 <- kmin - 1 + which.min(mse_a)
  # if(isTRUE(plot)){
  #  plot(x = kcand,
  #       y = mse_a,
  #       type = "l",
  #       bty = "l",
  #       ylab = "mean squared error",
  #       xlab = "k",
  #       panel.first = { abline(v = k0, lty = 2)})
  # }

  # Return the Hill estimator for the associated index
  gamma0 <- shapes[k0]

  # Estimator of the second order regular variation coefficient
  rho_mn <- function(m, k, shape_est) {
    log(abs(
      (shape_est[floor((m + k) / 4)] - shape_est[floor((m + k) / 2)]) /
        (shape_est[floor(m / 2)] - shape_est[m])
    )) /
      (log(2 * m) - log(m + k))
  }
  mse_rho_fun <- function(m, k, shape_est) {
    Rmn <- rho_mn(m = m, k = k, shape_est = shape_est)
    mean(sapply(k:(m - 1), function(j) {
      (log(abs(
        (shape_est[floor(j / 2)] - shape_est[j]) /
          (shape_est[floor(m / 2 + 0.5)] - shape_est[m + 1])
      )) -
        Rmn * (log(m + 1) - log(j)))^2
    }))
  }
  # Find index for second order regular variation coefficient
  mcand <- (k0 + 10):(n - 2)
  mse_rho <- sapply(mcand, function(m) {
    mse_rho_fun(m = m, k = k0, shape_est = shapes)
  })
  m0 <- mcand[which.min(mse_rho)]
  # if(isTRUE(plot)){
  #   plot(x = mcand,
  #        y = mse_rho,
  #        type = "l",
  #        bty = "l",
  #        ylab = "mean squared error",
  #        xlab = "k",
  #        panel.first = { abline(v = m0, lty = 2)})
  #   rho_est <- sapply(mcand, function(m){
  #     rho_mn(m = m, k = k0, shape_est = shapes)
  #   })
  #   plot(x = mcand, y = rho_est, bty = "l", type = "l",
  #        panel.first = {abline(v = m0, lty = 2)})
  # }
  # Initial nonpositive index of 2nd order RV
  rho0 <- rho_mn(m = m0, k = k0, shape_est = shapes)
  #if(rho0 < -2 | rho0 > 0){
  # if(gamma0 > 0){
  #   rho0 <- rho_fa(logdata = logx, k = k0, tau = 0)
  # }
  #}
  rho0 <- max(min(0, rho0), -1)
  # Weights
  # zetaf <- function(rho, k){
  #   ((((1:k)/(k+1))^(-rho) - 1)/rho)^2
  # }
  # d1k <- function(rho, k){
  #  (1-rho)^2 * mean(zetaf(rho = rho, k = k))
  # }
  # d2k <- function(rho, k){
  #   (1-rho)^2 / ((k+1) * k) * sum((1:k)*zetaf(rho = rho, k = k))
  # }
  # wopt <- function(k, rho){
  #   d2c <- d2k(rho = rho, k = k)
  #   d1c <- d1k(rho = rho, k = k)
  #   c2c <- c2k(k = k)
  #   c1c <- c1k(k = k)
  #   (d2c - c2c + (1:k)/k*(c1c-d1c))/(c1c*d2c-c2c*d1c)
  # }
  # wbias <- function(k, rho){
  #   d2c <- d2k(rho = rho, k = k)
  #   d1c <- d1k(rho = rho, k = k)
  #   c2c <- c2k(k)
  #   c1c <- c1k(k)
  #   (c2c - c1c*(1:k)/k)/(c1c*d2c-c2c*d1c)
  # }

  # Run for loop to determine the values
  kvals <- mvals <- integer(maxiter)
  rhovals <- gammavals <- numeric(maxiter)
  conv <- FALSE
  # Initialize values
  k_cur <- k0
  m_cur <- m0
  rho_cur <- rho0
  gamma_cur <- gamma0
  nstep <- 0
  while (!conv & (nstep < maxiter)) {
    # TODO figure out why the weights are increasing...
    mse_a_cur <- sapply((n - 1):kmin, function(k) {
      mse_a_fun(
        logdata = logUH,
        k = k,
        shape = gamma_cur,
        # Weights depend on k , but are increasing! nonsense?
        weights = wopt(k = k, shape = gamma_cur, rho = rho_cur)
      )
    })

    # plot(mse_a_cur, type = "l")
    k_new <- n - which.min(mse_a_cur)
    gamma_new <- shapes[k_new]
    m_new <- n -
      which.min(
        sapply((n - 1):(k_new + 1), function(m) {
          mse_rho_fun(m = m, k = k_new, shape_est = shapes)
        })
      )
    rho_new <- rho_mn(m = m_new, k = k_new, shape_est = shapes)
    # if(rho_new < -1 | rho_new > 0){
    #   if(gamma_new > 0){
    #    rho_new <- rho_fa(logdata = logx, k = k_new, tau = 1)
    #   } else{
    rho_new <- max(min(0, rho_new), -2)
    #   }
    # }
    nstep <- nstep + 1L
    if (abs(k_new - k_cur) < tol) {
      conv <- TRUE
    }
    # Store values
    kvals[nstep] <- k_cur
    mvals[nstep] <- m_cur
    rhovals[nstep] <- rho_cur
    gammavals[nstep] <- gamma_cur
    # Update current values
    k_cur <- k_new
    rho_cur <- rho_new
    gamma_cur <- gamma_new
    rho_cur <- rho_new
    # if(conv | nstep == maxiter){
    # At convergence
    # mse <-  mse_a_fun(logdata = logUH,
    #                   k = k_cur,
    #                   shape = gamma_cur,
    #                   weights = wopt(k = k_cur, rho = rho_cur, shape = gamma_cur))
    # sqbias <-  mse_a_fun(logdata = logx,
    #                     k = k_cur,
    #                     shape = gamma_cur,
    #                     weights = wbias(k = k_cur, rho = rho_cur))
  }
  # }
  # End of while-loop
  # if(mse < 0){
  #   warning("Invalid output for procedure: negative mean squared error")
  #   mse <- NA
  # }
  # if(rho_cur > 0){
  #   warning("Invalid output for procedure: second order regular variation index isn't nonpositive.")
  #   rho_cur <- NA
  # }
  # if(isTRUE(mse > sqbias)){
  #  se <- sqrt(mse - sqbias)
  # } else{
  #   se <- NA
  # }
  res <- list(
    k0 = k_cur,
    shape = gamma_cur,
    rho = rho_cur,
    thresh0 = exp(logx[k_cur + 1]),
    # mse = mse,
    # se = se,
    convergence = conv
  )
  class(res) <- "mev_thselect_expgqt"
  return(invisible(res))
}

#' @export
print.mev_thselect_expgqt <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  cat(
    "Threshold selection method: Beirlant, Vynckier and Teugels (1996)\nGeneralized quantile threshold selection\n\n"
  )
  if (isTRUE(x$convergence)) {
    cat("Selected threshold:", round(x$thresh0, digits), "\n")
    cat("Number of exceedances:", round(x$k0, digits), "\n")
    cat("Shape estimate:", round(x$shape, digits), "\n")
  } else {
    cat("Algorithm did not converge.\n")
  }
  return(invisible(NULL))
}
