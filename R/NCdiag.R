####################################################################################################### Functions to create plots and tests in article Northrop, Paul J and Claire L Coleman (2014), Extremes 'Improved threshold
####################################################################################################### diagnostic plots for extreme value analyses'
#------------------------------------------------------------------------------#
# Main function #
#------------------------------------------------------------------------------#
#' Score and likelihood ratio tests fit of equality of shape over multiple thresholds
#'
#' The function returns a P-value path for the score testand/or likelihood ratio
#' test for equality of the shape parameters over
#' multiple thresholds under the generalized Pareto model.
#'
#' @param xdat numeric vector of raw data
#' @param u \code{m}-vector of thresholds (sorted from smallest to largest)
#' @param GP.fit function used to optimize the generalized Pareto model.
#' @param do.LRT boolean indicating whether to perform the likelihood ratio test (in addition to the score test)
#' @param size level at which a horizontal line is drawn on multiple threshold plot
#' @param xi.tol numerical tolerance for threshold distance; if the absolute value of \code{xi1.hat} is less than \code{xi.tol} use linear interpolation
#'                to evaluate score vectors, expected Fisher information matrices, Hessians
#' @param plot logical; if \code{TRUE}, return a plot of p-values against threshold.
#' @param ... additional parameters passed to \code{plot}
#' @details The default method is \code{'Grimshaw'} using the reduction of the parameters to a one-dimensional
#' maximization. Other options are one-dimensional maximization of the profile the \code{nlm} function or \code{optim}.
#' Two-dimensional optimisation using 2D-optimization \code{\link[ismev]{ismev}} using the routine
#' from \code{gpd.fit} from the \code{ismev} library, with the addition of the algebraic gradient.
#' The choice of \code{GP.fit} should make no difference but the options were kept.
#' \bold{Warning}: the function will not recover from failure of the maximization routine, returning various error messages.
#'
#'
#' @references Grimshaw (1993). Computing Maximum Likelihood Estimates for the Generalized
#'  Pareto Distribution, \emph{Technometrics}, \bold{35}(2), 185--191.
#' @references Northrop & Coleman (2014). Improved threshold diagnostic plots for extreme value
#' analyses, \emph{Extremes}, \bold{17}(2), 289--303.
#' @references Wadsworth & Tawn (2012). Likelihood-based procedures for threshold
#' diagnostics and uncertainty in extreme value modelling, \emph{J. R. Statist. Soc. B}, \bold{74}(3), 543--567.
#' @export
#' @author Paul J. Northrop and Claire L. Coleman
#' @return a plot of P-values for the test at the different thresholds \code{u}
#' @examples
#' \dontrun{
#' data(nidd)
#' u <- seq(65,90, by = 1L)
#' NC.diag(nidd, u, size = 0.05)
#' }
#' @keywords internal
NC.diag <- function(
  xdat,
  u,
  GP.fit = c("Grimshaw", "nlm", "optim", "ismev"),
  do.LRT = FALSE,
  size = NULL,
  plot = TRUE,
  ...,
  xi.tol = 0.001
) {
  args <- list(...)
  if (any(diff(u) <= 0)) {
    warning("Thresholds supplied in u are not in increasing order")
  }
  x <- as.numeric(xdat[is.finite(xdat)])
  u <- sort(u)
  n_u <- length(u) # total number of thresholds
  #------------------------------------------------------------------------------#
  # 1. Fit GP distribution to excesses of u[i], i=1, ..., n_u #
  #------------------------------------------------------------------------------#
  GP.fit <- match.arg(GP.fit)
  z <- list() # list to store the results
  z$thresh <- u # all thresholds
  z$nexc <- unlist(lapply(u, function(y) {
    sum(x > y)
  })) # number of excesses of each threshold
  z$n.between <- c(-diff(z$nexc), z$nexc[n_u])
  # sample sizes between thresholds (and above the highest threshold)
  for (i in seq_len(n_u)) {
    # loop over all thresholds
    if (GP.fit == "nlm") {
      temp <- .gpd_1D_fit(x, u[i], show = FALSE, xi.tol = xi.tol, calc.se = F) # threshold u[j1], 1D max, algebraic Hessian
      phi.init <- temp$mle[2] / temp$mle[1] # better initial estimate of phi
      temp <- .GP_1D_fit_nlm(
        x,
        u[i],
        init.val = phi.init,
        gradtol = 1e-20,
        steptol = 1e-20,
        calc.se = F
      )
      # 1D max, use nlm to get gradients v close to zero
    } else if (GP.fit == "ismev") {
      temp <- .gpd_2D_fit(x, u[i], show = FALSE) # threshold u[i], ismev GP fitting function
      temp <- .gpd_2D_fit(
        x,
        u[i],
        show = FALSE,
        siginit = temp$mle[1],
        shinit = temp$mle[2],
        method = "BFGS",
        reltol = 1e-30,
        abstol = 1e-30
      )
    } else if (GP.fit == "optim") {
      temp <- .gpd_1D_fit(x, u[i], show = FALSE, xi.tol = xi.tol) # threshold u[i], 1D max, algebraic Hessian
      temp <- .gpd_1D_fit(
        x,
        u[i],
        show = FALSE,
        xi.tol = xi.tol,
        phi.input = temp$mle[2] / temp$mle[1],
        reltol = 1e-30,
        abstol = 1e-30
      )
      # threshold u[i], 1D max, algebraic Hessian
    } else if (GP.fit == "Grimshaw") {
      yy <- x[x > u[i]] - u[i] # thresholds excesses
      pjn <- .fit.gpd.grimshaw(yy) # Grimshaw (1993) function, note: k is -xi, a is sigma
      temp <- list()
      temp$mle <- c(pjn$a, -pjn$k) # mle for (sigma, xi)
      sc <- rep(temp$mle[1], length(yy))
      xi <- temp$mle[2]
      temp$nllh <- sum(log(sc)) + sum(log(1 + xi * yy / sc) * (1 / xi + 1))
    }
    z$xi.mle[i] <- temp$mle[2] # MLE of xi
    z$sigma.mle[i] <- temp$mle[1] # MLE of sigma1
    z$nllh[i] <- temp$nllh # negated log-likelihood at MLE
  }
  # .....................# end of loop over thresholds
  #------------------------------------------------------------------------------#
  # 2. For each threshold u[i], i=1, ..., m-1 calculate score-based p-value for # the test of H_0: xi[i]=...=xi[m], equivalently
  # H_0: phi_1=...=phi_m # ... and produce plot of the p-values vs threshold.  #
  #------------------------------------------------------------------------------#
  # Create functions score.test(), mult.u.gpd.fit() and mult.thresh.LR.test() ...
  #-------------------------- Start of function score.test() ------------------------------#
  ##########################################################################################
  score.test <- function(my.data, m, vi, wi, pars, sigma1, xi1) {
    #--------------------- if xi.hat is not very close to 0 ... ---------------------------#
    if (abs(xi1) >= xi.tol) {
      score <- .score_algebraic(my.data, pars, wi, vi, m) # score vector under H_0
      e.info <- n.exc * .exp_info_algebraic(pars, wi, vi, m) # expected information under H_0
    }
    #----------------------- if xi.hat is very close to 0 ... -----------------------------#
    if (abs(xi1) < xi.tol) {
      delta <- 2 * xi.tol # evaluate score/info at xi+delta and xi-delta
      phis <- 1 / (sigma1 / (xi1 + delta) + cumsum(c(0, wi[-m]))) # values of phi under H_0
      pars1 <- c(sigma1, phis) # estimates of sigma1, phi_1, ..., phi_m under H_0.
      phis <- 1 / (sigma1 / (xi1 - delta) + cumsum(c(0, wi[-m]))) # values of phi under H_0
      pars2 <- c(sigma1, phis) # estimates of sigma1, phi_1, ..., phi_m under H_0.
      score1 <- .score_algebraic(my.data, pars1, wi, vi, m) # score vector at xi+delta
      score2 <- .score_algebraic(my.data, pars2, wi, vi, m) # score vector at xi-delta
      score <- (score1 + score2) / 2 # average of the two scores
      info1 <- n.exc * .exp_info_algebraic(pars1, wi, vi, m) # expected information at xi+delta
      info2 <- n.exc * .exp_info_algebraic(pars2, wi, vi, m) # expected information at xi-delta
      e.info <- (info1 + info2) / 2 # average of the two information matrices
    }
    # .......... Start of test.stat.calc() ...........
    test.stat.calc <- function(score, info) {
      my.svd <- svd(info) # SVD of information matrix
      vec <- t(score) %*% my.svd$v # avoid matrix inversion
      stat <- vec %*% diag(1 / my.svd$d, nrow = m + 1) %*% t(vec) # score test statistic
      stat
    }
    # ........... End of test.stat.calc() ............

    e.stat <- test.stat.calc(score, e.info) # score test statistic
    e.p <- pchisq(e.stat, df = m - 1, lower.tail = FALSE) # p-value
    c(e.stat, e.p)
  }
  # .......................................................  end of function score.test

  ##########################################################################################
  #-------------------------- End of function score.test() --------------------------------#
  ##########################################################################################

  ##########################################################################################
  #----------------------- Start of function mult.u.gpd.fit() -----------------------------#
  ##########################################################################################

  mult.u.gpd.fit <- function(
    y,
    m,
    v,
    w,
    npy = 365,
    method = "Nelder-Mead",
    maxit = 10000,
    init.ests = NULL,
    ...
  ) {
    negated.mult.log.likelihood <- function(sig) {
      # sig: (sigma1, xi_1, ..., xi_m) y: excesses of threshold
      sigma1 <- sig[1] # sigma_1
      if (sigma1 <= 0) return(1e+30) # Need sigma_1 > 0
      xi <- sig[-1] # (xi_1, ..., xi_m)
      sigma <- sigma1 + cumsum(c(0, xi[-m] * w[-m])) # (sigma_1, ..., sigma_m)
      phi <- xi / sigma # (phi_1, ..., phi_m)
      if (any(1 + phi[-m] * w[-m] <= 0)) {
        return(1e+30) # Need all elements of 1+phi*w/sigma > 0
      }
      Ij <- unlist(lapply(y, function(sig) {
        sum(sig - v > 0)
      })) # interval indicators
      if (any(1 + phi[Ij] * (y - v[Ij]) <= 0)) return(1e+30) # Need all elements of 1+phi[Ij]*(y-v[Ij]) > 0
      aj <- c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m])) # -log(p_j), j=1, ..., m
      pj <- exp(-aj) # P(Y > v_j), j=1, ..., m
      bj <- log(sigma)
      dj <- log(1 + phi[Ij] * (y - v[Ij]))
      ej <- log(1 + phi[Ij] * (y - v[Ij])) / sigma[Ij] / phi[Ij]
      sum(aj[Ij] + bj[Ij] + dj + ej)
    }

    fscale <- negated.mult.log.likelihood(init.ests)
    temp <- optim(
      init.ests,
      negated.mult.log.likelihood,
      hessian = FALSE,
      method = method,
      control = list(maxit = maxit, fnscale = fscale, ...)
    )
    zz <- list()
    zz$mle <- temp$par
    zz$nllh <- temp$value
    zz$conv <- temp$convergence
    zz$counts <- temp$counts
    zz$message <- temp$message
    invisible(zz)
  }

  ##########################################################################################
  #------------------------- End of function mult.u.gpd.fit() -----------------------------#
  ##########################################################################################

  ##########################################################################################
  #--------------------- Start of function mult.thresh.LR.test() --------------------------#
  ##########################################################################################

  mult.thresh.LR.test <- function(
    my.data,
    m,
    vi,
    wi,
    init.ests = NULL,
    null.nllh = NULL
  ) {
    nllh2 <- mult.u.gpd.fit(my.data, m, vi, wi, init.ests = init.ests)$nllh
    LRT.stat <- ifelse(
      is.null(null.nllh),
      2 * (z$nllh[i] - nllh2),
      2 * (null.nllh - nllh2)
    )
    pvalue <- pchisq(LRT.stat, m - 1, lower.tail = FALSE)
    c(LRT.stat, pvalue)
  }

  ##########################################################################################
  #----------------------- End of function mult.thresh.LR.test() --------------------------#
  ##########################################################################################

  for (i in 1:(n_u - 1)) {
    # loop from u[1] to u[n_u-1]
    excess.data <- x[x > u[i]] - u[i] # excesses of threshold u[i]
    n.exc <- z$nexc[i] # number of excesses of current threshold
    # (calculated at top of function)
    m <- n_u - i + 1 # number of thresholds at or above u[i]
    z$df[i] <- m - 1 # df of null chi-squared distribution
    vi <- u[i:n_u] - u[i] # thresholds u[i], ..., u[n_u] relative to u[i]
    wi <- c(diff(vi), NA) # differences between thresholds (wi[m] not used)
    phis <- 1 / (z$sigma.mle[i] / z$xi.mle[i] + cumsum(c(0, wi[-m]))) # values of phi under H_0
    pars <- c(z$sigma.mle[i], phis) # estimates of sigma1, phi_1, ..., phi_m under H_0.
    # MLEs of sigma1 and xi1 under H_0, for passing to score.test() ...
    sigma1 <- z$sigma.mle[i]
    xi1 <- z$xi.mle[i]
    temp <- score.test(excess.data, m, vi, wi, pars, sigma1, xi1) # do score test
    z$e.test.stats[i] <- temp[1]
    z$e.p.values[i] <- temp[2]
    if (do.LRT) {
      null.init.ests <- c(sigma1, rep(xi1, m))
      temp <- mult.thresh.LR.test(
        excess.data,
        m,
        vi,
        wi,
        init.ests = null.init.ests
      ) # do LR test
      z$LRT.p.values[i] <- temp[2]
      z$LRT.test.stats[i] <- temp[1]
    }
  } #.......................# end of loop over thresholds
  z$u <- u[seq_len(n_u - 1)] # (lowest) thresholds for each test
  class(z) <- "mev_thdiag_northropcoleman"
  if (isTRUE(plot)) {
    plot(z, size = size, ...)
  }
  invisible(z)
}

#' @export
plot.mev_thdiag_northropcoleman <- function(x, size = 0.05, ...) {
  args <- list(...)
  n_u <- length(x$u) + 1L
  # Backward compatibility
  if (!is.null(args$my.xlab)) {
    args$xlab <- args$my.xlab
    args$my.xlab <- NULL
  }
  # Produce the plot ......
  if (is.null(args$ylab)) {
    args$ylab <- "p-value"
  }
  if (is.null(args$xlab)) {
    args$xlab <- "threshold"
  }
  if (is.null(args$type)) {
    args$type <- "b"
  }
  if (is.null(args$pch)) {
    args$pch <- 16
  }
  if (is.null(args$ylim)) {
    args$ylim <- c(0, 1)
  }
  args$x <- x$u
  args$y <- x$e.p.values
  do.call(what = "plot", args = args)
  axis(3, at = x$u, labels = x$nexc[-length(x$nexc)], cex.axis = 0.7)
  if (!is.null(x$LRT.p.values)) {
    lines(x$u, x$LRT.p.values, type = "b", lty = 4, pch = 2)
  }
  if (!is.null(size)) {
    if (isTRUE(all(is.numeric(size), size < 1, size > 0))) {
      abline(h = size, lty = 2)
    }
  }
  return(invisible(NULL))
}

#' @export
print.mev_thdiag_northropcoleman <-
  function(x, digits = min(3, getOption("digits") - 3), level = 0.05, ...) {
    cat("Threshold selection method: Northrop and Coleman penultimate model.\n")
    method <- "Score test"
    thselect <- x$thresh[which.max(which(x$e.p.values > level))]
    if (!is.null(x$LRT.p.values)) {
      method <- "Likelihood ratio test"
      thselect <- x$thresh[which.max(which(x$LRT.p.values > level))]
    }
    cat(method, "for piecewise generalized Pareto models.", "\n")
    cat(
      "Largest threshold above which we always fail to reject",
      "\n",
      "the null hypothesis of common generalized Pareto at level",
      level,
      "\n"
    )
    cat("Selected threshold:", round(thselect, digits), "\n")

    return(invisible(NULL))
  }
#------------------------------------------------------------------------------#
# Algebraic calculation of score vector #
#------------------------------------------------------------------------------#

# Algebraic score
#
# @param y vector of excesses of lowest threshold \code{u1}
# @param x parameter vector (\code{sigma1}, \code{phi_1}, \ldots, \code{phi_m})
# @param v thresholds relative to lowest threshold
# @param w differences between thresholds (\code{w[m]} not used)
# @param m number of thresholds
# @keywords internal
.score_algebraic <- function(y, x, w, v, m) {
  n <- length(y) # sample size
  sigma1 <- x[1] # sigma_1
  phi <- x[-1] # (phi_1, ..., phi_m)
  sigma <- sigma1 * c(1, cumprod(1 + phi[-m] * w[-m])) # (sigma_1, ..., sigma_m)
  xi <- sigma * phi # (xi_1, ..., xi_m)
  Ij <- unlist(lapply(y, function(x) {
    sum((x - v) > 0)
  })) # interval indicators
  h <- phi * sigma / sigma1 # (h_1, ..., h_m)
  ##################### Derivatives of bj ...
  db.s1 <- n / sigma1
  db.phi <- c(
    unlist(lapply(1:(m - 1), function(x) {
      sum(Ij > x) * w[x] / (1 + phi[x] * w[x])
    })),
    0
  )
  ##################### Derivatives of dj ...
  dd.s1 <- 0
  dd.phi <- unlist(lapply(1:m, function(x) {
    sum((y[Ij == x] - v[x]) / (1 + phi[x] * (y[Ij == x] - v[x])))
  }))
  ##################### Derivatives of ej ...
  de.s1 <- sum(-sigma1^(-2) * h[Ij]^(-1) * log(1 + phi[Ij] * (y - v[Ij])))
  de.phi.fn <- function(x) {
    temp1 <- h[x]^(-1) *
      (-phi[x]^(-1) *
        log(1 + phi[x] * (y[Ij == x] - v[x])) +
        (y[Ij == x] - v[x]) * (1 + phi[x] * (y[Ij == x] - v[x]))^(-1))
    if (x == m) {
      return(sigma1^(-1) * sum(temp1))
    }
    ind1 <- Ij %in% (x + 1):m
    ind2 <- Ij[ind1]
    temp2 <- -w[x] *
      (1 + phi[x] * w[x])^(-1) *
      (h[Ij][ind1]^(-1) * log(1 + phi[Ij][ind1] * (y[ind1] - v[Ij][ind1])))
    sigma1^(-1) * (sum(temp1) + sum(temp2))
  }
  de.phi <- unlist(lapply(1:m, de.phi.fn))
  ##################### Derivatives of aj ...
  aj <- c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m])) # -log(p_j), j=1, ..., m
  da.s1 <- -sigma1^(-1) * sum(aj[Ij])
  da.phi.fn <- function(x) {
    n.Ijs <- table(c(Ij, 1:m)) - 1 # frequencies of values of Ij
    n.Ijs <- n.Ijs[-1]
    temp1 <- h[x]^(-1) *
      (-phi[x]^(-1) * log(1 + phi[x] * w[x]) + w[x] * (1 + phi[x] * w[x])^(-1)) # B(x, x)
    my.zeros <- numeric(m - 1)
    my.zeros[x] <- temp1
    temp1 <- my.zeros
    n.Ijs <- as.numeric(n.Ijs)
    if (x == (m - 1)) {
      return(sigma1^(-1) * sum(temp1 * n.Ijs))
    }
    if (x < (m - 1)) {
      ind <- (x + 1):(m - 1)
      temp2 <- -w[x] *
        (1 + phi[x] * w[x])^(-1) *
        h[ind]^(-1) *
        log(1 + phi[ind] * w[ind])
    }
    temp1[(x + 1):(m - 1)] <- temp2
    sigma1^(-1) * sum(cumsum(temp1) * n.Ijs)
  }
  #
  da.phi <- c(unlist(lapply(1:(m - 1), da.phi.fn)), 0)
  ########### Return score vector ...
  -c(da.s1 + db.s1 + dd.s1 + de.s1, da.phi + db.phi + dd.phi + de.phi)
}

#------------------------------------------------------------------------------#
# Algebraic calculation of expected information #
#------------------------------------------------------------------------------#

# Algebraic calculation of the expected information
#
# @param x parameter vector: (\code{sigma1}, \code{phi_1}, \ldots,\code{phi_m})
# @param v thresholds relative to lowest threshold
# @param w differences between thresholds (\code{w[m]} not used)
# @param m number of thresholds
# @return a matrix with the expected information
# @keywords internal
.exp_info_algebraic <- function(x, w, v, m) {
  sigma1 <- x[1] # sigma_1
  phi <- x[-1] # (phi_1, ..., phi_m)
  sigma <- sigma1 * c(1, cumprod(1 + phi[-m] * w[-m])) # (sigma_1, ..., sigma_m)
  xi <- sigma * phi # (xi_1, ..., xi_m)
  aj <- c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m])) # -log(p_j), j=1, ..., m
  pj <- exp(-aj) # P(Y > v_j), j=1, ..., m
  qj <- c(pj[-m] * (1 - (1 + phi[-m] * w[-m])^(-1 / xi[-m])), pj[m]) # P(v_j < Y < v_{j+1}), j=1, ..., m
  h <- phi * sigma / sigma1 # (h_1, ..., h_m)
  # Various integrals ...
  I0b <- function(b, j) {
    t1 <- (1 + b * xi[j])^(-1)
    if (j == m) return(t1)
    t1 - (1 + phi[j] * w[j])^(-b - 1 / xi[j]) * t1
  }
  #
  I1b <- function(b, j) {
    t1 <- sigma[j] / (1 + (b - 1) * xi[j]) / (1 + b * xi[j])
    if (j == m) return(t1)
    t2 <- (1 + phi[j] * w[j])^(-b - 1 / xi[j]) / (1 + b * xi[j])
    t3 <- (1 + phi[j] * w[j])^(1 - b - 1 / xi[j]) / (1 + (b - 1) * xi[j])
    t1 + (t2 - t3) / phi[j]
  }
  #
  I2b <- function(b, j) {
    t1 <- 2 *
      sigma[j]^2 /
      (1 + (b - 2) * xi[j]) /
      (1 + (b - 1) * xi[j]) /
      (1 + b * xi[j])
    if (j == m) return(t1)
    t2 <- (1 + phi[j] * w[j])^(2 - b - 1 / xi[j]) / (1 + (b - 2) * xi[j])
    t3 <- (1 + phi[j] * w[j])^(1 - b - 1 / xi[j]) / (1 + (b - 1) * xi[j])
    t4 <- (1 + phi[j] * w[j])^(-b - 1 / xi[j]) / (1 + b * xi[j])
    t1 - (t2 - 2 * t3 + t4) / phi[j]^2
  }
  #
  J <- function(j) {
    t1 <- xi[j]
    if (j == m) return(t1)
    t2 <- -(1 + phi[j] * w[j])^(-1 / xi[j]) * log(1 + phi[j] * w[j])
    t3 <- -xi[j] * (1 + phi[j] * w[j])^(-1 / xi[j])
    t1 + t2 + t3
  }
  # Derivatives of bj (constant w.r.t. y) ...  Note: all contributions are zero apart from the diagonal elements for sigma1, phi_1,
  # ..., phi_{m-1}
  b.exp.info <- matrix(0, m + 1, m + 1) # matrix for expected information from b
  db.phi.phi <- c(
    unlist(lapply(1:(m - 1), function(x) {
      -w[x]^2 * (1 + phi[x] * w[x])^(-2) * sum(qj[(x + 1):m])
    })),
    0
  )
  diag(b.exp.info) <- c(-1 / sigma1^2, db.phi.phi)
  # Derivatives of dj ...  Note: all contributions are zero apart from the diagonal elements for phi_1, ..., phi_m.
  d.exp.info <- matrix(0, m + 1, m + 1) # matrix for expected information from b
  dd.phi.phi <- unlist(lapply(1:m, function(x) {
    -pj[x] * I2b(2, x)
  }))
  diag(d.exp.info) <- c(0, dd.phi.phi)
  ##################### Derivatives of ej ...
  e.exp.info <- matrix(0, m + 1, m + 1) # matrix for expected information from b
  de.s1.s1 <- 2 * sigma1^(-3) * sum(h^(-1) * pj * unlist(lapply(1:m, J)))
  #
  de.phi.fn <- function(x) {
    temp1 <- h[x]^(-1) * (-phi[x]^(-1) * J(x) + I1b(1, x))
    if (x == m) return(sigma1^(-1) * temp1 * pj[x])
    temp2 <- -w[x] *
      (1 + phi[x] * w[x])^(-1) *
      h[(x + 1):m]^(-1) *
      unlist(lapply((x + 1):m, J))
    sigma1^(-1) * (temp1 * pj[x] + sum(temp2 * pj[(x + 1):m]))
  }
  #
  de.phi <- unlist(lapply(1:m, de.phi.fn))
  #
  de.phi.phi.fn <- function(x) {
    temp1 <- h[x]^(-1) *
      (2 * phi[x]^(-2) * J(x) - 2 * phi[x]^(-1) * I1b(1, x) - I2b(2, x))
    if (x == m) return(sigma1^(-1) * temp1 * pj[x])
    temp2 <- 2 *
      w[x]^2 *
      (1 + phi[x] * w[x])^(-2) *
      h[(x + 1):m]^(-1) *
      unlist(lapply((x + 1):m, J))
    sigma1^(-1) * (temp1 * pj[x] + sum(temp2 * pj[(x + 1):m]))
  }
  #
  de.phi.phi <- unlist(lapply(1:m, de.phi.phi.fn))
  #
  de.phil.phik.fn <- function(x) {
    k <- x[1]
    l <- x[2]
    temp1 <- h[k]^(-1) *
      w[l] *
      (1 + w[l] * phi[l])^(-1) *
      (phi[k]^(-1) * J(k) - I1b(1, k))
    if (k == m) return(sigma1^(-1) * temp1 * pj[k])
    temp2 <- w[k] *
      (1 + phi[k] * w[k])^(-1) *
      w[l] *
      (1 + phi[l] * w[l])^(-1) *
      h[(k + 1):m]^(-1) *
      unlist(lapply((k + 1):m, J))
    sigma1^(-1) * (temp1 * pj[k] + sum(temp2 * pj[(k + 1):m]))
  }
  k.vals <- rep(2:m, times = 2:m - 1) # k > l
  l.vals <- unlist(lapply(2:m - 1, function(x) seq(from = 1, to = x)))
  kl.vals <- cbind(k.vals, l.vals)
  rev.kl.vals <- kl.vals[, 2:1]
  if (m == 2) {
    kl.vals <- matrix(kl.vals, nrow = 1, ncol = 2) # if m=2 make kl.vals a matrix
    rev.kl.vals <- matrix(kl.vals[, 2:1], nrow = 1, ncol = 2) # if m=2 make rev.kl.vals a matrix
  }
  de.phil.phik <- unlist(apply(kl.vals, 1, de.phil.phik.fn))
  #
  diag(e.exp.info) <- c(de.s1.s1, de.phi.phi)
  e.exp.info[1, 2:(m + 1)] <- e.exp.info[2:(m + 1), 1] <- -sigma1^(-1) * de.phi
  e.exp.info[kl.vals + 1] <- e.exp.info[rev.kl.vals + 1] <- de.phil.phik
  ##################### Derivatives of aj ...
  a.exp.info <- matrix(0, m + 1, m + 1) # matrix for expected information from b
  aj <- c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m])) # -log(p_j), j=1, ..., m
  da.s1.s1 <- 2 * sigma1^(-2) * sum(aj * qj)
  da.phi.fn <- function(x) {
    temp1 <- h[x]^(-1) *
      (-phi[x]^(-1) * log(1 + phi[x] * w[x]) + w[x] * (1 + phi[x] * w[x])^(-1)) # B(x, x)
    my.zeros <- numeric(m - 1)
    my.zeros[x] <- temp1
    temp1 <- my.zeros
    if (x == (m - 1)) return(sigma1^(-1) * sum(temp1 * qj[-1]))
    if (x < (m - 1)) {
      ind <- (x + 1):(m - 1)
      temp2 <- -w[x] *
        (1 + phi[x] * w[x])^(-1) *
        h[ind]^(-1) *
        log(1 + phi[ind] * w[ind])
    }
    temp1[(x + 1):(m - 1)] <- temp2
    sigma1^(-1) * sum(cumsum(temp1) * qj[-1])
  }
  #
  da.phi <- c(unlist(lapply(1:(m - 1), da.phi.fn)), 0)
  da.s1.phi <- -sigma1^(-1) * da.phi
  #
  da.phi.phi.fn <- function(x) {
    y.1v <- 1 + phi[x] * w[x]
    y.v <- w[x]
    temp1 <- h[x]^(-1) *
      (2 *
        phi[x]^(-2) *
        log(y.1v) -
        2 * phi[x]^(-1) * y.v * y.1v^(-1) -
        y.v^2 * y.1v^(-2))
    my.zeros <- numeric(m - 1)
    my.zeros[x] <- temp1
    temp1 <- my.zeros
    if (x == (m - 1)) return(sigma1^(-1) * sum(temp1 * qj[-1]))
    if (x < (m - 1)) {
      ind <- (x + 1):(m - 1)
      temp2 <- 2 *
        w[x]^2 *
        (1 + phi[x] * w[x])^(-2) *
        h[ind]^(-1) *
        log(1 + phi[ind] * w[ind])
    }
    temp1[(x + 1):(m - 1)] <- temp2
    sigma1^(-1) * sum(cumsum(temp1) * qj[-1])
  }
  #
  da.phi.phi <- c(unlist(lapply(1:(m - 1), da.phi.phi.fn)), 0)
  #
  da.phil.phik.fn <- function(x) {
    k <- x[1]
    l <- x[2]
    y.1v <- 1 + phi[k] * w[k]
    y.v <- w[k]
    temp1 <- h[k]^(-1) *
      w[l] *
      (1 + w[l] * phi[l])^(-1) *
      (phi[k]^(-1) * log(y.1v) - y.v * y.1v^(-1))
    my.zeros <- numeric(m - 1)
    my.zeros[k] <- temp1
    temp1 <- my.zeros
    if (k == (m - 1)) return(sigma1^(-1) * sum(temp1 * qj[-1]))
    if (k < (m - 1)) {
      ind <- (k + 1):(m - 1)
      temp2 <- w[k] *
        (1 + phi[k] * w[k])^(-1) *
        w[l] *
        (1 + phi[l] * w[l])^(-1) *
        h[ind]^(-1) *
        log(1 + phi[ind] * w[ind])
    }
    temp1[(k + 1):(m - 1)] <- temp2
    sigma1^(-1) * sum(cumsum(temp1) * qj[-1])
  }
  diag(a.exp.info) <- c(da.s1.s1, da.phi.phi)
  a.exp.info[1, 2:(m + 1)] <- a.exp.info[2:(m + 1), 1] <- da.s1.phi
  if (m > 2) {
    k.vals <- rep(2:(m - 1), times = 1:(m - 2)) # k > l
    l.vals <- unlist(lapply(1:(m - 2), function(x) {
      seq(from = 1, to = x)
    }))
    kl.vals <- cbind(k.vals, l.vals)
    rev.kl.vals <- kl.vals[, 2:1]
    if (m == 3) {
      kl.vals <- matrix(kl.vals, nrow = 1, ncol = 2) # if m=3 make kl.vals a matrix
      rev.kl.vals <- matrix(kl.vals[, 2:1], nrow = 1, ncol = 2) # if m=3 make rev.kl.vals a matrix
    }
    da.phil.phik <- unlist(apply(kl.vals, 1, da.phil.phik.fn))
    a.exp.info[kl.vals + 1] <- a.exp.info[rev.kl.vals + 1] <- da.phil.phik
  }
  # Return observed information ...
  exp.info <- a.exp.info + b.exp.info + d.exp.info + e.exp.info
  #
  exp.info
}

#' Northop and Coleman piecewise generalized Pareto threshold selection diagnostic
#'
#' The model tests the null hypothesis of a generalized Pareto above each threshold in \code{thresh} against the alternative of a piecewise generalized Pareto model with continuity constraints.
#' @param xdat [vector] observations
#' @param thresh [vector] candidate thresholds
#' @param test [string] indicating whether to perform \code{score} test or likelihood ratio (\code{lr}) test. The latter requires fitting the alternative model, and so is more computationally expensive.
#' @param plot [logical]; if \code{TRUE}, return a plot with the p-value path.
#' @param level [double] confidence level for confidence interval, defaults to 0.95
#' @param ... additional arguments, for backward compatibility purposes
#' @return an object of class \code{mev_thselect_ncpgp} containing the test statistic (\code{stat}), the p-values (\code{pval}), the threshold candidates (\code{thresh}) and the selected threshold (\code{thresh0}).
#' @export
thselect.ncpgp <- function(
  xdat,
  thresh,
  test = "score",
  plot = FALSE,
  level = 0.95,
  ...
) {
  test <- match.arg(
    test,
    choices = c("score", "lr"),
    several.ok = FALSE
  )
  ncdiag <- NC.diag(
    xdat = xdat,
    u = thresh,
    do.LRT = test == "lr",
    size = 1 - level,
    plot = FALSE
  )
  res <- list(
    thresh = ncdiag$thresh,
    shape = ncdiag$xi.mle,
    scale = ncdiag$sigma.mle,
    nexc = ncdiag$nexc,
    df = ncdiag$df,
    test = test,
    level = level
  )
  if (test == "score") {
    res$stat <- ncdiag$e.test.stats
    res$pval <- ncdiag$e.p.values
  } else {
    res$stat <- ncdiag$e.test.stats
    res$pval <- ncdiag$e.p.values
  }
  res$thresh0 <- ncdiag$thresh[which.max(which(res$pval > 1 - level))]
  if (length(res$thresh0) == 0L) {
    res$thresh0 <- tail(ncdiag$thresh, 1)
  }
  class(res) <- "mev_thselect_ncpgp"
  if (isTRUE(plot)) {
    plot(res, ...)
  }
  return(res)
}


#' @export
print.mev_thselect_ncpgp <-
  function(x, digits = min(3, getOption("digits") - 3), ...) {
    cat(
      "Threshold selection method: \nNorthrop and Coleman piecewise generalized Pareto model.\n"
    )
    method <- ifelse(x$test == "score", "score test", "likelihood ratio test")
    cat(method, "for piecewise generalized Pareto models.", "\n")
    cat(
      "Largest threshold above which we always fail to reject ",
      "\n",
      "the null hypothesis of common generalized Pareto at level",
      1 - x$level,
      "\n"
    )
    cat("Selected threshold:", round(x$thresh0, digits), "\n")

    return(invisible(NULL))
  }


#' @export
plot.mev_thselect_ncpgp <- function(x, ...) {
  args <- list(...)
  # Produce the plot ......
  if (is.null(args$ylab)) {
    args$ylab <- "p-value"
  }
  if (is.null(args$bty)) {
    args$bty <- "l"
  }
  if (is.null(args$xlab)) {
    args$xlab <- "threshold"
  }
  if (is.null(args$type)) {
    args$type <- "b"
  }
  if (is.null(args$pch)) {
    args$pch <- 20
  }
  if (is.null(args$ylim)) {
    args$ylim <- c(0, 1)
  }
  args$x <- x$thresh[-length(x$thresh)]
  args$y <- x$pval
  do.call(what = "plot", args = args)
  axis(
    side = 3,
    at = x$thresh[-length(x$thresh)],
    labels = x$nexc[-length(x$nexc)],
    cex.axis = 0.7
  )
  mtext(
    switch(x$test, lr = "likelihood ratio", score = "score"),
    side = 1,
    line = 2,
    adj = 1,
    outer = FALSE
  )
  abline(h = 1 - x$level, lty = 2)
  return(invisible(NULL))
}
