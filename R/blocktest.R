#' Generalized extreme value model by block
#'
#' Given a sample of block maxima of dimension \code{n} by \code{m}, split further into groups of size \code{m} (number of columns), fit a generalized extreme value distribution allowing the maximum of each subblock of size \code{m} to have different parameters under the alternative. Likelihood for left-censored and rounded (interval-censored) observations is
#' also supported.
#'
#' @param pars vector of location, scale and shape parameters for the GEV
#' @param xdat \code{n} by \code{m} matrix of observations assumed to arise from a GEV, ordered by row
#' @param null logical; if \code{TRUE}, the parameter vector for the largest observation of the block is obtained from the max stability relationship
#' @param alternative integer; \code{1} for a single shape parameter with parameters obtained from max-stability, \code{2} for a common shape parameter, but with free location and scale, and \code{3} for different parameters (location, scale and shape).
#' @param lb lower bound for left-censoring; default to none (\code{NULL})
#' @param rounding double, a positive number indicating the amount of censoring (e.g., 0.1 or 1)
#' @param logscale logical; if \code{TRUE}, scale parameter is provided on the log scale
#' @return log likelihood of sample
#' @keywords internal
#' @export
#' @examples
#' samp <- build.blocks(mev::rgev(50, scale = 10), m = 4)
#' pars <- c(0,10,0)
#' gevblock.ll(pars, xdat = round(samp,0), null = TRUE, rounding = 1, lb = -5)
#' gevblock.ll(pars, xdat = round(samp,0), null = TRUE, rounding = 1)
#' gevblock.ll(pars, xdat = samp, null = TRUE)
#' gevblock.ll(pars, xdat = samp, lb = -5, null = TRUE)
gevblock.ll <- function(
  pars,
  xdat,
  null = FALSE,
  lb = NULL,
  rounding = 0,
  alternative = 1:3,
  logscale = FALSE
) {
  if (!is.null(lb)) {
    stopifnot(length(lb) == 1L)
  }
  alternative <- as.integer(alternative[1])
  if (alternative == 4L) {
    # back compatibility
    alternative <- 1L
  }
  stopifnot(
    isTRUE(alternative %in% 1:3)
  )
  xdat <- as.matrix(xdat)
  n <- nrow(xdat)
  m <- ncol(xdat)
  if (isTRUE(logscale)) {
    pars[2] <- exp(pars[2])
  }
  if (isTRUE(null)) {
    pars2 <- mev::maxstable(pars[1:3], m = m)
  } else {
    if (alternative == 1L) {
      # Consider alternative with m * lambda and test for lambda=1
      pars2 <- mev::maxstable(pars[1:3], m = m * exp(pars[4]))
    } else if (alternative == 2L) {
      # shape is fixed, location varies
      pars2 <- pars[c(4:5, 3)]
    } else if (alternative == 3L) {
      # all parameters vary
      pars2 <- pars[c(4:6)]
    }
    if (isTRUE(logscale) & (alternative %in% 2:3)) {
      pars2[2] <- exp(pars2[2])
    }
  }
  # Sort row from smallest to largest
  xdat <- t(apply(xdat, 1, sort))
  delta <- abs(rounding[1]) / 2 #HALF ROUNDING since add/subtract
  rounded <- isTRUE(delta > 1e-6)
  leftcens <- !is.null(lb)
  if (leftcens & !rounded) {
    # Matrix of indicators for left-censoring
    lcens <- apply(xdat, 1:2, function(x) {
      x < lb
    })
    # GEV for smaller blocks, truncated above by block max
    out <- sum(mev::dgev(
      x = as.numeric(xdat[, -m][!lcens[, -m]]),
      loc = pars[1],
      scale = pars[2],
      shape = pars[3],
      log = TRUE
    )) +
      sum(lcens[, -m]) *
        mev::pgev(
          q = lb,
          loc = pars[1],
          scale = pars[2],
          shape = pars[3],
          log.p = TRUE
        ) -
      (m - 1) *
        sum(mev::pgev(
          q = as.numeric(pmax(lb, xdat[, m])),
          loc = pars[1],
          scale = pars[2],
          shape = pars[3],
          log = TRUE
        )) +
      # Null hypothesis model: maximum is GEV
      sum(mev::dgev(
        x = as.numeric(xdat[, m][!lcens[, m]]),
        loc = pars2[1],
        scale = pars2[2],
        shape = pars2[3],
        log = TRUE
      )) +
      sum(lcens[, m]) *
        mev::pgev(
          q = lb,
          loc = pars2[1],
          scale = pars2[2],
          shape = pars2[3],
          log.p = TRUE
        )
  } else if (leftcens & rounded) {
    # browser()
    weight_fn <- function(x, pars) {
      w <- rep(0, length(x))
      w[((x - delta) < lb) & ((x + delta) >= lb)]
      w[x - delta > lb] <- 1
      indet <- which(((x - delta) < lb) & ((x + delta) >= lb))
      if (length(indet) > 0) {
        plb <- pgev(lb, pars[1], pars[2], pars[3])
        pu <- pgev(
          x[indet] + delta,
          loc = pars[1],
          scale = pars[2],
          shape = pars[3]
        )
        pl <- pgev(
          x[indet] - delta,
          loc = pars[1],
          scale = pars[2],
          shape = pars[3]
        )
        w[indet] <- (pu - plb) / (pu - pl)
      }
      return(w)
    }
    wmat <- cbind(
      apply(xdat[, -m, drop = FALSE], 2, weight_fn, pars = pars),
      weight_fn(xdat[, m], pars = pars2)
    )
    # GEV for smaller blocks, truncated above by block max
    out <-
      sum(1 - wmat[, -m]) *
      mev::pgev(
        q = lb,
        loc = pars[1],
        scale = pars[2],
        shape = pars[3],
        log.p = TRUE
      ) +
      sum(
        c(wmat[, -m]) *
          log(pmax(
            1e-16, # to handle zero cases
            mev::pgev(
              q = pmax(lb, c(xdat[, -m]) + delta),
              loc = pars[1],
              scale = pars[2],
              shape = pars[3]
            ) -
              mev::pgev(
                q = pmax(lb, c(xdat[, -m]) - delta),
                loc = pars[1],
                scale = pars[2],
                shape = pars[3]
              )
          ))
      )
    -(m - 1) *
      sum(mev::pgev(
        q = pmax(lb, as.numeric(xdat[, m]) + delta),
        loc = pars[1],
        scale = pars[2],
        shape = pars[3],
        log = TRUE
      )) +
      # Null hypothesis model: maximum is GEV
      sum(1 - wmat[, m]) *
        mev::pgev(
          q = lb,
          loc = pars2[1],
          scale = pars2[2],
          shape = pars2[3],
          log.p = TRUE
        ) +
      sum(
        wmat[, m] *
          log(pmax(
            1e-16,
            mev::pgev(
              q = as.numeric(xdat[, m]) + delta,
              loc = pars2[1],
              scale = pars2[2],
              shape = pars2[3]
            ) -
              mev::pgev(
                q = pmax(lb, as.numeric(xdat[, m]) - delta),
                loc = pars2[1],
                scale = pars2[2],
                shape = pars2[3]
              )
          ))
      )
  } else if (rounded & !leftcens) {
    out <- sum(log(
      mev::pgev(
        q = as.numeric(xdat[, -m]) + delta,
        loc = pars[1],
        scale = pars[2],
        shape = pars[3]
      ) -
        mev::pgev(
          q = as.numeric(xdat[, -m]) - delta,
          loc = pars[1],
          scale = pars[2],
          shape = pars[3]
        )
    )) -
      (m - 1) *
        sum(mev::pgev(
          q = as.numeric(xdat[, m]) + delta,
          loc = pars[1],
          scale = pars[2],
          shape = pars[3],
          log = TRUE
        )) +
      # Null hypothesis model: maximum is GEV
      sum(log(
        mev::pgev(
          q = as.numeric(xdat[, m]) + delta,
          loc = pars2[1],
          scale = pars2[2],
          shape = pars2[3]
        ) -
          mev::pgev(
            q = as.numeric(xdat[, m]) - delta,
            loc = pars2[1],
            scale = pars2[2],
            shape = pars2[3]
          )
      ))
  } else if (!rounded & !leftcens) {
    # GEV for smaller blocks, truncated above by block max
    out <- sum(mev::dgev(
      x = as.numeric(xdat[, -m]),
      loc = pars[1],
      scale = pars[2],
      shape = pars[3],
      log = TRUE
    )) -
      (m - 1) *
        sum(mev::pgev(
          q = as.numeric(xdat[, m]),
          loc = pars[1],
          scale = pars[2],
          shape = pars[3],
          log = TRUE
        )) +
      # Null hypothesis model: maximum is GEV
      sum(mev::dgev(
        x = as.numeric(xdat[, m]),
        loc = pars2[1],
        scale = pars2[2],
        shape = pars2[3],
        log = TRUE
      ))
  }
  return(as.numeric(out))
}


#' Likelihood ratio test for max-stability
#'
#' Given a matrix of block maxima split into blocks of size \code{m},
#' calculate test statistics and return p-values based on the asymptotic
#' chi-square distribution.
#' @param xdat \code{n} by \code{m} matrix of observations assumed to arise from a GEV, ordered by row
#' @param alternative integer; \code{1} for a single shape parameter with parameters obtained from max-stability, \code{2} for a common shape parameter, but with free location and scale.
#' @param lb lower bound for left-censoring; default to none (\code{NULL})
#' @param rounding double, a positive number indicating the amount of censoring (e.g., 0.1 or 1)
#' @return a data frame containing likelihood ratio statistics (\code{stat}), the degrees of freedom, a vector of p-values \code{pval} and the name of the \code{alternative}.
#' @export
#' @examples
#' samp <- build.blocks(mev::rgev(50, scale = 10), m = 4)
#' test.blocksize(xdat = round(samp, 0), rounding = 1, lb = -5)
#' test.blocksize(xdat = round(samp, 0), rounding = 1)
#' test.blocksize(xdat = samp)
#' test.blocksize(xdat = samp, lb = -5, alternative = 1L)
test.blocksize <- function(
  xdat,
  rounding = 0,
  alternative = c(1L, 2L, 3L),
  lb = NULL
) {
  m <- ncol(xdat)
  # Get starting values from the null model fit using standard routines
  fitgev <- mev::fit.gev(xdat = as.numeric(xdat))
  start <- coef(fitgev)
  # Log likelihood of null model
  null <- gevblock.ll(
    pars = start,
    xdat = xdat,
    null = TRUE,
    rounding = rounding,
    lb = lb,
    logscale = FALSE
  )
  if (isTRUE(any(2:3 %in% alternative))) {
    fitgev_max <- mev::fit.gev(xdat = xdat[, m])
  }
  alt <- alternative[alternative %in% 1:3]
  if (length(alt) == 0) {
    stop("Invalid vector of alternative: must be integer between 1 and 3.")
  }
  maxstab <- mev::maxstable(pars = start, m = m)
  # Alternative where only shape varies
  alt_stat <- numeric(length(alternative))
  df <- numeric(length(alternative))
  for (i in seq_along(alternative)) {
    if (alternative[i] == 1L) {
      alt_opt1 <- optim(
        par = c(start, 0),
        fn = gevblock.ll,
        method = "Nelder",
        control = list(fnscale = -1, maxit = 2e4L),
        alternative = 1L,
        lb = lb,
        rounding = rounding,
        xdat = xdat
      )
      alt_stat[i] <- alt_opt1$value
      df[i] <- 1
    } else if (alt[i] == 2L) {
      # Alternative 2: fix shape, allow location and scale to vary
      alt_opt2 <- optim(
        par = c(start, coef(fitgev_max)[1:2]),
        fn = gevblock.ll,
        method = "Nelder",
        control = list(fnscale = -1, maxit = 2e4L),
        alternative = 2L,
        lb = lb,
        rounding = rounding,
        xdat = xdat
      )
      alt_stat[i] <- alt_opt2$value
      df[i] <- 2
    } else if (alt[i] == 3L) {
      # Alternative 3: fit model, but allowing all three parameters to vary
      alt_opt3 <- optim(
        par = c(start, coef(fitgev_max)),
        fn = gevblock.ll,
        method = "Nelder",
        control = list(fnscale = -1, maxit = 2e4L),
        alternative = 3L,
        lb = lb,
        rounding = rounding,
        xdat = xdat
      )
      alt_stat[i] <- alt_opt3$value
      df[i] <- 3
    } else {
      stop("Invalid value for alternative \"alt\".")
    }
  }
  # P-values
  stat <- 2 * (alt_stat - null)
  pvals <- pchisq(stat, df = df, lower.tail = FALSE)
  data.frame(alternative = paste0("A", alt), stat = stat, df = df, pval = pvals)
  # names(pvals) <- paste0("A", alt)
  # return(pvals)
}


#' Compute block maxima and order them by block
#'
#' Given a time series of observations in \code{xdat},
#' compute the maximum of blocks of size \code{block} (\eqn{b}),
#' and then order them by further blocks of size \code{m},
#' increasing by row from left to right. If the length of
#' \code{xdat} is not a multiple of \code{block},
#' the last observations are discarded without warning.
#' @export
#' @param xdat vector of length \code{n}
#' @param block integer, size of block over which to compute maxima
#' @param m number of columns for further sub-blocking
#' @return a matrix with \eqn{\lfloor n/b \rfloor} observations, ordered by row, with \code{m} columns.
build.blocks <- function(xdat, block = 1L, m = 2L) {
  data <- as.numeric(c(xdat))
  data <- data[is.finite(data)]
  block <- as.integer(block)
  nb <- floor(length(data) / (m * block))
  n <- nb * block * m
  # stopifnot(nb %% m == 0)
  if (block > 1) {
    out <- t(apply(
      matrix(apply(matrix(data[1:n], nrow = block), 2, max), nrow = m),
      2,
      sort
    ))
  } else {
    out <- t(apply(matrix(data[1:n], nrow = m), 2, sort))
  }
  if (m == 1L) {
    return(c(out))
  } else {
    return(out)
  }
}


#' Durbin's mapping of uniforms
#'
#' Given a set of uniform variates, obtained for example through application
#' of the probability integral  transform, consider alternative scaling which
#' are also uniform but may lead to more powerful tests.
#' @param xdat sample of uniform variates
#' @return a data frame with elements \code{stat} containing the test statistics and \code{pval} for the p-values from the null distributions, for the modified probability product, median and (one-sided) Kolmogorov--Smirnov tests.
#' @references Durbin, J. (1961). Some Methods of Constructing Exact Tests. \emph{Biometrika}, 48(\bold{1/2}), 41–55. \doi{10.2307/2333128}
#' @export
#' @keywords internal
durbin.unif <- function(xdat) {
  u_vec <- sort(c(xdat))
  n <- length(u_vec)
  stopifnot(u_vec[1] >= 0, u_vec[n] <= 1)
  c_vec <- sort(c(u_vec[1], diff(u_vec), 1 - u_vec[n]))
  g_vec <- (n + 2 - 1:(n + 1)) * diff(c(0, c_vec))
  # tinytest::expect_equal(sum(c_vec), 1)
  # tinytest::expect_equal(sum(g_vec), 1)
  w_vec <- cumsum(g_vec[-(n + 1)])
  # Modified probability product test
  stat_1 <- -2 * sum(log(w_vec))
  # More affected by rounding errors.
  pval_1 <- pchisq(q = stat_1, df = 2 * n, lower.tail = FALSE)
  # Modified median test
  r <- ifelse(n %% 2 == 0, round(n / 2), ceiling(n / 2))
  stat_2 <- r / (n + 1 - r) * (1 - w_vec[r]) / w_vec[r]
  pval_2 <- pf(
    q = stat_2,
    df1 = 2 * (n + 1 - r),
    df2 = 2 * r,
    lower.tail = FALSE
  )
  stat_3 <- ks.test(x = w_vec, y = "punif", alternative = "greater")
  # Formula in Durbin is max((1:n) / n - w_vec)
  return(data.frame(
    stat = c(stat_1, stat_2, stat_3$statistic),
    pval = c(pval_1, pval_2, stat_3$p.value)
  ))
}


#' Durbin's mapping of uniforms
#'
#' Given a set of uniform variates, obtained for example through application
#' of the probability integral  transform, consider alternative scaling which
#' are also uniform but may lead to more powerful tests.
#' @param xdat sample of uniform variates
#' @return a data frame with elements
#' \itemize{
#' \item \code{test} string indicating which of the modified probability product, median and (one-sided) Kolmogorov--Smirnov tests
#' \item \code{stat} the test statistics
#' \item \code{pval} for the p-values from the null distributions
#' }
#' @references Durbin, J. (1961). Some Methods of Constructing Exact Tests. \emph{Biometrika}, 48(\bold{1/2}), 41–55. \doi{10.2307/2333128}
#' @keywords internal
#' @export
#' @examples
#' test.unif(runif(1000))
#' test.unif(rbeta(1000, 1, 1.5))
test.unif <- function(xdat) {
  u_vec <- sort(c(xdat))
  n <- length(u_vec)
  stopifnot(u_vec[1] >= 0, u_vec[n] <= 1)
  c_vec <- sort(c(u_vec[1], diff(u_vec), 1 - u_vec[n]))
  g_vec <- (n + 2 - 1:(n + 1)) * diff(c(0, c_vec))
  # tinytest::expect_equal(sum(c_vec), 1)
  # tinytest::expect_equal(sum(g_vec), 1)
  w_vec <- cumsum(g_vec[-(n + 1)])
  # Modified probability product test
  stat_1 <- -2 * sum(log(w_vec))
  # More affected by rounding errors.
  pval_1 <- pchisq(q = stat_1, df = 2 * n, lower.tail = FALSE)
  # Modified median test
  r <- ifelse(n %% 2 == 0, round(n / 2), ceiling(n / 2))
  stat_2 <- r / (n + 1 - r) * (1 - w_vec[r]) / w_vec[r]
  pval_2 <- pf(
    q = stat_2,
    df1 = 2 * (n + 1 - r),
    df2 = 2 * r,
    lower.tail = FALSE
  )
  stat_3 <- ks.test(x = w_vec, y = "punif", alternative = "greater")
  # Formula in Durbin is max((1:n) / n - w_vec)
  return(data.frame(
    test = c("probability product", "median", "Kolmogorov-Smirnov"),
    stat = c(stat_1, stat_2, stat_3$statistic),
    pval = c(pval_1, pval_2, stat_3$p.value)
  ))
}


#' Optimization for the marginal likelihood for a sample of GEV order statistics
#'
#' Given a matrix of \code{n} ordered samples of \code{m} order statistics from a postulated GEV, fit the parameters of the latter based on the marginal likelihood of the first
#' \code{m-1} order statistics using maximum likelihood.
#'
#' Additionally, one can set \code{constraint} to \code{TRUE} to add a support constraint
#' to the optimization to ensure that all values of \code{xdat} are in the support of the resulting distribution.
#'
#' @param xdat matrix of observations of size \code{n} by \code{m}, ordered by rows
#' @param constraint logical; if \code{TRUE}, add support constraint
#' @param rounding double; indicate the amount of rounding around value; default to zero
#' @param start vector of length 3 for starting values for GEV; default to \code{NULL}
#' @param vcov logical; if \code{TRUE}, return as attribute the estimate of the covariance matrix of the parameters given by the inverse observed information matrix.
#' @param lb lower bound; any point below \code{lb} is left-censored
#' @return (constrained) maximum likelihood estimator of location, scale and shape parameters
#' @keywords internal
fit.gevblock.marginal <- function(
  xdat,
  constraint = TRUE,
  rounding = 0,
  lb = NULL,
  start = NULL,
  vcov = FALSE
) {
  m <- ncol(xdat)
  maxx <- max(xdat[, m])
  rounding <- rounding[1] / 2
  if (!is.null(lb)) {
    lcens <- xdat[, -m] < lb
    # We consider the actual measurement as lower bound
  }
  if (abs(rounding) < 1e-6) {
    if (is.null(lb)) {
      gev_nll_os <- function(pars, xdat, rounding, lb, ...) {
        m <- ncol(xdat)
        -sum(mev::dgev(
          x = c(xdat[, -m]),
          loc = pars[1],
          scale = pars[2],
          shape = pars[3],
          log = TRUE
        )) -
          sum(mev::pgev(
            q = c(xdat[, m - 1]),
            loc = pars[1],
            scale = pars[2],
            shape = pars[3],
            lower.tail = FALSE,
            log.p = TRUE
          ))
      }
    } else {
      gev_nll_os <- function(pars, xdat, rounding, lb, ...) {
        m <- ncol(xdat)
        -sum(mev::dgev(
          x = c(xdat[, -m])[!c(lcens)],
          loc = pars[1],
          scale = pars[2],
          shape = pars[3],
          log = TRUE
        )) +
          sum(lcens) *
            mev::pgev(
              q = lb,
              loc = pars[1],
              scale = pars[2],
              shape = pars[3],
              log.p = TRUE
            ) +
          -sum(mev::pgev(
            q = c(xdat[, m - 1])[!lcens[, m - 1]],
            loc = pars[1],
            scale = pars[2],
            shape = pars[3],
            lower.tail = FALSE,
            log.p = TRUE
          ))
      }
    }
  } else {
    if (is.null(lb)) {
      gev_nll_os <- function(pars, xdat, rounding, lb, ...) {
        m <- ncol(xdat)
        -sum(
          log(
            mev::pgev(
              q = c(xdat[, -m]) + rounding,
              loc = pars[1],
              scale = pars[2],
              shape = pars[3]
            ) -
              mev::pgev(
                q = c(xdat[, -m]) - rounding,
                loc = pars[1],
                scale = pars[2],
                shape = pars[3]
              )
          )
        ) -
          sum(mev::pgev(
            q = c(xdat[, m - 1]) - rounding,
            loc = pars[1],
            scale = pars[2],
            shape = pars[3],
            lower.tail = FALSE,
            log.p = TRUE
          ))
      }
    } else {
      #rounding and truncation
      gev_nll_os <- function(pars, xdat, rounding, lb, ...) {
        m <- ncol(xdat)
        -sum(
          log(
            mev::pgev(
              q = pmax(lb, c(xdat[, -m]) + rounding),
              loc = pars[1],
              scale = pars[2],
              shape = pars[3]
            ) -
              ifelse(
                c(lcens[, -m]),
                0,
                mev::pgev(
                  q = as.numeric(xdat[, -m]) - rounding,
                  loc = pars[1],
                  scale = pars[2],
                  shape = pars[3]
                )
              )
          )
        ) - # The largest only contributes if the (m-1) order statistic is not left-censored
          sum(mev::pgev(
            q = c(xdat[, m - 1] - rounding)[xdat[, m - 1] > lb],
            loc = pars[1],
            scale = pars[2],
            shape = pars[3],
            lower.tail = FALSE,
            log.p = TRUE
          ))
      }
    }
  }
  if (is.null(start)) {
    start <- coef(mev::fit.gev(c(xdat)))
  }
  if (isTRUE(constraint)) {
    obj_fn <- function(xp, xdat, maxx, rounding, lb, ...) {
      out <- try(
        gev_nll_os(
          pars = c(xp[1], exp(xp[2]), xp[3]),
          xdat = xdat,
          rounding = rounding,
          lb = lb
        ),
        silent = TRUE
      )
      if (inherits(out, "try-error")) {
        out <- 1e10
      } else {
        if (!is.finite(out)) {
          out <- 1e10
        }
      }
      return(out)
    }
    start[2] <- log(start[2])
    nll_st <- obj_fn(
      start,
      xdat = xdat,
      maxx = maxx,
      rounding = rounding,
      lb = lb
    )
    if (isTRUE(nll_st == 1e10)) {
      start <- coef(mev::fit.gev(c(xdat)))
      start[2] <- log(start[2])
    }
    opt <- alabama::auglag(
      par = start,
      control.outer = list(trace = 0),
      fn = obj_fn,
      hin = function(xp, xdat, maxx, rounding, lb, ...) {
        c(exp(xp[2]) + xp[3] * (maxx + rounding - xp[1]), xp[3] + 1)
      },
      xdat = xdat,
      maxx = maxx,
      rounding = rounding,
      lb = lb
    )
    #if (opt$kkt1 & opt$kkt2) {
    if (opt$convergence == 0 | opt$value == 1e10) {
      out <- opt$par
      out[2] <- exp(out[2])
      # because the gradient isn't zero
      # so the approximation should be centered at opt - gradient *
      # not evaluated at the mode where gradient is zero
      if (isTRUE(vcov)) {
        opt_fn2 <- function(pars) {
          gev_nll_os(
            pars = pars,
            xdat = xdat,
            rounding = rounding,
            lb = lb
          )
        }
        grad_est <- try(
          numDeriv::grad(func = opt_fn2, x = out),
          silent = TRUE
        )
        vcov_est <- try(
          solve(numDeriv::hessian(func = opt_fn2, x = out)),
          silent = TRUE
        )
        if (!inherits(vcov_est, "try-error")) {
          attr(out, "mean") <- c(out - vcov_est %*% grad_est)
          attr(out, "vcov") <- vcov_est
        }
      }
      return(out)
    } else {
      stop("Optimization did not converge")
    }
  } else {
    opt <- optim(
      par = start,
      fn = gev_nll_os,
      xdat = xdat,
      rounding = rounding,
      hessian = isTRUE(vcov)
    )
    if (isTRUE(opt$convergence == 0)) {
      out <- opt$par
      if (isTRUE(vcov)) {
        vcov_est <- try(solve(opt$hessian), silent = TRUE)
        if (!inherits(vcov_est, "try-error")) {
          attr(out, "vcov") <- vcov_est
        }
      }
      return(out)
    } else {
      stop("Optimization did not converge")
    }
  }
}

#' Generate exponential
#'
#' Given integers \code{n} and \code{m}, simulate a matrix of
#' independent unit exponential order statistics,
#' but ordered by row from largest to smallest
#' @param n number of rows
#' @param m number of columns
#' @return an \code{n} by \code{m} matrix of order statistics from unit exponential
#' @keywords internal
generate.exponential <- function(n, m) {
  n <- as.integer(n)
  m <- as.integer(m)
  stopifnot(m >= 2L, n >= 1)
  # esamp <- matrix(rexp(n * m), nrow = n, ncol = m)
  # xdat <- matrix(NA, nrow = n, ncol = m)
  # xdat[, m] <- rexp(n = n, rate = m)
  # for (j in (m - 1):1) {
  #   xdat[, j] <- xdat[, j + 1] + rexp(n, rate = j)
  # }
  # return(xdat)
  t(apply(matrix(rexp(n = n * m), ncol = m), 1, sort, decreasing = TRUE))
}

#' Transformation of matrix from exponential to GEV order statistics
#' Given a matrix of exponential draws (ordered from largest to smallest),
#' apply a transformation to map them to the generalized extreme value distribution (GEV)
#' where now the rows have increasing order statistics from left to right
#' @param edat matrix of exponential order statistics (ordered by row from largest to smallest)
#' @param par vector of location, scale and shape parameter of the GEV
#' @return a matrix of the same dimension as \code{edat}
#' @keywords internal
ordexp2gev <- function(edat, par) {
  stopifnot(length(par) == 3L)
  if (abs(par[3]) > 1e-8) {
    par[1] + par[2] * (edat^(-par[3]) - 1) / par[3]
  } else {
    par[1] - par[2] * log(edat)
  }
}


#' Transform order statistics of GEV to exponential
#' @param xdat matrix of observations (ordered by row) from GEV
#' @param par vector of location, scale and shape parameters of the GEV
#' @return a matrix of the same dimension as \code{xdat} containing the
#' ordered unit exponential quantiles
#' @keywords internal
gev2exp <- function(xdat, par) {
  stopifnot(length(par) == 3L)
  m <- ncol(xdat)
  if (abs(par[3]) > 1e-8) {
    txdat <- (1 + par[3] * (xdat - par[1]) / par[2])^(-1 / par[3])
  } else {
    txdat <- exp(-(xdat - par[1]) / par[2])
  }
  texp <- matrix(nrow = nrow(xdat), ncol = ncol(xdat)) # Copy the same length
  texp[, 1] <- txdat[, m] * m
  for (j in 1:(m - 1)) {
    texp[, j + 1] <- (m - j) * (txdat[, m - j] - txdat[, m - j + 1])
  }
  return(texp)
}


#' Optimization for the GEV likelihood for blocks
#'
#' Given a matrix of \code{n} ordered samples of \code{m} order statistics
#' from a postulated GEV, fit the parameters of the latter based on the
#'  marginal likelihood of the first \code{m-1} order statistics using
#'  maximum likelihood.
#'
#' One can set \code{constraint} to \code{TRUE} to add a support constraint
#' to the optimization to ensure that all values of \code{xdat} are in the
#' support of the resulting distribution (only for the marginal likelihood).
#'
#' @param xdat matrix of observations of size \code{n} by \code{m}, ordered by rows
#' @param constraint logical; if \code{TRUE}, add support constraint
#' @param marginal logical; if \code{TRUE}, use marginal likelihood of lower order statistics
#' @param rounding double; indicate the amount of rounding around value; default to zero
#' @param start vector of length 3 for starting values for GEV; default to \code{NULL}
#' @param lb lower bound; any point below \code{lb} is left-censored
#' @param vcov logical; if \code{TRUE}, return as attribute the estimate of the covariance matrix of the parameters given by the inverse observed information matrix.
#' @return (constrained) maximum likelihood estimator of location, scale and shape parameters
#' @export
#' @examples
#' set.seed(2026)
#' xdat <- build.blocks(mev::rgev(n = 200, shape = 0.1), m = 4)
#' fit.gevblock(xdat, marginal = TRUE)
#' fit.gevblock(round(xdat, 1), marginal = TRUE, lb = NULL, rounding = 0.1)
#' fit.gevblock(round(xdat, 1), marginal = TRUE, lb = -2, rounding = 0.1)
#' fit.gevblock(xdat, marginal = TRUE, lb = -2)
#' fit.gevblock(xdat)
#' fit.gevblock(round(xdat, 1), lb = NULL, rounding = 0.1)
#' fit.gevblock(round(xdat, 1), lb = -2, rounding = 0.1)
#' fit.gevblock(xdat, lb = -2)
fit.gevblock <- function(
  xdat,
  marginal = FALSE,
  constraint = TRUE,
  rounding = 0,
  lb = NULL,
  start = NULL,
  vcov = FALSE
) {
  if (is.null(start)) {
    start <- coef(mev::fit.gev(c(xdat)))
  }
  if (isTRUE(marginal)) {
    stopifnot(is.matrix(xdat), ncol(xdat) >= 1)
    return(fit.gevblock.marginal(
      xdat = xdat,
      constraint = constraint,
      rounding = rounding,
      start = start,
      vcov = vcov
    ))
  } else {
    opt <- optim(
      par = start,
      fn = function(pars, xdat, lb, rounding, ...) {
        out <- gevblock.ll(
          pars = pars,
          xdat = xdat,
          null = TRUE,
          lb = lb,
          rounding = rounding
        )
        if (!is.finite(out)) {
          out <- -1e10
        }
        return(out)
      },
      rounding = rounding,
      lb = lb,
      null = TRUE,
      lower = c(-Inf, 1e-8, -1 + 1e-8),
      upper = c(rep(Inf, 2), 10),
      control = list(fnscale = -1),
      method = "L-BFGS",
      hessian = isTRUE(vcov),
      xdat = xdat,
    )
    out <- opt$par
    if (isTRUE(vcov)) {
      vcov_est <- try(solve(-opt$hessian), silent = TRUE)
      if (!inherits(vcov_est, "try-error")) {
        attr(out, "vcov") <- vcov_est
      }
    }
    return(out)
  }
}


#' @export
plot.mev_plot_blocksize <- function(
  x,
  type = c("max", "range", "all"),
  scale = "unif",
  ...
) {
  type <- unique(match.arg(type, several.ok = TRUE))
  type_plot <- lapply(x$plots, function(x) {
    x$type
  })
  lcens <- isTRUE(all(is.numeric(x$lb), is.finite(x$lb)))

  scale <- match.arg(scale, c("unif", "gev", "exp"), several.ok = TRUE)
  type <- type[type %in% type_plot]
  scale <- rep(scale, length.out = length(type))
  # The following should never be executed...
  if (lcens & "range" %in% type) {
    scale <- scale[type != "range"]
    type <- type[type != "range"]
  }
  if (length(type) == 0L) {
    stop("Invalid type: cannot produce plot of range with left-censored data")
  }

  for (t in seq_along(type)) {
    ind <- which(type_plot == type[t])
    lb_unif <- 0
    xx <- x$plots[[ind]]
    if (scale[t] == "gev") {
      gevpars <- switch(
        type[t],
        "max" = mev::maxstable(x$mle, m = x$m),
        "all" = x$mle,
        "range" = c(0, 1, 0)
      ) # gumbel for lack of better scale
    }
    # if (!lcens) {
    xpos <- switch(
      scale[t],
      "unif" = xx$x,
      "exp" = qexp(xx$x),
      "gev" = mev::qgev(
        xx$x,
        loc = gevpars[1],
        scale = gevpars[2],
        shape = gevpars[3]
      )
    )

    ypos <- switch(
      scale[t],
      "unif" = xx$y,
      "exp" = qexp(xx$y),
      "gev" = mev::qgev(
        xx$y,
        loc = gevpars[1],
        scale = gevpars[2],
        shape = gevpars[3]
      )
    )
    lower_confint <- switch(
      scale[t],
      "unif" = xx$confint$pointwise[, "lower"],
      "exp" = qexp(xx$confint$pointwise[, "lower"]),
      "gev" = mev::qgev(
        xx$confint$pointwise[, "lower"],
        loc = gevpars[1],
        scale = gevpars[2],
        shape = gevpars[3]
      )
    )
    upper_confint <- switch(
      scale[t],
      "unif" = xx$confint$pointwise[, "upper"],
      "exp" = qexp(xx$confint$pointwise[, "upper"]),
      "gev" = mev::qgev(
        xx$confint$pointwise[, "upper"],
        loc = gevpars[1],
        scale = gevpars[2],
        shape = gevpars[3]
      )
    )
    lower_simult <- switch(
      scale[t],
      "unif" = xx$confint$simultaneous[, "lower"],
      "exp" = qexp(xx$confint$simultaneous[, "lower"]),
      "gev" = mev::qgev(
        xx$confint$simultaneous[, "lower"],
        loc = gevpars[1],
        scale = gevpars[2],
        shape = gevpars[3]
      )
    )
    upper_simult <- switch(
      scale[t],
      "unif" = xx$confint$simultaneous[, "upper"],
      "exp" = qexp(xx$confint$simultaneous[, "upper"]),
      "gev" = mev::qgev(
        xx$confint$simultaneous[, "upper"],
        loc = gevpars[1],
        scale = gevpars[2],
        shape = gevpars[3]
      )
    )
    np <- length(xpos)
    ran <- switch(
      scale[t],
      "unif" = c(lb_unif, 1),
      "exp" = c(
        qexp(lb_unif),
        max(
          tail(xpos[is.finite(xpos)], 1),
          tail(ypos[is.finite(ypos)], 1),
          tail(upper_confint[is.finite(upper_confint)], 1),
          tail(upper_simult[is.finite(upper_simult)], 1)
        )
      ),
      "gev" = c(
        min(
          xpos[is.finite(xpos)][1],
          ypos[is.finite(ypos)][1],
          lower_confint[is.finite(lower_confint)][1],
          lower_simult[is.finite(lower_simult)][1]
        ),
        max(
          tail(xpos[is.finite(xpos)], 1),
          tail(ypos[is.finite(ypos)], 1),
          tail(upper_confint[is.finite(upper_confint)], 1),
          tail(upper_simult[is.finite(upper_simult)], 1)
        )
      )
    )
    plot(
      x = NULL,
      ylim = ran,
      xlim = ran,
      xaxs = "i",
      yaxs = "i",
      xlab = "theoretical quantiles",
      ylab = "empirical quantiles",
      bty = "l",
    )
    polygon(
      x = c(xpos, rev(xpos)),
      c(pmax(ran[1], lower_simult), rev(pmin(ran[2], upper_simult))),
      col = "grey95",
    )
    segments(
      x0 = xpos,
      x1 = xpos,
      y0 = pmax(ran[1], lower_confint),
      y1 = pmin(ran[2], upper_confint),
      col = "grey70"
    )
    abline(a = 0, b = 1)
    points(
      x = xpos,
      y = ypos,
      pch = 20,
      col = ifelse((ypos < lower_simult) | (ypos > upper_simult), 2, 1)
    )
  }
}

#' Diagnostic plots for max-stability based on blocks of GEV samples
#'
#' Given a sample of ordered GEV draws, calculate the ingredients of diagnostic
#' quantile-quantile plots using the bootstrap
#'
#' @param xdat \code{n} by \code{m} matrix of GEV observations, ordered by row from smallest to largest
#' @param type string; the statistic to return. Either the maximum of each row (\code{max}), the standardized difference between the penultimate and largest value (\code{spacing}), the ratio of maximum to spacing (\code{ratio}) or the whole sample (\code{all})
#' @param B number of bootstrap samples
#' @param marginal logical; if \code{TRUE}, estimates are based on the marginal likelihood of the \eqn{m-1} smallest order statistics of the sample
#' @param rounding amount of rounding
#' @param lb lower bound for left-censoring, default to \code{NULL} in absence
#' @param plot logical; if \code{TRUE} (default), returns a quantile-quantile plot
#' @param level confidence level for confidence and tolerance intervals
#' @param np number of points at which to evaluate quantile-quantile plots. Must be either \code{NULL}, or a vector of integer of the same length as \code{type} (otherwise it is recycled).
#' @export
#' @return a list with elements for building quantile-quantile plots, including
#' \itemize{
#' \item \code{plots} list of plots with elements \code{x}, \code{y}, a list \code{confint} with matrices \code{simultaneous} and \code{pointwise}, \code{type} of value and \code{distribution} (currently only uniform)
#' \item \code{mle}: maximum likelihood estimate of the location, scale, and shape
#' \item \code{param} \code{B} by 3 matrix of bootstrap parameter estimates
#' \item \code{type} vector of string with statistics
#' \item \code{bootstrap} type of bootstrap, only \code{parametric} for now
#' \item \code{n} number of rows of \code{xdat}
#' \item \code{m} number of columns of \code{xdat} for comparison
#' \item \code{marginal} logical; if \code{TRUE}, uses the marginal likelihood of the \eqn{m-1} smallest order statistics per block for estimation
#' \item \code{icens} logical; if \code{TRUE}, data treated as rounded (interval-censored)
#' \item \code{lcens} logical; if \code{TRUE}, data are left-censored below \code{lb}
#' \item \code{lb} lower bound for left-censoring
#' \item \code{rounding} double \eqn{\delta} indicating the amount of rounding, assuming \eqn{\delta/2} on either size of the reported value
#' \item \code{xdat} matrix of original observations
#' }
#' @examples
#' xdat <- build.blocks(mev::rgev(n = 50), m = 2)
#' \dontrun{
#' qqplot.blocksize(xdat, type = "max", marginal = TRUE, B = 100)
#' }
qqplot.blocksize <- function(
  xdat,
  type = c("max", "range", "all"),
  # scale = "unif",
  B = 1e3L,
  marginal = FALSE,
  rounding = 0,
  lb = NULL,
  plot = TRUE,
  level = 0.95,
  np = NULL
) {
  scale <- "unif"
  if (is.null(lb)) {
    lcens <- FALSE
  }
  if (is.null(rounding)) {
    icens <- FALSE
  } else {
    icens <- abs(rounding[1]) > 1e-12
  }
  lcens <- !is.null(lb)
  if (is.numeric(lb)) {
    lcens <- is.finite(lb[1])
  }
  scale <- match.arg(
    scale,
    choices = c("unif", "gev", "exp"),
    several.ok = TRUE
  )
  type <- unique(match.arg(type, several.ok = TRUE))
  scale <- rep(scale, length.out = length(type))
  if ("range" %in% type & lcens) {
    warning("May not be able to determine \"range\" with left-censoring.")
    scale <- scale[type != "range"]
    type <- type[type != "range"]
  }
  if (isTRUE(any(icens, lcens))) {
    qqplot.blocksize.rounded(
      xdat = xdat,
      type = type,
      B = B,
      marginal = marginal,
      rounding = rounding,
      lb = lb,
      scale = scale,
      plot = plot,
      level = level,
      np = np
    )
  } else {
    qqplot.blocksize.parametric(
      xdat = xdat,
      type = type,
      B = B,
      marginal = marginal,
      scale = scale,
      plot = plot,
      level = level,
      np = np
    )
  }
}

# Dispatch to different methods depending on fully observed / rounded / left-censored
qqplot.blocksize.parametric <- function(
  xdat,
  type = c("max", "range", "all"),
  B = 1e3L,
  marginal = FALSE,
  scale = "unif",
  plot = TRUE,
  level = 0.95,
  np = 500L
) {
  scale <- "unif"
  n <- nrow(xdat)
  m <- ncol(xdat)
  B <- as.integer(B)
  stopifnot(B > 0)
  level <- rep(level, length.out = 2)
  stopifnot(is.numeric(level), isTRUE(all(level > 0, level < 1)))
  alpha <- 1 - rep(level, length.out = 2L)
  type <- match.arg(type, several.ok = TRUE)
  # Create containers
  boot_out <- list()
  xdat_out <- list()
  pp <- list()
  if (!isTRUE(is.finite(np))) {
    np <- n * m
  }
  nobs <- integer(length(type))
  np <- rep(as.integer(np), length.out = length(type))
  for (t in seq_along(type)) {
    pp[[t]] <- (1:(np[t] - 1)) / np[t]
    boot_out[[t]] <- matrix(
      nrow = np[t] - 1,
      ncol = B
    )
  }
  mle_coefs_boot <- matrix(nrow = B, ncol = 3)
  # 0. obtain MLE
  mle <- fit.gevblock(
    xdat = xdat,
    marginal = marginal,
    constraint = TRUE
  )
  unif_xdat <- apply(
    xdat,
    2,
    mev::pgev,
    loc = mle[1],
    scale = mle[2],
    shape = mle[3]
  )
  for (t in seq_along(type)) {
    pivots_xdat <- switch(
      type[t],
      "max" = pbeta(q = unif_xdat[, m], shape1 = m, shape2 = 1),
      "range" = pbeta(
        q = unif_xdat[, m] - unif_xdat[, 1],
        shape1 = m - 1,
        shape2 = 2
      ),
      "all" = c(unif_xdat)
    )
    nobs[t] <- length(pivots_xdat)
    xdat_out[[t]] <- ecdf(pivots_xdat)(pp[[t]])
    #quantile(pivots_xdat, pp[[t]])
  }

  gamma <- matrix(nrow = B, ncol = length(type))
  # Bootstrap loop
  for (b in seq_len(B)) {
    # 1. Generate bootstrap sample
    xdat_boot <- build.blocks(
      mev::rgev(
        n = n * m,
        loc = mle[1],
        scale = mle[2],
        shape = mle[3]
      ),
      m = m
    )
    # 2. Fit MLE to bootstrap sample, map to uniform
    mle_boot <- try(
      fit.gevblock(
        xdat = xdat_boot,
        marginal = marginal,
        start = mle
      ),
      silent = TRUE
    )
    if (!inherits(mle_boot, "try-error")) {
      # 3. Map to uniform, calculate pivot and map to uniform again
      mle_coefs_boot[b, ] <- mle_boot
      unif_boot <- matrix(
        pgev(
          xdat_boot,
          loc = mle_boot[1],
          scale = mle_boot[2],
          shape = mle_boot[3]
        ),
        ncol = m
      )
      for (t in seq_along(type)) {
        pivots_boot <- switch(
          type[t],
          "max" = pbeta(q = unif_boot[, m], shape1 = m, shape2 = 1),
          "range" = pbeta(
            q = unif_boot[, m] - unif_boot[, 1],
            shape1 = m - 1,
            shape2 = 2
          ),
          "all" = c(unif_boot)
        )

        Fz_boot <- ecdf(pivots_boot)(pp[[t]])
        gamma[b, t] <- 2 *
          min(
            pbinom(
              length(pivots_boot) * Fz_boot,
              size = length(pivots_boot),
              prob = pp[[t]]
            ),
            1 -
              pbinom(
                length(pivots_boot) * Fz_boot - 1,
                size = length(pivots_boot),
                prob = pp[[t]]
              )
          )
        boot_out[[t]][, b] <- Fz_boot
      }
    }
  }

  distribution <- rep("unif", length.out = length(type))
  plots <- list()
  for (t in seq_along(type)) {
    ptwise_conf <- cbind(
      qbinom(alpha[1] / 2, prob = pp[[t]], size = nobs[t]) / nobs[t],
      qbinom(1 - alpha[1] / 2, prob = pp[[t]], size = nobs[t]) / nobs[t]
    )
    alpha_star <- quantile(gamma[, t], probs = alpha[2], na.rm = TRUE)
    simult_conf <- cbind(
      qbinom(alpha_star / 2, prob = pp[[t]], size = nobs[t]) / nobs[t],
      qbinom(1 - alpha_star / 2, prob = pp[[t]], size = nobs[t]) / nobs[t]
    )
    # colnames(conf) <-
    colnames(simult_conf) <- colnames(ptwise_conf) <- c(
      "lower",
      "upper"
    )

    plots[[t]] <- list(
      x = pp[[t]],
      y = as.numeric(xdat_out[[t]]),
      # confidence = conf,
      confint = list(pointwise = ptwise_conf, simultaneous = simult_conf),
      type = type[t],
      dist = distribution[t],
      alpha = c(alpha[1], alpha_star)
    )
  }
  out <- list(
    plots = plots,
    mle = mle,
    param = mle_coefs_boot,
    type = type,
    bootstrap = "parametric",
    B = B,
    n = n,
    m = m,
    marginal = marginal,
    icens = FALSE,
    lcens = FALSE,
    lb = NULL,
    rounding = 0,
    xdat = xdat
  )
  class(out) <- "mev_plot_blocksize"
  if (isTRUE(plot)) {
    plot(out, type = type, scale = scale)
  }
  return(invisible(out))
}


# Dispatch to different methods depending on fully observed / rounded / left-censored
qqplot.blocksize.rounded <- function(
  xdat,
  type = c("max", "range", "all"),
  B = 1e3L,
  marginal = FALSE,
  rounding = 0,
  lb = NULL,
  scale = "unif",
  plot = TRUE,
  level = 0.95,
  np = 500L
) {
  lcens <- !is.null(lb)
  if (lcens) {
    stopifnot(length(lb) == 1L, is.numeric(lb))
  }
  icens <- !is.null(rounding)
  if (icens) {
    if (abs(rounding) < 1e-12) {
      icens <- FALSE
    } else {
      stopifnot(length(rounding) == 1L, is.numeric(rounding))
      delta <- rounding / 2
    }
  }

  n <- nrow(xdat) # this is the number of replications
  # the number of exceedances of lb may be lower than n*m
  m <- ncol(xdat)
  B <- as.integer(B)
  stopifnot(B > 0)
  level <- rep(level, length.out = 2)
  stopifnot(is.numeric(level), isTRUE(all(level > 0, level < 1)))
  alpha <- 1 - rep(level, length.out = 2L)
  type <- unique(match.arg(type, several.ok = TRUE))
  # Create containers
  boot_out <- list()
  xdat_out <- list()
  if (lcens & is.null(np)) {
    np <- sapply(seq_along(type), function(t) {
      ifelse(type[t] == "all", length(xdat), nrow(xdat))
    })
  }
  pp <- list()
  if (!isTRUE(is.finite(np))) {
    np <- n * m
  }
  nobs <- integer(length(type))
  np <- rep(as.integer(np), length.out = length(type))
  for (t in seq_along(type)) {
    pp[[t]] <- (1:(np[t] - 1)) / np[t]
    boot_out[[t]] <- matrix(
      nrow = np[t] - 1,
      ncol = B
    )
  }
  # Utility function to simulate continuous observations from rounded records
  impute_rounded <- function(x, rounding, pars) {
    p_lb <- pgev(
      x - rounding / 2,
      loc = pars[1],
      scale = pars[2],
      shape = pars[3]
    )
    p_ub <- pgev(
      x + rounding / 2,
      loc = pars[1],
      scale = pars[2],
      shape = pars[3]
    )
    qgev(
      p_lb + runif(length(x)) * (p_ub - p_lb),
      loc = pars[1],
      scale = pars[2],
      shape = pars[3]
    )
  }
  mle_coefs_boot <- matrix(nrow = B, ncol = 3)
  # 0. obtain MLE
  mle <- fit.gevblock(
    xdat = xdat,
    marginal = marginal,
    rounding = rounding,
    lb = lb,
    constraint = TRUE
  )
  xdat_new <- c(t(xdat))
  if (icens) {
    xdat_new <- impute_rounded(x = xdat_new, rounding = rounding, pars = mle)
  }

  # Calculate statistics on the original sample
  unif_xdat <- build.blocks(
    mev::pgev(
      xdat_new,
      loc = mle[1],
      scale = mle[2],
      shape = mle[3]
    ),
    m = m
  )
  if (lcens) {
    lb_unif <- pgev(lb, loc = mle[1], scale = mle[2], shape = mle[3])
    mmle <- maxstable(mle, m = m)
    lb_max_unif <- pgev(lb, loc = mmle[1], scale = mmle[2], shape = mmle[3])
    lcens_mat <- unif_xdat < lb_unif
  }
  for (t in seq_along(type)) {
    if (lcens) {
      pivots_xdat <- switch(
        type[t],
        # Focus on either spacing X_{(m)} - X_{(m-1)} or X_{(m)}
        "max" = ifelse(
          lcens_mat[, m],
          NA,
          pbeta(q = unif_xdat[, m], shape1 = m, shape2 = 1)
        ),
        "range" = ifelse(
          lcens_mat[, 1],
          NA,
          pbeta(
            q = unif_xdat[, m] - unif_xdat[, 1],
            shape1 = m - 1,
            shape2 = 2
          )
        ),
        "all" = ifelse(c(lcens_mat), NA, c(unif_xdat))
      )
      if (type[t] == "max") {
        lb_unif_pivot <- pbeta(q = lb_unif, shape1 = m, shape2 = 1)
      } else {
        lb_unif_pivot <- lb_unif
      }
      pivots_xdat <- (pivots_xdat[!is.na(pivots_xdat)] - lb_unif_pivot) /
        (1 - lb_unif_pivot)
    } else {
      pivots_xdat <- switch(
        type[t],
        # Focus on either spacing X_{(m)} - X_{(m-1)} or X_{(m)}
        "max" = pbeta(q = unif_xdat[, m], shape1 = m, shape2 = 1),
        "range" = pbeta(
          q = unif_xdat[, m] - unif_xdat[, 1],
          shape1 = m - 1,
          shape2 = 2
        ),
        "all" = c(unif_xdat)
      )
    }
    nobs[t] <- length(pivots_xdat)
    xdat_out[[t]] <- ecdf(pivots_xdat)(pp[[t]])
  }
  gamma <- matrix(nrow = B, ncol = length(type))
  # Bootstrap loop
  for (b in seq_len(B)) {
    # 1. Generate bootstrap sample
    xdat_boot <- build.blocks(
      mev::rgev(
        n = n * m,
        loc = mle[1],
        scale = mle[2],
        shape = mle[3]
      ),
      m = m
    )
    if (icens) {
      # Round to same precision as observation
      xdat_boot <- round(xdat_boot / rounding, 0) * rounding
    }
    # 2. Fit MLE to bootstrap sample,
    mle_boot <- try(
      fit.gevblock(
        xdat = xdat_boot,
        marginal = marginal,
        start = mle,
        constraint = TRUE,
        rounding = rounding,
        lb = lb
      ),
      silent = TRUE
    )
    if (!inherits(mle_boot, "try-error")) {
      xdat_boot_new <- c(t(xdat_boot))
      if (icens) {
        # make continuous
        xdat_boot_new <- impute_rounded(
          x = xdat_boot_new,
          rounding = rounding,
          pars = mle_boot
        )
      }
      # 3. Map to uniform, calculate pivot and map to uniform again
      mle_coefs_boot[b, ] <- mle_boot
      unif_boot <- build.blocks(
        pgev(
          xdat_boot_new,
          loc = mle_boot[1],
          scale = mle_boot[2],
          shape = mle_boot[3]
        ),
        m = m
      )
      if (lcens) {
        lb_unif_boot <- pgev(
          lb,
          loc = mle_boot[1],
          scale = mle_boot[2],
          shape = mle_boot[3]
        )
        lcens_mat_boot <- unif_boot < lb_unif_boot
      }
      for (t in seq_along(type)) {
        if (lcens) {
          pivots_boot <- switch(
            type[t],
            # Focus on either spacing X_{(m)} - X_{(m-1)} or X_{(m)}
            "max" = pbeta(q = unif_boot[, m], shape1 = m, shape2 = 1)[
              !lcens_mat_boot[, m]
            ],
            "all" = c(unif_boot)[c(!lcens_mat_boot)]
          )
          if (type[t] == "max") {
            lb_unif_boot_pivot <- pbeta(
              q = lb_unif_boot,
              shape1 = m,
              shape2 = 1
            )
          } else if (type[t] == "all") {
            lb_unif_boot_pivot <- lb_unif_boot
          }
          pivots_boot <- (pivots_boot - lb_unif_boot_pivot) /
            (1 - lb_unif_boot_pivot)
        } else {
          pivots_boot <- sort(switch(
            type[t],
            "max" = pbeta(q = unif_boot[, m], shape1 = m, shape2 = 1),
            "range" = pbeta(
              q = unif_boot[, m] - unif_boot[, 1],
              shape1 = m - 1,
              shape2 = 2
            ),
            "all" = c(unif_boot)
          ))
        }
        Fz_boot <- ecdf(pivots_boot)(pp[[t]])
        gamma[b, t] <- 2 *
          min(
            pbinom(
              length(pivots_boot) * Fz_boot,
              size = length(pivots_boot),
              prob = pp[[t]]
            ),
            1 -
              pbinom(
                length(pivots_boot) * Fz_boot - 1,
                size = length(pivots_boot),
                prob = pp[[t]]
              )
          )
        boot_out[[t]][, b] <- Fz_boot # quantile(pivots_boot, probs = pp[[t]])
      }
    }
  }

  distribution <- rep("unif", length.out = length(type))
  plots <- list()
  for (t in seq_along(type)) {
    ptwise_conf <- cbind(
      qbinom(alpha[1] / 2, prob = pp[[t]], size = nobs[t]) / nobs[t],
      qbinom(1 - alpha[1] / 2, prob = pp[[t]], size = nobs[t]) / nobs[t]
    )
    alpha_star <- quantile(gamma[, t], probs = alpha[2], na.rm = TRUE)
    simult_conf <- cbind(
      qbinom(alpha_star / 2, prob = pp[[t]], size = nobs[t]) / nobs[t],
      qbinom(1 - alpha_star / 2, prob = pp[[t]], size = nobs[t]) / nobs[t]
    )
    colnames(simult_conf) <- colnames(ptwise_conf) <- c(
      "lower",
      "upper"
    )
    plots[[t]] <- list(
      x = pp[[t]],
      y = as.numeric(sort(xdat_out[[t]])),
      confint = list(pointwise = ptwise_conf, simultaneous = simult_conf),
      type = type[t],
      dist = distribution[t],
      alpha = c(alpha[1], alpha_star)
    )
  }
  out <- list(
    plots = plots,
    xdat = xdat,
    mle = mle,
    param = mle_coefs_boot,
    type = type,
    bootstrap = "parametric",
    B = B,
    n = n,
    m = m,
    marginal = marginal,
    icens = icens,
    lcens = lcens,
    lb = lb,
    rounding = rounding
  )
  class(out) <- "mev_plot_blocksize"
  if (isTRUE(plot)) {
    plot(out, type = type, scale = scale)
  }
  return(invisible(out))
}


#' Pointwise and simultaneous binomial confidence intervals for uniform via simulation
#'
#' Given a vector of draws transformed using the probability integral transform scale
#' to what should be uniform positions, produce plots with pointwise and simultaneous confidence intervals for uniformity.
#'
#' @param xdat vector of \code{N} postulated uniform samples, obtained by applying the ECDF
#' @param K number of evaluation points for the plotting positions
#' @param B number of Monte Carlo samples
#' @param level vector of pointwise and simultaneous confidence levels, recycled if necessary
#' @param plot logical; if \code{TRUE}, produce a plot of the empirical distribution function
#' @references Sailynoja, T., Burkner, P.C. and Vehtari, A. (2022). Graphical test for discrete uniformity and its applications in goodness-of-fit evaluation and multiple sample comparison, \emph{Statistics and Computing}, 32, \doi{10.1007/s11222-022-10090-6}
#' @export
#' @examples
#' xdat <- runif(200)
#' qqplot.unif(xdat)
qqplot.unif <- function(
  xdat,
  K = 100,
  B = 1000,
  level = 0.95,
  plot = TRUE
) {
  level <- rep(level, length = 2L)
  stopifnot(isTRUE(all(is.finite(level), level > 0, level < 1)))
  alpha <- (1 - level)
  u <- sort(xdat)
  n <- length(u)
  stopifnot(u[1] >= 0, u[n] <= 1)
  # Plotting positions, equally spaced (do not evaluate 0 and N
  z <- 1:(K - 1) / K
  ecdf_u <- ecdf(u)
  Fz <- ecdf_u(z)
  # plot(z, Fz)
  ## Pointwise confidence intervals
  ptwise_conf <- cbind(
    qbinom(alpha[1] / 2, prob = z, size = n) / n,
    qbinom(1 - alpha[1] / 2, prob = z, size = n) / n
  )
  # Simultaneous confidence intervals via simulation
  gamma <- numeric(B)
  for (b in seq_len(B)) {
    Fz_boot <- ecdf(runif(n))(z)
    gamma[b] <- 2 *
      min(
        pbinom(n * Fz_boot, size = n, prob = z),
        1 - pbinom(n * Fz_boot - 1, size = n, prob = z)
      )
  }
  alpha_star <- quantile(gamma, alpha[2])
  simult_conf <- cbind(
    qbinom(alpha_star / 2, prob = z, size = n) / n,
    qbinom(1 - alpha_star / 2, prob = z, size = n) / n
  )
  out <- list(x = z, y = Fz, ptwise = ptwise_conf, simult = simult_conf)
  class(out) <- "mev_ecdf_unif"
  if (isTRUE(plot)) {
    plot(out)
  }
  return(invisible(out))
}
#' @export
plot.mev_ecdf_unif <- function(x, ...) {
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    bty = "n",
    xlab = "",
    ylab = ""
  )
  segments(
    x0 = x$x,
    y0 = x$simult[, 1],
    y1 = x$simult[, 2],
    lwd = 1,
    col = "grey"
  )
  segments(x0 = x$x, y0 = x$ptwise[, 1], y1 = x$ptwise[, 2], lwd = 0.5)
  points(
    x = x$x,
    y = x$y,
    pch = 20,
    col = ifelse((x$y > x$simult[, 2] | x$y < x$simult[, 1]), "red", "black")
  )
}
