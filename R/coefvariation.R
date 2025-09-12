#' Threshold selection via coefficient of variation
#'
#' This function computes the empirical coefficient of variation and
#' computes a weighted statistic comparing the squared distance with
#' the theoretical coefficient variation corresponding to a specific
#' shape parameter (estimated from the data using a moment estimator
#' as the value minimizing the test statistic, or using maximum likelihood).
#' The procedure stops if there are no more than 10 exceedances above the
#' highest threshold.
#'
#' @note The authors suggest transformation of \deqn{Y = -1/(X + c) + 1/c,}
#' where \eqn{X} are exceedances and \eqn{c=\sigma/\xi} is the ratio of estimated scale and shape parameters. For heavy-tailed distributions with \eqn{\xi > 0.25}, this may be preferable, but must be conducted outside of the function.
#'
#' @param xdat [vector] vector of observations
#' @param thresh [vector] vector of threshold. If missing, set to \eqn{p^k} for \eqn{k=0} to \eqn{k=}\code{nthresh}
#' @param method [string], either moment estimator for the (weighted) coefficient of variation (\code{wcv} and \code{cv}) or maximum likelihood (\code{mle})
#' @param nsim [integer] number of bootstrap replications
#' @param nthresh [integer] number of thresholds, if \code{thresh} is not supplied by the user
#' @param level [numeric] probability level for sequential testing procedure
#' @param lazy [logical] compute the bootstrap p-value until the test stops rejecting at level \code{level}? Default to \code{FALSE}
#' @param plot [logical] if \code{TRUE}, returns a plot of the p-value path
#' @return a list with elements
#' \itemize{
#' \item \code{thresh}: value of threshold returned by the procedure, \code{NA} if the hypothesis is rejected at all thresholds
#' \item \code{thresh0}: sorted vector of candidate thresholds
#' \item \code{cindex}: index of selected threshold among \code{thresh0} or \code{NA} if none returned
#' \item \code{pval}: bootstrap p-values, with \code{NA} if \code{lazy} and the p-value exceeds level at lower thresholds
#' \item \code{shape}: shape parameter estimates
#' \item \code{nexc}: number of exceedances of each threshold \code{thresh0}
#' \item \code{method}: estimation method for the shape parameter
#' }
#' @export
#' @references del Castillo, J. and M. Padilla (2016). \emph{Modelling extreme values by the residual coefficient of variation}, SORT, 40(\bold{2}), pp. 303--320.
#' @examples
#' thselect.cv(
#'  xdat = rgp(1000),
#'  thresh = qgp(seq(0,0.9, by = 0.1)),
#'  nsim = 99,
#'  lazy = TRUE,
#'  plot = TRUE)
#'
thselect.cv <- function(
  xdat,
  thresh,
  method = c("mle", "wcv", "cv"),
  nsim = 999L,
  nthresh = 10L,
  level = 0.05,
  lazy = FALSE,
  plot = FALSE
) {
  method <- match.arg(method)
  xdat <- as.numeric(xdat[is.finite(xdat)])
  stopifnot(
    length(level) == 1L,
    is.finite(level),
    level >= 0,
    level <= 1,
    is.logical(lazy),
    length(lazy) == 1L
  )
  if (isTRUE(plot)) {
    lazy <- FALSE
  }
  # Set grid of thresholds or order them and keep exceedances
  if (!missing(thresh)) {
    thresh <- thresh[is.finite(thresh)]
    thresh <- sort(thresh, decreasing = TRUE)
    nthresh <- length(thresh)
    xdat <- xdat[xdat >= min(thresh)]
    # Weight vector is survival probability of thresholds
    n <- length(xdat)
  } else {
    n <- length(xdat)
    stopifnot(
      length(nthresh) == 1L,
      is.finite(nthresh),
      nthresh > 1
    )
    nthresh <- as.integer(nthresh)
    # Set thresholds at 1-p^k for k=0, ..., nthresh quantiles
    # Ensure there are enough observations
    p <- (10 / n)^(1 / nthresh)
    # Survival probability of the thresholds
    survprob <- p^((nthresh - 1):0)
    thresh <- quantile(xdat, 1 - survprob)
  }
  # Number of exceedances above thresholds
  nexc <- sapply(thresh, function(x) {
    sum(xdat > x)
  })
  if (nexc[1] < 10) {
    stop(
      "Threshold is too small: the procedure requires at least 10 exceedances"
    )
  }
  survprob <- nexc / n
  # Coefficient of variation
  coefvar <- function(x, na.rm = TRUE) {
    sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
  }
  # Coefficient of variation for exceedances
  cv <- sapply(thresh, function(xth) {
    coefvar(xdat[xdat > xth] - xth)
  })
  if (method == "wcv") {
    cv0 <- cumsum(cv * survprob) / cumsum(survprob)
    shape <- 0.5 * (1 - 1 / cv0^2)
  } else if (method == "mle") {
    shape <- as.numeric(sapply(thresh, function(th) {
      coef(mev::fit.gpd(xdat = xdat, threshold = th))[2]
    }))
    shape2cv <- function(shape) {
      ifelse(shape < 0.5, (1 - 2 * shape)^(-0.5), NA)
    }
    cv0 <- shape2cv(shape)
  } else if (method == "cv") {
    cv0 <- cv
    shape <- 0.5 * (1 - 1 / cv0^2)
  }
  if (isTRUE(any(shape > 0.25))) {
    warning(
      "Estimated shape parameter larger than 1/4. Variance of coefficient of variation may be infinite."
    )
  }
  stat <- sapply(seq_along(shape), function(i) {
    sum(nexc[1:i] * (cv[1:i] - cv0[i])^2)
  })
  # Function to loop over
  compute_cvstat <- function(xdat, survprob, method = "mle") {
    # n*plevel is number of exceedances above threshold
    survprob <- sort(survprob, decreasing = TRUE)
    n <- length(xdat)
    # Thresholds are ordered from largest to smallest
    thresh <- quantile(xdat, 1 - survprob)
    # Compute coefficient of variation
    cv <- sapply(thresh, function(x) {
      coefvar(xdat[xdat > x] - x)
    })
    # Logical: compute by minimizing stat
    if (method == "mle") {
      shape <- as.numeric(
        coef(
          mev::fit.gpd(xdat = xdat, threshold = 0)
        )[2]
      )
      if (!isTRUE(shape < 0.5)) {
        return(NA)
      }
      # Null value
      cv0 <- (1 - 2 * shape)^(-0.5)
    } else if (method == "cv") {
      cv0 <- cv[1]
    } else if (method == "wcv") {
      cv0 <- cumsum(cv * survprob) / cumsum(survprob)
    }
    n * sum(survprob * (cv - cv0)^2)
  }
  # Bootstrap function
  bootfun <- function(n, shape, survprob, method = "mle") {
    dat <- rgp(n = n, loc = 0, scale = 1, shape = shape)
    compute_cvstat(
      xdat = dat,
      method = method,
      survprob = survprob
    )
  }
  boot_pval <- rep(NaN, nthresh)
  thselect <- NaN
  cindex <- NULL
  # Bootstrap loop
  for (r in seq_along(thresh)) {
    # Number of exceedances, null shape and value of the statistic
    # return the bootstrap p-value
    rk <- nthresh - r + 1
    boot_stat <- replicate(
      n = nsim,
      expr = bootfun(
        n = nexc[rk],
        survprob = nexc[1:rk] / nexc[rk],
        shape = shape[rk]
      )
    )
    boot_pval[r] <- mean(c(boot_stat >= stat[rk], TRUE), na.rm = TRUE)
    if (boot_pval[r] >= level[1] & is.na(thselect)) {
      thselect <- thresh[rk]
      cindex <- r
      if (lazy) {
        break
      }
    }
  }
  ret <- list(
    thresh = rev(as.numeric(thresh)),
    thresh0 = as.numeric(thselect),
    cindex = cindex,
    pval = as.numeric(boot_pval),
    shape = as.numeric(shape),
    nexc = as.numeric(nexc),
    method = method,
    level = level
  )
  class(ret) <- "mev_thselect_cv"
  if (isTRUE(plot)) {
    plot(ret)
  }
  return(invisible(ret))
}
#' @export
print.mev_thselect_cv <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  cat("Threshold selection method: coefficient of variation\n")
  cat(
    "Method:",
    switch(
      x$method,
      wcv = "weighted moment estimator",
      cv = "moment estimator",
      mle = "maximum likelihood estimator"
    ),
    "bootstrap p-value:",
    x$pval[x$cindex],
    "\n"
  )
  cat("Selected threshold:", round(x$thresh0, digits), "\n")
  return(invisible(NULL))
}

#' @export
plot.mev_thselect_cv <- function(x, ...) {
  plot(
    y = x$pval,
    x = x$thresh,
    pch = 20,
    bty = "l",
    ylim = c(0, 1),
    yaxs = "i",
    xlab = "threshold",
    ylab = "p-value",
    panel.first = {
      abline(h = x$level, col = "grey")
    }
  )
  abline(v = x$thresh0, lty = 2)
}

#' Coefficient of variation threshold stability plot
#'
#' This function calculates parametric estimates of the coefficient of variation with pointwise Wald confidence intervals along with empirical estimates and returns a threshold stability plot.
#' @inheritParams thselect.cv
#' @param ... additional parameters, notably for package \code{boot}, for the \code{type} of confidence intervals.
#' @export
#' @examples
#' tstab.cv(
#' xdat = rgp(1000),
#' thresh = qgp(seq(0,0.9, by = 0.1)),
#' method = "cv")
#'  tstab.cv(
#' xdat = rgp(1000),
#' thresh = qgp(seq(0,0.9, by = 0.1)),
#' method = "empirical")
tstab.cv <- function(
  xdat,
  thresh,
  method = c("empirical", "mle", "wcv", "cv"),
  nthresh = 10L,
  nsim = 99L,
  plot = TRUE,
  level = 0.95,
  ...
) {
  args <- list(...)
  method <- match.arg(method)
  if (method == "empirical") {
    if (!is.null(args$type)) {
      type <- match.arg(
        arg = args$type,
        choices = c("norm", "basic", "stud", "perc", "bca"),
        several.ok = FALSE
      )
    } else {
      type <- "basic"
    }
  }
  method <- match.arg(method)
  xdat <- as.numeric(xdat[is.finite(xdat) & xdat > 0])
  # Set grid of thresholds or order them and keep exceedances
  if (!missing(thresh)) {
    thresh <- thresh[is.finite(thresh)]
    thresh <- sort(thresh, decreasing = TRUE)
    nthresh <- length(thresh)
    xdat <- xdat[xdat >= thresh[length(thresh)]]
    # Weight vector is survival probability of thresholds
    n <- length(xdat)
  } else {
    # xdat <- xdat[xdat > 0]
    n <- length(xdat)
    stopifnot(length(nthresh) == 1L, is.finite(nthresh), nthresh > 1)
    nthresh <- as.integer(nthresh)
    # Set thresholds at 1-p^k for k=0, ..., nthresh quantiles
    # Ensure there are enough observations
    p <- (10 / n)^(1 / nthresh)
    # Survival probability of the thresholds
    survprob <- p^((nthresh - 1):0)
    thresh <- quantile(xdat, 1 - survprob)
  }
  # Number of exceedances above thresholds
  nexc <- sapply(thresh, function(x) {
    sum(xdat > x)
  })
  if (nexc[1] < 10) {
    stop(
      "Invalid threshold: the procedure requires at least 10 exceedances for all candidate thresholds."
    )
  }
  survprob <- nexc / n
  # Coefficient of variation
  # Syntax is there to match the boot function call
  coefvar <- function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  }
  # Coefficient of variation for exceedances
  cv <- sapply(thresh, function(x) {
    coefvar(xdat[xdat > x] - x)
  })

  if (method == "wcv") {
    cv0 <- cumsum(cv * survprob) / cumsum(survprob)
    shape <- 0.5 * (1 - 1 / cv0^2)
  } else if (method == "mle") {
    shape <- as.numeric(sapply(thresh, function(th) {
      coef(mev::fit.gpd(xdat = xdat, threshold = th))[2]
    }))
    shape2cv <- function(shape) {
      ifelse(shape < 0.5, (1 - 2 * shape)^(-0.5), NA)
    }
    cv0 <- shape2cv(shape)
  } else if (method == "cv") {
    cv0 <- cv
    shape <- 0.5 * (1 - 1 / cv0^2)
  } else if (method == "empirical") {
    cv0 <- cv
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop(
        "Package \"boot\" must be installed to use this function.",
        call. = FALSE
      )
    }
    lower <- upper <- numeric(nthresh)
    for (i in seq_along(thresh)) {
      coefvar_boot <- function(data, ind, ...) {
        sd(data[ind], na.rm = TRUE) / mean(data[ind], na.rm = TRUE)
      }
      boot.out <- boot::boot(
        data = xdat[xdat > thresh[i]] - thresh[i],
        statistic = coefvar_boot,
        R = nsim
      )
      ci <- boot::boot.ci(
        boot.out = boot.out,
        conf = level,
        type = type
      )
      lower[i] <- ci[[type]][1, 4]
      upper[i] <- ci[[type]][1, 5]
    }
  }
  if (method != "empirical") {
    type <- NULL
    alpha <- 1 - level
    cv.sd <- function(shape, nexc) {
      stopifnot(length(shape) == length(nexc))
      ifelse(
        shape >= 0.25,
        NA,
        sqrt(
          (1 - shape)^2 *
            (6 * shape^2 - shape + 1) /
            ((1 - 2 * shape)^2 * (1 - 3 * shape) * (1 - 4 * shape)) /
            nexc
        )
      )
    }
    sd_cv <- cv.sd(shape, nexc)
    lower <- cv0 + qnorm(alpha / 2) * sd_cv
    upper <- cv0 + qnorm(1 - alpha / 2) * sd_cv
  }
  res <- list(
    thresh = thresh,
    cv = cv0,
    lower = lower,
    upper = upper,
    method = method,
    type = type
  )
  class(res) <- "mev_tstab_cv"
  if (isTRUE(plot)) {
    plot(res)
  }
  return(invisible(res))
}

#' @export
plot.mev_tstab_cv <- function(x, level = 0.05, legend = TRUE, ...) {
  args <- list(...)
  if (is.null(args$xlab)) {
    args$xlab <- 'threshold'
  }
  if (is.null(args$ylab)) {
    args$ylab <- 'coefficient of variation'
  }
  if (is.null(args$pch)) {
    args$pch <- 20
  }
  stopifnot(
    is.numeric(level),
    length(level) == 1L,
    level < 1,
    level > 0
  )
  ylims <- c(
    min(x$lower, na.rm = TRUE),
    max(x$upper, na.rm = TRUE)
  )
  plot(
    x = x$thresh,
    y = x$cv,
    ylim = ylims,
    xlab = args$xlab,
    ylab = args$ylab,
    pch = args$pch,
    bty = "l"
  )
  for (i in seq_along(x$thresh)) {
    segments(
      x0 = x$thresh[i],
      x1 = x$thresh[i],
      y0 = x$lower[i],
      y1 = x$upper[i]
    )
  }
  mtext(
    side = 1,
    adj = 1,
    line = 3,
    cex = 0.8,
    outer = FALSE,
    text = paste("method:", x$method)
  )
  return(invisible(x))
}

#' Threshold selection via coefficient of variation
#'
#' This function computes the empirical coefficient of variation and
#' computes a weighted statistic comparing the squared distance with
#' the theoretical coefficient variation corresponding to a specific
#' shape parameter (estimated from the data using a moment estimator
#' as the value minimizing the test statistic, or using maximum likelihood).
#' The procedure stops if there are no more than 10 exceedances above the
#' highest threshold
#'
#' @param xdat [vector] vector of observations
#' @param thresh [vector] vector of threshold. If missing, set to \eqn{p^k} for \eqn{k=0} to \eqn{k=}\code{nthresh}
#' @param method [string], either moment estimator for the (weighted) coefficient of variation (\code{wcv} and \code{cv}) or maximum likelihood (\code{mle})
#' @param nsim [integer] number of bootstrap replications
#' @param nthresh [integer] number of thresholds, if \code{thresh} is not supplied by the user
#' @param level [numeric] probability level for sequential testing procedure
#' @param lazy [logical] compute the bootstrap p-value until the test stops rejecting at level \code{level}? Default to \code{FALSE}
#' @return a list with elements
#' \itemize{
#' \item \code{thresh0}: value of threshold returned by the procedure, \code{NA} if the hypothesis is rejected at all thresholds
#' \item \code{thresh}: sorted vector of candidate thresholds
#' \item \code{cindex}: index of selected threshold among \code{thresh} or \code{NA} if none returned
#' \item \code{pval}: bootstrap p-values, with \code{NA} if \code{lazy} and the p-value exceeds level at lower thresholds
#' \item \code{shape}: shape parameter estimates
#' \item \code{nexc}: number of exceedances of each threshold \code{thresh}
#' \item \code{method}: estimation method for the shape parameter
#' }
#' @export
#' @references del Castillo, J. and M. Padilla (2016). \emph{Modelling extreme values by the residual coefficient of variation}, SORT, 40(\bold{2}), pp. 303--320.
#' @keywords internal
cvselect <- function(
  xdat,
  thresh,
  method = c("mle", "wcv", "cv"),
  nsim = 999L,
  nthresh = 10L,
  level = 0.05,
  lazy = FALSE
) {
  .Deprecated(new = "thselect.cv", package = "mev")
  thselect.cv(
    xdat = xdat,
    thresh = thresh,
    method = method,
    nsim = nsim,
    nthresh = nthresh,
    level = level,
    lazy = lazy,
    plot = FALSE
  )
}
