#' Extended generalised Pareto families
#'
#' @description This function provides the log-likelihood and quantiles for the three different families presented
#' in Papastathopoulos and Tawn (2013) and the two proposals of Gamet and Jalbert (2022), plus exponential tilting. All of the models contain an additional parameter, \eqn{\kappa \ge 0}.
#' All families share the same tail index as the generalized Pareto distribution, while allowing for lower thresholds.
#' For most models, the distribution reduce to the generalised Pareto when \eqn{\kappa=1} (for models \code{gj-tnorm} and \code{logist}, on the boundary of the parameter space when \eqn{kappa \to 0}.
#'
#' @references Papastathopoulos, I. and J. Tawn (2013). Extended generalised Pareto models for tail estimation, \emph{Journal of Statistical Planning and Inference} \bold{143}(3), 131--143, <doi:10.1016/j.jspi.2012.07.001>.
#' @references Gamet, P. and Jalbert, J. (2022). A flexible extended generalized Pareto distribution for tail estimation. \emph{Environmetrics}, \bold{33}(6), <doi:10.1002/env.2744>.
#' @section Usage:
#' \code{egp.ll(xdat, thresh, model, par)}
#' @section Usage:
#' \code{egp.retlev(xdat, thresh, par, model, p, plot=TRUE)}
#' @description \code{egp.retlev} gives the return levels for the extended generalised Pareto distributions
#' @name egp
#' @param xdat vector of observations, greater than the threshold
#' @param thresh threshold value
#' @param par parameter vector (\eqn{\kappa}, \eqn{\sigma}, \eqn{\xi}).
#' @param model a string indicating which extended family to fit
#' @param show logical; if \code{TRUE}, print the results of the optimization
#' @param p extreme event probability; \code{p} must be greater than the rate of exceedance for the calculation to make sense. See \bold{Details}.
#' @param plot logical; if \code{TRUE}, a plot of the return levels
#' @importFrom grDevices rainbow
#'
#' @details
#'
#' For return levels, the \code{p} argument can be related to \eqn{T} year exceedances as follows:
#' if there are \eqn{n_y} observations per year, than take \code{p}
#' to equal \eqn{1/(Tn_y)} to obtain the \eqn{T}-years return level.
#' @author Leo Belzile
#' @return \code{egp.ll} returns the log-likelihood value.
#' @return \code{egp.retlev} returns a plot of the return levels if \code{plot=TRUE} and a matrix of return levels.
#' @examples
#' set.seed(123)
#' xdat <- mev::rgp(1000, loc = 0, scale = 2, shape = 0.5)
#' par <- fit.egp(xdat, thresh = 0, model = 'gj-beta')$par
#' p <- c(1/1000, 1/1500, 1/2000)
#' # With multiple thresholds
#' th <- c(0, 0.1, 0.2, 1)
#' opt <- tstab.egp(xdat, th, model = 'gj-beta')
#' egp.retlev(xdat, opt$thresh, opt$par, 'gj-beta', p = p)
#' opt <- tstab.egp(xdat, th, model = 'pt-power', plots = NA)
#' egp.retlev(xdat, opt$thresh, opt$par, 'pt-power', p = p)
NULL

#' Fit of extended GP models and parameter stability plots
#'
#' This function is an alias of \code{\link{fit.egp}}.
#'
#' Supported for backward compatibility
#'
#' @export
#' @keywords internal
egp.fit <- function(
  xdat,
  thresh,
  model = c(
    "pt-beta",
    "pt-gamma",
    "pt-power",
    "gj-tnorm",
    "gj-beta",
    "exptilt",
    "logist"
  ),
  init,
  show = FALSE
) {
  fit.egp(
    xdat = xdat,
    thresh = thresh,
    model = model,
    init = init,
    show = show
  )
}

#' Extended generalised Pareto log likelihood
#'
#' This function provides the log-likelihood and quantiles for the different extended generalized Pareto distributions. All families share the same tail index as the generalized Pareto, and reduce to the latter when \code{kappa=1}.
#' @export
#' @inheritParams egp
#' @name egp-function
#' @keywords internal
egp.ll <- function(
  xdat,
  thresh,
  par,
  model = c(
    "pt-beta",
    "pt-gamma",
    "pt-power",
    "gj-tnorm",
    "gj-beta",
    "exptilt",
    "logist"
  )
) {
  if (model %in% c("egp1", "egp2", "egp3") & length(model) == 1) {
    model <- switch(
      model,
      egp1 = "pt-beta",
      egp2 = "pt-gamma",
      egp3 = "pt-power"
    )
  }
  par <- as.numeric(par) # strip names
  model <- match.arg(model)
  if (isTRUE(any(xdat < thresh))) {
    xdat = xdat[xdat > thresh] - thresh
  }
  kappa = par[1]
  sigma = par[2]
  xi = par[3]
  if (sigma < 0 || kappa < 0) {
    return(-Inf)
  }
  args = pmax(0, (1 + xi * xdat / sigma))
  if (abs(xi) > 1e-08) {
    sum(degp(
      x = xdat,
      scale = sigma,
      shape = xi,
      kappa = kappa,
      model = model,
      log = TRUE
    ))
  } else {
    # if xi=0
    if (model %in% c("pt-beta", "pt-gamma")) {
      sum(dgamma(x = xdat, shape = kappa, scale = sigma, log = TRUE))
    } else {
      sum(degp(
        x = xdat,
        scale = sigma,
        shape = xi,
        kappa = kappa,
        model = model,
        log = TRUE
      ))
    }
  }
}

#' @name egp-function
#' @export
#' @inheritParams egp
#' @keywords internal
#' @importFrom grDevices hcl.colors
egp.retlev <- function(
  xdat,
  thresh,
  par,
  model = c(
    "pt-beta",
    "pt-gamma",
    "pt-power",
    "gj-tnorm",
    "gj-beta",
    "exptilt",
    "logist"
  ),
  p,
  plot = TRUE
) {
  if (model %in% c("egp1", "egp2", "egp3") & length(model) == 1) {
    model <- switch(
      model,
      egp1 = "pt-beta",
      egp2 = "pt-gamma",
      egp3 = "pt-power"
    )
  }
  model <- match.arg(model)
  if (length(par) %% 3 != 0) {
    stop("Invalid parameter input")
  }
  if (!inherits(par, "matrix")) {
    par = matrix(c(par), ncol = 3)
  }
  rate <- sapply(thresh, function(u) {
    length(xdat[xdat > u]) / length(xdat)
  })
  if (!isTRUE(all.equal(length(rate), length(thresh), nrow(par)))) {
    stop("Input dimension does not match")
  }
  retlev <- matrix(0, nrow = length(thresh), ncol = length(p))
  if (
    any(sapply(rate, function(zeta) {
      zeta < p
    }))
  ) {
    warning(
      "Some probabilities \"p\" are higher than the exceedance rate. Evaluate those empirically."
    )
  }
  p <- sort(p)
  pq <- rev(1 / p)
  np <- length(p)
  for (i in 1:length(thresh)) {
    for (j in 1:np) {
      pl = 1 - p[j] / rate[i]
      retlev[i, np - j + 1] <- thresh[i] +
        qegp(p = pl, scale = par[i, 2], shape = par[i, 3], kappa = par[i, 1])
    }
  }
  if (plot) {
    if (length(thresh) > 1L) {
      cols <- hcl.colors(
        n = length(thresh),
        palette = "viridis"
      )
    } else {
      cols <- "black"
    }
    matplot(
      pq,
      t(retlev),
      pch = 20,
      type = "b",
      lty = rep(1, length(thresh)),
      col = cols,
      xlab = "return period",
      ylab = "return level",
      bty = "l"
    )
    text <- switch(
      model,
      "pt-beta" = "Papastathopoulos-Tawn's EGP 1",
      "pt-gamma" = "Papastathopoulos-Tawn's EGP 2",
      "pt-power" = "Papastathopoulos-Tawn's EGP 3 (power)",
      "gj-tnorm" = "Gamet-Jonathan's truncated normal EGP",
      "gj-beta" = "Gamet-Jonathan's beta EGP",
      "logist" = "logistic EGP",
      "exptilt" = "exponential tilting EGP"
    )
    mtext(text, side = 3, cex = 0.9, line = 0.5, adj = 0)
    if (length(thresh) > 1) {
      legend(
        x = "bottomright",
        inset = c(0, 1),
        cex = 0.8,
        xpd = TRUE,
        horiz = TRUE,
        bty = "n",
        legend = thresh,
        col = cols,
        pch = 20
      )
    }
  }
  res <- retlev
  colnames(res) <- pq
  rownames(res) <- thresh
  return(invisible(res))
}

#' Parameter stability plot and maximum likelihood routine for extended GP models
#'
#' \code{fit.egp} is a numerical optimization routine to fit the extended generalised Pareto models of Papastathopoulos and Tawn (2013),
#' using maximum likelihood estimation.
#'
#' @references Papastathopoulos, I. and J. Tawn (2013). Extended generalised Pareto models for tail estimation, \emph{Journal of Statistical Planning and Inference} \bold{143}(3), 131--143.
#' @inheritParams egp
#' @author Leo Belzile
#' @param init vector of initial values, with \eqn{\log(\kappa)}{log(\kappa)} and \eqn{\log(\sigma)}{log(\sigma)}; can be omitted.
#' @return \code{fit.egp} outputs the list returned by \link[stats]{optim}, which contains the parameter values, the hessian and in addition the standard errors
#' @name fit.egp
#' @description The function \code{tstab.egp} provides classical threshold stability plot for (\eqn{\kappa}, \eqn{\sigma}, \eqn{\xi}).
#' The fitted parameter values are displayed with pointwise normal 95\% confidence intervals.
#' The function returns an invisible list with parameter estimates and standard errors, and p-values for the Wald test that \eqn{\kappa=1}.
#'  The plot is for the modified scale (as in the generalised Pareto model) and as such it is possible that the modified scale be negative.
#' \code{tstab.egp} can also be used to fit the model to multiple thresholds.
#' @param plots vector of integers specifying which parameter stability to plot (if any); passing \code{NA} results in no plots
#' @inheritParams egp
#' @param umin optional minimum value considered for threshold (if \code{thresh} is not provided)
#' @param umax optional maximum value considered for threshold (if \code{thresh} is not provided)
#' @param nint optional integer number specifying the number of thresholds to test.
#' @param changepar logical; if \code{TRUE}, the graphical parameters (via a call to \code{par}) are modified.
#' @return \code{tstab.egp} returns a plot(s) of the parameters fit over the range of provided thresholds, with pointwise normal confidence intervals; the function also returns an invisible list containing notably the matrix of point estimates (\code{par}) and standard errors (\code{se}).
#' @importFrom graphics arrows points polygon title
#' @export
#' @examples
#' xdat <- mev::rgp(
#'   n = 100,
#'   loc = 0,
#'   scale = 1,
#'   shape = 0.5)
#' fitted <- fit.egp(
#'   xdat = xdat,
#'   thresh = 1,
#'   model = "egp2",
#'   show = TRUE)
#' thresh <- mev::qgp(seq(0.1, 0.5, by = 0.05), 0, 1, 0.5)
#' tstab.egp(
#'    xdat = xdat,
#'    thresh = thresh,
#'    model = "egp2",
#'    plots = 1:3)
fit.egp <- function(
  xdat,
  thresh,
  model = c(
    "pt-beta",
    "pt-gamma",
    "pt-power",
    "gj-tnorm",
    "gj-beta",
    "exptilt",
    "logist"
  ),
  init,
  show = FALSE
) {
  if (model %in% c("egp1", "egp2", "egp3") & length(model) == 1) {
    model <- switch(
      model,
      egp1 = "pt-beta",
      egp2 = "pt-gamma",
      egp3 = "pt-power"
    )
  }
  model <- match.arg(model)
  if (length(thresh) > 1) {
    warning(
      "Length of threshold vector greater than one. Selecting first component."
    )
    thresh <- thresh[1]
  }
  if (sum(xdat > thresh[1]) < 4L) {
    stop("Not enough observations to fit an extended generalized Pareto model.")
  }
  # If no initial values are provided, fit a GP distribution to obtain them
  changinit <- missing(init)
  if (!changinit) {
    if (!isTRUE(any(init[1:2] < 0))) {
      changinit <- TRUE
    }
  }
  gpfit <- fit.gpd(xdat, threshold = thresh[1], show = FALSE)
  if (changinit) {
    init <- c(
      kappa = ifelse(model %in% c("gj-tnorm", "logist"), 0.01, 1.01),
      suppressWarnings(gpfit$est)
    )
    # Change starting values for boundary cases, otherwise the optimization stalls
    if (init[3] < -0.99) {
      init[3] <- -0.9
    }
  }
  if (
    !is.finite(egp.ll(par = init, xdat = xdat, thresh = thresh, model = model))
  ) {
    stop("Invalid starting parameters.")
  }
  # Keep exceedances only
  xdata <- xdat[xdat > thresh] - thresh
  xmax <- max(xdata)
  mle <- alabama::auglag(
    par = init,
    fn = function(par, xdat, thresh, model) {
      -egp.ll(par = par, xdat = xdata, thresh = thresh, model = model)
    },
    hin = function(par, ...) {
      c(
        par[1] - 1e-10,
        par[2] - 1e-10,
        par[3] + 1,
        ifelse(par[3] < 0, thresh - par[2] / par[3] - xmax, 1)
      )
    },
    xdat = xdata,
    thresh = 0,
    model = model,
    control.outer = list(trace = FALSE, method = "BFGS"),
    control.optim = list(maxit = 500, reltol = 1e-10)
  )
  browser()
  boundary <- FALSE
  if (model %in% c("gj-tnorm", "logist")) {
    if (isTRUE(mle$value <= gpfit$nllh)) {
      boundary <- TRUE
      mle$par <- c(0, coef(gpfit))
      mle$value <- gpfit$nllh
      # Compute standard errors by hand!
      if (model == "gj-tnorm") {
        info_kappa <- function(par) {
          sigma <- par[1]
          xi <- par[2]
          c(
            length(xdata) / 45,
            sum(xdata * (xdata * xi / sigma + 1)^(-2 / xi - 1) / sigma^2),
            -sum(
              (-log(xdata * xi / sigma + 1) /
                xi^2 +
                xdata / (sigma * (xdata * xi / sigma + 1) * xi)) /
                (xdata * xi / sigma + 1)^(2 / xi)
            )
          )
        }
        info_coef <- info_kappa(coef(gpfit))
        mle$hessian <- rbind(
          kappa = info_coef,
          cbind(info_coef[-1], solve(gpfit$vcov))
        )
      }
    }
  }
  fitted <- list()
  fitted$estimate <- fitted$param <- mle$par
  fitted$deviance <- 2 * mle$value
  fitted$nllh <- mle$value
  if (mle$convergence == 0) {
    fitted$convergence <- "successful"
    fitted$vcov <- try(solve(mle$hessian))
    fitted$std.err <- try(sqrt(diag(fitted$vcov)))
    if (inherits(fitted$std.err, what = "try-error") || mle$par[3] < -0.5) {
      fitted$vcov <- NULL
      fitted$se <- rep(NA, 3)
    }
  } else {
    fitted$convergence <- mle$convergence
    warning(
      "Maximization routine may have failed; check output and try providing better starting values"
    )
  }
  names(fitted$estimate) <- names(fitted$std.err) <- c(
    "kappa",
    "scale",
    "shape"
  )
  if (!is.null(fitted$vcov)) {
    colnames(fitted$vcov) <- rownames(fitted$vcov) <- c(
      "kappa",
      "scale",
      "shape"
    )
  }
  fitted$counts <- mle$counts
  fitted$threshold <- thresh
  fitted$nat <- length(xdata)
  fitted$pat <- length(xdata) / length(xdat)
  fitted$exceedances <- xdata
  fitted$hessian <- mle$hessian
  fitted$method <- "copt"
  fitted$model <- model
  class(fitted) <- c("mev_egp")
  if (show) {
    print(fitted)
  }
  return(invisible(fitted))
}

# @param x A fitted object of class \code{mev_gpd}.
# @param digits Number of digits to display in \code{print} call.
# @param ... Additional argument passed to \code{print}.
#' @export
print.mev_egp <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  text <- switch(
    x$model,
    "pt-beta" = "Papastathopoulos-Tawn's EGP 1",
    "pt-gamma" = "Papastathopoulos-Tawn's EGP 2",
    "pt-power" = "Papastathopoulos-Tawn's EGP 3 (power)",
    "gj-tnorm" = "Gamet-Jonathan's truncated normal EGP",
    "gj-beta" = "Gamet-Jonathan's beta EGP",
    "logist" = "logistic EGP",
    "exptilt" = "exponential tilting EGP"
  )
  cat("Model:", text, "\n")
  cat("Deviance:", round(x$deviance, digits), "\n")

  cat("\nThreshold:", round(x$threshold, digits), "\n")
  cat("Number Above:", x$nat, "\n")
  cat("Proportion Above:", round(x$pat, digits), "\n")

  cat("\nEstimates\n")
  print.default(
    format(x$estimate, digits = digits),
    print.gap = 2,
    quote = FALSE,
    ...
  )
  if (!is.na(x$std.err[1]) && x$estimate[3] > -0.5) {
    cat("\nStandard Errors\n")
    print.default(
      format(x$std.err, digits = digits),
      print.gap = 2,
      quote = FALSE,
      ...
    )
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  cat("  Function Evaluations:", x$counts["function"], "\n")
  cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
  cat("\n")
  invisible(x)
}

#' Deprecated function for parameter stability plots
#'
#' @export
#' @keywords internal
egp.fitrange <-
  function(
    xdat,
    thresh,
    model = c("egp1", "egp2", "egp3"),
    plots = 1:3,
    umin,
    umax,
    nint
  ) {
    tstab.egp(
      xdat = xdat,
      thresh = thresh,
      model = model,
      plots = plots,
      umin = umin,
      umax = umax,
      nint = nint
    )
  }

#' @inheritParams egp
#' @param ... additional arguments for the plot function, currently ignored
#' @rdname fit.egp
#' @export
tstab.egp <- function(
  xdat,
  thresh,
  model = c(
    "pt-beta",
    "pt-gamma",
    "pt-power",
    "gj-tnorm",
    "gj-beta",
    "exptilt",
    "logist"
  ),
  plots = 1:3,
  type = c("wald", "profile"),
  umin,
  umax,
  nint,
  changepar = TRUE,
  ...
) {
  if (model %in% c("egp1", "egp2", "egp3") & length(model) == 1) {
    model <- switch(
      model,
      egp1 = "pt-beta",
      egp2 = "pt-gamma",
      egp3 = "pt-power"
    )
  }
  model <- match.arg(model)
  type <- match.arg(type)
  if (missing(thresh) && isTRUE(any(c(missing(umin), missing(umax))))) {
    stop(
      "Must provide either minimum and maximum threshold values, or a vector of threshold \"thresh\"."
    )
  } else if (missing(thresh)) {
    stopifnot(
      inherits(umin, c("integer", "numeric")),
      inherits(umax, c("integer", "numeric")),
      length(umin) == 1,
      length(umax) == 1,
      umin < umax
    )
    thresh <- seq(umin, umax, length = nint)
  } else if (length(thresh) <= 1) {
    stop(
      "Invalid argument\"thresh\" provided;\n please use a vector of threshold candidates of length at least 2"
    )
  }
  pe <- se <- matrix(NA, ncol = 4, nrow = length(thresh))
  conv <- rep(0, length(thresh))
  fit <- try(suppressWarnings(fit.egp(
    xdat = xdat,
    thresh = thresh[1],
    model = model
  )))
  if (!inherits(fit, "try-error")) {
    pe[1, -4] <- fit$param
    colnames(pe) <- colnames(se) <- c(names(fit$param), "modif scale")
    se[1, -4] <- fit$std.err
    conv[1] <- ifelse(is.character(fit$convergence), 0, fit$convergence)
    se[1, 4] <- sqrt(
      cbind(1, -thresh[1]) %*%
        solve(fit$hessian[-1, -1]) %*%
        rbind(1, -thresh[1])
    )[1]
  }
  for (i in 2:length(thresh)) {
    fit <- try(
      suppressWarnings(
        fit.egp(
          xdat = xdat,
          thresh = thresh[i],
          model = model,
          init = pe[i - 1, -4]
        )
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      fit <- try(
        suppressWarnings(
          fit.egp(
            xdat = xdat,
            thresh = thresh[i],
            model = model
          )
        ),
        silent = TRUE
      )
    }
    if (!inherits(fit, "try-error")) {
      pe[i, -4] <- fit$param
      se[i, -4] <- fit$std.err
      conv[i] <- ifelse(is.character(fit$convergence), 0, fit$convergence)
      # Standard error for the modified scale via the delta-method
      se[i, 4] <- sqrt(
        cbind(1, -thresh[i]) %*%
          solve(fit$hessian[-1, -1]) %*%
          rbind(1, -thresh[i])
      )[1]
      # Modify point estimates for the modif scale (all at once)
      pe[, 4] <- pe[, 2] - pe[, 3] * thresh
    }
  }
  # Graphics
  plots <- plots[is.finite(plots)]
  if (!isTRUE(all(conv == 0))) {
    warning(paste("Convergence failed for", sum(conv != 0), "thresholds"))
    plots <- NULL
  }
  if (length(plots) > 0 & !isTRUE(all(is.na(pe)))) {
    plots <- sort(unique(plots))
    if (!isTRUE(all(plots %in% 1:3))) {
      stop(
        "Invalid plot selection. Must be a vector of integers containing indices 1, 2 or 3."
      )
    }
    if (changepar) {
      old.par <- par(no.readonly = TRUE)
      on.exit(par(old.par))
      par(mfrow = c(1, length(plots)), mar = c(4.5, 4.5, 3.1, 0.1))
    }
    for (i in plots) {
      if (i == 2) {
        i <- 4
      } #Get modified scale
      # Plotting devices limits
      ylims = c(
        min(pe[, i]) - qnorm(0.975) * max(se[, i]),
        max(pe[, i]) + qnorm(0.975) * max(se[, i])
      )
      plot(
        x = thresh,
        y = pe[, i],
        pch = 20,
        xlab = "threshold",
        bty = "l",
        ylab = switch(
          i,
          expression(kappa),
          expression(sigma),
          expression(xi),
          expression(tilde(sigma))
        ),
        ylim = ylims,
        type = "n"
      )
      polygon(
        x = c(thresh, rev(thresh)),
        y = c(
          pe[, i] - qnorm(0.975) * se[, i],
          rev(pe[, i] + qnorm(0.975) * se[, i])
        ),
        col = "gray95",
        border = FALSE
      )
      # if (i == min(plots)) {
      # title(
      #   paste0("Parameter stability plots for EGP", substr(model, 4, 4), ""),
      #   outer = FALSE
      # )
      # }
      if (i == max(plots)) {
        text <- switch(
          model,
          "pt-beta" = "Papastathopoulos-Tawn's EGP 1",
          "pt-gamma" = "Papastathopoulos-Tawn's EGP 2",
          "pt-power" = "Papastathopoulos-Tawn's EGP 3 (power)",
          "gj-tnorm" = "Gamet-Jonathan's truncated normal EGP",
          "gj-beta" = "Gamet-Jonathan's beta EGP",
          "logist" = "logistic EGP",
          "exptilt" = "exponential tilting EGP"
        )
        mtext(text, side = 3, cex = 0.9, line = 0.5, adj = 1)
      }
      if (i == 1) {
        abline(h = 1, lwd = 0.5, col = "gray20", lty = 2)
      }
      arrows(
        x0 = thresh,
        y0 = pe[, i] - qnorm(0.975) * se[, i],
        y1 = pe[, i] + qnorm(0.975) * se[, i],
        length = 0.05,
        angle = 90,
        code = 3
      )
      points(x = thresh, y = pe[, i], type = "p", pch = 20)
    }
  }
  pval <- 2 * pnorm(abs(pe[, 1] - 1) / se[, 1], lower.tail = FALSE)
  return(invisible(list(
    par = pe[, -4],
    se = se[, -4],
    model = model,
    pval = pval,
    conv = conv,
    thresh = thresh
  )))
}


pegp.G1 <- function(x, kappa, shape) {
  x <- pmin(1, pmax(0, x))
  stopifnot(length(kappa) == 1L, kappa > 0)
  pbeta(1 - (1 - x)^(abs(shape)), kappa, 1 / abs(shape))
}

pegp.G2 <- function(x, kappa) {
  x <- pmin(1, pmax(0, x))
  stopifnot(length(kappa) == 1L, kappa > 0)
  pgamma(-log(1 - x), shape = kappa) # INCORRECT FORMULA
}

pegp.G3 <- function(x, kappa) {
  x <- pmin(1, pmax(0, x))
  stopifnot(length(kappa) == 1L, kappa > 0)
  x^kappa
}

pegp.G4 <- function(x, kappa, a = 1 / 32) {
  x <- pmin(1, pmax(0, x))
  stopifnot(length(a) == 1L, a > 0, a < 0.5, length(kappa) == 1L, kappa > 0)
  cst <- (pbeta(0.5, kappa, kappa) - pbeta(a, kappa, kappa))
  (pbeta((0.5 - a) * x + a, shape1 = kappa, shape2 = kappa) -
    pbeta(a, kappa, kappa)) /
    cst
}

pegp.G5 <- function(x, kappa) {
  x <- pmin(1, pmax(0, x))
  if (kappa > 0) {
    (pnorm(q = x, mean = 1, sd = 1 / sqrt(kappa)) -
      pnorm(q = 0, mean = 1, sd = 1 / sqrt(kappa))) /
      (0.5 - pnorm(-sqrt(kappa)))
  } else {
    x
  }
}

pegp.G6 <- function(x, kappa) {
  x <- pmin(1, pmax(0, x))
  stopifnot(length(kappa) == 1L, kappa > 0)
  if (!isTRUE(all.equal(as.numeric(kappa), 1))) {
    (kappa^x - 1) / (kappa - 1)
  } else {
    x
  }
}

pegp.G7 <- function(x, kappa) {
  x <- pmin(1, pmax(0, x))
  stopifnot(length(kappa) == 1L, kappa >= 0)
  if (kappa == 0) {
    x
  } else {
    cst <- plogis(0, 1, 1 / kappa)
    (plogis(x, location = 1, scale = 1 / kappa) - cst) / (0.5 - cst)
  }
}

degp.G1 <- function(x, kappa, shape, log = FALSE) {
  # problems with zero and 1
  logdens <- rep(-Inf, length(x))
  logdens[is.na(x)] <- NA
  xval <- which(x < 1 & x > 0)
  if (kappa != 1) {
    logdens[xval] <- log(abs(shape)) +
      (abs(shape) - 1) * log(1 - x[xval]) +
      dbeta(1 - (1 - x[xval])^(abs(shape)), kappa, 1 / abs(shape), log = TRUE)
  } else {
    logdens[xval] <- 0
  }
  if (isTRUE(log)) {
    return(logdens)
  } else {
    return(exp(logdens))
  } # TODO limit when shape = 0
}

degp.G2 <- function(x, kappa, log = FALSE) {
  logdens <- rep(-Inf, length(x))
  logdens[is.na(x)] <- NA
  xval <- which(x < 1 & x >= 0)
  logdens[xval] <- dgamma(-log(1 - x[xval]), shape = kappa, log = TRUE) -
    log(1 - x[xval])
  if (isTRUE(log)) {
    return(logdens)
  } else {
    return(exp(logdens))
  }
}

degp.G3 <- function(x, kappa, log = FALSE) {
  logdens <- rep(-Inf, length(x))
  logdens[is.na(x)] <- NA
  xval <- which(x <= 1 & x >= 0)
  if (kappa != 1) {
    logdens[xval] <- log(kappa) + (kappa - 1) * log(x)
  } else {
    logdens[xval] <- 0
  }
  if (isTRUE(log)) {
    return(logdens)
  } else {
    return(exp(logdens))
  }
}

degp.G4 <- function(x, kappa, a = 1 / 32, log = FALSE) {
  logdens <- rep(-Inf, length(x))
  logdens[is.na(x)] <- NA
  xval <- which(x <= 1 & x >= 0)
  stopifnot(
    length(a) == 1L,
    a > 0,
    a < 0.5,
    length(kappa) == 1L,
    kappa > 0
  )
  logdens[xval] <- log(0.5 - a) +
    dbeta((0.5 - a) * x[xval] + a, shape1 = kappa, shape2 = kappa, log = TRUE) -
    log(pbeta(0.5, kappa, kappa) - pbeta(a, kappa, kappa))
  if (isTRUE(log)) {
    return(logdens)
  } else {
    return(exp(logdens))
  }
}

degp.G5 <- function(x, kappa, log = FALSE) {
  logdens <- rep(-Inf, length(x))
  logdens[is.na(x)] <- NA
  xval <- which(x <= 1 & x >= 0)
  stopifnot(kappa >= 0)
  if (kappa > 0) {
    logdens[xval] <- dnorm(
      x = x[xval],
      mean = 1,
      sd = 1 / sqrt(kappa),
      log = TRUE
    ) -
      log(0.5 - pnorm(-sqrt(kappa)))
  } else if (kappa == 0) {
    logdens[xval] <- 0
  }
  if (isTRUE(log)) {
    return(logdens)
  } else {
    return(exp(logdens))
  }
}

degp.G6 <- function(x, kappa, log = FALSE) {
  logdens <- rep(-Inf, length(x))
  logdens[is.na(x)] <- NA
  xval <- which(x <= 1 & x >= 0)
  stopifnot(length(kappa) == 1L, kappa > 0)
  if (!isTRUE(all.equal(as.numeric(kappa), 1))) {
    logdens[xval] <- x[xval] *
      log(kappa) +
      log(abs(log(kappa))) -
      log(abs(kappa - 1))
  } else {
    logdens[xval] <- 0
  }
  if (isTRUE(log)) {
    return(logdens)
  } else {
    return(exp(logdens))
  }
}

degp.G7 <- function(x, kappa, log = FALSE) {
  logdens <- rep(-Inf, length(x))
  logdens[is.na(x)] <- NA
  xval <- which(x <= 1 & x >= 0)
  if (kappa > 0) {
    logdens[xval] <- dlogis(x, location = 1, scale = 1 / kappa, log = TRUE) -
      log(0.5 - plogis(0, location = 1, scale = 1 / kappa))
  } else if (kappa == 0) {
    logdens[xval] <- 0
  }
  if (isTRUE(log)) {
    return(logdens)
  } else {
    return(exp(logdens))
  }
}


qegp.G1 <- function(x, kappa, shape) {
  x <- pmin(1, pmax(0, x))
  1 - (1 - qbeta(x, kappa, 1 / abs(shape)))^(1 / (abs(shape)))
}

qegp.G2 <- function(x, kappa) {
  x <- pmin(1, pmax(0, x))
  1 - exp(-qgamma(x, shape = kappa))
}

qegp.G3 <- function(x, kappa) {
  x <- pmin(1, pmax(0, x))
  x^(1 / kappa)
}

qegp.G4 <- function(x, kappa, a = 1 / 32) {
  x <- pmin(1, pmax(0, x))
  cst <- pbeta(0.5, kappa, kappa) - pbeta(a, kappa, kappa)
  (qbeta(x * cst + pbeta(a, kappa, kappa), kappa, kappa) - a) / (0.5 - a)
}

qegp.G5 <- function(x, kappa, a = 1 / 32) {
  x <- pmin(1, pmax(0, x))
  cst <- 0.5 - pnorm(-sqrt(kappa))
  qnorm(
    x * cst + pnorm(q = 0, mean = 1, sd = 1 / sqrt(kappa)),
    mean = 1,
    sd = 1 / sqrt(kappa)
  )
}

qegp.G6 <- function(x, kappa) {
  x <- pmin(1, pmax(0, x))
  if (kappa != 1) {
    log1p(x * (kappa - 1)) / log(kappa)
  } else {
    x
  }
}

qegp.G7 <- function(x, kappa, a = 1 / 32) {
  x <- pmin(1, pmax(0, x))
  stopifnot(length(kappa) == 1L, kappa >= 0)
  if (kappa > 0) {
    cst <- plogis(0, 1, scale = 1 / kappa)
    qlogis(x * (0.5 - cst) + cst, location = 1, scale = 1 / kappa)
  } else if (kappa == 0) {
    x
  }
}


pegp <- function(
  q,
  scale,
  shape,
  kappa,
  model = c(
    "pt-beta",
    "pt-gamma",
    "pt-power",
    "gj-tnorm",
    "gj-beta",
    "exptilt",
    "logist"
  ),
  lower.tail = TRUE,
  log.p = FALSE
) {
  model <- match.arg(model)
  stopifnot(
    length(kappa) == 1L,
    length(scale) == 1L,
    length(shape) == 1L,
    kappa > 0,
    scale > 0
  )
  if (model %in% c("pt-beta", "pt-gamma") & abs(shape) < 1e-8) {
    return(pgamma(
      q = q,
      scale = scale,
      shape = kappa,
      lower.tail = lower.tail,
      log.p = log.p
    ))
  }
  p <- pgp(q, loc = 0, scale = scale, shape = shape, lower.tail = lower.tail)
  pg <- switch(
    model,
    "pt-beta" = pegp.G1(pg, kappa = kappa, shape = shape),
    "pt-gamma" = pegp.G2(pg, kappa = kappa),
    "pt-power" = pegp.G3(pg, kappa = kappa),
    "gj-tnorm" = pegp.G5(pg, kappa = kappa),
    "gj-beta" = pegp.G4(pg, kappa = kappa),
    "exptilt" = pegp.G6(pg, kappa = kappa),
    "logist" = pegp.G7(pg, kappa = kappa)
  )
  if (!isTRUE(log.p)) {
    return(pg)
  } else {
    return(log(pg))
  }
}

degp <- function(
  x,
  scale,
  shape,
  kappa,
  model = c(
    "pt-beta",
    "pt-gamma",
    "pt-power",
    "gj-tnorm",
    "gj-beta",
    "exptilt",
    "logist"
  ),
  log = FALSE
) {
  stopifnot(
    length(kappa) == 1L,
    length(scale) == 1L,
    length(shape) == 1L,
    kappa >= 0,
    scale > 0
  )
  model <- match.arg(model)
  pg <- pgp(q = x, scale = scale, shape = shape)
  logdens <- switch(
    model,
    "pt-beta" = degp.G1(pg, kappa = kappa, shape = shape, log = TRUE),
    "pt-gamma" = degp.G2(pg, kappa = kappa, log = TRUE),
    "pt-power" = degp.G3(pg, kappa = kappa, log = TRUE),
    "gj-tnorm" = degp.G5(pg, kappa = kappa, log = TRUE),
    "gj-beta" = degp.G4(pg, kappa = kappa, log = TRUE),
    "exptilt" = degp.G6(pg, kappa = kappa, log = TRUE),
    "logist" = degp.G7(pg, kappa = kappa, log = TRUE)
  ) +
    dgp(x = x, scale = scale, shape = shape, log = TRUE)
  if (isTRUE(log)) {
    return(logdens)
  } else {
    return(exp(logdens))
  }
}

qegp <- function(
  p,
  scale,
  shape,
  kappa,
  model = c(
    "pt-beta",
    "pt-gamma",
    "pt-power",
    "gj-tnorm",
    "gj-beta",
    "exptilt",
    "logist"
  ),
  lower.tail = TRUE
) {
  stopifnot(
    length(kappa) == 1L,
    length(scale) == 1L,
    length(shape) == 1L,
    kappa >= 0,
    scale > 0
  )
  model <- match.arg(model)
  if (!isTRUE(lower.tail)) {
    p <- 1 - p
  }
  qu <- rep(NA, length(p))
  vals <- is.finite(p) & p >= 0 & p <= 1
  if (model %in% c("pt-beta", "pt-gamma") & abs(shape) < 1e-8) {
    qu[vals] <- qgamma(
      p = p[vals],
      scale = scale,
      shape = kappa,
      lower.tail = lower.tail
    )
  } else {
    qu[vals] <- qgp(
      switch(
        model,
        "pt-beta" = qegp.G1(p[vals], kappa = kappa, shape = shape),
        "pt-gamma" = qegp.G2(p[vals], kappa = kappa),
        "pt-power" = qegp.G3(p[vals], kappa = kappa),
        "gj-tnorm" = qegp.G5(p[vals], kappa = kappa),
        "gj-beta" = qegp.G4(p[vals], kappa = kappa),
        "exptilt" = qegp.G6(p[vals], kappa = kappa),
        "logist" = qegp.G7(p[vals], kappa = kappa)
      ),
      scale = scale,
      shape = shape
    )
  }
  return(qu)
}


regp <- function(
  n,
  scale,
  shape,
  kappa,
  model = c(
    "pt-beta",
    "pt-gamma",
    "pt-power",
    "gj-tnorm",
    "gj-beta",
    "exptilt",
    "logist"
  )
) {
  stopifnot(
    length(kappa) == 1L,
    length(scale) == 1L,
    length(shape) == 1L,
    kappa >= 0,
    scale > 0
  )
  model <- match.arg(model)
  pg <- switch(
    model,
    "pt-beta" = qegp.G1(runif(n), kappa = kappa, shape = shape),
    "pt-gamma" = qegp.G2(runif(n), kappa = kappa),
    "pt-power" = qegp.G3(runif(n), kappa = kappa),
    "gj-tnorm" = qegp.G5(runif(n), kappa = kappa),
    "gj-beta" = qegp.G4(runif(n), kappa = kappa),
    "exptilt" = qegp.G6(runif(n), kappa = kappa),
    "logist" = qegp.G7(runif(n), kappa = kappa)
  )
}


egp.pll <- function(
  psi,
  model = c(
    "pt-beta",
    "pt-gamma",
    "pt-power",
    "gj-tnorm",
    "gj-beta",
    "exptilt",
    "logist"
  ),
  param = c("shape", "kappa"),
  mle = NULL,
  xdat,
  threshold = NULL,
  plot = FALSE
) {
  param <- match.arg(param)
  if (is.null(threshold)) {
    threshold <- 0
  } else {
    stopifnot(
      is.numeric(threshold),
      length(threshold) == 1L
    )
    stopifnot(threshold >= 0)
    xdat <- xdat[xdat > threshold] - threshold
  }
  if (is.null(mle)) {
    mle <- try(
      fit.egp(xdat = xdat, thresh = threshold, model = model),
      silent = TRUE
    )
    if (inherits(mle, "try-error")) {
      stop("Could not find maximum likelihood estimator.")
    }
  } else {
    stopifnot(inherits(mle, what = "mev_egp"))
  }
  ind <- switch(param, kappa = 1L, shape = 3L)
  coef_mle <- coef(mle)[ind]
  se_mle <- sqrt(diag(vcov(mle))[ind])
  if (missing(psi)) {
    if (is.numeric(se_mle) & is.finite(se_mle)) {
      psi <- seq(-3 * se_mle, 3 * se_mle, length.out = 55) + coef_mle
    } else {
      stop("Could not determine a suitable sequence of values for profiling.")
    }
  } else {
    psi <- sort(unique(c(psi, coef_mle)))
  }
  if (param == "kappa") {
    psi <- psi[psi > 0]
  } else if (param == "shape") {
    psi <- psi[psi > -1]
  }
  mid <- which(psi == coef_mle)
  # should be 28 for seq of psi if no zero values
  pars <- matrix(nrow = length(psi), ncol = 3L)
  pll <- numeric(length(psi))
  pll[mid] <- -mle$nllh
  pars[mid, ] <- coef(mle)
  xmax <- max(xdat)
  if (param == "kappa") {
    nll_egp_kappa <- function(eta, psi, xdat, model) {
      -egp.ll(
        xdat = xdat,
        thresh = 0,
        par = c(psi, eta),
        model = model
      )
    }
    for (i in seq_len(mid - 1)) {
      fit_const <- alabama::auglag(
        par = pars[mid - i + 1, -1],
        fn = nll_egp_kappa,
        hin = function(par, ...) {
          c(
            par[1] - 1e-10,
            par[2] + 1,
            ifelse(par[2] < 0, -par[1] / par[2] - xmax, 1)
          )
        },
        xdat = xdat,
        psi = psi[mid - i],
        model = model,
        control.outer = list(trace = FALSE, method = "nlminb")
      )
      pars[mid - i, ] <- c(psi[mid - i], fit_const$par)
      pll[mid - i] <- -fit_const$value
    }
    for (i in seq_len(length(psi) - mid)) {
      fit_const <- alabama::auglag(
        par = pars[mid + i - 1, -1],
        fn = nll_egp_kappa,
        hin = function(par, ...) {
          c(
            par[1] - 1e-10,
            par[2] + 1,
            ifelse(par[2] < 0, -par[1] / par[2] - xmax, 1)
          )
        },
        xdat = xdat,
        psi = psi[mid + i],
        model = model,
        control.outer = list(trace = FALSE, method = "nlminb")
      )
      pars[mid + i, ] <- c(psi[mid + i], fit_const$par)
      pll[mid + i] <- -fit_const$value
    }
  } else if (param == "shape") {
    nll_egp_shape <- function(eta, psi, xdat, model) {
      -egp.ll(
        xdat = xdat,
        thresh = 0,
        par = c(exp(eta), psi),
        model = model
      )
    }
    for (i in seq_len(mid - 1)) {
      fit_const <- alabama::auglag(
        par = log(pars[mid - i + 1, -3]),
        fn = nll_egp_shape,
        hin = function(par, psi, ...) {
          # TODO check if this is correct
          c(ifelse(psi < 0, -exp(par[2]) / psi - xmax, 1))
        },
        xdat = xdat,
        psi = psi[mid - i],
        model = model,
        control.outer = list(trace = FALSE, method = "nlminb")
      )
      pars[mid - i, ] <- c(exp(fit_const$par), psi[mid - i])
      pll[mid - i] <- -fit_const$value
    }
    for (i in seq_len(length(psi) - mid)) {
      fit_const <- alabama::auglag(
        par = log(pars[mid + i - 1, -3]),
        fn = nll_egp_shape,
        hin = function(par, psi, ...) {
          c(ifelse(psi < 0, -exp(par[2]) / psi - xmax, 1))
        },
        xdat = xdat,
        psi = psi[mid + i],
        model = model,
        control.outer = list(trace = FALSE, method = "nlminb")
      )
      pars[mid + i, ] <- c(exp(fit_const$par), psi[mid + i])
      pll[mid + i] <- -fit_const$value
    }
  }
  ans <- list(
    mle = coef(mle),
    psi.max = coef_mle,
    param = param,
    std.error = se_mle,
    psi = psi,
    pll = pll,
    maxpll = -mle$nllh,
    family = "egp",
    threshold = threshold
  )
  ans$r <- sign(coef_mle - psi) * sqrt(2 * (ans$maxpll - ans$pll))
  ans$normal <- c(ans$psi.max, ans$std.error)
  class(ans) <- "eprof"
  if (isTRUE(plot)) {
    plot(ans)
  }
  return(invisible(ans))
}
