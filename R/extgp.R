#' Extended generalised Pareto families
#'
#' @description This function provides the log-likelihood and quantiles for the three different families presented
#' in Papastathopoulos and Tawn (2013) and the two proposals of Gamet and Jalbert (2022), plus exponential tilting. All of the models contain an additional parameter, \eqn{\kappa \ge 0}.
#' All families share the same tail index as the generalized Pareto distribution, while allowing for lower thresholds.
#' For most models, the distribution reduce to the generalised Pareto when \eqn{\kappa=1} (for models \code{gj-tnorm} and \code{logist}, on the boundary of the parameter space when \eqn{\kappa \to 0}).
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
#' @importFrom grDevices hcl.colors
#'
#' @details For return levels, the \code{p} argument can be related to \eqn{T} year exceedances as follows:
#' if there are \eqn{n_y} observations per year, than take \code{p}
#' to equal \eqn{1/(Tn_y)} to obtain the \eqn{T}-years return level.
#' @author Leo Belzile
#' @return \code{egp.ll} returns the log-likelihood value, while \code{egp.retlev} returns a plot of the return levels if \code{plot=TRUE} and a list with tail probabilities \code{p}, return levels \code{retlev}, thresholds \code{thresh} and model name \code{model}.
#' @examples
#' set.seed(123)
#' xdat <- rgp(1000, loc = 0, scale = 2, shape = 0.5)
#' par <- fit.egp(xdat, thresh = 0, model = 'gj-beta')$par
#' p <- c(1/1000, 1/1500, 1/2000)
#' # With multiple thresholds
#' th <- c(0, 0.1, 0.2, 1)
#' opt <- tstab.egp(xdat, thresh = th, model = 'gj-beta')
#' egp.retlev(xdat = xdat, thresh = th, model = 'gj-beta', p = p)
#' opt <- tstab.egp(xdat, th, model = 'pt-power', plots = NA)
#' egp.retlev(xdat = xdat, thresh = th, model = 'pt-power', p = p)
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
  if (!missing(thresh)) {
    if (!isTRUE(all.equal(thresh, 0))) {
      xdat = xdat[xdat > thresh] - thresh
    }
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
  confint = FALSE,
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
  if (missing(par)) {
    par <- matrix(ncol = 3, nrow = length(thresh))
    for (i in seq_along(thresh)) {
      par[i, ] <- coef(
        fit.egp(
          xdat = xdat,
          thresh = thresh[i],
          model = model
        )
      )
    }
  } else {
    if (length(par) %% 3 != 0) {
      stop("Invalid parameter input")
    }
    if (!inherits(par, "matrix")) {
      par <- matrix(data = as.numeric(par), ncol = 3)
    }
    stopifnot(length(par) == (length(thresh) * 3L))
  }
  rate <- sapply(thresh, function(u) {
    length(xdat[xdat > u]) / length(xdat)
  })
  if (!isTRUE(all.equal(length(rate), length(thresh), nrow(par)))) {
    stop("Input dimension does not match")
  }
  retlev <- matrix(
    data = 0,
    nrow = length(thresh),
    ncol = length(p)
  )
  if (
    any(sapply(rate, function(zeta) {
      zeta < p
    }))
  ) {
    warning(
      "Some probabilities \"p\" are higher than the exceedance rate. Evaluate those empirically."
    )
  }
  p <- sort(p, decreasing = TRUE)
  pq <- rev(1 / p)
  np <- length(p)
  for (i in 1:length(thresh)) {
    for (j in 1:np) {
      pl = 1 - p[j] / rate[i]
      retlev[i, j] <- thresh[i] +
        qegp(
          p = pl,
          scale = par[i, 2],
          shape = par[i, 3],
          kappa = par[i, 1],
          model = model
        )
    }
  }
  if (length(thresh) == 1L) {
    retlev <- c(retlev)
    names(retlev) <- round(1 - p, 4)
  } else {
    rownames(retlev) <- thresh
    colnames(retlev) <- round(1 - p, 4)
  }
  obj <- list(
    p = p,
    thresh = thresh,
    retlev = retlev,
    model = model
  )
  class(obj) <- c("mev_egp_retlev")
  if (plot) {
    plot(obj)
  }
  return(invisible(obj))
}

#' @export
plot.mev_egp_retlev <- function(x, ...) {
  p <- x$p
  thresh <- x$thresh
  retlev <- x$retlev

  if (length(thresh) > 1L) {
    cols <- hcl.colors(
      n = length(thresh),
      palette = "viridis"
    )
  } else {
    cols <- "black"
  }
  matplot(
    rev(1 / p),
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
    x$model,
    "pt-beta" = "Papastathopoulos-Tawn's EGP 1",
    "pt-gamma" = "Papastathopoulos-Tawn's EGP 2",
    "pt-power" = "Papastathopoulos-Tawn's EGP 3 (power)",
    "gj-tnorm" = "Gamet-Jalbert's truncated normal EGP",
    "gj-beta" = "Gamet-Jalbert's beta EGP",
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

#' Parameter stability plot and maximum likelihood routine for extended GP models
#'
#' \code{fit.egp} is a numerical optimization routine to fit the extended generalised Pareto models of Papastathopoulos and Tawn (2013),
#' using maximum likelihood estimation.
#'
#' @references Papastathopoulos, I. and J. Tawn (2013). Extended generalised Pareto models for tail estimation, \emph{Journal of Statistical Planning and Inference} \bold{143}(3), 131--143.
#' @inheritParams egp
#' @author Leo Belzile
#' @param start optional named list of initial values, with \eqn{\kappa}{\kappa}, \eqn{sigma}{\sigma} or \eqn{xi}{\xi}.
#' @return \code{fit.egp} outputs the list returned by \link[stats]{optim}, which contains the parameter values, the hessian and in addition the standard errors
#' @name fit.egp
#' @description The function \code{tstab.egp} provides classical threshold stability plot for (\eqn{\kappa}, \eqn{\sigma}, \eqn{\xi}).
#' The fitted parameter values are displayed with pointwise normal 95\% confidence intervals.
#' The function returns an invisible list with parameter estimates and standard errors, and p-values for the Wald test that \eqn{\kappa=1}.
#'  The plot is for the modified scale (as in the generalised Pareto model) and as such it is possible that the modified scale be negative.
#' \code{tstab.egp} can also be used to fit the model to multiple thresholds.
#' @inheritParams egp
#' @inheritParams fit.gpd
#' @param ... additional parameters, for backward compatibility purposes
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
#'   model = "pt-gamma",
#'   show = TRUE)
#' thresh <- mev::qgp(seq(0.1, 0.5, by = 0.05), 0, 1, 0.5)
#' tstab.egp(
#'    xdat = xdat,
#'    thresh = thresh,
#'    model = "pt-gamma")
#' xdat <- regp(
#'   n = 100,
#'   scale = 1,
#'   shape = 0.1,
#'   kappa = 0.5,
#'   model = "pt-power"
#' )
#' fit.egp(
#'  xdat = xdat,
#'  model = "pt-power",
#'  show = TRUE,
#'  fpar = list(kappa = 1),
#'  method = "Nelder"
#' )
fit.egp <- function(
  xdat,
  thresh = 0,
  model = c(
    "pt-beta",
    "pt-gamma",
    "pt-power",
    "gj-tnorm",
    "gj-beta",
    "exptilt",
    "logist"
  ),
  start = NULL,
  method = c("Nelder", "nlminb", "BFGS"),
  fpar = NULL,
  show = FALSE,
  ...
) {
  # Include backward compatibility for init
  args <- list(...)
  if (!is.null(args$init)) {
    stopifnot(length(args$init) == 2L)
    if (is.null(start)) {
      start <- list(
        kappa = exp(args$init[1]),
        scale = exp(args$init[2]),
        shape = 0.1
      )
    }
  }
  # Backward compatibility for model names
  if (model %in% c("egp1", "egp2", "egp3") & length(model) == 1) {
    model <- switch(
      model,
      egp1 = "pt-beta",
      egp2 = "pt-gamma",
      egp3 = "pt-power"
    )
  }
  model <- match.arg(
    model,
    choices = c(
      "pt-beta",
      "pt-gamma",
      "pt-power",
      "gj-tnorm",
      "gj-beta",
      "exptilt",
      "logist"
    ),
    several.ok = FALSE
  )
  if (length(thresh) > 1) {
    warning(
      "Length of threshold vector greater than one. Selecting first component."
    )
    thresh <- thresh[1]
  }
  stopifnot(thresh >= 0)
  if (sum(xdat > thresh) < 5L) {
    stop("Not enough observations to fit an extended generalized Pareto model.")
  }
  xdat <- as.numeric(xdat[is.finite(xdat)])
  # Keep exceedances only
  xdata <- xdat[xdat > thresh] - thresh
  xmax <- max(xdata)
  # Fit submodel to check convergence afterwards
  gpfit <- fit.gpd(xdata, threshold = 0, show = FALSE)
  # Can also be used for initial values in case
  # there were not provided by the user
  param_names <- c("kappa", "scale", "shape")
  stopifnot(is.null(fpar) | is.list(fpar))
  wf <- (param_names %in% names(fpar))
  if (sum(wf) == 3L) {
    # TODO turn this into a warning and evaluate nllh?
    stop("Invalid input: all of the model parameters are fixed.")
  }
  if (is.list(fpar) && (length(fpar) >= 1L)) {
    #NULL has length zero
    if (is.null(names(fpar))) {
      stop("\"fpar\" must be a named list")
    }
    if (!isTRUE(all(names(fpar) %in% param_names))) {
      stop(
        "Unknown fixed parameter: must be one of \"kappa\",\"scale\" or \"shape\". "
      )
    }
    if (!isTRUE(all(unlist(lapply(fpar, length)) == rep(1L, sum(wf))))) {
      stop("Each fixed parameter must be of length one.")
    }
  }
  method <- match.arg(method)
  if (is.null(start)) {
    spar <- c(
      kappa = ifelse(model %in% c("gj-tnorm", "logist"), 0.01, 1.01),
      suppressWarnings(gpfit$est)
    )
    # Change starting values for boundary cases, otherwise the optimization stalls
    if (spar[3] < -0.99) {
      spar[3] <- -0.9
    }
    names(spar) <- param_names
  } else {
    stopifnot(length(start) == (3L - sum(wf)))
    spar <- vector(mode = "numeric", length = 3L)
    names(spar) <- param_names
    if (is.null(names(start))) {
      spar[!wf] <- unlist(start) # assume order, for better or worse
    } else {
      stopifnot(isTRUE(all(names(start) %in% param_names)))
      for (name in names(start)) {
        spar[name] <- unlist(start[name])
      }
    }
  }
  for (i in seq_along(fpar)) {
    spar[names(fpar[i])] <- unlist(fpar[i])[1]
  }
  stopifnot(
    spar[1] > 0,
    spar[2] > 0,
    spar[3] > -1 - 1e-8
  )
  if (
    !is.finite(
      egp.ll(
        par = spar,
        xdat = xdata,
        thresh = 0,
        model = model
      )
    )
  ) {
    stop("Invalid starting parameters.")
  }
  start_vals <- spar[!wf]
  fixed_vals <- spar[wf] #when empty, a num vector of length zero
  wfo <- order(c(which(!wf), which(wf)))
  if (method != "Nelder") {
    mle <- try(suppressWarnings(
      alabama::auglag(
        par = start_vals,
        fpar = fixed_vals,
        wfixed = wf,
        wfo = wfo,
        fn = function(par, fpar, wfixed, wfo) {
          params <- c(par, fpar)[wfo]
          nll <- -egp.ll(params, xdat = xdata, thresh = 0, model = model)
          ifelse(is.finite(nll), nll, 1e10)
        },
        hin = function(par, fpar, wfixed, wfo) {
          params <- c(par, fpar)[wfo]
          # TODO check whether parameters are
          # constrained to be positive
          c(
            params[1],
            params[2],
            params[3] + 1,
            params[2] + params[3] * xmax
          )
        },
        control.outer = list(method = method, trace = FALSE),
        control.optim = switch(
          method,
          nlminb = list(
            iter.max = 500L,
            rel.tol = 1e-10,
            step.min = 1e-10
          ),
          list(maxit = 1000L, reltol = 1e-10)
        )
      )
    ))
  } else {
    mle <- try(suppressWarnings(
      optim(
        par = start_vals,
        fpar = fixed_vals,
        method = "Nelder",
        hessian = TRUE,
        wfixed = wf,
        wfo = wfo,
        fn = function(par, fpar, wfixed, wfo) {
          params <- c(par, fpar)[wfo]
          nll <- -egp.ll(params, xdat = xdata, thresh = 0, model = model)
          ifelse(is.finite(nll), nll, 1e10)
        }
      )
    ))
  }
  boundary <- FALSE
  if (model %in% c("gj-tnorm", "logist")) {
    if (isTRUE(mle$value >= gpfit$nllh) & sum(wf) == 0L) {
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
  fitted$estimate <- mle$par
  fitted$param <- c(mle$par, spar[wf])[wfo]
  fitted$deviance <- 2 * mle$value
  fitted$nllh <- mle$value
  if (
    mle$convergence == 0 |
      isTRUE(mle$kkt1 & mle$kkt2)
  ) {
    fitted$convergence <- "successful"
    fitted$vcov <- try(
      expr = solve(mle$hessian),
      silent = TRUE
    )
    fitted$std.err <- try(
      expr = sqrt(diag(fitted$vcov)),
      silent = TRUE
    )
    if (
      inherits(fitted$std.err, what = "try-error") || fitted$param[3] < -0.5
    ) {
      fitted$vcov <- NULL
      fitted$std.err <- rep(NA, 3)
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
  )[!wf]
  if (!is.null(fitted$vcov)) {
    colnames(fitted$vcov) <- rownames(fitted$vcov) <- c(
      "kappa",
      "scale",
      "shape"
    )[!wf]
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
print.mev_egp <- function(x, digits = min(3, getOption("digits") - 3), ...) {
  text <- switch(
    x$model,
    "pt-beta" = "Papastathopoulos-Tawn's EGP 1",
    "pt-gamma" = "Papastathopoulos-Tawn's EGP 2",
    "pt-power" = "Papastathopoulos-Tawn's EGP 3 (power)",
    "gj-tnorm" = "Gamet-Jalbert's truncated normal EGP",
    "gj-beta" = "Gamet-Jalbert's beta EGP",
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
  if (!is.na(x$std.err[1])) {
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

#' Threshold stability plots for extended generalized Pareto models
#' @inheritParams fit.egp
#' @param transform logical; if \code{TRUE} and \code{type="wald"}, intervals for \code{kappa} are computed on the log-scale and back-transformed.
#' @param ... additional arguments for the plot function, currently ignored.
#' @param level [double] confidence interval level, default to 0.95.
#' @param type [string] confidence interval type, either \code{wald} or \code{profile}.
#' @param param [string] parameter, either \code{shape} or additional parameter \code{kappa}
#' @param plot [logical] if \code{TRUE} (default), return a threshold stability plot
#' @export
#' @examples
#' xdat <- rgp(n = 1000)
#' tstab.egp(
#'  xdat = xdat,
#'  thresh = c(0, quantile(xdat, 0.5)),
#'  model = "gj-tnorm",
#'  param = "kappa",
#'  transform = TRUE)
#'
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
  param = c("shape", "kappa"),
  type = c("wald", "profile"),
  transform = FALSE,
  level = 0.95,
  plot = TRUE,
  ...
) {
  args <- list(...)
  plots <- args$plots
  changepar <- args$changepar
  if (is.null(changepar)) {
    changepar <- TRUE
  } else {
    changepar <- isTRUE(args$changepar)
  }
  if (missing(thresh)) {
    if (!is.null(args$umin) & !is.null(args$umax) & !is.null(args$nint)) {
      umin <- args$umin
      umax <- args$umax
      nint <- args$nint
    } else {
      stop("Threshold vector not provided")
    }
  }
  model <- match.arg(
    model,
    choices = c(
      "pt-beta",
      "pt-gamma",
      "pt-power",
      "gj-tnorm",
      "gj-beta",
      "exptilt",
      "logist",
      "egp1",
      "egp2",
      "egp3"
    ),
    several.ok = FALSE
  )

  if (model %in% c("egp1", "egp2", "egp3")) {
    model <- switch(
      model,
      egp1 = "pt-beta",
      egp2 = "pt-gamma",
      egp3 = "pt-power"
    )
  }
  type <- match.arg(type)
  if (is.numeric(plots)) {
    plots <- as.integer(plots)
  } else {
    plots <- c()
    plots <- c(1L, 3L)[c("kappa", "shape") %in% param]
  }
  transform <- isTRUE(as.logical(transform)[1])
  # Option for backward compatibility
  if (2L %in% plots) {
    warning("Modified scale not available for EGPD models.")
  }
  plots <- sort(plots[plots %in% c(1L, 3L)])
  if (length(plots) == 0) {
    stop("No option for plots")
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
  if (1L %in% plots) {
    kappa_pars <- matrix(
      data = NA,
      nrow = length(thresh),
      ncol = 3
    )
  } else {
    kappa_pars <- NULL
  }
  if (3L %in% plots) {
    shape_pars <- matrix(
      data = NA,
      nrow = length(thresh),
      ncol = 3
    )
  } else {
    shape_pars <- NULL
  }
  for (i in seq_along(thresh)) {
    mle <- try(
      fit.egp(
        xdat = xdat,
        thresh = thresh[i],
        model = model,
        method = "BFGS"
      ),
      silent = TRUE
    )
    if (inherits(mle, "try-error")) {
      next
    } else {
      est <- coef(mle)
      if (type == "wald") {
        if (1L %in% plots) {
          if (transform) {
            est[1] <- log(est[1] + 1e-8)
            fn_numderiv <- function(par) {
              -egp.ll(
                par = c(exp(par[1]), par[2], par[3]),
                xdat = xdat,
                thresh = thresh[i],
                model = model
              )
            }
            opt <- optim(
              par = est,
              fn = fn_numderiv,
              method = "Nelder",
              hessian = TRUE
            )
            se <- try(
              expr = suppressWarnings(
                sqrt(diag(solve(opt$hessian)))[1]
              ),
              silent = TRUE
            )
            if (inherits(se, "try-error")) {
              se <- NA
            }
          } else {
            se <- mle$std.err[1]
          }
          kappa_pars[i, 1] <- coef(mle)[1]
          if (
            isTRUE(all.equal(kappa_pars[i, 1], 0, check.attributes = FALSE))
          ) {
            boundary <- TRUE
            crit <- c(-1, 1) * sqrt(0.5 * qchisq(level, df = 1))
          } else {
            boundary <- FALSE
            crit <- qnorm(c((1 - level) / 2, 1 - (1 - level) / 2))
          }
          kappa_pars[i, 2:3] <- est[1] + crit * se
          if (transform) {
            kappa_pars[i, 2:3] <- exp(kappa_pars[i, 2:3])
          }
        }
        if (3L %in% plots) {
          shape_pars[i, 1] <- coef(mle)[3]
          shape_pars[i, 2:3] <- coef(mle)[3] +
            qnorm(c((1 - level) / 2, 1 - (1 - level) / 2)) * mle$std.err[3]
        }
      } else if (type == "profile") {
        if (1L %in% plots) {
          prof_kappa <- try(
            egp.pll(
              model = model,
              param = "kappa",
              mle = mle,
              plot = FALSE,
              thresh = thresh[i],
              xdat = xdat
            ),
            silent = TRUE
          )
          if (!inherits(prof_kappa, "try-error")) {
            boundary <- ifelse(
              isTRUE(all.equal(
                pmax(0, coef(mle)[1]),
                0,
                check.attributes = FALSE
              )),
              TRUE,
              FALSE
            )
            kappa_pars[i, ] <- confint(
              prof_kappa,
              level = level,
              boundary = boundary
            )
            if (boundary) {
              kappa_pars[i, 1:2] <- 0
            }
          }
        }
        if (3L %in% plots) {
          prof_shape <- try(
            egp.pll(
              model = model,
              param = "shape",
              mle = mle,
              plot = FALSE,
              thresh = thresh[i],
              xdat = xdat
            ),
            silent = TRUE
          )
          if (!inherits(prof_shape, "try-error")) {
            shape_pars[i, ] <- confint(prof_kappa, level = level)
          }
        }
      }
    }
  }
  if (!is.null(kappa_pars)) {
    colnames(kappa_pars) <- c("estimate", "lower", "upper")
    kappa_pars[, 2] <- pmax(0, kappa_pars[, 2])
  }
  if (!is.null(shape_pars)) {
    colnames(shape_pars) <- c("estimate", "lower", "upper")
  }
  obj <- list(
    thresh = thresh,
    kappa = kappa_pars,
    shape = shape_pars,
    level = level,
    model = model,
    type = type
  )
  class(obj) <- c("mev_egp_tstab")
  if (isTRUE(plot)) {
    plot(obj, changepar)
  }
  return(invisible(obj))
}
#' @export
plot.mev_egp_tstab <- function(x, ...) {
  args <- list(...)
  thresh <- x$thresh
  kappa <- x$kappa
  shape <- x$shape
  if (!is.null(args$param)) {
    param <- match.arg(
      arg = args$param,
      choices = c("shape", "kappa"),
      several.ok = TRUE
    )
    if ("shape" %in% param & is.null(shape)) {
      stop("Invalid plot choice: \"shape\" is missing from object.")
    }
    if ("kappa" %in% param & is.null(kappa)) {
      stop("Invalid plot choice: \"kappa\" is missing from object")
      ng <- length(param)
    }
  } else {
    param <- c("kappa", "shape")[
      !c(
        is.null(kappa),
        is.null(shape)
      )
    ]
    ng <- length(param)
  }
  changepar <- args$changepar
  if (is.null(changepar)) {
    changepar <- TRUE
  }
  if (!(ng > 0)) {
    stop("Invalid inputs")
  }
  # Graphics
  if (isTRUE(changepar)) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(
      mfrow = c(1, ng),
      mar = c(4.5, 4.5, 3.1, 0.1)
    )
  }
  for (i in seq_along(param)) {
    pars <- get(param[i])
    ylims <- range(pars, na.rm = TRUE)
    plot(
      x = thresh,
      y = pars[, 1],
      pch = 20,
      xlab = "threshold",
      bty = "l",
      ylab = c(expression(kappa), expression(xi))[i],
      ylim = ylims,
      type = "n"
    )
    if (i == ng) {
      text <- switch(
        x$model,
        "pt-beta" = "Papastathopoulos-Tawn's EGP 1",
        "pt-gamma" = "Papastathopoulos-Tawn's EGP 2",
        "pt-power" = "Papastathopoulos-Tawn's EGP 3 (power)",
        "gj-tnorm" = "Gamet-Jalbert's truncated normal EGP",
        "gj-beta" = "Gamet-Jalbert's beta EGP",
        "logist" = "logistic EGP",
        "exptilt" = "exponential tilting EGP"
      )
      mtext(text, side = 3, cex = 0.9, line = 0.5, adj = 1)
    }
    if (param[i] == "kappa") {
      abline(h = 1, lwd = 0.5, col = "gray20", lty = 2)
    }
    arrows(
      x0 = thresh,
      y0 = pars[, 2],
      y1 = pars[, 3],
      length = 0.05,
      angle = 90,
      code = 3
    )
    points(
      x = thresh,
      y = pars[, 1],
      type = "p",
      pch = 20
    )
  }
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


#' Extended generalized Pareto distribution
#'
#' Density function, distribution function, quantile function and
#' random number generation for various extended generalized Pareto
#' distributions
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n scalar number of observations
#' @param kappa shape parameter for the tilting distribution.
#' @param scale scale parameter, strictly positive.
#' @param shape shape parameter.
#' @param model string giving the distribution of the model
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log.p,log logical; if \code{FALSE} (default), values are returned on the probability scale.
#' @references Papastathopoulos, I. and J. Tawn (2013). Extended generalised Pareto models for tail estimation, \emph{Journal of Statistical Planning and Inference} \bold{143}(3), 131--143, <doi:10.1016/j.jspi.2012.07.001>.
#' @references Gamet, P. and Jalbert, J. (2022). A flexible extended generalized Pareto distribution for tail estimation. \emph{Environmetrics}, \bold{33}(6), <doi:10.1002/env.2744>.
#' @name egpdist

#' @rdname egpdist
#' @export
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

#' @rdname egpdist
#' @export
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

#' @rdname egpdist
#' @export
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


#' @rdname egpdist
#' @export
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
  qgp(pg, loc = 0, scale = scale, shape = shape)
}

#' Profile log likelihood for extended generalized Pareto models
#'
#' Computes the profile log likelihood over a grid of values of \eqn{\psi} for various parameters, including return levels.
#'
#' @export
#'
#' @param psi grid of values for the parameter to profile
#' @param model string; choice of extended eneralized Pareto model.
#' @param param string; parameter to profile
#' @param mle a vector or matrix with maximum likelihood estimates of \code{kappa}, \code{scale}, \code{shape}. This can be a matrix if there are multiple threshold
#' @param xdat vector of observations
#' @param thresh vector of positive thresholds. If \code{NULL}, defaults to zero.
#' @param method string giving the optimization method for the outer optimization in the augmented Lagrangian routine; one of \code{nlminb} or \code{BFGS}
#' @param plot logical; if \code{TRUE}, returns a plot of the profile log likelihood
#' @param p tail probability for return level if \code{param="retlev"}.
#' @param ... additional arguments, currently ignored
#' @return an object of class \code{eprof}
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
  param = c("kappa", "scale", "shape", "retlev"),
  mle = NULL,
  xdat,
  thresh = NULL,
  plot = FALSE,
  method = c("Nelder", "nlminb", "BFGS"),
  p,
  ...
) {
  param <- match.arg(param)
  method <- match.arg(method)
  if (param == "retlev") {
    if (missing(p)) {
      stop("Tail probability \"p\" missing.")
    }
  }
  if (is.null(thresh)) {
    thresh <- 0
    rate <- 1
  } else {
    stopifnot(
      is.numeric(thresh),
      length(thresh) == 1L
    )
    stopifnot(thresh >= 0)
    xdata <- xdat[xdat > thresh] - thresh
    rate <- length(xdata) / length(xdat)
  }
  if (is.null(mle)) {
    mle <- try(
      fit.egp(
        xdat = xdata,
        thresh = 0,
        model = model
      ),
      silent = TRUE
    )
    if (inherits(mle, "try-error")) {
      stop("Could not find maximum likelihood estimator.")
    }
  } else {
    stopifnot(inherits(mle, what = "mev_egp"))
  }
  if (param != "retlev") {
    ind <- switch(
      param,
      kappa = 1L,
      scale = 2L,
      shape = 3L
    )
    coef_mle <- coef(mle)[ind]
    se_mle <- mle$std.err[ind]
    if (missing(psi)) {
      if (is.numeric(se_mle) & is.finite(se_mle)) {
        psi <- coef_mle +
          seq(
            from = -3 * se_mle,
            to = 3 * se_mle,
            length.out = 55
          )
      } else {
        stop("Could not determine a suitable sequence of values for profiling.")
      }
    }
    if (param %in% c("kappa", "scale")) {
      psi <- psi[psi > 0]
    } else if (param == "shape") {
      psi <- psi[psi >= -1]
    }
    psi <- sort(unique(c(psi, coef_mle)))
    mid <- which(psi == coef_mle)
    if (mid == 1L) {
      psi <- c(
        seq(0.01, coef_mle, length.out = 10L)[-10],
        psi
      )
      mid <- which(psi == coef_mle)
    }
    # should be 28 for seq of psi if no zero values
    pars <- matrix(nrow = length(psi), ncol = 3L)
    colnames(pars) <- c("kappa", "scale", "shape")
    pll <- numeric(length(psi))
    pll[mid] <- -mle$nllh
    pars[mid, ] <- coef(mle)
    pars[, ind] <- psi
    if (mid > 1) {
      for (i in (mid - 1):1) {
        fpar <- list(as.numeric(pars[i, ind]))
        names(fpar) <- param
        fit_sub <- fit.egp(
          xdat = xdata,
          thresh = 0,
          model = model,
          fpar = fpar,
          start = c(sapply(
            pars[i + 1, -ind],
            function(x) {
              list(x)
            }
          ))
        )
        pars[i, -ind] <- coef(fit_sub)
        pll[i] <- -fit_sub$nllh
      }
    }
    if (mid < length(psi)) {
      for (i in (mid + 1):length(psi)) {
        fpar <- list(as.numeric(pars[i, ind]))
        names(fpar) <- param
        fit_sub <- fit.egp(
          xdat = xdata,
          thresh = 0,
          model = model,
          fpar = fpar,
          method = method,
          start = c(sapply(
            pars[i - 1, -ind],
            function(x) {
              list(x)
            }
          ))
        )
        pars[i, -ind] <- coef(fit_sub)
        pll[i] <- -fit_sub$nllh
      }
    }

    ans <- list(
      mle = coef_mle,
      psi.max = coef_mle,
      param = param,
      std.error = se_mle,
      psi = psi,
      pll = pll,
      maxpll = -mle$nllh,
      family = "egp",
      threshold = thresh
    )
  } else if (param == "retlev") {
    # Get point estimate and
    # standard errors (via Delta-method)
    mle_retlev <- qegp(
      p = 1 - p / rate,
      scale = coef(mle)[2],
      kappa = coef(mle)[1],
      shape = coef(mle)[3],
      model = model
    )
    grad_g <- numDeriv::grad(
      func = function(pars) {
        qegp(
          p = p / rate,
          scale = pars[2],
          shape = pars[3],
          kappa = pars[1],
          model = model,
          lower.tail = FALSE
        )
      },
      x = coef(mle)
    )
    se_retlev <- c(t(grad_g) %*% vcov(mle) %*% grad_g)
    if (missing(psi)) {
      psi <- mle_retlev +
        seq(
          from = -1.5 * se_retlev,
          to = 4 * se_retlev,
          length.out = 55L
        )
      psi <- psi[
        psi >
          quantile(
            xdata,
            probs = min(0.75, 1 - 2 * p)
          )
      ]
    } else {
      psi <- psi - thresh
      psi <- psi[psi > 0]
    }
    psi <- sort(unique(c(psi, mle_retlev)))
    mid <- which(psi == mle_retlev)
    pars <- matrix(nrow = length(psi), ncol = 3L)
    colnames(pars) <- c("kappa", "retlev", "shape")
    pll <- numeric(length(psi))
    egp_pll_retlev <- function(
      par,
      retlev,
      xdat,
      model
    ) {
      scale <- retlev /
        qegp(
          p = 1 - p / rate,
          scale = 1,
          kappa = par[1],
          shape = par[2],
          model = model
        )
      nll <- -egp.ll(
        xdat = xdat,
        thresh = 0,
        par = c(par[1], scale, par[2]),
        model = model
      )
      ifelse(is.finite(nll), nll, 1e10)
    }
    xmax <- max(xdata)
    egp_retlev_hin <- function(
      par,
      retlev,
      xdat,
      model
    ) {
      scale <- retlev /
        qegp(
          p = 1 - p / rate,
          scale = 1,
          kappa = par[1],
          shape = par[2],
          model = model
        )
      params <- c(par[1], scale, par[2])
      c(
        params[1],
        params[2],
        params[3] + 1,
        params[2] + params[3] * xmax
      )
    }
    mid <- which(psi == mle_retlev)
    ind <- 2L
    pars[mid, -ind] <- coef(mle)[-2]
    pars[, ind] <- psi
    pll[mid] <- -mle$nllh
    if (method == "Nelder") {
      opt_fun <- function(pars, retlev) {
        optim(
          par = pars,
          fn = egp_pll_retlev,
          retlev = retlev,
          method = "Nelder",
          xdat = xdata,
          model = model
        )
      }
    } else {
      opt_fun <- function(pars, retlev) {
        alabama::auglag(
          par = pars,
          fn = egp_pll_retlev,
          retlev = retlev,
          hin = egp_retlev_hin,
          control.outer = list(method = method, trace = FALSE),
          control.optim = switch(
            method,
            nlminb = list(
              iter.max = 500L,
              rel.tol = 1e-10,
              step.min = 1e-10
            ),
            list(maxit = 1000L, reltol = 1e-10)
          ),
          xdat = xdata,
          model = model
        )
      }
    }
    if (mid > 1) {
      for (i in (mid - 1):1) {
        opt <- opt_fun(pars = pars[i + 1, -ind], retlev = psi[i])
        pars[i, -ind] <- opt$par
        pll[i] <- -opt$value
      }
    }
    for (i in (mid + 1):length(psi)) {
      opt <- opt_fun(pars = pars[i - 1, -ind], retlev = psi[i])
      pars[i, -ind] <- opt$par
      pll[i] <- -opt$value
    }
    coef_mle <- mle_retlev
    ans <- list(
      mle = coef_mle + thresh,
      psi.max = coef_mle + thresh,
      param = param,
      std.error = se_retlev,
      psi = psi + thresh,
      pll = pll,
      maxpll = -mle$nllh,
      family = "egp",
      threshold = thresh
    )
  }

  ans$r <- sign(coef_mle - psi) *
    sqrt(2 * (ans$maxpll - ans$pll))
  ans$normal <- c(
    coef = as.numeric(ans$psi.max),
    "std.error" = as.numeric(ans$std.error)
  )
  class(ans) <- "eprof"
  if (isTRUE(plot)) {
    plot(ans)
  }
  return(invisible(ans))
}
