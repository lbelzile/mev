#' Metric-based threshold selection
#'
#' Adaptation of Varty et al.'s metric-based threshold
#' automated diagnostic for the  independent and identically distributed case with no rounding.
#'
#' The algorithm proceeds by first computing the maximum
#' likelihood algorithm and then simulating replication datasets
#' using either a parametric or nonparametric bootstrap.
#' For each bootstrap sample, we refit the
#'  model and convert the quantiles to
#'  exponential or uniform variates depending on \code{type}, or else if \code{eqd} by calculating the expected plotting positions for the simulated sample.
#'
#' If \code{uq = TRUE} and we specify \code{bootstrap = "parametric"}, the
#' estimation uncertainty is taken into consideration and each sample
#' is drawn from a generalized Pareto distribution, but with different parameters
#' reflecting the sampling distribution.
#'
#'
#' The mean absolute or mean squared distance is calculated on
#' each bootstrap sample at each threshold, and then aggregated into
#' a single average at each \code{thresh} value. The threshold
#' returned is the one with the lowest average value of the metric.
#'
#' Collings et al. (2025) recommend to use quantile-quantile plot, but with
#' \code{pp} starting from some minimal threshold and going no further than the \eqn{1-10/n} probability level. This can be supplied via \code{pp}. When choosin  \code{type = "tails"}, only probability points exceeding the threshold level are kept, so the metric is evaluated at the same levels, but with fewer points, as we increase the threshold level.
#'
#' @param xdat vector of observations
#' @param thresh vector of thresholds
#' @param type string indicating scale, either \code{exp} for exponential quantile-quantile plot  as in Varty et al. (2021) or \code{pp} for probability-probability plot (uniform). The method \code{qq} or expected quantile discrepancy (\code{eqd}) corresponds to Murphy et al. (2024) with the quantiles on the generalized Pareto scale.
#' @param dist string indicating norm, either \code{l1} for absolute error or \code{l2} for quadratic error
#' @param B number of bootstrap replications
#' @param uq logical; if \code{TRUE}, generate bootstrap samples accounting for the sampling distribution of parameters. Only valid when \code{bootstrap = "parametric"}.
#' @param pp plotting positions for the uniform. If \code{type = "tails"}, only the values exceeding the threshold probability level are kept. Default to 250 uniform plotting positions on the unit interval.
#' @param level level of symmetric confidence interval. Default to 0.95
#' @param bootstrap string, one of \code{nonparametric} (sampling with replacement from exceedances) or \code{parametric} (sampling from generalized Pareto).
#' @param plot logical; if \code{TRUE}, returns a plot
#' @param ... additional arguments for backward compatibility
#' @return an invisible list with components
#' \itemize{
#' \item \code{thresh}: scalar threshold minimizing criterion
#' \item \code{thresh0}: vector of candidate thresholds
#' \item \code{metric}: value of the metric criterion evaluated at each threshold
#' \item \code{type}: argument \code{type}
#' \item \code{dist}: argument \code{dist},
#' \item \code{level}: level of confidence interval, from \code{level}
#' \item \code{bootstrap}: type of bootstrap, either \code{parametric} or \code{nonparametric}.
#' }
#' @export
#' @references Varty, Z. and J.A. Tawn and P.M. Atkinson and S. Bierman (2021+), Inference for extreme earthquake magnitudes accounting for a time-varying measurement process.
#' @references Murphy, C., Tawn, J. A., & Varty, Z. (2024). \emph{Automated Threshold Selection and Associated Inference Uncertainty for Univariate Extremes}. Technometrics, 67(\bold{2}), 215–224. <doi:10.1080/00401706.2024.2421744>
#' @references Collings, T.P., C. Murphy-Barltrop, C. Murphy, I.D. Haigh, P.D. Bates, and N.D. Quinn (2025). \emph{Automated tail-informed threshold selection for extreme coastal sea levels}, Natural Hazards and Earth System Sciences, 25(\bold{11}), 4545–4562, <doi:10.5194/nhess-25-4545-2025>.
#'
#' @examples
#' xdat <- rexp(1000, rate = 1/2)
#' thresh <- quantile(xdat, prob = c(0.25,0.5, 0.75))
#' # Method of Murphy, Tawn and Varty (2024) - EQD
#' thv <- thselect.vmetric(xdat, thresh, B = 99)
#' plot(thv)
#' plot(thv, type = "metric")
#' print(thv)
#' # TAILS method
#' tails <- thselect.vmetric(
#'   xdat,
#'   thresh = thresh,
#'   type = "tails",
#'   pp = seq(0.8, 1-10/length(xdat), length.out = 250))
thselect.vmetric <- function(
  xdat,
  thresh,
  B = 199L,
  type = c("eqd", "exp", "qq", "pp", "tails"),
  dist = c("l1", "l2"),
  uq = FALSE,
  bootstrap = c("nonparametric", "parametric"),
  pp = ppoints(250),
  level = 0.95,
  plot = FALSE,
  ...
) {
  args <- list(...)
  if (!is.null(args$ci)) {
    level <- args$ci
  }
  if (!is.null(args$neval)) {
    pp <- ppoints(neval)
    warning(
      "Overwriting \"pp\" supplied values since argument \"neval\" was provided"
    )
  }
  # sort pp and remove missing values
  pp <- sort(pp, na.last = NA)
  neval <- length(pp)
  stopifnot(pp[1] > 0, pp[neval] < 1)
  ci <- level
  stopifnot(ci > 0, ci < 1, length(ci) == 1L, is.finite(ci))
  dist <- match.arg(dist)
  type <- match.arg(type)
  bootstrap <- match.arg(bootstrap)
  type0 <- type
  tails <- type0 == "tails"
  if (type %in% c("eqd", "tails")) {
    type <- "qq"
  }
  B <- as.integer(B)
  stopifnot(is.finite(B), isTRUE(B > 1), is.numeric(xdat))
  xdat <- xdat[is.finite(xdat)]
  thresh <- sort(thresh[is.finite(thresh)])
  nt <- length(thresh)

  # Compute metric - xdat should be exponential/gp/uniform
  compute_metric <- function(
    xdat,
    ppoints,
    type = c("exp", "qq", "pp"),
    dist = c("l1", "l2"),
    pars = c(1, 0)
  ) {
    dist <- match.arg(dist)
    type <- match.arg(type)
    if (type == "exp" & dist == "l1") {
      m <- mean(abs(-log(1 - ppoints) - quantile(xdat, ppoints)))
    } else if (type == "exp" & dist == "l2") {
      m <- mean((-log(1 - ppoints) - quantile(xdat, ppoints))^2)
    } else if (type == "qq" & dist == "l1") {
      m <- mean(abs(
        mev::qgp(ppoints, scale = pars[1], shape = pars[2]) -
          quantile(xdat, ppoints)
      ))
    } else if (type == "qq" & dist == "l2") {
      m <- mean(
        (mev::qgp(ppoints, scale = pars[1], shape = pars[2]) -
          quantile(xdat, ppoints))^2
      )
    } else if (type == "pp" & dist == "l1") {
      # Points are uniform from the bootstrap
      # this avoids having the transformation to
      #  unit exponential / backtransformation,
      #  which is unnecessary
      m <- mean(
        (ppoints * (1 - ppoints) / sqrt(length(xdat)))^(-1 / 2) *
          abs(ppoints - ecdf(xdat)(ppoints))
      )
    } else if (type == "pp" & dist == "l2") {
      m <- mean(
        (ppoints * (1 - ppoints) / sqrt(length(xdat)))^(-1 / 2) *
          (ppoints - ecdf(xdat)(ppoints))^2
      )
    }
    return(m)
  }
  # Container for results at each threshold
  metric <- se_metric <- numeric(nt)
  scale <- numeric(nt)
  shape <- numeric(nt)
  boot_tolerance <- list()
  for (i in seq_along(thresh)) {
    # 1) Fit model
    exc <- xdat[xdat > thresh[i]] - thresh[i]
    mle <- mev::fit.gpd(
      xdat = exc,
      threshold = 0
    )
    # Probability of exceedance and levels (for TAILS)
    if (tails) {
      tau_u <- length(exc) / length(xdat)
      pps <- pp > (1 - tau_u)
      pp_u <- 1 - (1 - pp[pps]) / tau_u
    } else {
      pp_u <- pp
    }
    if (isTRUE(uq)) {
      if (bootstrap == "nonparametric") {
        warning(
          "Invalid option: cannot get uncertainty quantification with nonparametric bootstrap. Setting \"uq = FALSE\"."
        )
        uq <- FALSE
      } else {
        # default is false
        boot_par <- mev::gpd.boot(
          object = mle,
          B = B,
          method = "post"
        )
      }
    }
    scale[i] <- coef(mle)[1]
    shape[i] <- coef(mle)[2]
    spars <- as.numeric(coef(mle))
    # Create GP sample, refit model
    stat <- rep(NA, B)
    boot_samp <- matrix(0, nrow = nobs(mle), ncol = B)
    for (j in seq_len(B)) {
      if (bootstrap == "parametric") {
        if (isTRUE(uq)) {
          spars <- boot_par[j, ]
        }
        boot_samp[, j] <- mev::rgp(
          n = nobs(mle),
          scale = spars[1],
          shape = spars[2]
        )
      } else if (bootstrap == "nonparametric") {
        boot_samp[, j] <- sample(
          exc,
          size = length(exc),
          replace = TRUE
        )
      }
      boot_mle <- try(mev::fit.gpd(boot_samp[, j], thresh = 0))
      if (!inherits(boot_mle, "try-error")) {
        if (type %in% c("exp", "pp")) {
          boot_samp[, j] <- pgp(
            q = boot_samp[, j],
            scale = coef(boot_mle)[1],
            shape = coef(boot_mle)[2]
          )
          if (type == "exp") {
            boot_samp[, j] <- qexp(boot_samp[, j])
          }
        }
        stat[j] <- compute_metric(
          xdat = boot_samp[, j],
          ppoints = pp_u,
          type = type,
          dist = dist,
          pars = coef(boot_mle)
        )
      } else {
        # if the fit failed, try again!
        # generates a new sample
        j = j - 1
      }
    }
    stat <- stat[is.finite(stat)]
    metric[i] <- mean(stat)
    if (type %in% c("exp", "pp")) {
      exc <- pgp(
        q = exc,
        scale = coef(mle)[1],
        shape = coef(mle)[2]
      )
      if (type == "exp") {
        exc <- qexp(exc)
      }
    }
    se_metric[i] <- sd(stat) / sqrt(length(stat))
    boot_tolerance[[i]] <- apply(
      cbind(boot_samp, exc),
      2,
      quantile,
      probs = pp_u
    )
  }
  cindex <- which.min(metric)
  res <- list(
    thresh = thresh,
    thresh0 = thresh[cindex],
    metric = metric,
    type = type,
    dist = dist,
    scale = scale,
    shape = shape,
    xdat = xdat,
    level = level,
    bootstrap = bootstrap,
    tolerance = t(apply(
      boot_tolerance[[cindex]],
      1,
      quantile,
      c((1 - ci) / 2, 1 - (1 - ci) / 2)
    ))
  )
  class(res) <- "mev_thselect_vmetric"
  if (isTRUE(plot)) {
    plot(res)
  }
  return(invisible(res))
}

#' @export
print.mev_thselect_vmetric <-
  function(x, digits = min(3, getOption("digits") - 3), ...) {
    cat("Threshold selection method: metric-based assessment\n")
    cat(
      switch(
        x$type,
        "qq" = "expected quantile discrepancy,",
        "exp" = "exponential,",
        "pp" = "uniform,"
      ),
      paste(
        switch(x$type, "qq" = "", "pp" = "weighted"),
        switch(
          x$dist,
          "l1" = "mean absolute error",
          "l2" = "mean squared error"
        ),
        "\n"
      )
    )
    cat(paste(x$bootstrap, "bootstrap \n"))
    cat("Selected threshold:", round(x$thresh0, digits), "\n")
    return(invisible(NULL))
  }


#' Plots for Varty and al. metric-based threshold selection
#'
#' This S3 method produces quantile-quantile plots with confidence and tolerance bands on various scale (uniform, exponential, generalized Pareto), or a plot of the metric as a function of the threshold.
#'
#' @export
#' @param x an object of class \code{mev_thselect_vmetric} produced by a call to \code{thselect.vmetric}
#' @param type string; a single string indicating the choice of plot
#' @param B number of simulations for variability of estimation
#' @param probs quantile levels for intervals.
#' @param ... additional arguments, currently ignored
#' @return NULL; the function is used to produce a plot
plot.mev_thselect_vmetric <-
  function(
    x,
    type = c("qq", "pp", "exp", "metric"),
    B = 1000L,
    probs = c(0.025, 0.975),
    ...
  ) {
    type <- match.arg(type)
    if (type == "metric") {
      plot(
        x = x$thresh,
        y = x$metric,
        bty = "l",
        pch = 20,
        type = "b",
        xlab = "threshold",
        ylab = "metric"
      )
      mtext(
        side = 3,
        cex = 0.75,
        adj = 1,
        text = paste0(
          switch(x$type, "qq" = "", "pp" = "weighted"),
          switch(
            x$dist,
            "l1" = "mean absolute error",
            "l2" = "mean squared error"
          ),
          " (",
          switch(
            x$type,
            "exp" = "exponential",
            "pp" = "uniform",
            "qq" = "expected quantile discrepancy"
          ),
          ")"
        )
      )
    } else {
      thresh <- x$thresh0
      thid <- which(x$thresh == thresh)
      probs <- sort(probs)[1:2]
      stopifnot(probs[1] > 0, probs[2] < 1)
      if (x$type == "pp") {
        obs_quant_sim <- switch(
          type,
          pp = x$tolerance, # uniform
          exp = qexp(x$tolerance),
          qq = apply(
            x$tolerance,
            2,
            qgp,
            scale = x$scale[thid],
            shape = x$shape[thid]
          )
        )
      } else if (x$type == "exp") {
        obs_quant_sim <- switch(
          type,
          pp = pexp(x$tolerance), # uniform
          exp = x$tolerance,
          qq = apply(
            pexp(x$tolerance),
            2,
            qgp,
            scale = x$scale[thid],
            shape = x$shape[thid]
          )
        )
      } else if (x$type == "qq") {
        obs_quant_sim <- switch(
          type,
          pp = apply(
            x$tolerance,
            2,
            pgp,
            scale = x$scale[thid],
            shape = x$shape[thid]
          ),
          exp = qexp(apply(
            x$tolerance,
            2,
            pgp,
            scale = x$scale[thid],
            shape = x$shape[thid]
          )),
          qq = x$tolerance
        )
      }
      np <- nrow(x$tolerance)
      if (type == "exp") {
        tdat <- qexp(mev::pgp(
          q = sort(x$xdat[x$xdat > thresh]),
          loc = thresh,
          scale = x$scale[thid],
          shape = x$shape[thid]
        ))
        theo_quant_sim <-
          apply(
            apply(matrix(rexp(n = B * np), ncol = np), 1, sort),
            1,
            quantile,
            probs = probs
          )
        xp <- qexp(ppoints(np))
      } else if (type == "qq") {
        tdat <- sort(x$xdat[x$xdat > thresh]) - thresh
        xp <- qgp(
          p = ppoints(np),
          loc = 0,
          scale = x$scale[thid],
          shape = x$shape[thid]
        )
        theo_quant_sim <-
          apply(
            apply(
              matrix(
                rgp(n = B * np, scale = x$scale[thid], shape = x$shape[thid]),
                ncol = np
              ),
              1,
              sort
            ),
            1,
            quantile,
            probs = probs
          )
      } else if (type == "pp") {
        tdat <- mev::pgp(
          q = sort(x$xdat[x$xdat > thresh]),
          loc = thresh,
          scale = x$scale[thid],
          shape = x$shape[thid]
        )
        xp <- ppoints(np)
        theo_quant_sim <-
          apply(
            apply(matrix(runif(n = B * np), ncol = np), 1, sort),
            1,
            quantile,
            probs = probs
          )
      }
      n <- length(tdat)
      above_sim_int <- obs_quant_sim[, 1] > theo_quant_sim[2, ]
      below_sim_int <- obs_quant_sim[, 2] < theo_quant_sim[1, ]
      point_colours <- rep(1, np) #black
      point_colours[above_sim_int] <- 2 #red
      point_colours[below_sim_int] <- 4 #blue
      lims <- c(
        0,
        pmax(tdat[n], obs_quant_sim[np, 2], theo_quant_sim[2, np])
      )
      plot(
        x = 1,
        y = 1,
        ylab = 'sample quantiles',
        xlab = 'theoretical quantiles',
        type = 'n', # don't plot
        bty = "l",
        ylim = lims,
        xlim = c(0, xp[np] + 0.1),
        xaxs = "i",
        yaxs = "i"
      )
      polygon(
        x = c(xp, rev(xp)),
        y = c(theo_quant_sim[1, ], rev(theo_quant_sim[2, ])),
        col = "grey90",
        border = NA
      )

      segments(
        x0 = xp,
        x1 = xp,
        y0 = obs_quant_sim[, 1],
        y1 = obs_quant_sim[, 2],
        col = point_colours,
        cex = 0.5
      )
    }
    return(invisible(NULL))
  }

#' Bootstrap approximation for generalized Pareto parameters
#'
#' Given an object of class \code{mev_gpd},
#' returns a matrix of parameter values to mimic
#' the estimation uncertainty.
#'
#' Two options are available: a normal approximation to
#' the scale and shape based on the maximum likelihood
#' estimates and the observed information matrix.
#' This method uses forward sampling to simulate
#' from a bivariate normal distribution that satisfies
#' the support and positivity constraints
#'
#' The second approximation uses the ratio-of-uniforms
#' method to obtain samples from the posterior
#' distribution with uninformative priors, thus
#' mimicking the joint distribution of maximum likelihood.
#' The benefit of the latter is that it is more reliable
#' in small samples and when the shape is negative.
#'
#' @param object object of class \code{mev_gpd}
#' @param B number of pairs to sample
#' @param method string; one of \code{'norm'} for the
#' normal approximation or \code{'post'} (default) for posterior sampling
#' @return a matrix of size B by 2 whose columns contain scale and shape parameters
#' @export
#' @examples
#' set.seed(2025)
#' xdat <- rgev(100, loc = 0, scale = 2, shape = -0.1)
#' fgp <- fit.gpd(xdat)
#' plot(
#'  gpd.boot(fgp, method = "post")
#' )
#' points(
#'  gpd.boot(fgp, method = "norm"),
#'  col = 2,
#'  pch = 20
#' )
gpd.boot <- function(object, B = 1000L, method = c("post", "norm")) {
  method <- match.arg(method)
  B <- as.integer(B)
  stopifnot(B > 1L, inherits(object, "mev_gpd"))
  if (is.null(object$exceedances)) {
    stop("Exported object does not contain exceedances.")
  }
  if (method == "post") {
    if (!requireNamespace("revdbayes", quietly = TRUE)) {
      stop(
        "Package \"revdbayes\" must be installed to use this function.",
        call. = FALSE
      )
    }
    rpostsamp <- suppressWarnings(
      try(
        revdbayes::rpost(
          n = B,
          model = "gp",
          prior = revdbayes::set_prior(
            prior = "flat",
            model = "gp",
            min_xi = -1
          ),
          thresh = 0,
          data = object$exceedances,
          init_ests = coef(object),
          trans = "BC"
        )$sim_vals
      )
    )
    if (inherits(rpostsamp, "try-error")) {
      stop("Ratio-of-uniform method failed.")
    }
    boot_par <- rpostsamp
  } else if (method == "norm") {
    if (isTRUE(coef(object)[2] < -0.5)) {
      stop("Observed information undefined:\ncannot use normal approximation")
    }
    boot_par <- matrix(NA, ncol = 2, nrow = B)
    stopifnot(isTRUE(all(eigen(object$vcov, only.values = TRUE)$values > 0)))
    boot_par[, 2] <- rnorm(
      n = B,
      mean = coef(object)[2],
      sd = object$std.err[2]
    )
    vmat <- vcov(object)
    cmean <- coef(object)[1] +
      vmat[1, 2] / vmat[2, 2] * (boot_par[, 2] - coef(object)[2])
    csd <- sqrt(vmat[1, 1] - vmat[1, 2]^2 / vmat[2, 2])
    # This breaks down if the mean is too small,
    # below -8.3 lower bound, once standardised
    #  but the cases we consider here will have
    #  positive mean
    maxexc <- max(object$exceedances)
    # Sample one-sided truncated normal
    rltnorm <- function(n, mean, sd, lb) {
      stopifnot(
        isTRUE(length(lb) %in% c(1L, n)),
        isTRUE(length(mean) %in% c(1L, n)),
        isTRUE(length(sd) %in% c(1L, n))
      )
      lbs <- (lb - mean) / sd
      mean +
        sd *
          qnorm(
            pnorm(lbs) +
              pnorm(lbs, lower.tail = FALSE) * runif(n)
          )
    }
    if (requireNamespace("TruncatedNormal", quietly = TRUE)) {
      boot_par[, 1] <-
        TruncatedNormal::rtnorm(
          n = 1,
          mu = cmean,
          sd = csd,
          lb = ifelse(boot_par[, 2] < 0, -boot_par[, 2] * maxexc, 0),
          ub = rep(Inf, B)
        )
    } else {
      # This works most of the time, but try-catch
      boot_par[, 1] <-
        rltnorm(
          n = B,
          mean = cmean,
          sd = csd,
          lb = ifelse(boot_par[, 2] < 0, -boot_par[, 2] * maxexc, 0)
        )
    }
    colnames(boot_par) <- c("scale", "shape")
    return(boot_par)
  }
}
