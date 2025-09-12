#' Metric-based threshold selection
#'
#' Adaptation of Varty et al.'s metric-based threshold
#' automated diagnostic for the  independent and identically distributed case with no rounding.
#'
#' The algorithm proceeds by first computing the maximum
#' likelihood algorithm and then simulating datasets from
#' replication with parameters drawn from a bivariate normal
#' approximation to the maximum likelihood estimator distribution.
#'
#' For each bootstrap sample, we refit the
#'  model and convert the quantiles to
#'  exponential or uniform variates.
#' The mean absolute or mean squared distance
#' is calculated on these. The threshold
#' returned is the one with the lowest value
#' of the metric.
#'
#' @param xdat vector of observations
#' @param thresh vector of thresholds
#' @param type string indicating scale, either \code{exp} for exponential quantile-quantile plot  as in Varty et al. (2021) or \code{pp} for probability-probability plot (uniform). The method \code{qq} or expected quantile discrepancy (\code{eqd}) corresponds to Murphy et al. (2024) with the quantiles on the generalized Pareto scale.
#' @param dist string indicating norm, either \code{l1} for absolute error or \code{l2} for quadratic error
#' @param B number of bootstrap replications
#' @param uq logical; if \code{TRUE}, generate bootstrap samples accounting for the sampling distribution of parameters
#' @param neval number of points at which to estimate the metric. Default to 1000L
#' @param ci level of symmetric confidence interval. Default to 0.95
#' @return an invisible list with components
#' \itemize{
#' \item \code{thresh}: scalar threshold minimizing criterion
#' \item \code{thresh0}: vector of candidate thresholds
#' \item \code{metric}: value of the metric criterion evaluated at each threshold
#' \item \code{type}: argument \code{type}
#' \item \code{dist}: argument \code{dist}
#' }
#' @export
#' @references Varty, Z. and J.A. Tawn and P.M. Atkinson and S. Bierman (2021+), Inference for extreme earthquake magnitudes accounting for a time-varying measurement process.
#' @references Murphy, C., Tawn, J. A., & Varty, Z. (2024). \emph{Automated Threshold Selection and Associated Inference Uncertainty for Univariate Extremes}. Technometrics, 67(\bold{2}), 215–224. <doi:10.1080/00401706.2024.2421744>
#'
#' @examples
#' xdat <- rexp(1000, rate = 1/2)
#' thresh <- quantile(xdat, prob = c(0.25,0.5, 0.75))
#'
thselect.vmetric <- function(
  xdat,
  thresh,
  B = 199L,
  type = c("eqd", "exp", "qq", "pp"),
  dist = c("l1", "l2"),
  uq = FALSE,
  neval = 1000L,
  ci = 0.95
) {
  stopifnot(ci > 0, ci < 1, length(ci) == 1L, is.finite(ci))
  dist <- match.arg(dist)
  type <- match.arg(type)
  if (type == "eqd") {
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
  ps <- ppoints(n = neval)
  boot_tolerance <- list()
  for (i in seq_along(thresh)) {
    # 1) Fit model
    exc <- xdat[xdat > thresh[i]] - thresh[i]
    mle <- mev::fit.gpd(
      xdat = exc,
      threshold = 0
    )
    if (isTRUE(uq)) {
      boot_par <- mev::gpd.boot(
        object = mle,
        B = B,
        method = "post"
      )
    }
    scale[i] <- coef(mle)[1]
    shape[i] <- coef(mle)[2]
    spars <- as.numeric(coef(mle))
    # Create GP sample, refit model
    stat <- rep(NA, B)
    boot_samp <- matrix(0, nrow = nobs(mle), ncol = B)
    for (j in seq_len(B)) {
      if (isTRUE(uq)) {
        spars <- boot_par[j, ]
      }
      boot_samp[, j] <- mev::rgp(
        n = nobs(mle),
        scale = spars[1],
        shape = spars[2]
      )
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
          ppoints = ps,
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
      probs = ps
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
    tolerance = t(apply(
      boot_tolerance[[cindex]],
      1,
      quantile,
      c((1 - ci) / 2, 1 - (1 - ci) / 2)
    ))
  )
  class(res) <- "mev_thselect_vmetric"
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
      probs <- sort(probs)
      stopifnot(probs[1] > 0, probs[2] < 1)
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
