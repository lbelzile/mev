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
#' @param type string indicating scale, either \code{qq} for exponential quantile-quantile plot or \code{pp} for probability-probability plot (uniform)
#' @param dist string indicating norm, either \code{l1} for absolute error or \code{l2} for quadratic error
#' @param B number of bootstrap replications
#' @param neval number of points at which to estimate the metric. Default to 1000L
#' @param ci level of symmetric confidence interval. Default to 0.95
#' @return an invisible list with components
#' \itemize{
#' \item \code{thresh}: scalar threshold minimizing criterion
#' \item \code{cthresh}: vector of candidate thresholds
#' \item \code{metric}: value of the metric criterion evaluated at each threshold
#' \item \code{type}: argument \code{type}
#' \item \code{dist}: argument \code{dist}
#' }
#' @export
#' @references Varty, Z. and J.A. Tawn and P.M. Atkinson and S. Bierman (2021+), Inference for extreme earthquake magnitudes accounting for a time-varying measurement process
vmetric.diag <- function(
      xdat,
      thresh,
      B = 199L,
      type = c("qq", "pp"),
      dist = c("l1", "l2"),
      neval = 1000L,
      ci = 0.95){
  dist <- match.arg(dist)
  type <- match.arg(type)
  B <- as.integer(B)
  stopifnot(is.finite(B),
            isTRUE(B > 1),
            is.numeric(xdat))
  xdat <- xdat[is.finite(xdat)]
  thresh <- sort(thresh[is.finite(thresh)])
  nt <- length(thresh)

  # Compute metric from exponential data
  compute_metric <- function(expdat,
                             ppoints,
                             type = c("qq", "pp"),
                             dist = c("l1", "l2")){
    dist <- match.arg(dist)
    type <- match.arg(type)
    if(type == "qq" & dist == "l1"){
      m <- mean(abs(-log(1-ppoints) -
                      quantile(expdat, ppoints)))
    } else if(type == "qq" & dist == "l2"){
      m <- mean((-log(1-ppoints) - quantile(expdat,ppoints))^2)
    } else if(type == "pp" & dist == "l1"){
      m <- mean((ppoints*(1-ppoints)/sqrt(length(expdat)))^(-1/2)*abs(ppoints - ecdf(expdat)(-log(1-ppoints))))
    } else if(type == "pp" & dist == "l2"){
      m <- mean((ppoints*(1-ppoints)/sqrt(length(expdat)))^(-1/2)*(ppoints - ecdf(expdat)(-log(1-ppoints)))^2)
    }
    return(m)
  }
  # Container for results at each threshold
  metric <- se_metric <- numeric(nt)
  scale <- numeric(nt)
  shape <- numeric(nt)
  ps <- ppoints(n = neval)
  # Container for all bootstrap tolerance intervals
  boot_tolerance <- list()
  stopifnot(length(ci) == 1,
            ci > 0,
            ci < 1)
  halfalpha <- (1-ci)/2
  probs <- c(halfalpha, 1-halfalpha)
  for(i in seq_along(thresh)){
   # 1) Fit model
   mle <- mev::fit.gpd(xdat = xdat,
                 threshold = thresh[i])
   scale[i] <- coef(mle)[1]
   shape[i] <- coef(mle)[2]
   boot_par <- mev::gpd.boot(mle,
                        B = B,
                        method = "post")
   # Create GP sample, refit model
   stat <- rep(NA, B)
   boot_samp_mat <- matrix(0, nrow = nobs(mle), ncol = B)
   for(j in seq_len(B)){
     # boot_samp <- mev::rgp(n = nobs(mle),
     #                        scale = boot_par[j,1],
     #                        shape = boot_par[j,2])
     # boot_mle <- try(mev::fit.gpd(boot_samp,
     #                          threshold = 0))
     # if(!inherits(boot_mle, "try-error")){
     boot_samp_mat[, j] <-
       sort(qexp(mev::pgp(q = mle$exceedances,
                      scale = boot_par[j,1],
                      shape = boot_par[j,2])))
       # qexp(mev::pgp(q = mle$exceedances,
       #           scale = coef(boot_mle)[1],
       #           shape = coef(boot_mle)[2]))
    stat[j] <- compute_metric(expdat = boot_samp_mat[, j],
                              ppoints = ps,
                              type = type,
                              dist = dist)
     # } else {
     #   j = j - 1
     # }
   }
   stat <- stat[is.finite(stat)]
   metric[i] <- mean(stat)

   se_metric[i] <- sd(stat)/sqrt(length(stat))
   boot_tolerance[[i]] <- t(apply(boot_samp_mat, 1,
                                  quantile,
                                  probs = probs))
  }
 result <- list(
   thresh = thresh[which.min(metric)],
   cthresh = thresh,
   metric = metric,
   type = type,
   dist = dist,
   scale = scale,
   shape = shape,
   xdat = xdat,
   tolerance = boot_tolerance)
 class(result) <- "mev_thdiag_varty"
 result
}

#' @export
print.mev_thdiag_varty <-
  function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Threshold selection method: metric-based assessment\n")
  cat(switch(x$type,
             "qq" = "exponential (quantile-quantile plot)",
             "pp" = "uniform (probability-probability plot)"), "\n")
  cat(paste(
    switch(x$type,
           "qq" = "",
           "pp" = "weighted"),
    switch(x$dist,
             "l1" = "mean absolute error",
             "l2" = "mean squared error"), "\n"))
  cat("Selected threshold:", round(x$thresh, digits), "\n")
  return(invisible(NULL))
}

#' @export
plot.mev_thdiag_varty <-
  function(x,
           type = c("qqplot", "ppplot", "metric"),
           thid = NULL,
           B = 1000L,
           probs = c(0.025,0.975),
           ...) {
    type <- match.arg(type)
    if(type == "metric"){
   plot(x = x$cthresh,
        y = x$metric,
        bty = "l",
        pch = 20,
        type = "b",
        xlab = "threshold",
        ylab = "metric")
      mtext(side = 3,
            cex = 0.75,
            adj = 1,
            text = paste0(
        switch(x$type,
               "qq" = "",
               "pp" = "weighted"),
        switch(x$dist,
               "l1" = "mean absolute error",
               "l2" = "mean squared error"),
        " (",
        switch(x$type,
               "qq" = "exponential",
               "pp" = "uniform"),
        ")"
        ))
    } else if(type == "qqplot"){
      if(is.null(thid)){
        thresh <- x$thresh
        thid <- which(x$cthresh == thresh)
      } else{
        thid <- as.integer(thid)
        stopifnot(length(thid) == 1L,
                  is.finite(thid),
                  thid <= length(x$cthresh),
                  thid >= 1)
        thresh <- x$cthresh[thid]
      }
      expdat <- qexp(mev::pgp(q = sort(x$xdat[x$xdat > thresh]),
                          loc = thresh,
                          scale = x$scale[thid],
                          shape = x$shape[thid]))
      n <- length(expdat)
      xp <- qexp(ppoints(n))
      theo_quant_sim <-
        apply(
          apply(matrix(rexp(n = B*n), ncol = n),
                1, sort),
                1, quantile, probs = probs)
      obs_quant_sim <- x$tolerance[[thid]]
      above_sim_int <- obs_quant_sim[,1] > theo_quant_sim[2,]
      below_sim_int <- obs_quant_sim[,2] < theo_quant_sim[1,]
      point_colours <- rep(1, n) #black
      point_colours[above_sim_int] <- 2 #red
      point_colours[below_sim_int] <- 4 #blue
      lims <- c(0, pmax(expdat[n], obs_quant_sim[n,2], theo_quant_sim[2,n]))
      plot(
        x = 1,
        y = 1,
        ylab = 'sample quantiles',
        xlab = 'theoretical exponential quantiles',
        type = 'n', # don't plot
        bty = "l",
        ylim = lims,
        xlim = c(0, xp[n] + 0.1),
        xaxs = "i",
        yaxs = "i")
      polygon(
        x = c(xp, rev(xp)),
        y = c(theo_quant_sim[1,], rev(theo_quant_sim[2,])),
        col = "grey90",
        border = NA
      )
      for(j in seq_len(n)){
        lines(
          x = rep(xp[j], 2),
          y = obs_quant_sim[j,],
          col = point_colours[j])
      }
      #points(xp, expdat, pch = 20)
    }  else if(type == "ppplot"){
      if(is.null(thid)){
        thresh <- x$thresh
        thid <- which(x$cthresh == thresh)
      } else{
        thid <- as.integer(thid)
        stopifnot(length(thid) == 1L,
                  is.finite(thid),
                  thid <= length(x$cthresh),
                  thid >= 1)
        thresh <- x$cthresh[thid]
      }
      unifdat <- mev::pgp(q = sort(x$xdat[x$xdat > thresh]),
                               loc = thresh,
                               scale = x$scale[thid],
                               shape = x$shape[thid])
      n <- length(unifdat)
      xp <- ppoints(n)
      theo_quant_sim <-
        apply(
          apply(matrix(runif(n = B*n), ncol = n),
                1, sort),
          1, quantile, probs = probs)
      obs_quant_sim <- pexp(x$tolerance[[thid]])
      above_sim_int <- obs_quant_sim[,1] > theo_quant_sim[2,]
      below_sim_int <- obs_quant_sim[,2] < theo_quant_sim[1,]
      point_colours <- rep(1, n) #black
      point_colours[above_sim_int] <- 2 #red
      point_colours[below_sim_int] <- 4 #blue
      plot(
        x = 1,
        y = 1,
        ylab = 'sample quantiles',
        xlab = 'theoretical uniform quantiles',
        type = 'n', # don't plot
        bty = "l",
        ylim = c(0,1),
        xlim = c(0,1),
        xaxs = "i",
        yaxs = "i")
      polygon(
        x = c(xp, rev(xp)),
        y = c(theo_quant_sim[1,], rev(theo_quant_sim[2,])),
        col = "grey90",
        border = NA
      )
      for(j in seq_len(n)){
        lines(
          x = rep(xp[j], 2),
          y = obs_quant_sim[j,],
          col = point_colours[j])
      }
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
gpd.boot <- function(object,
                        B = 1000L,
                        method = c("post","norm")){
  method <- match.arg(method)
  B <- as.integer(B)
  stopifnot(B > 1L,
            inherits(object, "mev_gpd"))
  if(is.null(object$exceedances)){
    stop("Exported object does not contain exceedances.")
  }
  if(method == "post"){
    if (!requireNamespace("revdbayes", quietly = TRUE)) {
    stop(
      "Package \"revdbayes\" must be installed to use this function.",
      call. = FALSE
    )
    }
    rpostsamp <- suppressWarnings(
      try(revdbayes::rpost(n = B,
                     model = "gp",
                     prior = revdbayes::set_prior(
                       prior = "flat",
                       model = "gp",
                       min_xi = -1),
                     thresh = 0,
                     data = object$exceedances,
                     init_ests = coef(object),
                     trans = "BC")$sim_vals))
    if(inherits(rpostsamp, "try-error")){
      stop("Ratio-of-uniform method failed.")
    }
    boot_par <- rpostsamp
  } else if (method == "norm"){
    if(isTRUE(coef(object)[2] < -0.5)){
      stop("Observed information undefined:\ncannot use normal approximation")
    }
    boot_par <- matrix(NA, ncol = 2, nrow = B)
    stopifnot(isTRUE(all(eigen(object$vcov,
              only.values = TRUE)$values > 0)))
    boot_par[,2] <- rnorm(n = B,
                          mean = coef(object)[2],
                          sd = object$std.err[2])
    vmat <- vcov(object)
    cmean <- coef(object)[1] +
      vmat[1,2]/vmat[2,2]*
      (boot_par[,2] - coef(object)[2])
    csd <- sqrt(vmat[1,1] - vmat[1,2]^2/vmat[2,2])
    # This breaks down if the mean is too small,
    # below -8.3 lower bound, once standardised
    #  but the cases we consider here will have
    #  positive mean
    maxexc <- max(object$exceedances)
    # Sample one-sided truncated normal
    rltnorm <- function(n, mean, sd, lb){
      stopifnot(isTRUE(length(lb) %in% c(1L, n)),
                isTRUE(length(mean) %in% c(1L, n)),
                isTRUE(length(sd) %in% c(1L, n)))
      lbs <- (lb - mean)/sd
      mean + sd*qnorm(
        pnorm(lbs) +
          pnorm(lbs, lower.tail = FALSE)*runif(n))
    }
    if (requireNamespace("TruncatedNormal", quietly = TRUE)) {
      boot_par[,1] <-
        TruncatedNormal::rtnorm(
          n = 1,
          mu = cmean,
          sd = csd,
          lb = ifelse(boot_par[,2]< 0,
                            -boot_par[,2]*maxexc,
                            0),
          ub = rep(Inf, B))
    } else{
    # This works most of the time, but try-catch
    boot_par[,1] <-
      rltnorm(n = B,
              mean = cmean,
              sd = csd,
              lb = ifelse(boot_par[,2]< 0,
                          -boot_par[,2]*maxexc,
                          0))
}
  colnames(boot_par) <- c("scale", "shape")
  return(boot_par)
 }
}
