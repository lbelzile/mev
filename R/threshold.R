#' Parameter stability plots for peaks-over-threshold
#'
#' This function computes the maximum likelihood estimate
#' at each provided threshold and plots the estimates (pointwise),
#' along with 95% confidence/credible intervals obtained using Wald or profile confidence intervals,
#' or else from 1000 independent draws from the posterior distribution under
#' vague independent normal prior on the log-scale and shape.
#' The latter two methods better reflect the asymmetry of the estimates than the Wald confidence intervals.
#'
#'
#' @param xdat a vector of observations
#' @param thresh a vector of candidate thresholds at which to compute the estimates.
#' @param method string indicating the method for computing confidence or credible intervals.
#' Must be one of \code{"wald"}, \code{"profile"} or \code{"post"}.
#' @param level confidence level of the intervals. Default to 0.95.
#' @param plot logical; should parameter stability plots be displayed? Default to \code{TRUE}.
#' @param which character vector with elements \code{scale} or \code{shape}
#' @param changepar logical; if \code{TRUE}, changes the graphical parameters.
#' @param ... additional arguments passed to \code{plot}.
#' @return a list with components
#' \itemize{
#' \item{\code{threshold}:} vector of numerical threshold values.
#' \item{\code{mle}:} matrix of modified scale and shape maximum likelihood estimates.
#' \item{\code{lower}:} matrix of lower bounds for the confidence or credible intervals.
#' \item{\code{upper}:} matrix of lower bounds for the confidence or credible intervals.
#' \item{\code{method}:} method for the confidence or coverage intervals.
#' }
#' @note The function is hard coded to prevent fitting a generalized Pareto distribution to samples of size less than 10. If the estimated shape parameters are all on the boundary of the parameter space (meaning \eqn{\hat{\xi}=-1}), then the plots return one-sided confidence intervals for both the modified scale and shape parameters: these typically suggest that the chosen thresholds are too high for estimation to be reliable.
#' @author Leo Belzile
#' @seealso \code{\link[ismev]{gpd.fitrange}}
#' @return plots of the modified scale and shape parameters, with pointwise confidence/credible intervals
#' and an invisible data frame containing the threshold \code{thresh} and the modified scale and shape parameters.
#' @export
#' @examples
#' dat <- abs(rnorm(10000))
#' u <- qnorm(seq(0.9,0.99, by= 0.01))
#' par(mfrow = c(1,2))
#' tstab.gpd(xdat = dat, thresh = u, changepar = FALSE)
#' \dontrun{
#' tstab.gpd(xdat = dat, thresh = u, method = "profile")
#' tstab.gpd(xdat = dat, thresh = u, method = "post")
#' }
tstab.gpd <- function(xdat,
                      thresh,
                      method = c("wald", "profile", "post"),
                      level = 0.95,
                      plot = TRUE,
                      which = c("scale", "shape"),
                      changepar = TRUE,
                      ...) {
  args <- list(...)
  if (isTRUE(all(which %in% 1:2))) {
    which <- c("scale", "shape")[which]
  }
  which <-
    match.arg(which,
              choices = c("scale", "shape"),
              several.ok = TRUE)
  if (missing(xdat) && !is.null(args$dat)) {
    dat <- args$dat
    args$dat <- NULL
  } else{
    dat <- xdat
  }
  dat <- as.vector(dat)
  thresh <- unique(sort(thresh))
  # Strip threshold vectors for which the fit will likely fail
  thresh <- thresh[sapply(thresh, function(u) {
    sum(dat > u)
  }) > 10]
  stopifnot(length(level) == 1, level > 0, level < 1)
  method <- match.arg(method)
  if (method == "post") {
    if (!requireNamespace("revdbayes", quietly = TRUE)) {
      stop(
        "Package \"revdbayes\" needed for this function to work. Please install it.",
        call. = FALSE
      )
    }
  }

  alpha <- (1 - level)
  if (length(thresh) < 2) {
    stop("Invalid threshold: should be of length greater than 2")
  }
  # if(!is.vector(dat) || inherits(dat, "numeric")){
  #   stop("Invalid data")
  # }
  nt <- length(thresh)
  np <- 2
  parmat <- matrix(NA, ncol = np, nrow = nt)
  confintmat <- matrix(NA, ncol = 2 * np, nrow = nt)
  #Some quantities for profile
  if (method == "profile") {
    xmax <- max(dat)
    pllsigmainv <- function(par, sigmat, dat, thresh) {
      sigma <- sigmat + par * thresh
      nll <- -gpd.ll(par = c(sigma, par), dat = dat)
      if (!is.finite(nll)) {
        return(-1e10)
      } else{
        return(nll)
      }
    }
  }
  maxmodifsc <- max(dat) - thresh[1]
  for (i in seq_len(nt)) {
    gpdu <- try(fit.gpd(xdat = dat, threshold = thresh[i]))
    if (!inherits(gpdu, "try-error")) {
      parmat[i, ] <- gpdu$estimate
      parmat[i, 1] <- parmat[i, 1] - parmat[i, 2] * (thresh[i] - thresh[1])
      if (is.matrix(gpdu$vcov)) {
        stderr.transfo <-
          try(sqrt(t(c(1, -thresh[i] + thresh[1])) %*% gpdu$vcov %*% c(1, -thresh[i] +
                                                                         thresh[1]))[1, 1], silent = TRUE)
      }
      if (!is.matrix(gpdu$vcov) |
          inherits(stderr.transfo, "try-error")) {
        stderr.transfo <- NA
      }
      if (method == "wald") {
        confintmat[i, 3] <-
          pmax(-1, gpdu$estimate['shape'] - gpdu$std.err['shape'] * qnorm(1 - alpha / 2))
        confintmat[i, 4] <-
          gpdu$estimate['shape'] + gpdu$std.err['shape'] * qnorm(1 - alpha / 2)
        confintmat[i, 1] <-
          parmat[i, 1] - stderr.transfo * qnorm(1 - alpha / 2)
        confintmat[i, 2] <-
          pmin(maxmodifsc, parmat[i, 1] + stderr.transfo * qnorm(1 - alpha / 2))
      } else if (method == "profile") {
        if ("shape" %in% which) {
          if (gpdu$estimate['shape'] == -1) {
            # Specify grid of psi values (only one-sided)
            profxi <- gpd.pll(
              psi = seq(-1, 0, by = 0.01),
              param = "shape",
              mod = "profile",
              mle = gpdu$estimate,
              dat = gpdu$exceedances,
              plot = FALSE
            )
            confintmat[i, 3:4] <-
              confint(profxi, level = level, print = FALSE)[2:3]
            confintmat[i, 3] <- -1
          } else{
            profxi <- gpd.pll(
              param = "shape",
              mod = "profile",
              mle = gpdu$estimate,
              dat = gpdu$exceedances,
              plot = FALSE
            )

            confintmat[i, 3:4] <-
             pmax(-1, confint(profxi, level = level, print = FALSE)[2:3])
          }
        }
        if ("scale" %in% which) {
          k <- 30L
          prof_vals <- rep(NA_real_, k)
          xi_sigma_vals <- rep(NA_real_, k)
          if (!is.na(stderr.transfo)) {
            grid_psi <-
              parmat[i, 1] + seq(-3 * stderr.transfo,
                                 3.5 * stderr.transfo,
                                 length = k)
          } else{
            grid_psi <-
              seq(-3 * parmat[i, 1] / sqrt(gpdu$nat),
                  3.5 * parmat[i, 1] / sqrt(gpdu$nat),
                  length = k)
          }
          #Profile for scale := sigma_u - xi (u - u_0)
          for (j in seq_len(k)) {
            opt_prof <- try(optimize(
              f = pllsigmainv,
              upper = 1.5,
              lower = pmin(0, pmax(
                -1, -grid_psi[j] / (xmax - thresh[1])
              ) + 1e-10),
              sigmat = grid_psi[j],
              dat = gpdu$exceedances,
              thresh = thresh[i] - thresh[1]
            ),
            silent = TRUE)
            if (!inherits(opt_prof, "try-error")) {
              xi_sigma_vals[j] <- opt_prof$minimum
              prof_vals[j] <- opt_prof$objective
            }
          }
          prof <- structure(
            list(
              psi = grid_psi[-prof_vals != 1e10],
              psi.max = parmat[i, 1],
              pll = -prof_vals[-prof_vals != 1e10],
              maxpll = -gpdu$nllh,
              std.err = stderr.transfo
            ),
            class = "eprof"
          )
          if(isTRUE(parmat[i,2] == -1)){
            spl <- smooth.spline(y = c(prof$psi),
                                 x = c(-2*(prof$pll - prof$maxpll)))
            conf <- sort(c(predict(spl, qchisq(0.95,1))$y, parmat[i,1]))
          } else{
          conf <-  try(confint(
              prof,
              level = level,
              print = FALSE)[2:3],
            silent = TRUE)
          }
          if (!inherits(conf, "try-error")) {
            confintmat[i, 1:2] <- pmin(maxmodifsc, conf)
          }
        }
      } else if (method == "post") {
        postsim <- #suppressWarnings(
          revdbayes::rpost_rcpp(
            n = 1000,
            thresh = 0,
            model = "gp",
            init_ests = gpdu$estimate,
            data = c(-1, gpdu$exceedances),
            trans = "BC",
            prior = revdbayes::set_prior(
              prior = "norm",
              model = "gp",
              mean = rep(0, 2),
              cov = c(100 *
                        gpdu$estimate[1], 1) * diag(np)
            )
          )$sim_vals
        # )
        postsim[, 1] <- postsim[, 1] - postsim[, 2] * (thresh[i] - thresh[1])
        confintmat[i, ] <-
          apply(postsim, 2, quantile, c(alpha / 2, 1 - alpha / 2))
      }

    }
  }
  lower <- confintmat[, c(1, 3)]
  upper <- confintmat[, c(2, 4)]
  colnames(parmat) <-
    colnames(lower) <- colnames(upper) <-  c("modif. scale", "shape")
  ret <-
    structure(
      list(
        threshold = thresh,
        mle = parmat,
        lower = lower,
        upper = upper,
        method = method,
        level = level
      ),
      class = "mev_tstab.gpd"
    )
  if (plot) {
    if (isTRUE(all(ret$mle[, 'shape'] == -1))) {
      warning(
        "The estimated parameters are constant for all thresholds, with shape parameter estimates on the boundary of the parameter space."
      )
    }
    plot(ret,
           which = (1:2)[c("scale", "shape") %in% which],
           changepar = changepar, ...)
  }
  return(invisible(ret))
}

#'@export
plot.mev_tstab.gpd <-
  function(x,
           which = 1:2,
           changepar = TRUE,
           ...) {
    ellipsis <- list(...)
    names_ell <- names(ellipsis)
    if (!is.logical(changepar)) {
      changepar <- FALSE
    }
    ellipsis$changepar <- NULL
    if (!changepar) {
      oldpar <- par(no.readonly = TRUE)
      if (length(which) == 2) {
        par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
      }
      on.exit(par(oldpar))
    }
    if ("main" %in% names_ell) {
      main <- ellipsis$main
      ellipsis$main <- NULL
    } else{
      main <- c("Parameter stability plot", "")

    }
    if (!"sub" %in% names_ell) {
      sub <- c(switch(
        x$method,
        wald = paste("Wald", x$level * 100, "% pointwise confidence intervals"),
        profile = paste(
          "profile likelihood",
          x$level * 100,
          "% pointwise confidence intervals"
        ),
        post = paste(x$level * 100, "% pointwise credible intervals")
      ), "")
    } else{
      sub <- ellipsis$sub
      ellipsis$sub <- NULL
    }
    if (length(which) == 2L) {
      main[2] <- main[1]
      sub[2] <- sub[1]
    }
    if (!"ylab" %in% names_ell) {
      ylab <- c("modif. scale", "shape")
    } else{
      ylab <- ellipsis$ylab
      if (length(ylab) != length(which)) {
        stop(paste(
          "Invalid \"ylab argument: must be a vector of size",
          length(which)
        ))
      }
      ellipsis$ylab <- NULL
    }
    if (!"xlab" %in% names_ell) {
      xlab <- "threshold"
    } else{
      xlab <- ellipsis$xlab
      ellipsis$xlab <- NULL
    }
    if (!"bty" %in% names_ell) {
      bty <- "l"
    } else{
      bty <- ellipsis$bty
      ellipsis$bty <- NULL
    }
    if (!"pch" %in% names_ell) {
      pch <- 20
    } else{
      pch <- ellipsis$pch
      ellipsis$pch <- NULL
    }
    #modified scale
    if (1 %in% which) {
      ylim <- c(min(x$lower[, 1], na.rm = TRUE),
                x$mle[,1] + c(-0.1, 0.1),
                max(x$upper[, 1], na.rm = TRUE))
      # If there are only missing values for upper/lower, this
      # returns a vector of length 0 and +/- Inf
      ylim <- range(ylim[is.finite(ylim)])
      pars <-
        list(
          x = x$threshold,
          y =  x$mle[, 1],
          pch = pch,
          ylab = ylab[1],
          xlab = xlab,
          ylim = ylim,
          bty = bty,
          main = main[1]
        )
      do.call(plot, pars)
      mtext(
        side = 3,
        line = -0.2,
        adj = 0.5,
        text = sub[1]
      )
      for (i in seq_along(x$threshold)) {
        lines(c(x$threshold[i], x$threshold[i]), c(x$lower[i, 1], x$upper[i, 1]))
      }
    }
    if (2 %in% which) {
      ylim <- c(min(x$lower[, 2], na.rm = TRUE),
                x$mle[,2] + c(-0.1, 0.1),
                max(x$upper[, 2], na.rm = TRUE))
      # If there are only missing values for upper/lower, this
      # returns a vector of length 0 and +/- Inf
      ylim <- range(ylim[is.finite(ylim)])
      #shape
      pars <-
        list(
          x = x$threshold,
          y = x$mle[, 2],
          pch = pch,
          ylab = ylab[2],
          xlab = xlab,
          ylim = ylim,
          bty = bty,
          main = main[2]
        )
      do.call(plot,  pars)
      mtext(
        side = 3,
        line = -0.2,
        adj = 0.5,
        text = sub[2]
      )
      for (i in 1:length(x$threshold)) {
        lines(c(x$threshold[i], x$threshold[i]), c(x$lower[i, 2], x$upper[i, 2]))
      }
    }
  }
