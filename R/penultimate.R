### This file contains penultimate models
### (1) the extended GP families of Papastathopoulos and Tawn (2013)
### (2) the penultimate approximations of Smith (1987)


#' Extended generalised Pareto families
#'
#' @description This function provides the log-likelihood and quantiles for the three different families presented
#' in Papastathopoulos and Tawn (2013). The latter include an additional parameter, \eqn{\kappa}.
#' All three families share the same tail index as the generalized Pareto distribution, while allowing for lower thresholds.
#' In the case \eqn{\kappa=1}, the models reduce to the generalised Pareto.
#'
#' @references Papastathopoulos, I. and J. Tawn (2013). Extended generalised Pareto models for tail estimation, \emph{Journal of Statistical Planning and Inference} \bold{143}(3), 131--143.
#' @section Usage:
#' \code{egp.ll(xdat, thresh, par, model=c('egp1','egp2','egp3'))}
#' @section Usage:
#' \code{egp.retlev(xdat, thresh, par, model=c('egp1','egp2','egp3'), p, plot=TRUE)}
#' @description \code{egp.retlev} gives the return levels for the extended generalised Pareto distributions
#' @name egp
#' @param xdat vector of observations, greater than the threshold
#' @param thresh threshold value
#' @param par parameter vector (\eqn{\kappa}, \eqn{\sigma},\eqn{\xi}).
#' @param model a string indicating which extended family to fit
#' @param show logical; if \code{TRUE}, print the results of the optimization
#' @param p extreme event probability; \code{p} must be greater than the rate of exceedance for the calculation to make sense. See \bold{Details}.
#' @param plot boolean indicating whether or not to plot the return levels
#' @importFrom grDevices rainbow
#'
#' @details
#'
#'For return levels, the \code{p} argument can be related to \eqn{T} year exceedances as follows:
#'if there are \eqn{n_y} observations per year, than take \code{p}
#'to equal \eqn{1/(Tn_y)} to obtain the \eqn{T}-years return level.
#' @author Leo Belzile
#' @return \code{egp.ll} returns the log-likelihood value.
#' @return \code{egp.retlev} returns a plot of the return levels if \code{plot=TRUE} and a matrix of return levels.
#' @examples
#' set.seed(123)
#' xdat <- evd::rgpd(1000, loc = 0, scale = 2, shape = 0.5)
#' par <- fit.egp(xdat, thresh = 0, model = 'egp3')$par
#' p <- c(1/1000, 1/1500, 1/2000)
#' #With multiple thresholds
#' th <- c(0, 0.1, 0.2, 1)
#' opt <- tstab.egp(xdat, th, model = 'egp1')
#' egp.retlev(xdat, opt$thresh, opt$par, 'egp1', p = p)
#' opt <- tstab.egp(xdat, th, model = 'egp2', plots = NA)
#' egp.retlev(xdat, opt$thresh, opt$par, 'egp2', p = p)
#' opt <- tstab.egp(xdat, th, model = 'egp3', plots = NA)
#' egp.retlev(xdat, opt$thresh, opt$par, 'egp3', p = p)
NULL

#' Fit of extended GP models and parameter stability plots
#'
#' This function is an alias of \code{\link{fit.egp}}.
#'
#' Supported for backward compatibility
#'
#' @export
#' @keywords internal
egp.fit <- function(xdat, thresh, model = c("egp1", "egp2", "egp3"), init, show = FALSE){
  fit.egp(xdat = xdat, thresh = thresh, model = model, init = init, show = show)
}

#' Extended generalised Pareto families of Papastathopoulos and Tawn (functions)
#'
#' This function provides the log-likelihood and quantiles for the three different families presented
#' in Papastathopoulos and Tawn (2013). The latter include an additional parameter, \eqn{\kappa}.
#' All three families share the same tail index than the GP model, while allowing for lower thresholds.
#' @export
#' @inheritParams egp
#' @name egp-function
#' @keywords internal
egp.ll <- function(xdat, thresh, par, model = c("egp1", "egp2", "egp3")) {
    if (!(model %in% c("egp1", "egp2", "egp3")) || length(model) != 1) {
        stop("Invalid model selection")
    }
    if (isTRUE(any(xdat < thresh))) {
        xdat = xdat[xdat > thresh]
    }
    kappa = par[1]
    sigma = par[2]
    xi = par[3]
    if (sigma < 0 || kappa < 0) {
        return(-Inf)
    }
    args = pmax(0, (1 + xi * (xdat - thresh)/sigma))
    if (abs(xi) > 1e-08) {
        switch(model, egp1 = length(xdat) * (log(abs(xi)) - log(sigma) - lbeta(kappa, 1/abs(xi))) + (kappa - 1) * sum(log(pmax(0,
            1 - args^(-sign(xi))))), egp2 = length(xdat) * (-log(sigma) - lgamma(kappa)) + (kappa - 1) * sum(log(log(args)/xi)), egp3 = length(xdat) *
            (log(kappa) - log(sigma)) + (kappa - 1) * sum(log(1 - args^(-1/xi)))) - (1/xi + 1) * sum(log(args))
    } else {
        # if xi=0
        switch(model, egp1 = length(xdat) * (-log(sigma) - lgamma(kappa)) + (kappa - 1) * sum(log(xdat - thresh)), egp2 = length(xdat) *
            (-log(sigma) - lgamma(kappa)) + (kappa - 1) * sum(log(xdat - thresh)), egp3 = length(xdat) * (-log(sigma) + log(kappa)) +
            (kappa - 1) * sum(log(1 - exp(-(xdat - thresh)/sigma)))) - sum(xdat - thresh)/sigma
    }
}

#' @name egp-function
#' @export
#' @inheritParams egp
#' @keywords internal
egp.retlev <- function(xdat, thresh, par, model = c("egp1", "egp2", "egp3"), p, plot = TRUE) {
    if (!(model %in% c("egp1", "egp2", "egp3")) || length(model) != 1) {
        stop("Invalid model selection")
    }
    if (length(par)%%3 != 0) {
        stop("Invalid parameter input")
    }
    if (class(par) != "matrix") {
        par = matrix(c(par), ncol = 3)
    }
    rate <- sapply(thresh, function(u) {
        length(xdat[xdat > u])/length(xdat)
    })
    if (!isTRUE(all.equal(length(rate), length(thresh), nrow(par)))) {
        stop("Input dimension does not match")
    }
    retlev <- matrix(0, nrow = length(thresh), ncol = length(p))
    if (any(sapply(rate, function(zeta) {
        zeta < p
    }))) {
        warning("Some probabilities `p` are higher than the exceedance rate. Evaluate those empirically")
    }
    p <- sort(p)
    for (i in 1:length(thresh)) {
        for (j in 1:length(p)) {
            pl = 1 - p[j]/rate[i]
            if (par[i, 3] == 0) {
                retlev[i, j] <- thresh[i] - par[i, 2] * switch(model, egp1 = -qgamma(pl, scale = 1, shape = par[i, 1]), egp2 = -qgamma(pl,
                  scale = 1, shape = par[i, 1]), egp3 = log(1 - pl^(1/par[i, 1])))
            } else {
                retlev[i, j] <- thresh[i] + par[i, 2]/par[i, 3] * (switch(model, egp1 = (1 - qbeta(pl, par[i, 1], 1/abs(par[i, 3])))^(-sign(par[i,
                  3])), egp2 = exp(par[i, 3] * qgamma(pl, scale = 1, shape = par[i, 1])), egp3 = (1 - pl^(1/par[i, 1]))^(-par[i, 3])) -
                  1)
            }
        }
    }
    if (plot) {
        matplot(1/p, t(retlev), , type = "b", lty = rep(1, length(thresh)), col = rainbow(n = length(thresh), start = 2/6, end = 0),
            xlab = "Return period", ylab = "Estimated return level", main = paste0("Return level plot for EGP", substr(model, 4, 4),
                ""))
    }
    return(retlev)
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
#' @description The function \code{tstab.egp} provides classical threshold stability plot for (\eqn{\kappa}, \eqn{\sigma}, \eqn{\xi}). The fitted parameter values are displayed with pointwise normal 95\% confidence intervals.
#'  The plot is for the modified scale (as in the generalised Pareto model) and as such it is possible that the modified scale be negative.
#' \code{tstab.egp} can also be used to fit the model to multiple thresholds.
#' @param plots vector of integers specifying which parameter stability to plot (if any); passing \code{NA} results in no plots
#' @inheritParams egp
#' @param umin optional minimum value considered for threshold (if \code{thresh} is not provided)
#' @param umax optional maximum value considered for threshold (if \code{thresh} is not provided)
#' @param nint optional integer number specifying the number of thresholds to test.
#' @return \code{tstab.egp} returns a plot(s) of the parameters fit over the range of provided thresholds, with pointwise normal confidence intervals; the function also returns an invisible list containing notably the matrix of point estimates (\code{par}) and standard errors (\code{se}).
#' @importFrom graphics arrows points polygon title
#' @export
#' @examples
#' xdat <- evd::rgpd(n = 100, loc = 0, scale = 1, shape = 0.5)
#' fitted <- fit.egp(xdat = xdat, thresh = 1, model = "egp2", show = TRUE)
#' thresh <- evd::qgpd(seq(0.1, 0.5, by = 0.05), 0, 1, 0.5)
#' tstab.egp(xdat = xdat, thresh = thresh, model = "egp2", plots = 1:3)
fit.egp <- function(xdat, thresh, model = c("egp1", "egp2", "egp3"), init, show = FALSE) {
    if (!(model %in% c("egp1", "egp2", "egp3")) || length(model) != 1) {
        stop("Invalid model  argument: must be one of `egp1', `egp2' or `egp3'.")
    }
    if (length(thresh) > 1) {
        warning("Length of threshold vector greater than one. Selecting first component.")
      thresh <- thresh[1]
    }
    # If no initial values are provided, fit a GP distribution to obtain them
    changinit <- missing(init)
    if(!changinit){
     if(any(init[1:2] < 0)){
       changinit <- TRUE
     }
    }
    if (changinit) {
        init <- c(kappa = 1.01, suppressWarnings(fit.gpd(xdat, threshold = thresh[1], show = FALSE)$est))
    }

    # Keep exceedances only
    xdata = xdat[xdat > thresh]
    xmax <- max(xdata);
    mle <- alabama::auglag(par = init,
                           fn = function(par, xdat, thresh, model){-egp.ll(par = par, xdat = xdat, thresh = thresh, model = model)},
                           hin = function(par, ...){c(par[1]-1e-10, par[2]-1e-10, par[3]+1,
                                                      ifelse(par[3] < 0, thresh-par[2]/par[3] - xmax, 1))},
                           xdat = xdata, thresh = thresh, model = model,
                           control.outer = list(trace = FALSE, method = "BFGS"),
                           control.optim = list(maxit = 500, reltol = 1e-10))
    fitted <- list()
    fitted$estimate <- fitted$param <- mle$par
    fitted$deviance <- 2*mle$value
    if (mle$convergence == 0) {
      fitted$convergence <- "successful"
      fitted$vcov <- try(solve(mle$hessian))
      fitted$std.err <- try(sqrt(diag(fitted$vcov)))
      if(is.character(mle$se) || mle$par[3] < -0.5){
        fitted$vcov <- NULL
        fitted$se <- rep(NA, 3)
      }
    } else{
      fitted$convergence <- mle$convergence
      warning("Maximization routine may have failed; check output and try providing better starting values")
    }
    names(fitted$estimate) <- names(fitted$std.err) <-  c("kappa", "scale", "shape")
    fitted$counts <- mle$counts
    fitted$threshold <- thresh
    fitted$nat <- length(xdata)
    fitted$pat <- length(xdata)/length(xdat)
    fitted$exceedances <- xdata
    fitted$hessian <- mle$hessian
    fitted$method <- "copt"
    fitted$model <- model
    class(fitted) <- c("mev_egp")
    if(show){
      print(fitted)
    }
    return(invisible(fitted))
}

# @param x A fitted object of class \code{mev_gpd}.
# @param digits Number of digits to display in \code{print} call.
# @param ... Additional argument passed to \code{print}.
#' @export
print.mev_egp <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Model:", x$model, "\n")
  cat("Deviance:", round(x$deviance, digits), "\n")

  cat("\nThreshold:", round(x$threshold, digits), "\n")
  cat("Number Above:", x$nat, "\n")
  cat("Proportion Above:", round(x$pat, digits), "\n")

  cat("\nEstimates\n")
  print.default(format(x$estimate, digits = digits), print.gap = 2, quote = FALSE, ...)
  if (!is.na(x$std.err[1]) && x$estimate[3] > -0.5) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, quote = FALSE, ...)
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
egp.fitrange <- function(xdat, thresh, model = c("egp1", "egp2", "egp3"), plots = 1:3, umin, umax, nint){
  tstab.egp(xdat = xdat, thresh = thresh, model = model, plots = plots, umin = umin, umax = umax, nint = nint)
}

#' @inheritParams egp
#' @rdname fit.egp
#' @export
tstab.egp <- function(xdat, thresh, model = c("egp1", "egp2", "egp3"), plots = 1:3, umin, umax, nint) {
    if (!(model %in% c("egp1", "egp2", "egp3")) || length(model) != 1) {
        stop("Invalid model selection")
    }
    if (missing(thresh) && isTRUE(any(c(missing(umin), missing(umax))))) {
        stop("Must provide either minimum and maximum threshold values, or a vector of threshold `thresh'")
    } else if (missing(thresh)) {
        stopifnot(class(umin) %in% c("numeric", "integer"),
                  class(umax) %in% c("numeric", "integer"),
                  length(umin) == 1, length(umax) == 1, umin < umax)
        thresh <- seq(umin, umax, length = nint)
    } else if (length(thresh) <= 1) {
        stop("Invalid `thresh' provided; please use a vector of threshold candidates of length at least 2")
    }
    pe <- se <- matrix(0, ncol = 4, nrow = length(thresh))
    conv <- rep(0, length(thresh))
    fit <- suppressWarnings(fit.egp(xdat = xdat, thresh = thresh[1], model = model))
    pe[1, -4] <- fit$param
    se[1, -4] <- fit$std.err
    conv[1] <- ifelse(is.character(fit$convergence), 0, fit$convergence)
    se[1, 4] <- sqrt(cbind(1, -thresh[1]) %*% solve(fit$hessian[-1, -1]) %*% rbind(1, -thresh[1]))[1]
    for (i in 2:length(thresh)) {
        fit <- suppressWarnings(fit.egp(xdat = xdat, thresh = thresh[i], model = model, init = pe[i - 1, -4]))
        pe[i, -4] <- fit$param
        se[i, -4] <- fit$std.err
        conv[i] <- ifelse(is.character(fit$convergence), 0, fit$convergence)
        # Standard error for the modified scale via the delta-method
        se[i, 4] <- sqrt(cbind(1, -thresh[i]) %*% solve(fit$hessian[-1, -1]) %*% rbind(1, -thresh[i]))[1]
    }
    # Modify point estimates for the scale
    pe[, 4] <- pe[, 2] - pe[, 3] * thresh

    # Graphics
    if (!(length(plots) == 1 && is.na(plots))) {
        plots <- sort(unique(plots))
        if (!isTRUE(all(plots %in% 1:3))) {
            stop("Invalid plot selection. Must be a vector of integers containing indices 1, 2 or 3.")
        }

        old.par <- par(no.readonly = TRUE)
        on.exit(par(old.par))
        par(mfrow = c(length(plots), 1), mar = c(4.5, 4.5, 3.1, 0.1))
        for (i in plots) {
            if (i == 2) {  i <- 4  }  #Get modified scale
            # Plotting devices limits
            ylims = c(min(pe[, i]) - qnorm(0.975) * max(se[, i]), max(pe[, i]) + qnorm(0.975) * max(se[, i]))
            plot(x = thresh, y = pe[, i], pch = 20, xlab = "Threshold", bty = "l", ylab = switch(i, expression(kappa), expression(sigma),
                expression(xi), expression(tilde(sigma))), ylim = ylims, type = "n")  #,cex.lab=1.25)
            polygon(c(thresh, rev(thresh)), c(pe[, i] - qnorm(0.975) * se[, i], rev(pe[, i] + qnorm(0.975) * se[, i])), col = "gray95",
                border = FALSE)
            if (i == min(plots)) {
                title(paste0("Parameter stability plots for EGP", substr(model, 4, 4), ""), outer = FALSE)
            }
            if (i == 1) {
                abline(h = 1, lwd = 0.5, col = "gray20", lty = 2)
            }
            arrows(x0 = thresh, y0 = pe[, i] - qnorm(0.975) * se[, i], y1 = pe[, i] + qnorm(0.975) * se[, i], length = 0.05, angle = 90,
                code = 3)
            points(x = thresh, y = pe[, i], type = "b", pch = 20)
        }
      }
    return(invisible(list(par = pe[, -4], se = se[, -4], model = model, conv = conv, thresh = thresh)))
}


## (2) Smith penultimate approximations

#' Smith's penultimate approximations
#'
#' The function takes as arguments the distribution and density functions. There are two options:
#' \code{method='bm'} yields block maxima and the user should provide in such case the block sizes via the
#' argument \code{m}. If instead \code{method='pot'} is provided, a vector of threshold values must be
#' provided. The other argument (\code{u} or \code{m} depending on the method) is ignored.
#'
#' Alternatively, the user can provide functions \code{densF}, \code{quantF} and \code{distF} for the density,
#' quantile function and distribution functions, respectively. The user can also supply the derivative
#' of the density function, \code{ddensF}. If the latter is missing, it will be approximated using finite-differences.
#'
#' @param family the name of the parametric family. Will be used to obtain \code{dfamily}, \code{pfamily}, \code{qfamily}
#' @param method either block maxima (\code{'bm'}) or peaks-over-threshold (\code{'pot'}) are supported
#' @param u vector of thresholds for method \code{'pot'}
#' @param qu vector of quantiles for method \code{'pot'}. Ignored if argument \code{u} is provided.
#' @param m vector of block sizes for method \code{'bm'}
#' @param returnList logical; should the arguments be returned as a list or as a matrix of parameter
#' @param ... additional arguments passed to \code{densF} and \code{distF}
#' @author Leo Belzile
#' @importFrom methods formalArgs
#' @import stats
#' @return either a vector, a matrix if either \code{length(m)>1} or \code{length(u)>1} or a list (if \code{returnList}) containing
#' \itemize{
#' \item{\code{loc}:} {location parameters (\code{method='bm'})}
#' \item{\code{scale}:} {scale parameters}
#' \item{\code{shape}:} {shape parameters}
#' \item{\code{u:}} {thresholds (if \code{method='pot'})}
#' \item{\code{u:}} {percentile corresponding to threshold (if \code{method='pot'})}
#' \item{\code{m:}} {block sizes (if \code{method='bm'})}
#' }
#' @references Smith, R.L. (1987). Approximations in extreme value theory. \emph{Technical report 205}, Center for Stochastic Process, University of North Carolina, 1--34.
#' @examples
#' #Threshold exceedance for Normal variables
#' qu <- seq(1,5,by=0.02)
#' penult <- smith.penult(family = "norm", ddensF=function(x){-x*dnorm(x)},
#'    method = 'pot', u = qu)
#' plot(qu, penult$shape, type='l', xlab='Quantile',
#'    ylab='Penultimate shape', ylim=c(-0.5,0))
#' #Block maxima for Gamma variables -
#' #User must provide arguments for shape (or rate)
#' m <- seq(30, 3650, by=30)
#' penult <- smith.penult(family = 'gamma', method = 'bm', m=m, shape=0.1)
#' plot(m, penult$shape, type='l', xlab='Quantile', ylab='Penultimate shape')
#'
#' #Comparing density of GEV approximation with true density of maxima
#' m <- 100 #block of size 100
#' p <- smith.penult(family='norm',
#'    ddensF=function(x){-x*dnorm(x)}, method='bm', m=m, returnList=FALSE)
#' x <- seq(1, 5, by = 0.01)
#' plot(x, m*dnorm(x)*exp((m-1)*pnorm(x,log.p=TRUE)),type='l', ylab='Density',
#' main='Distribution of the maxima of\n 100 standard normal variates')
#' lines(x, evd::dgev(x,loc=p[1], scale=p[2], shape=0),col=2)
#' lines(x, evd::dgev(x,loc=p[1], scale=p[2], shape=p[3]),col=3)
#' legend(x = 'topright',lty = c(1,1,1,1), col = c(1,2,3,4),
#'    legend = c('exact', 'ultimate', 'penultimate'), bty = 'n')
#' @export
smith.penult <- function(family, method = c("bm", "pot"), u, qu, m, returnList = TRUE, ...) {
  ellips <- list(...)
  #Compatibility condition with version 1.12 and before, via ellipsis argument
  if(!is.null(ellips$model) && length(method) == 2){
    method <- ellips$model
  }
    # Redefine density, quantile and distribution functions from family
    if (!missing(family)) {
        densF <- paste0("d", family)
        distF <- paste0("p", family)
        quantF <- paste0("q", family)
        computeQuant <- TRUE
    } else{ #compatibility - copy from previous
      if(any(c(is.null(ellips$densF),
               is.null(ellips$distF)))){
        stop("Argument `family` missing.")
      } else{
       densF <- ellips$densF
       distF <- ellips$distF
       if(!is.null(ellips$quantF)){
        quantF <-  ellips$quantF
        computeQuant <- TRUE
       } else{
        computeQuant <- FALSE
       }
      }
    }
    # Matching extra arguments with additional ones passed via ellipsis
    # Which are formals of the function
    indf <- names(ellips) %in% formalArgs(densF)
    indF <- names(ellips) %in% formalArgs(distF)
    if (!is.null(quantF)) {
        indQ <- names(ellips) %in% formalArgs(quantF)
    }
    fn.arg <- ellips[which(indf * (indf == indF) == 1)]
    method <- match.arg(method)
    # Distribution function, density and density derivative
    densFn <- function(x) {
        do.call(densF, c(x = x, fn.arg))
    }
    distFn <- function(x) {
        do.call(distF, c(q = x, fn.arg))
    }
    if (is.null(ellips$ddensF)) {
        ddensFn <- function(x) {
          tol <- 6e-6
            #numDeriv::grad(densFn, x = x, method = "Richardson")
            (densFn(x+tol) - densFn(x-tol))/(2*tol)
        }
    } else {
      ddensF <- ellips$ddensF
        if (class(ddensF) != "function") {
            stop("Invalid arguments. Please provide valid functions.")
        }
        ddensFn <- function(x) {
            do.call(ddensF, c(x, fn.arg))
        }
    }
    # not checking for concordance via numerical derivatives, but could be added Block maxima
    if (method == "bm") {
        if (missing(m)) {
            stop("Sequence of block size must be provided.")
        }
        # Normalizing sequence
        bm <- sapply(m, function(n) {
            if (!is.null(quantF)) {
                bmroot <- do.call(quantF, c(p = exp(-1/n), fn.arg))
            } else {
                bmroot <- uniroot(f = function(bn) {
                  do.call(distF, c(bn, fn.arg)) - exp(-1/n)
                }, lower = -1e+15, upper = 1e+15, f.lower = -exp(-1/n), f.upper = 1 - exp(-1/n), tol = 1e-08, check.conv = TRUE)
                if (abs(bmroot$f.root) < 1e-05) {
                  return(bmroot$root)
                } else {
                  warning("Could not find `bm` using numerical root finder.")
                  return(NA)
                }
            }
        })
        # Penultimate scale and shape functions
        phi <- function(x) {
            -sapply(x, function(xval) {
                distFn(xval) * log(distFn(xval))/densFn(xval)
            })
        }
        dphi <- function(x) {
            sapply(x, function(xval) {
                -(1 + log(distFn(xval))) + distFn(xval) * log(distFn(xval)) * ddensFn(xval)/(densFn(xval)^2)
            })
        }
        if (returnList) {
            params <- list(loc = bm, scale = phi(bm), shape = dphi(bm), m = m)
        } else {
            params <- cbind(loc = bm, scale = phi(bm), shape = dphi(bm), m = m)
            if(nrow(params) == 1L){
              params <- params[1,]
            }
        }
        return(params)

    } else if (method == "pot") {
        if (missing(u) && missing(qu)) {
            stop("Sequence of thresholds must be provided.")
        } else if(missing(u) && !missing(qu)){
          if(computeQuant){
            u <-  sapply(qu, function(p){do.call(quantF, c(p = p, fn.arg))})
          } else{
            u <- rep(NA, length(qu))
          }
        }  else if(!missing(u) && missing(qu)){
         qu <-  sapply(u, function(q){distFn(x = q)})
        }
        phi <- function(x) {
            sapply(x, function(xval) {
                (1 - distFn(xval))/densFn(xval)
            })
        }
        dphi <- function(x) {
            sapply(x, function(xval) {
                -1 - (1 - distFn(xval)) * ddensFn(xval)/(densFn(xval)^2)
            })
        }
        if (returnList) {
            params <- list(u = u, scale = phi(u), shape = dphi(u),  qu = qu)
        } else {
            params <- cbind(u = u, scale = phi(u), shape = dphi(u),  qu = qu)
            if(nrow(params) == 1L){
             params <- params[1,]
            }
        }
        return(params)
    }
    ## The approximations are for \eqn{F^m(x)} for the GEV distribution and for \eqn{1-F(u+x)}{1-F(u)} for the GP distribution.

}


#' Smith's third penultimate approximation
#'
#' This function returns the density and distribution functions
#' of the 3rd penultimate approximation for extremes of Smith (1987). It requires
#' knowledge of the exact constants \eqn{\epsilon} and \eqn{\rho} described in the paper.
#'
#' Let \eqn{F}, \eqn{f} denote respectively the distribution and density functions and define the function \eqn{\phi(x)}  as
#' \deqn{\phi(x)=-\frac{F(x)\log F(x)}{f(x)}}{\phi(x)=-F(x)log F(x)/f(x)}
#' for block maxima.
#' The sequence \code{loc} corresponds to \eqn{b_n} otherwise, defined as the solution of \eqn{F(b_n)=\exp(-1/n)}{F(b_n)=exp(-1/n)}.
#'
#' The \code{scale} is given by \eqn{a_n=\phi(b_n)}, the \code{shape} as \eqn{\gamma_n=\phi'(b_n)}. These are returned by a call to \link{smith.penult}.
#'
#' For threshold exceedances, \eqn{b_n} is replaced by the sequence of thresholds \eqn{u} and we
#' take instead \eqn{\phi(x)} to be the reciprocal hazard function \eqn{\phi(x)=(1-F(x))/f(x)}{\phi(x)=(1-F(x))/f(x)}.
#'
#' In cases where the distribution function is in the maximum domain of
#' attraction of the Gumbel distribution, \eqn{\rho} is possibly undetermined and
#' \eqn{\epsilon} can be equal to \eqn{\phi(b_n)\phi''(b_n)}.
#'
#' For distributions in the maximum domain of
#' attraction of the Gumbel distribution and that are class N, it is also possible to abstract from the \eqn{\rho} parameter by substituting the function \eqn{H_{\rho}}{H(x;\rho)} by \eqn{x^3/6} without affecting the rate of convergence. This can be done by setting \code{mdaGumbel=TRUE} in the function call.
#'
#' @section Warning:
#' The third penultimate approximation does not yield a valid distribution function over the whole range of the original distribution, but is rather valid in a neighborhood of the true support of the distribution of maxima/threshold exceedance.
#' The function handles the most standard failure (decreasing distribution function and negative densities), but any oscillatory behaviour will not necessarily be captured.
#' This is inherent to the method and can be resolved by `not' evaluating the functions \eqn{F} and \eqn{f} at the faulty points.
#' @param loc location parameter returned by \code{\link{smith.penult}} or threshold vector
#' @param scale scale parameter returned by \code{\link{smith.penult}}
#' @param shape shape parameter returned by \code{\link{smith.penult}}
#' @param eps parameter vector, see \strong{Details}.
#' @param rho second-order parameter, model dependent
#' @param method one of \code{pot} for the generalized Pareto or \code{bm} for the generalized extreme value distribution
#' @param mdaGumbel logical indicating whether the function \eqn{H_{\rho}}{H(x;\rho)} should be replaced by \eqn{x^3/6}; see \strong{Details}.
#' @param ... additional parameters, currently ignored. These are used for backward compatibility due to a change in the names of the arguments.
#' @references Smith, R.L. (1987). Approximations in extreme value theory. \emph{Technical report 205}, Center for Stochastic Process, University of North Carolina, 1--34.
#' @examples
#' #Normal maxima example from Smith (1987)
#' m <- 100 #block of size 100
#' p <- smith.penult(family='norm',
#'    ddensF=function(x){-x*dnorm(x)}, method='bm', m=m, returnList=FALSE)
#' approx <- smith.penult.fn(loc=p[1], scale=p[2], shape=p[3],
#'    eps=p[3]^2+p[3]+p[2]^2, mdaGumbel=TRUE, method='bm')
#' x <- seq(0.5,6,by=0.001)
#' #First penultimate approximation
#' plot(x, exp(m*pnorm(x, log.p=TRUE)),type='l', ylab='CDF',
#' main='Distribution of the maxima of\n 100 standard normal variates')
#' lines(x, evd::pgev(x,loc=p[1], scale=p[2], shape=0),col=2)
#' lines(x, evd::pgev(x,loc=p[1], scale=p[2], shape=p[3]),col=3)
#' lines(x, approx$F(x),col=4)
#' legend(x='bottomright',lty=c(1,1,1,1),col=c(1,2,3,4),
#'    legend=c('Exact','1st approx.','2nd approx.','3rd approx'),bty='n')
#' plot(x, m*dnorm(x)*exp((m-1)*pnorm(x,log.p=TRUE)),type='l', ylab='Density',
#' main='Distribution of the maxima of\n 100 standard normal variates')
#' lines(x, evd::dgev(x,loc=p[1], scale=p[2], shape=0),col=2)
#' lines(x, evd::dgev(x,loc=p[1], scale=p[2], shape=p[3]),col=3)
#' lines(x, approx$f(x),col=4)
#' legend(x='topright',lty=c(1,1,1,1),col=c(1,2,3,4),
#'  legend=c('Exact','1st approx.','2nd approx.','3rd approx'),bty='n')
#'
#' #Threshold exceedances
#' par <- smith.penult(family = "norm", ddensF=function(x){-x*dnorm(x)},
#' method='pot', u=4, returnList=FALSE)
#' approx <- smith.penult.fn(loc=par[1], scale=par[2], shape=par[3],
#'  eps=par[3]^2+par[3]+par[2]^2, mdaGumbel=TRUE, method='pot')
#' x <- seq(4.01,7,by=0.01)
#' #Distribution function
#' plot(x, 1-(1-pnorm(x))/(1-pnorm(par[1])),type='l', ylab='Conditional CDF',
#' main='Exceedances over 4\n for standard normal variates')
#' lines(x, evd::pgpd(x, loc=par[1], scale=par[2], shape=0),col=2)
#' lines(x, evd::pgpd(x, loc=par[1], scale=par[2], shape=par[3]),col=3)
#' lines(x, approx$F(x),col=4)
#' #Density
#' plot(x, dnorm(x)/(1-pnorm(par[1])),type='l', ylab='Conditional density',
#' main='Exceedances over 4\n for standard normal variates')
#' lines(x, evd::dgpd(x, loc=par[1], scale=par[2], shape=0),col=2)
#' lines(x, evd::dgpd(x, loc=par[1], scale=par[2], shape=par[3]),col=3)
#' lines(x, approx$f(x),col=4)
#' @export
smith.penult.fn <- function(loc, scale, shape, eps, rho = NULL, method = c("bm", "pot"), mdaGumbel = FALSE, ...) {
  ellips <- list(...)
  #Compatibility condition with version 1.12 and before, via ellipsis argument
  if(!is.null(ellips$model) && length(method) == 2){
    method <- ellips$model
  }
    bn = loc
    an = scale
    gamma = shape
    # Functions appearing in the theory of slow variation
    hrho <- function(x, rho) {
        if (rho == 0) {
            log(x)
        } else {
            (x^rho - 1)/rho
        }
    }
    H <- function(x, eta, rho) {
        H0 <- function(x, eta) (0.5 * (log(pmax(1 + x * eta, 0))^2) - log(pmax(1 + x * eta, 0)) + 1 - 1/(1 + x * eta))/eta^3
        Hm1 <- function(x, eta) (log(pmax(1 + x * eta, 0))/(1 + x * eta) + log(pmax(1 + x * eta, 0)) - 2 * (1 - 1/(1 + x * eta)))/eta^3
        Hrho <- function(x, eta, rho) ifelse((1 + x * eta) > 0, hrho(1 + x * eta, rho) + rho * hrho(1 + x * eta, -1) - (rho + 1) *
            log(1 + x * eta), 0)
        if (rho == 0) {
            H0(x, eta)
        } else if (rho == -1) {
            Hm1(x, eta)
        } else {
            Hrho(x, eta, rho)
        }
    }
    # derivative of dH
    dH <- function(x, eta, rho) {
        ifelse((1 + x * eta) > 0, ((1 + x * eta)^(rho - 1) + rho * (1 + x * eta)^(-2) - (rho + 1)/(1 + x * eta))/(rho * (rho + 1) *
            eta^2), 0)
    }

    ## Block maxima - GEV-like distribution functions and densities Distribution function of third penultimate approximation
    if (method == "bm") {
        if (!mdaGumbel) {
            if (is.null(rho)) {
                stop("Invalid `rho' parameter")
            }
            GEV3rd <- function(x) {
                # , an, bn, gamma, eps, rho
                q <- (x - bn)/an
                if (gamma == 0)
                  p <- exp(-exp(-q) * (1 + eps * H(q, gamma, rho))) else p <- exp(-pmax(1 + gamma * q, 0)^(-1/gamma) * (1 + eps * H(q, gamma, rho)))
                invalid <- which(diff(p) < 0)
                p[invalid] <- 0
                p <- pmin(1, pmax(0, p))
                return(p)
            }
            # Density of third penultimate approximation , an, bn, gamma, eps, rho
            dGEV3rd <- function(x) {
                q <- (x - bn)/an
                if (gamma == 0)
                  stop("Not yet implemented")
                pmax(0, GEV3rd(x)/an * ((1 + eps * H(q, gamma, rho)) * exp((-1/gamma - 1) * log(pmax(1 + q * gamma, 0))) - exp((-1/gamma) *
                  log(pmax(1 + q * gamma, 0))) * eps * dH(q, gamma, rho)))
            }
            return(list(F = GEV3rd, f = dGEV3rd))
        } else {
            # mdaGumbel Distribution function of the third penultimate approximation, replacing H_rho(x, eta) by x^3/6 , an, bn, gamma, eps
            GEV3rda <- function(x) {
                q <- (x - bn)/an
                if (gamma == 0)
                  p <- exp(-exp(-q) * (1 + eps * q^3/6)) else p <- exp(-pmax(1 + gamma * q, 0)^(-1/gamma) * (1 + eps * q^3/6))
                invalid <- which(diff(p) < 0)
                p[invalid] <- 0
                p <- pmin(1, pmax(0, p))
                return(p)
            }
            # Density of third penultimate approximation, replacing H_rho(x, eta) by x^3/6 , an, bn, gamma, eps
            dGEV3rda <- function(x) {
                q <- (x - bn)/an
                if (gamma == 0)
                  stop("Not yet implemented")
                pmax(0, GEV3rda(x)/an * ((1 + eps * q^3/6) * exp((-1/gamma - 1) * log(pmax(1 + q * gamma, 0))) - exp((-1/gamma) *
                  log(pmax(1 + q * gamma, 0))) * eps * q^2/2))
            }
            return(list(F = GEV3rda, f = dGEV3rda))
        }
    } else if (method == "pot") {
        if (!mdaGumbel) {
            if (is.null(rho)) {
                stop("Invalid `rho' parameter")
            }
            # Peaks-over-threshold - GP-like distribution functions and densities , an, bn, gamma, eps, rho
            GP3rd <- function(x) {
                q <- (x - bn)/an
                p <- 1 - pmax(0, (1 + gamma * q))^(-1/gamma) * (1 + eps * H(q, gamma, rho))
                p[which(diff(p) < 0)] <- 0
                p <- pmin(1, pmax(0, p))
                p
            }
            dGP3rd <- function(x) {
                # , an, bn, gamma, eps, rho
                q <- (x - bn)/an
                pmax(0, pmax(0, (1 + gamma * q))^(-1/gamma - 1)/an * ((1 + eps * H(q, gamma, rho)) - pmax(0, (1 + gamma * q)) * eps *
                  dH(q, gamma, rho)))
            }
            return(list(F = GP3rd, f = dGP3rd))
        } else {

            GP3rda <- function(x) {
                # , an, bn, gamma, eps
                q <- (x - bn)/an
                p <- 1 - pmax(0, (1 + gamma * q))^(-1/gamma) * (1 + eps * q^3/6)
                p[which(diff(p) < 0)] <- 0
                p <- pmin(1, pmax(0, p))
                p
            }

            dGP3rda <- function(x) {
                # , an, bn, gamma, eps
                q <- (x - bn)/an
                pmax(0, pmax(0, (1 + gamma * q))^(-1/gamma - 1)/an * ((1 + eps * q^3/6) - pmax(0, (1 + gamma * q)) * eps * q^2/2))
            }
            return(list(F = GP3rda, f = dGP3rda))
        }
    }


}
