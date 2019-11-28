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
#' @param dat a vector of observations
#' @param thresh a vector of candidate thresholds at which to compute the estimates.
#' @param method string indicating the method for computing confidence or credible intervals.
#' Must be one of \code{"wald"}, \code{"profile"} or \code{"post"}.
#' @param level confidence level of the intervals. Default to 0.95.
#' @param plot logical; should parameter stability plots be displayed? Default to \code{TRUE}.
#' @param ... additional arguments passed to \code{plot}.
#' @return a list with components
#' \itemize{
#' \item{\code{threshold}:} vector of numerical threshold values.
#' \item{\code{mle}:} matrix of modified scale and shape maximum likelihood estimates.
#' \item{\code{lower}:} matrix of lower bounds for the confidence or credible intervals.
#' \item{\code{upper}:} matrix of lower bounds for the confidence or credible intervals.
#' \item{\code{method}:} method for the confidence or coverage intervals.
#' }
#' @author Leo Belzile
#' @seealso \code{\link[ismev]{gpd.fitrange}}
#' @return plots of the modified scale and shape parameters, with pointwise confidence/credible intervals
#' and an invisible data frame containing the threshold \code{thresh} and the modified scale and shape parameters.
#' @export
#' @examples
#' dat <- abs(rnorm(10000))
#' u <- qnorm(seq(0.9,0.99, by= 0.01))
#' tstab.gpd(dat = dat, thresh = u)
#' \dontrun{
#' tstab.gpd(dat = dat, thresh = u, method = "profile")
#' tstab.gpd(dat = dat, thresh = u, method = "post")
#' }
tstab.gpd <- function(dat, thresh, method = c("wald", "profile", "post"), level = 0.95, plot = TRUE, ...){
  args <- list(...)
  if(missing(dat) && !is.null(args$xdat)){
    dat <- args$xdat
  }
  dat <- as.vector(dat)
  thresh <- unique(sort(thresh))
  stopifnot(length(level) == 1, level > 0, level < 1)
  method <- match.arg(method)
  if (method == "post") {
    if (!requireNamespace("revdbayes", quietly = TRUE)) {
      stop("Package \"revdbayes\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }

  alpha <- (1-level)
  if(length(thresh) < 2){
    stop("Invalid threshold: should be of length greater than 2")
  }
  # if(!is.vector(dat) || inherits(dat, "numeric")){
  #   stop("Invalid data")
  # }
  nt <- length(thresh)
  np <- 2
  parmat <- matrix(0, ncol = np, nrow = nt)
  confintmat <- matrix(0, ncol = 2 * np, nrow = nt)
  #Some quantities for profile
  if(method == "profile"){
    xmax <- max(dat)
    pllsigmainv <- function(par, sigmat, dat, thresh){
      sigma <- sigmat + par*thresh
      nll <- -gpd.ll(par = c(sigma, par), dat = dat)
      if(!is.finite(nll)){
        return(-1e10)
      } else{
        return(nll)
      }
    }
  }
  for(i in 1:nt){
    gpdu <- fit.gpd(xdat = dat, threshold = thresh[i])
    parmat[i,] <- gpdu$estimate
    parmat[i,1] <- parmat[i,1] - parmat[i,2]*(thresh[i]-thresh[1])
    stderr.transfo <- sqrt(t(c(1, -thresh[i]+thresh[1])) %*% gpdu$vcov %*% c(1, -thresh[i]+thresh[1]))[1,1]
   if(method == "wald"){
    confintmat[i,3] <- gpdu$estimate['shape'] - gpdu$std.err['shape']*qnorm(1-alpha/2)
    confintmat[i,4] <- gpdu$estimate['shape'] + gpdu$std.err['shape']*qnorm(1-alpha/2)
    confintmat[i,1] <- parmat[i,1] - stderr.transfo * qnorm(1-alpha/2)
    confintmat[i,2] <- parmat[i,1] + stderr.transfo * qnorm(1-alpha/2)
   } else if(method == "profile"){
     profxi <- gpd.pll(param = "shape", mod = "profile", mle = gpdu$estimate, dat = gpdu$exceedances, plot = FALSE)
     confintmat[i,3:4] <- confint(profxi, level = level, print = FALSE)[2:3]
     k <- 30L
     prof_vals <- rep(0, k)
     xi_sigma_vals <- rep(0, k)
     grid_psi <- seq(parmat[i,1] - 3 * stderr.transfo, parmat[i,1] + 3.5 * stderr.transfo, length = k)
     xmaxui <- xmax - thresh[i]
     #Profile for scale := sigma_u - xi (u - u_0)
      for(j in 1:k){
        opt_prof <- optimize(f = pllsigmainv, upper = 1.5,
                             lower = max(-grid_psi[j]/(thresh[i]-thresh[1]), -grid_psi[j]/(xmaxui+thresh[i]-thresh[1]))+1e-10,
                             sigmat = grid_psi[j], dat = gpdu$exceedances, thresh = thresh[i]-thresh[1])
        xi_sigma_vals[j] <- opt_prof$minimum
        prof_vals[j] <- opt_prof$objective
      }
     prof <- structure(list(psi = grid_psi, psi.max = parmat[i,1], pll = -prof_vals,
                            maxpll = -gpdu$nllh, std.err = stderr.transfo), class = "eprof")
     confintmat[i, 1:2] <- confint(prof, level = level, print = FALSE)[2:3]
   } else if(method == "post"){
      postsim <- suppressWarnings(revdbayes::rpost_rcpp(n = 1000, model = "gp", init_ests = gpdu$estimate,
                                     data = c(-1, gpdu$exceedances), trans = "BC",
                                     prior = revdbayes::set_prior(prior = "norm",
                                                                  model = "gp",
                                                                  mean = rep(0, 2),
                                                                  cov = c(100*gpdu$estimate[1], 1)*diag(np)))$sim_vals)
    postsim[,1] <- postsim[,1] - postsim[,2]*(thresh[i]-thresh[1])
    confintmat[i,] <- apply(postsim, 2, quantile, c(alpha/2, 1-alpha/2))
    }

  }
  lower <- confintmat[,c(1,3)]
  upper <- confintmat[,c(2,4)]
  colnames(parmat) <- colnames(lower) <- colnames(upper) <-  c("modif. scale", "shape")
  ret <- structure(list(threshold = thresh, mle = parmat, lower = lower,
                                    upper = upper, method = method, level = level), class = "mev_tstab.gpd")
  if(plot){
    plot(ret, ...)
  }
  return(invisible(ret))
}

#'@export
plot.mev_tstab.gpd <- function(x, which = 1:2, ...){
  oldpar <- par(no.readonly = TRUE)
  if(length(which) == 2){
    par(mfrow = c(2,1), mar = c(4,4.5,3,1))
  }
  on.exit(par(oldpar))
  ellipsis <- list(...)
  names_ell <- names(ellipsis)
  if("main" %in% names_ell){
    main <- ellipsis$main
    ellipsis$main <- NULL
  } else{
    main <- c("Parameter stability plot","");
  }
  if(!"sub" %in% names_ell){
   sub <- c(switch(x$method,
                 wald = paste("Wald", x$level*100, "% pointwise confidence intervals"),
                 profile = paste("profile likelihood", x$level*100, "% pointwise confidence intervals"),
                 post = paste(x$level*100, "% pointwise credible intervals")), "")
  } else{
    sub <- ellipsis$sub
    ellipsis$sub <- NULL
  }
  if(isTRUE(all.equal(which, 2, check.attributes = FALSE))){
    main[2] <- main[1]
    sub[2] <- sub[1]
  }
  if(!"ylab" %in% names_ell){
   ylab <- c("modif. scale", "shape")
  } else{
   ylab <- ellipsis$ylab
   if(length(ylab)!= length(which)){
     stop(paste("Invalid `ylab` argument: must be a vector of size", length(which)))
   }
   ellipsis$ylab <- NULL
   }
   if(!"xlab" %in% names_ell){
    xlab <- "threshold"
   } else{
    xlab <- ellipsis$xlab
    ellipsis$xlab <- NULL
   }
   if(!"bty" %in% names_ell){
    bty <- "l"
   } else{
    bty <- ellipsis$bty
    ellipsis$bty <- NULL
   }
   if(!"pch" %in% names_ell){
     pch <- 20
   } else{
     pch <- ellipsis$pch
     ellipsis$pch <- NULL
   }
   #modified scale
  if(1 %in% which){
    ylim <- c(min(x$lower[,1]), max(x$upper[,1]))
    pars <- list(x = x$threshold, y =  x$mle[,1], pch = pch, ylab = ylab[1], xlab = xlab,
                 ylim = ylim, bty = bty, main = main[1])
    do.call(plot, c(pars, ellipsis))
    mtext(side = 3, line = -0.2, adj = 0.5, text = sub[1])
    for(i in 1:length(x$threshold)){
      lines(c(x$threshold[i], x$threshold[i]), c(x$lower[i,1], x$upper[i,1]))
    }
    }
  if(2 %in% which){
  ylim <- c(min(x$lower[,2]), max(x$upper[,2]))
  #shape
  pars <- list(x = x$threshold, y = x$mle[,2], pch = pch, ylab = ylab[2], xlab = xlab,
  ylim = ylim, bty = bty, main = main[2])
  do.call(plot,  c(pars,ellipsis))
  mtext(side = 3, line = -0.2, adj = 0.5, text = sub[2])
  for(i in 1:length(x$threshold)){
    lines(c(x$threshold[i], x$threshold[i]), c(x$lower[i,2], x$upper[i,2]))
  }
  }
}
