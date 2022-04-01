#' Automated mean residual life plots
#'
#' This function implements the automated proposal from
#' Section 2.2 of Langousis et al. (2016)
#' for mean residual life plots. It returns the threshold
#' that minimize the weighted mean square error and
#' moment estimators for the scale and shape parameter
#' based on weighted least squares.
#'
#' The procedure consists in estimating the usual
#'
#' @references Langousis, A., A. Mamalakis, M. Puliga and R. Deidda (2016).
#' \emph{Threshold detection for the generalized Pareto distribution:
#' Review of representative methods and application to the
#' NOAA NCDC daily rainfall database}, Water Resources Research, \strong{52}, 2659--2681.
#'
#' @return a list containing
#' \itemize{
#' \item{thresh}{selected threshold}
#' \item{scale}{scale parameter estimate}
#' \item{shape}{shape parameter estimate}
#' }
#'
#' @param xdat [numeric] vector of observations
#' @param kmax [integer] maximum number of order statistics
#' @param thresh [numeric] vector of thresholds; if missing, uses all order statistics from the 20th largest until \code{kmax} as candidates
#' @param plot [logical] if \code{TRUE}, return a plot of the mean residual life plot with the fitted slope
#' and the chosen threshold
automrl <- function(
      xdat,
      kmax,
      thresh,
      plot = FALSE){
   k <- 9
   # expected conditional exceedances based
   # on at least k+1 obs (default to 10, hardcoded)
   xdat <- sort(xdat[is.finite(xdat)], decreasing = TRUE)
   if(!missing(thresh)){
     stopifnot(isTRUE(all(is.finite(thresh))))
     thresh <- sort(thresh)
     nt <- length(thresh)
     xdat <- xdat[xdat > thresh[1]]
     if(isTRUE(findInterval(x = xdat[10], vec = thresh) < nt)){
       stop("Not enough observations for reliable estimation:\n users must provide at least 10 exceedances over largest threshold.")
     }
     user_thresh <- TRUE
   } else{
     if(!missing(kmax)){
        kmax <- as.integer(kmax)
        stopifnot(kmax <= length(xdat))
        xdat <- xdat[seq_len(kmax)]
     }
     # Use order statistics
     thresh <- xdat[-(1:20)]
     nt <- length(thresh)
     user_thresh <- FALSE
   }
   n <- length(xdat)
   # At least 20 observations
   stopifnot(n >= 20)
   # Compute running mean and variance efficiently
   meanX <- cumsum(xdat) / seq_along(xdat)
   # Running variance
   varX <- cumsum((xdat[-1] - meanX[-1])*
                    (xdat[-1] - meanX[-n])) /
      seq.int(1L, to = n-1, by = 1L)
   # Exclude 10 largest observations to ensure that
   # we have not too much variability
   excu <- (meanX[-n] - xdat[-1])[-(1:k)]
   weights <- (seq_len(n-1) / varX)[-(1:k)]
   xk <- xdat[-(1:(k+1))]
   # Containers
   mse <- numeric(nt)
   coefs <- matrix(0, nrow = nt, ncol = 2)
   # Could use rank-one update, but is improvement worth it?
   for(i in seq_along(thresh)){
      fit <- lm(excu ~ xk,
                weights = weights,
                subset = xk > thresh[i])
      mse[i] <- weighted.mean(x = fit$residuals^2,
                              w = weights[seq_len(nobs(fit))])
      coefs[i,] <- coef(fit)
   }
   # plot(x = thresh,
   #      y = sqrt(mse),
   #      type = "l",
   #      bty = "l")
   index <- which.min(mse)
   cthresh <- thresh[index]
   shape <- coefs[index,2]/(1 + coefs[index,2])
   scale <- coefs[index, 1]*(1-shape)+shape*xdat[k + index]
   if(plot){
      plot(x = xk,
           y = excu,
           ylab = "mean excess value",
           xlab = "observations",
           pch = 20,
           bty = "l")
     if(user_thresh){
      rug(thresh)
     }
      abline(v = cthresh, lty = 2)
      abline(a = coefs[index, 1],
             b = coefs[index, 2])
   }
   return(list(thresh = cthresh,
               scale = scale,
               shape = shape))
}
