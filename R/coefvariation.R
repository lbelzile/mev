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
#' \item{thresh}{value of threshold returned by the procedure, \code{NA} if the hypothesis is rejected at all thresholds}
#' \item{cthresh}{sorted vector of candidate thresholds}
#' \item{cindex}{index of selected threshold among \code{cthresh} or \code{NA} if none returned}
#' \item{pval}{bootstrap p-values, with \code{NA} if \code{lazy} and the p-value exceeds level at lower thresholds}
#' \item{shape}{shape parameter estimates}
#' \item{nexc}{number of exceedances of each threshold \code{cthresh}}
#' \item{method}{estimation method for the shape parameter}
#' }
#' @export
#' @references del Castillo, J. and M. Padilla (2016). \emph{Modelling extreme values by the residual coefficient of variation}, SORT, 40(\bold{2}), pp. 303--320.
cvselect <- function(
    xdat,
    thresh,
    method = c("mle", "wcv", "cv"),
    nsim = 999L,
    nthresh = 10L,
    level =  0.05,
    lazy = FALSE){
  method <- match.arg(method)
  xdat <- as.numeric(xdat[is.finite(xdat)])
  stopifnot(length(level) == 1L,
            is.finite(level),
            level >= 0,
            level <= 1,
            is.logical(lazy),
            length(lazy) == 1L
            )
  # Set grid of thresholds or order them and keep exceedances
  if(!missing(thresh)){
     thresh <- sort(thresh, decreasing = TRUE)
     nthresh <- length(thresh)
     xdat <- xdat[xdat >= min(thresh)]
     # Weight vector is survival probability of thresholds
     n <- length(xdat)
  } else{
     xdat <- xdat[xdat > 0]
     n <- length(xdat)
     stopifnot(length(nthresh) == 1L,
               is.finite(nthresh),
               nthresh > 1)
     nthresh <- as.integer(nthresh)
     # Set thresholds at 1-p^k for k=0, ..., nthresh quantiles
     # Ensure there are enough observations
     p <- (10/n)^(1/nthresh)
     # Survival probability of the thresholds
     survprob <- p^((nthresh-1):0)
     thresh <- quantile(xdat, 1 - survprob)
  }
  # Number of exceedances above thresholds
  nexc <- sapply(thresh, function(x){sum(xdat > x)})
  if(nexc[1] < 10){
     stop("Threshold is too small: the procedure requires at least 10 exceedances")
  }
  survprob <- nexc / n
  # Coefficient of variation
  coefvar <- function(x, na.rm = TRUE){
     sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm)
  }
  # Coefficient of variation for exceedances
  cv <- sapply(thresh, function(x){
     coefvar(xdat[xdat > x] - x)})

  if(method == "wcv"){
     cv0 <- cumsum(cv*survprob)/cumsum(survprob)
     shape <- 0.5*(1-1/cv0^2)
  } else if(method == "mle"){
     shape <- as.numeric(sapply(thresh, function(th){
        coef(mev::fit.gpd(xdat = xdat, threshold = th))[2]
     }))
     shape2cv <- function(shape){
        ifelse(shape < 0.5, (1-2*shape)^(-0.5), NA)
     }
     cv0 <- shape2cv(shape)
  } else{
     cv0 <- cv
     shape <- 0.5*(1-1/cv0^2)
  }
  if(isTRUE(any(shape > 0.25))){
     warning("Estimated shape parameter larger than 1/4.")
  }
  stat <- sapply(seq_along(shape), function(i){
     sum(nexc[1:i]*(cv[1:i] - cv0[i])^2)
  })
  # Function to loop over
  compute_cvstat <- function(xdat, survprob, method = "mle"){
     # n*plevel is number of exceedances above threshold
     survprob <- sort(survprob, decreasing = TRUE)
     n <- length(xdat)
     # Thresholds are ordered from largest to smallest
     thresh <- quantile(xdat, 1 - survprob)
     # Compute coefficient of variation
     cv <- sapply(thresh, function(x){
        coefvar(xdat[xdat > x] - x)})
     # Logical: compute by minimizing stat
     if(method == "mle"){
        shape <- as.numeric(coef(mev::fit.gpd(xdat = xdat, threshold = 0))[2])
      if(!isTRUE(shape < 0.5)){return(NA)}
        # Null value
        cv0 <- (1-2*shape)^(-0.5)
     } else if (method == "cv"){
        cv0 <- cv[1]
     } else if(method == "wcv"){
        cv0 <- cumsum(cv*survprob)/cumsum(survprob)
     }
     n*sum(survprob*(cv - cv0)^2)
  }
  # Bootstrap function
  bootfun <- function(n, shape, survprob, method = "mle"){
     dat <- revdbayes::rgp(n = n, loc = 0, scale = 1, shape = shape)
     compute_cvstat(xdat = dat, method = method, survprob = survprob)
  }
  boot_pval <- rep(NaN, nthresh)
  thselect <- NaN
  cindex <- NULL
  # Bootstrap loop
  for(r in seq_along(thresh)){
     # Number of exceedances, null shape and value of the statistic
     # return the bootstrap p-value
     rk <- nthresh - r + 1
     boot_stat <- replicate(n = nsim,
           bootfun(n = nexc[rk],
                   survprob = nexc[1:rk]/nexc[rk],
                   shape = shape[rk]))
     boot_pval[r] <- mean(c(boot_stat >= stat[rk], TRUE), na.rm = TRUE)
   if(boot_pval[r] >= level[1] & is.na(thselect)){
      thselect <- thresh[rk]
      cindex <- r
      if(lazy){break}
   }
  }
  list(thresh = as.numeric(thselect),
       cthresh = as.numeric(rev(thresh)),
       cindex = cindex,
       pval = as.numeric(boot_pval),
       shape = as.numeric(rev(shape)),
       nexc = as.numeric(rev(nexc)),
       method = method)
 # TODO add print method
}
