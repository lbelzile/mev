
#' Generalized Pareto maximum likelihood estimates for various quantities of interest
#'
#' This function calls the \code{fit.gpd} routine on the sample of excesses and returns maximum likelihood
#' estimates for all quantities of interest, including scale and shape parameters, quantiles and value-at-risk,
#' expected shortfall and mean and quantiles of maxima of \code{N} threshold exceedances
#'
#' @param xdat sample vector of excesses
#' @param args vector of strings indicating which arguments to return the maximum likelihood values for
#' @param m number of observations of interest for return levels. Required only for \code{args} values \code{'VaR'} or \code{'ES'}
#' @param N size of block over which to take maxima. Required only for \code{args} \code{Nmean} and \code{Nquant}.
#' @param p tail probability, equivalent to \eqn{1/m}. Required only for \code{args} \code{quant}.
#' @param q level of quantile for N-block maxima. Required only for \code{args} \code{Nquant}.
#' @return named vector with maximum likelihood values for arguments \code{args}
#' @export
#' @examples
#' xdat <- evd::rgpd(n = 30, shape = 0.2)
#' gpd.mle(xdat = xdat, N = 100, p = 0.01, q = 0.5, m = 100)
gpd.mle <- function(xdat,
                    args = c("scale",
                             "shape",
                             "quant",
                             "VaR",
                             "ES",
                             "Nmean",
                             "Nquant"),
                    m,
                    N,
                    p,
                    q) {
  args <- match.arg(args, c("scale", "shape", "quant", "VaR", "ES", "Nmean", "Nquant"), several.ok = TRUE)
  fitted <- try(gp.fit(xdat = na.omit(as.vector(xdat)), threshold = 0, method = "Grimshaw"))
  sigma <- fitted$estimate[1]
  xi <- fitted$estimate[2]
  # Does not handle the case xi=0 because the optimizer does not return this value!
  a <- sapply(args, switch, scale = sigma, shape = xi, quant = sigma/xi * (p^(-xi) - 1),
              Nquant = sigma/xi * ((1 - q^(1/N))^(-xi) - 1),
              Nmean = (exp(lgamma(N + 1) + lgamma(1 - xi) - lgamma(N + 1 - xi)) - 1) * sigma/xi,
              VaR = sigma/xi * (m^xi - 1),
              ES = ifelse(xi < 1, (sigma/xi * (m^xi - 1) + sigma)/(1 - xi), Inf))
  a <- as.vector(unlist(a))
  names(a) = args
  return(a)
}




#' Generalized extreme value maximum likelihood estimates for various quantities of interest
#'
#' This function calls the \code{fit.gev} routine on the sample of block maxima and returns maximum likelihood
#' estimates for all quantities of interest, including location, scale and shape parameters, quantiles and mean and
#' quantiles of maxima of \code{N} blocks.
#' @export
#' @param xdat sample vector of maxima
#' @param args vector of strings indicating which arguments to return the maximum likelihood values for.
#' @param N size of block over which to take maxima. Required only for \code{args} \code{Nmean} and \code{Nquant}.
#' @param p tail probability. Required only for \code{arg} \code{quant}.
#' @param q level of quantile for maxima of \code{N} exceedances. Required only for \code{args} \code{Nquant}.
#' @return named vector with maximum likelihood estimated parameter values for arguments \code{args}
#' @examples
#' dat <- evd::rgev(n = 100, shape = 0.2)
#' gev.mle(xdat = dat, N = 100, p = 0.01, q = 0.5)
#'
gev.mle <- function(xdat, args = c("loc", "scale", "shape", "quant", "Nmean", "Nquant"), N, p, q) {
  args <- match.arg(args, c("loc", "scale", "shape", "quant", "Nmean", "Nquant"), several.ok = TRUE)

  if(missing(q) && "Nquant" %in% args){
    stop("Argument \"q\" missing for \"Nquant\"")
  }
  if(missing(p) && "quant" %in% args){
    stop("Argument \"p\" missing for \"quant\"")
  }
  if(missing(N) && any(c("Nmean", "Nquant") %in% args)){
    stop("Argument \"N\" missing for \"Nquant\" or \"Nmean\"")
  }
  fitted <- suppressWarnings(fit.gev(xdat = xdat))
  mu <- fitted$estimate[1]
  sigma <- fitted$estimate[2]
  xi <- fitted$estimate[3]
  # Does not handle the case xi=0 because the optimizer does not return this value!
  a <- sapply(args, switch, loc = mu, scale = sigma, shape = xi,
              quant = evd::qgev(p = 1 -p, loc = mu, scale = sigma, shape = xi),
              Nquant = ifelse(xi != 0,
                              mu - sigma/xi * (1 - (N/log(1/q))^xi),
                              mu + sigma * (log(N) - log(log(1/q)))),
              Nmean = ifelse(xi != 0,
                             mu - sigma/xi * (1 - N^xi * gamma(1 - xi)),
                             mu + sigma * (log(N) - psigamma(1))))
  a <- as.vector(unlist(a))
  names(a) <- args
  a
}

#' Maximum likelihood estimation for the generalized Pareto distribution
#'
#' Numerical optimization of the generalized Pareto distribution for
#' data exceeding \code{threshold}.
#' This function returns an object of class \code{mev_gpd}, with default methods for printing and quantile-quantile plots.
#'
#' @param xdat a numeric vector of data to be fitted.
#' @param threshold the chosen threshold.
#' @param show logical; if \code{TRUE} (the default), print details of the fit.
#' @param method the method to be used. See \bold{Details}. Can be abbreviated.
#' @param MCMC \code{NULL} for frequentist estimates, otherwise a boolean or a list with parameters passed. If \code{TRUE}, runs a Metropolis-Hastings sampler to get posterior mean estimates. Can be used to pass arguments \code{niter}, \code{burnin} and \code{thin} to the sampler as a list.
#' @param k bound on the influence function (\code{method = "obre"}); the constant \code{k} is a robustness parameter
#' (higher bounds are more efficient, low bounds are more robust). Default to 4, must be larger than \eqn{\sqrt{2}}.
#' @param tol numerical tolerance for OBRE weights iterations (\code{method = "obre"}). Default to \code{1e-8}.
#' @param fpar a named list with fixed parameters, either \code{scale} or \code{shape}
#' @param warnSE logical; if \code{TRUE}, a warning is printed if the standard errors cannot be returned from the observed information matrix when the shape is less than -0.5.
#' @seealso \code{\link[evd]{fpot}} and \code{\link[ismev]{gpd.fit}}
#'
#' @details The default method is \code{'Grimshaw'}, which maximizes the profile likelihood for the ratio scale/shape.  Other options include \code{'obre'} for optimal \eqn{B}-robust estimator of the parameter of Dupuis (1998), vanilla maximization of the log-likelihood using constrained optimization routine \code{'auglag'}, 1-dimensional optimization of the profile likelihood using \code{\link[stats]{nlm}} and \code{\link[stats]{optim}}. Method \code{'ismev'} performs the two-dimensional optimization routine \code{\link[ismev]{gpd.fit}} from the \code{\link[ismev]{ismev}} library, with in addition the algebraic gradient.
#' The approximate Bayesian methods (\code{'zs'} and \code{'zhang'}) are extracted respectively from Zhang and Stephens (2009) and Zhang (2010) and consists of a approximate posterior mean calculated via importance
#' sampling assuming a GPD prior is placed on the parameter of the profile likelihood.
#' @note Some of the internal functions (which are hidden from the user) allow for modelling of the parameters using covariates. This is not currently implemented within \code{gp.fit}, but users can call internal functions should they wish to use these features.
#' @author Scott D. Grimshaw for the \code{Grimshaw} option. Paul J. Northrop and Claire L. Coleman for the methods \code{optim}, \code{nlm} and \code{ismev}.
#' J. Zhang and Michael A. Stephens (2009) and Zhang (2010) for the \code{zs} and \code{zhang} approximate methods and L. Belzile for methods \code{auglag} and \code{obre}, the wrapper and MCMC samplers.
#'
#' If \code{show = TRUE}, the optimal \eqn{B} robust estimated weights for the largest observations are printed alongside with the
#' \eqn{p}-value of the latter, obtained from the empirical distribution of the weights. This diagnostic can be used to guide threshold selection:
#' small weights for the \eqn{r}-largest order statistics indicate that the robust fit is driven by the lower tail
#' and that the threshold should perhaps be increased.
#'
#' @references Davison, A.C. (1984). Modelling excesses over high thresholds, with an application, in
#' \emph{Statistical extremes and applications}, J. Tiago de Oliveira (editor), D. Reidel Publishing Co., 461--482.
#' @references Grimshaw, S.D. (1993). Computing Maximum Likelihood Estimates for the Generalized
#'  Pareto Distribution, \emph{Technometrics}, \bold{35}(2), 185--191.
#' @references Northrop, P.J. and C. L. Coleman (2014). Improved threshold diagnostic plots for extreme value
#' analyses, \emph{Extremes}, \bold{17}(2), 289--303.
#' @references Zhang, J. (2010). Improving on estimation for the generalized Pareto distribution, \emph{Technometrics} \bold{52}(3), 335--339.
#' @references Zhang, J.  and M. A. Stephens (2009). A new and efficient estimation method for the generalized Pareto distribution.
#' \emph{Technometrics} \bold{51}(3), 316--325.
#' @references Dupuis, D.J. (1998). Exceedances over High Thresholds: A Guide to Threshold Selection,
#' \emph{Extremes}, \bold{1}(3), 251--261.
#'
#' @return If \code{method} is neither \code{'zs'} nor \code{'zhang'}, a list containing the following components:
#' \itemize{
#' \item \code{estimate} a vector containing the \code{scale} and \code{shape} parameters (optimized and fixed).
#' \item \code{std.err} a vector containing the standard errors. For \code{method = "obre"}, these are Huber's robust standard errors.
#' \item \code{vcov} the variance covariance matrix, obtained as the numerical inverse of the observed information matrix. For \code{method = "obre"},
#' this is the sandwich Godambe matrix inverse.
#' \item \code{threshold} the threshold.
#' \item \code{method} the method used to fit the parameter. See details.
#' \item \code{nllh} the negative log-likelihood evaluated at the parameter \code{estimate}.
#' \item \code{nat} number of points lying above the threshold.
#' \item \code{pat} proportion of points lying above the threshold.
#' \item \code{convergence} components taken from the list returned by \code{\link[stats]{optim}}.
#' Values other than \code{0} indicate that the algorithm likely did not converge (in particular 1 and 50).
#' \item \code{counts} components taken from the list returned by \code{\link[stats]{optim}}.
#' \item \code{exceedances} excess over the threshold.
#' }
#' Additionally, if \code{method = "obre"}, a vector of OBRE \code{weights}.
#'
#' Otherwise, a list containing
#' \itemize{
#' \item \code{threshold} the threshold.
#' \item \code{method} the method used to fit the parameter. See \bold{Details}.
#' \item \code{nat} number of points lying above the threshold.
#' \item \code{pat} proportion of points lying above the threshold.
#' \item \code{approx.mean} a vector containing containing the approximate posterior mean estimates.
#' }
#' and in addition if MCMC is neither \code{FALSE}, nor \code{NULL}
#' \itemize{
#' \item \code{post.mean} a vector containing the posterior mean estimates.
#' \item \code{post.se} a vector containing the posterior standard error estimates.
#' \item \code{accept.rate} proportion of points lying above the threshold.
#' \item \code{niter} length of resulting Markov Chain
#' \item \code{burnin} amount of discarded iterations at start, capped at 10000.
#' \item \code{thin} thinning integer parameter describing
#' }
#'
#' @export
#'
#' @examples
#' data(eskrain)
#' fit.gpd(eskrain, threshold = 35, method = 'Grimshaw', show = TRUE)
#' fit.gpd(eskrain, threshold = 30, method = 'zs', show = TRUE)
fit.gpd <- function(xdat,
                    threshold = 0,
                    method = "Grimshaw",
                    show = FALSE,
                    MCMC = NULL,
                    k = 4,
                    tol = 1e-8,
                    fpar = NULL,
                    warnSE = FALSE){
 if(!method == "obre"){
   gp.fit(xdat = na.omit(as.vector(xdat)),
          threshold = threshold,
          method = method,
          show = show,
          MCMC = MCMC,
          fpar = fpar,
          warnSE = warnSE)
 } else{
   if(!is.null(fpar)){
     warning("\"fpar\" argument ignored for OBRE method.")
   }
   .fit.gpd.rob(dat = na.omit(as.vector(xdat)),
                thresh = threshold,
                show = show,
                k = k,
                tol = tol)
 }
}

#' Maximum likelihood estimation of the point process of extremes
#'
#' Data above \code{threshold} is modelled using the limiting point process
#' of extremes.
#' @inheritParams fit.gpd
#' @param start named list of starting values
#' @param npp number of observation per period. See \bold{Details}
#' @param fpar a named list with optional fixed components \code{loc}, \code{scale} and \code{shape}
#' @param np number of periods of data, if \code{xdat} only contains exceedances.
#' @details The parameter \code{npp} controls the frequency of observations.
#' If data are recorded on a daily basis, using a value of \code{npp = 365.25}
#' yields location and scale parameters that correspond to those of the
#'  generalized extreme value distribution fitted to block maxima.
#'
#' @references Coles, S. (2001), An introduction to statistical modelling of extreme values. Springer : London, 208p.
#'
#' @return a list containing the following components:
#' \itemize{
#' \item \code{estimate} a vector containing all parameters (optimized and fixed).
#' \item \code{std.err} a vector containing the standard errors.
#' \item \code{vcov} the variance covariance matrix, obtained as the numerical inverse of the observed information matrix.
#' \item \code{threshold} the threshold.
#' \item \code{method} the method used to fit the parameter. See details.
#' \item \code{nllh} the negative log-likelihood evaluated at the parameter \code{estimate}.
#' \item \code{nat} number of points lying above the threshold.
#' \item \code{pat} proportion of points lying above the threshold.
#' \item \code{convergence} components taken from the list returned by \code{\link[stats]{optim}}.
#' Values other than \code{0} indicate that the algorithm likely did not converge (in particular 1 and 50).
#' \item \code{counts} components taken from the list returned by \code{\link[stats]{optim}}.
#' }
#' @export
#' @examples
#' data(eskrain)
#' pp_mle <- fit.pp(eskrain, threshold = 30, np = 6201)
#' plot(pp_mle)
fit.pp <- function(xdat,
                   threshold = 0,
                   npp = 1,
                   np = NULL,
                   method = c("nlminb", "BFGS"),
                   start = NULL,
                   show = FALSE,
                   fpar = NULL,
                   warnSE = FALSE){
  xdat <- as.vector(xdat)
  method <- match.arg(method)
  xdat <- xdat[is.finite(xdat)]
  n <- length(xdat)
  if (length(threshold) != 1 || mode(threshold) !=  "numeric")
    stop("\"threshold\" must be a numeric value")
  u <- as.double(threshold)
  xdatu <- xdat[xdat > u] #keep data above
  nu <- length(xdatu) #number above
  if(is.null(np)){
    np <- n/npp
  }
  # Fixed parameters
  param_names <- c("loc", "scale", "shape")
  stopifnot(is.null(fpar) | is.list(fpar))
  wf <- (param_names %in% names(fpar))
  if(sum(wf) == 3L){
    stop("Invalid input: all of the model parameters are fixed.")
  }
  if(is.list(fpar) && (length(fpar) >= 1L)){ #NULL has length zero
    if(is.null(names(fpar))){
      stop("\"fpar\" must be a named list")
    }
    if(!isTRUE(all(names(fpar) %in% param_names))){
      stop("Unknown fixed parameter: must be one of \"loc\",\"scale\" or \"shape\". ")
    }
    if(!isTRUE(all(unlist(lapply(fpar, length)) == rep(1L, sum(wf))))){
      stop("Each fixed parameter must be of length one.")
    }
  }
  spar <- vector(mode = "numeric", length = 3L)
  names(spar) <- param_names
  for(i in seq_along(fpar)){
    spar[names(fpar[i])] <- unlist(fpar[i])[1]
  }
  # Without covariates, we have (almost) exactly np*(1+xi/sigma*(u-mu))^(-1/xi)=nu
  # this follows from the Poisson approximation to the binomial
  # the mle for the latter is known (sample proportion of exceedance)
  # So we can effectively reduce the dimension of the optimization from 3 to 2 parameters
 pp.nll <- function(par, fpar, wf, xdat, u, np){
    param <- numeric(length = 3L)
    param[!wf] <- par
    param[wf] <- fpar
    nll <- -pp.ll(par = param, dat = xdat, u = u, np = np)
    ifelse(is.finite(nll), nll, 1e10)
 }
 pp.ngr <- function(par, fpar, wf, xdat, u, np){
   param <- numeric(length = 3L)
   param[!wf] <- par
   param[wf] <- fpar
   grad <- -pp.score(par = param, dat = xdat, u = u, np = np)[!wf]
   ifelse(is.finite(grad), grad, 1e10)
 }
 xmax <- max(xdatu); xmin <- min(xdatu)
 #hin Inequalities for Augmented Lagrangian
 pp.hin <- function(par, fpar, wf, xdat, u, np){
   param <- numeric(length = 3L)
   param[!wf] <- par
   param[wf] <- fpar
   c(param[2], param[3] + 1,
     (u - param[1])*param[3] + param[2],
     (xmin - param[1])*param[3] + param[2],
     (xmax - param[1])*param[3] + param[2])
 }

 # Starting values
 if(!is.null(start) && is.list(start)){
   if(!isTRUE(all(param_names[!wf] %in% names(start)))){
     stop(paste("Invalid starting value: named list must have components",
                paste(param_names[!wf], collapse = ", ")))
   }
   for(i in seq_along(start)){
     if(names(start[i]) %in% param_names[!wf]){
      spar[names(start[i])] <- unlist(start[i])[1]
     }
   }
   if(!isTRUE(all(pp.hin(par = spar[!wf], fpar = spar[wf], wf = wf, u = u, xdat = xdatu, np = np)>0))){
     stop("Starting values do not satisfy the inequality constraints.")
   }
 } else{
   gppars <- suppressWarnings(fit.gpd(xdat = xdatu, threshold = threshold)$estimate)
   sigma_init <- gppars['scale']*(length(xdatu)/np)^(gppars['shape'])
   mu_init <-  threshold - sigma_init*(((length(xdatu)/np)^(-gppars['shape']))-1)/gppars['shape']
   spar0 <- c(mu_init, sigma_init, gppars['shape'])
   if(all(!wf)){ # no missing values
     spar <- spar0
   } else{
    # Normal approximation around MLE (quadratic function)
   umle <- c(spar0[1], log(spar0[2]), spar0[3])
   # Compute precision matrix at MLE
   prec <- diag(c(1, spar0[2], 1)) %*% pp.infomat(par = spar0, dat = xdat, u = u, np = np) %*% diag(c(1, spar0[2], 1))
   # Best linear prediction
   spar[!wf] <- as.vector(c(umle[!wf] - solve(prec[which(!wf), which(!wf), drop = FALSE]) %*% prec[which(!wf), which(wf), drop = FALSE] %*% (spar[wf] - spar0[wf])))
   if(!wf[2]){
     # Backtransform scale if not fixed
     spar[2] <- exp(spar[2])
   }
   }
   # Check the starting values are feasible
   if(!isTRUE(all(pp.hin(par = spar[!wf], fpar = spar[wf], wf = wf, u = u, xdat = xdatu, np = np)>0))){
     stop("Starting values do not satisfy the inequality constraints.")
   }
 }
   # check_init_ll <- try(pp.ll(par = spar, dat = xdat, u = u, np = np))
   # if(inherits(check_init_ll, "try-error")){
   #   stop("Invalid starting values")
   #   # Invalid starting values
   # }
  # Optimization - basically started at MLE
 mle <- suppressWarnings(
   alabama::auglag(par = spar[!wf],
                   fpar = spar[wf],
                   wf = wf,
                   fn = pp.nll,
                   gr = pp.ngr,
                   hin = pp.hin,
                   u = u,
                   xdat = xdatu,
                   np = np,
                   control.outer = list(trace = FALSE,
                                        method = method)))
 if((mle$convergence == 0  || isTRUE(all.equal(mle$gradient, rep(0, sum(wf)), tolerance = 1e-3))) && isTRUE(all(mle$kkt1, mle$kkt2)) ){
   mle$convergence <- "successful"
 } else if(!wf[3] && isTRUE(all.equal(mle$par['shape'], -1, check.attributes = FALSE, tolerance = 1e-6))){
   mle$convergence <- "successful"
 } else {
  warning("Optimization routine may not have succeeded.")
 }
 wfo <- order(c(which(!wf), which(wf)))
 notf <- sum(!wf)
 fitted <- list()
 fitted$estimate <- mle$par
 fitted$param <- c(mle$par, spar[wf])[wfo]
 if(fitted$param[3] > -0.5){
   fitted$vcov <- try(solve(pp.infomat(par = fitted$param, u = u, np = np, dat = xdatu, method = "obs")[!wf,!wf]))
   fitted$std.err <- try(sqrt(diag(fitted$vcov)))
   if(inherits(fitted$std.err, what = "try-error")){
     fitted$vcov <- NULL
     fitted$std.err <- rep(NA, notf)
     if(warnSE){
       warning("Cannot calculate standard error based on observed information")
     }
   }
 } else{
   if(warnSE){
     warning("Cannot calculate standard error based on observed information")
   }
   fitted$vcov <- NULL
   fitted$std.err <- rep(NA, notf)
 }

 fitted$nllh <- mle$value
 names(fitted$param)<-  param_names
 names(fitted$estimate) <- names(fitted$std.err) <-  param_names[!wf]
 fitted$convergence <- mle$convergence
 fitted$counts <- mle$counts
 fitted$threshold <- u
 fitted$np <- np
 fitted$npp <- npp
 fitted$nat <- nu
 fitted$pat <- nu/n
 fitted$xdat <- xdat
 fitted$exceedances <- xdatu
 fitted$start <- spar
 fitted$wfixed <- wf
 class(fitted) <- c("mev_pp")
 if(show){
   print(fitted)
 }
 return(invisible(fitted))
}


#' Maximum likelihood estimation for the generalized extreme value distribution
#'
#' This function returns an object of class \code{mev_gev}, with default methods for printing and quantile-quantile plots.
#' The default starting values are the solution of the probability weighted moments.
#' @inheritParams gp.fit
#' @export
#' @importFrom alabama auglag
#' @param fpar a named list with optional fixed components \code{loc}, \code{scale} and \code{shape}
#' @param start named list of starting values
#' @param method string indicating the outer optimization routine for the augmented Lagrangian. One of \code{nlminb} or \code{BFGS}.
#' @return a list containing the following components:
#' \itemize{
#' \item \code{estimate} a vector containing the maximum likelihood estimates.
#' \item \code{std.err} a vector containing the standard errors.
#' \item \code{vcov} the variance covariance matrix, obtained as the numerical inverse of the observed information matrix.
#' \item \code{method} the method used to fit the parameter.
#' \item \code{nllh} the negative log-likelihood evaluated at the parameter \code{estimate}.
#' \item \code{convergence} components taken from the list returned by \code{\link[alabama]{auglag}}.
#' Values other than \code{0} indicate that the algorithm likely did not converge.
#' \item \code{counts} components taken from the list returned by \code{\link[alabama]{auglag}}.
#' \item \code{xdat} vector of data
#' }
#' @examples
#' xdat <- evd::rgev(n = 100)
#' fit.gev(xdat, show = TRUE)
#' # Example with fixed parameter
#' fit.gev(xdat, show = TRUE, fpar = list(shape = 0))
fit.gev <- function(xdat,
                    start = NULL,
                    method = c("nlminb","BFGS"),
                    show = FALSE,
                    fpar = NULL,
                    warnSE = FALSE){
  fitted <- list() # container
  param_names <- c("loc", "scale", "shape")
  stopifnot(is.null(fpar) | is.list(fpar))
  wf <- (param_names %in% names(fpar))
  if(sum(wf) == 3L){
    stop("Invalid input: all of the model parameters are fixed.")
  }
  if(is.list(fpar) && (length(fpar) >= 1L)){ #NULL has length zero
    if(is.null(names(fpar))){
      stop("\"fpar\" must be a named list")
    }
    if(!isTRUE(all(names(fpar) %in% param_names))){
      stop("Unknown fixed parameter: must be one of \"loc\",\"scale\" or \"shape\". ")
    }
    if(!isTRUE(all(unlist(lapply(fpar, length)) == rep(1L, sum(wf))))){
      stop("Each fixed parameter must be of length one.")
    }
  }
  method <- match.arg(method)
  xdat <- as.double(xdat[is.finite(xdat)])
  n <- length(xdat)
  xmean <- mean(xdat)
  if(is.null(start)){
    xdat <- sort(xdat)
    xmax <- xdat[n] #sorted data
    xmin <- xdat[1]
    #Optimization routine, with PWM as default starting values
    pwm <- function(dat, r){
      # Data must be sorted!
      r <- as.integer(r)
      n <- length(dat)
      stopifnot(n > r, r > 0)
      sum(exp(lgamma((r+1):n) - lgamma(((r+1):n)-r))*dat[-(1:r)])/ exp(lgamma(n+1) - lgamma(n-r))
    }
    bpwm <- c(xmean, pwm(xdat,1), pwm(xdat,2))
    kst <- (2*bpwm[2]-bpwm[1])/(3*bpwm[3]-bpwm[1])-log(2)/log(3)
    xi_start <- -(7.859*kst + 2.9554*kst^2)
    sigma_start <- -(2*bpwm[2]-bpwm[1])*xi_start/(gamma(1-xi_start)*(1-2^(xi_start)))
    mu_start <- bpwm[1]-sigma_start*(gamma(1-xi_start)-1)/xi_start
    spar <- c(mu_start, sigma_start, xi_start)
    names(spar) <- param_names
  } else{ #start is provided by user
    # Check if list or vector + sanity checks for vectors/lists
    xmax <- max(xdat)
    xmin <- min(xdat)
    stopifnot(length(start) == (3L - sum(wf)))
    spar <- vector(mode = "numeric", length = 3L)
    names(spar) <- param_names
    if(is.null(names(start))){
      spar[!wf] <- unlist(start) # assume order, for better or worse
    } else {
      stopifnot(isTRUE(all(names(start) %in% param_names)))
      for(name in names(start)){
        spar[name] <- unlist(start[name])
      }
    }
  }
  for(i in seq_along(fpar)){
    spar[names(fpar[i])] <- unlist(fpar[i])[1]
  }
  stopifnot(spar[2] > 0, spar[3] > -1-1e-8)
  # check if initial value satisfy inequality constraints
  # multiple clauses because we cannot modify a fixed parameter
  if((spar[3] < 0)&((spar[2] + spar[3] * (xmax - spar[1])) < 0)){
    if(!wf[1]){
      spar[1] <- spar[2]/spar[3] + xmax + 0.1
    } else if(!wf[2]){
      spar[2] <- 1.1*(spar[1]-xmax)*spar[3]
    } else if(!wf[3]){
      spar[3] <- -0.9*spar[2]/(xmax-spar[1])
    }
  } else if((spar[3] > 0)&((spar[2] + spar[3] * (xmin - spar[1])) < 0)){
    if(!wf[1]){
      spar[1] <- spar[2]/spar[3] + xmin - 0.1
    } else if(!wf[2]){
      spar[2] <- 1.1*(spar[1]-xmin)*spar[3]
    } else if(!wf[3]){
      spar[3] <- 0.9*spar[2]/(spar[1]-xmin)
    }
  }
  start_vals <- spar[!wf]
  fixed_vals <- spar[wf] #when empty, a num vector of length zero
  wfo <- order(c(which(!wf), which(wf)))
  mle <- try(suppressWarnings(
    alabama::auglag(par = start_vals, fpar = fixed_vals, wfixed = wf, wfo = wfo,
                    fn = function(par, fpar, wfixed, wfo){
                      params <- c(par, fpar)[wfo]
                      nll <- -gev.ll(params, dat = xdat)
                      ifelse(is.finite(nll), nll, 1e10)
                    }, gr = function(par, fpar, wfixed, wfo) {
                      params <- c(par, fpar)[wfo]
                      grad <- -gev.score(params, dat = xdat)[!wfixed]
                      ifelse(is.finite(grad), grad, 1e6)
                    }, hin = function(par, fpar, wfixed, wfo) {
                      params <- c(par, fpar)[wfo]
                      c(params[2] + params[3] * (xmax - params[1]),
                        params[2] + params[3] * (xmin - params[1]),
                        params[2],
                        params[3] + 1)
                    }, control.outer = list(method = method, trace = FALSE),
                    control.optim = switch(method,
                                           nlminb = list(iter.max = 500L, rel.tol = 1e-10, step.min = 1e-10),
                                           list(maxit = 1000L, reltol = 1e-10)

                    ))))

  # Special case of MLE on the boundary xi = -1
  if(inherits(mle, what = "try-error")){
    stop("Optimization routine for the GEV did not converge.")
  }
  fitted$nllh <- mle$value
  fitted$estimate <- mle$par
  fitted$param <- c(mle$par, spar[wf])[wfo]
  if(!any(wf) | all(isTRUE(all.equal(wf, c(FALSE, FALSE, TRUE))), isTRUE(all.equal(fixed_vals, -1, check.attributes = FALSE)))){
    par_boundary <- c(xmean, xmax-xmean, -1)
    nll_boundary <- n*(1+log(par_boundary[2]))
    #Extract information and store
    if(!((nll_boundary > mle$value)&(mle$par[3] >= -1))){
      fitted$nllh <- nll_boundary
      fitted$estimate <- par_boundary[!wf]
      fitted$param <- par_boundary
      fitted$conv <- 0
    }
  }
  #Observed information matrix and standard errors

  fitted$vcov <- matrix(NA, ncol = length(mle$par), nrow = length(mle$par))
  fitted$std.err <- rep(NA, length(mle$par))
  if(fitted$param[3] > -0.5){
    vcovt <- try(solve(gev.infomat(par = fitted$param, dat = xdat, method = "obs")[!wf,!wf]))
    if(!inherits(vcovt, what = "try-error")){
      fitted$vcov <- vcovt
      fitted$std.err <- sqrt(diag(fitted$vcov))
      if(warnSE){
        warning("Cannot calculate standard error based on observed information")
      }
    }
  } else{
    if(warnSE){
      warning("Cannot calculate standard error based on observed information")
    }
  }
  names(fitted$param) <- names(wf) <- c("loc","scale","shape")
  names(fitted$std.err) <- names(fitted$estimate) <- c("loc","scale","shape")[!wf]
  fitted$method <- "auglag"
  fitted$nobs <- length(xdat)
  if(isTRUE(all(mle$kkt1, mle$kkt2))){
    fitted$convergence <- "successful"
  } else{
    fitted$convergence <- "converge dubious"
  }
  fitted$counts <- mle$counts
  fitted$xdat <- xdat
  fitted$start <- spar #if start not provided, this is PWM
  fitted$wfixed <- wf
  class(fitted) <- "mev_gev"
  if(show){
    print(fitted)
  }
  invisible(fitted)
}

#' Maximum likelihood estimates of point process for the r-largest observations
#'
#' This uses a constrained optimization routine to return the maximum likelihood estimate
#' based on an \code{n} by \code{r} matrix of observations. Observations should be ordered, i.e.,
#' the \code{r}-largest should be in the last column.
#'
#' @export
#' @inheritParams fit.gpd
#' @inheritParams fit.gev
#' @return a list containing the following components:
#' \itemize{
#' \item \code{estimate} a vector containing all the maximum likelihood estimates.
#' \item \code{std.err} a vector containing the standard errors.
#' \item \code{vcov} the variance covariance matrix, obtained as the numerical inverse of the observed information matrix.
#' \item \code{method} the method used to fit the parameter.
#' \item \code{nllh} the negative log-likelihood evaluated at the parameter \code{estimate}.
#' \item \code{convergence} components taken from the list returned by \code{\link[alabama]{auglag}}.
#' Values other than \code{0} indicate that the algorithm likely did not converge.
#' \item \code{counts} components taken from the list returned by \code{\link[alabama]{auglag}}.
#' \item \code{xdat} an \code{n} by \code{r} matrix of data
#' }
#' @examples
#' xdat <- rrlarg(n = 10, loc = 0, scale = 1, shape = 0.1, r = 4)
#' fit.rlarg(xdat)
fit.rlarg <- function(xdat,
                      start = NULL,
                      method = c("nlminb","BFGS"),
                      show = FALSE,
                      fpar = NULL,
                      warnSE = FALSE){
  param_names <- c("loc", "scale", "shape")
  stopifnot(is.null(fpar) | is.list(fpar))
  wf <- (param_names %in% names(fpar))
  if(sum(wf) == 3L){
    stop("Invalid input: all of the model parameters are fixed.")
  }
  if(is.list(fpar) && (length(fpar) >= 1L)){ #NULL has length zero
    if(is.null(names(fpar))){
      stop("\"fpar\" must be a named list")
    }
    if(!isTRUE(all(names(fpar) %in% param_names))){
      stop("Unknown fixed parameter: must be one of \"loc\",\"scale\" or \"shape\". ")
    }
    if(!isTRUE(all(unlist(lapply(fpar, length)) == rep(1L, sum(wf))))){
      stop("Each fixed parameter must be of length one.")
    }
  }
  xdat <- na.omit(as.matrix(xdat))
  method <- match.arg(method)
  r <- ncol(xdat)
  if(which.max(xdat[1,]) != 1){ #only check first row
    stop("Input should be ordered from largest to smallest in each row")
  }
#Optimization routine, with default starting values
xmax <- max(xdat)
xmin <- min(xdat)
spar <- vector(mode = "numeric", length = 3L)
if(is.null(start)){
    if(nrow(xdat) > 15L){ # Fit a generalized extreme value distribution to largest
    in2 <- sqrt(6 * var(xdat[,1]))/pi
    in1 <- mean(xdat[,1]) - 0.57722 * in2
    shape <- suppressWarnings(fit.gev(xdat[,1])$estimate[3])
    spar[1:3] <- c(in1, in2, shape + 0.2)
    if(spar[3] > 0 && (spar[2] + spar[3] * (xmin - spar[1]) <= 0)){
     spar[2] <- abs(spar[3]*(xmin - spar[1]))*1.1
    } else if(spar[3] < 0 && (spar[2] + spar[3] * (xmax - spar[1]) <= 0)){
      spar[2] <- abs(spar[3]*(xmax - spar[1]))*1.1
    }
  } else {
    spar <- suppressWarnings(fit.pp(as.vector(xdat), threshold = xmin, np = 1)$estimate)
  }
} else{ # start is provided
  stopifnot(length(start) == (3L - sum(wf)))
  if(is.null(names(start))){
    spar[!wf] <- unlist(start) # assume order, for better or worse
  } else {
    stopifnot(isTRUE(all(names(start) %in% param_names)))
    for(name in names(start)){
      spar[name] <- unlist(start[name])
    }
  }
}

names(spar) <- param_names
# end of start - attempt to find initial values.
for(i in seq_along(fpar)){
  spar[names(fpar[i])] <- unlist(fpar[i])[1]
}
if(isTRUE(any(
  (spar[3] > 0) && (spar[2] + spar[3] * (xmin - spar[1]) <= 0),
  (spar[3] < 0) && (spar[2] + spar[3] * (xmax - spar[1]) <= 0),
  spar[2] < 0,
  spar[3] < -1-1e-8))){
  stop("Invalid starting values in \"start\"")
}

# check if initial value satisfy inequality constraints
# multiple clauses because we cannot modify a fixed parameter
if((spar[3] < 0)&((spar[2] + spar[3] * (xmax - spar[1])) < 0)){
  if(!wf[1]){
    spar[1] <- spar[2]/spar[3] + xmax + 0.1
  } else if(!wf[2]){
    spar[2] <- 1.1*(spar[1]-xmax)*spar[3]
  } else if(!wf[3]){
    spar[3] <- -0.9*spar[2]/(xmax-spar[1])
  }
} else if((spar[3] > 0)&((spar[2] + spar[3] * (xmin - spar[1])) < 0)){
  if(!wf[1]){
    spar[1] <- spar[2]/spar[3] + xmin - 0.1
  } else if(!wf[2]){
    spar[2] <- 1.1*(spar[1]-xmin)*spar[3]
  } else if(!wf[3]){
    spar[3] <- 0.9*spar[2]/(spar[1]-xmin)
  }
}

start_vals <- spar[!wf]
fixed_vals <- spar[wf] #when empty, a num vector of length zero
wfo <- order(c(which(!wf), which(wf)))
mle <- try(suppressWarnings(
  alabama::auglag(par = start_vals,
                  fpar = fixed_vals,
                  wfixed = wf,
                  wfo = wfo,
                  fn = function(par, fpar, wfixed, wfo){
                    params <- c(par, fpar)[wfo]
                    nll <-  -rlarg.ll(params, dat = xdat)
                    ifelse(is.finite(nll), nll, 1e10)
                  }, gr = function(par, fpar, wfixed, wfo) {
                    params <- c(par, fpar)[wfo]
                    grad <- -rlarg.score(params, dat = xdat)[!wfixed]
                    ifelse(is.finite(grad), grad, 1e10)
                  }, hin = function(par, fpar, wfixed, wfo) {
                    params <- c(par, fpar)[wfo]
                    c(params[2] + params[3] * (xmax - params[1]),
                      params[2] + params[3] * (xmin - params[1]),
                      params[2],
                      params[3] + 1-1e-8)
                  }, control.outer = list(method = method, trace = FALSE, NMinit = TRUE),
                  control.optim = switch(method,
                                         nlminb = list(iter.max = 1000L, rel.tol = 1e-10, step.min = 1e-10),
                                         list(maxit = 1000L, reltol = 1e-10)

                  ))))

if(inherits(mle, what = "try-error")){
  stop("Optimization routine for r-largest observations did not converge")
}
#Extract information and store
fitted <- list()
fitted$convergence <- mle$convergence
if(isTRUE(all(mle$kkt1, mle$kkt2, mle$convergence == 0))){
  fitted$convergence <- "successful"
} else{
  fitted$convergence <- "dubious convergence"
}
fitted$nllh <- mle$value
fitted$estimate <- mle$par
fitted$param <- c(mle$par, spar[wf])[wfo]
#Point estimate
if(sum(wf) == 0){
  xmeanr <- mean(xdat[,r])
  # check MLE at xi=-1
  par_boundary <- c(xmax - (xmax-xmeanr)/r, (xmax-xmeanr)/r, -1)
  nll_boundary <- -rlarg.ll(par_boundary, dat = xdat)
 if(!isTRUE(all(nll_boundary > mle$value, mle$par[3] > -1))){
   fitted$nllh <- nll_boundary
   fitted$estimate <- par_boundary
   fitted$convergence <- "successful"
 }
}
#Observed information matrix and standard errors
fitted$vcov <- matrix(NA, ncol = length(mle$par), nrow = length(mle$par))
fitted$std.err <- rep(NA, length(mle$par))
if(fitted$param[3] > -0.5){
  vcovt <- try(solve(rlarg.infomat(par = fitted$param, dat = xdat, method = "obs")[!wf,!wf]))
  if(!inherits(vcovt, what = "try-error")){

    fitted$vcov <- vcovt
    fitted$std.err <- sqrt(diag(fitted$vcov))
  }
} else{
  if(warnSE){
    warning("Cannot calculate standard error based on observed information")
  }
}
names(fitted$param) <- names(wf) <- c("loc","scale","shape")
names(fitted$std.err) <- names(fitted$estimate) <- c("loc","scale","shape")[!wf]
fitted$method <- "auglag"

fitted$counts <- mle$counts
fitted$xdat <- xdat
fitted$nobs <- length(xdat)
fitted$start <- spar #if start not provided, this is PWM
fitted$wfixed <- wf
class(fitted) <- c("mev_rlarg", "mev_gev")
if(show){
  print(fitted)
}
invisible(fitted)
}


# @param x A fitted object of class \code{gpd}.
# @param main title for the Q-Q plot #' @param xlab x-axis label
# @param ylab y-axis label
# @param ... additional argument passed to \code{matplot}.
#' @importFrom evd qgpd
#' @export
plot.mev_gpd <- function(x, which = 1:2, main, xlab = "Theoretical quantiles", ylab = "Sample quantiles", add = TRUE, ...) {
  if (!is.vector(x$exceedances)) {
    stop("Object \"x\" does not contain \"exceedances\", or else the latter is not a vector")
  }
  if (!is.numeric(which) || any(which < 1) || any(which > 2)){
    stop("\"which\" must be in 1:2")
  }
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  if(!add){
    old.par <- par(no.readonly = TRUE)
    if(sum(show) == 2){
      par(mfrow = c(1,2), mar = c(5,5,4,1))
    }
    on.exit(par(old.par))
  }
  if (missing(main)) {
    main <- c("Probability-probability plot", "Quantile-quantile plot")
  } else{
    if(length(main) != sum(show)){
     stop("Invalid input: \"main\" must be of the same length as \"which\"")
    }
  }

  dat <- sort(x$exceedances)
  n <- length(dat)
  pp_confint_lim <- t(sapply(1:n, function(i) {qbeta(c(0.025, 0.975), i, n - i + 1) }))
  qq_confint_lim <- apply(pp_confint_lim, 2, function(y){ evd::qgpd(y, loc = 0, scale = x[["param"]][1], shape = x[["param"]][2])})
  pobs <- (1:n)/(n + 1)
  quant <- evd::qgpd(pobs, loc = 0, scale = x[["param"]][1], shape = x[["param"]][2])
  if(show[1]){
    matplot(pobs, cbind(pp_confint_lim,
                        evd::pgpd(dat, loc = 0, scale = x[["param"]][1], shape = x[["param"]][2])
                        ), main = main[1], xlab = xlab, ylab = ylab, type = "llp",
            pch = 20, col = c("grey", "grey", 1), ylim = c(0,1), xlim = c(0,1),
            lty = c(2, 2, 1), bty = "l", pty = "s", first={abline(0, 1, col="grey")}, ..., add = FALSE)
  }

  if(show[2]){
    limqq <- c(0, max(c(quant[n], dat[n])))
    matplot(quant, cbind(qq_confint_lim, dat), main = main[2], xlab = xlab, ylab = ylab, type = "llp",
            pch = 20, col = c("grey", "grey", 1), ylim = limqq, xlim = limqq,
          lty = c(2, 2, 1), bty = "l", pty = "s", first={abline(0, 1, col="grey")}, ..., add = FALSE)
  }

  matlim <- cbind(quant, qq_confint_lim)
  colnames(matlim) <- c("quantile", "lower","upper")
  invisible(matlim)
}


# @param x A fitted object of class \code{mev_gpdbayes}.
# @param main title for the Q-Q plot #' @param xlab x-axis label
# @param ylab y-axis label
# @param ... additional argument passed to \code{matplot}.
#' @importFrom evd qgpd
#' @export
plot.mev_gpdbayes <- function(x, which = 1:2, main, xlab = "Theoretical quantiles", ylab = "Sample quantiles", add = TRUE, ...) {
  if (!is.vector(x$exceedances)) {
    stop("Object \"x\" does not contain \"exceedances\", or else the latter is not a vector")
  }
  if (!is.numeric(which) || any(which < 1) || any(which > 2)){
    stop("\"which\" must be in 1:2")
  }
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  if(!add){
    old.par <- par(no.readonly = TRUE)
    if(sum(show) == 2){
      par(mfrow = c(1,2), mar = c(5,5,4,1))
    }
    on.exit(par(old.par))
  }
  if (missing(main)) {
    main <- c("Probability-probability plot", "Quantile-quantile plot")
  } else{
    if(length(main) != sum(show)){
      stop("Invalid input: \"main\" must be of the same length as \"which\"")
    }
  }

  dat <- sort(x$exceedances)
  n <- length(dat)
  pp_confint_lim <- t(sapply(1:n, function(i) {qbeta(c(0.025, 0.975), i, n - i + 1) }))
  qq_confint_lim <- apply(pp_confint_lim, 2, function(y){ evd::qgpd(y, loc = 0, scale = x[["estimate"]][1], shape = x[["estimate"]][2])})
  pobs <- (1:n)/(n + 1)
  quant <- evd::qgpd(pobs, loc = 0, scale = x[["estimate"]][1], shape = x[["estimate"]][2])
  if(show[1]){
    matplot(pobs, cbind(pp_confint_lim,
                        evd::pgpd(dat, loc = 0, scale = x[["estimate"]][1], shape = x[["estimate"]][2])
    ), main = main[1], xlab = xlab, ylab = ylab, type = "llp",
    pch = 20, col = c("grey", "grey", 1), ylim = c(0,1), xlim = c(0,1),
    lty = c(2, 2, 1), bty = "l", pty = "s", first={abline(0, 1, col="grey")}, ..., add = FALSE)
  }

  if(show[2]){
    limqq <- c(0, max(c(quant[n], dat[n])))
    matplot(quant, cbind(qq_confint_lim, dat), main = main[2], xlab = xlab, ylab = ylab, type = "llp",
            pch = 20, col = c("grey", "grey", 1), ylim = limqq, xlim = limqq,
            lty = c(2, 2, 1), bty = "l", pty = "s", first={abline(0, 1, col="grey")}, ..., add = FALSE)
  }

  matlim <- cbind(quant, qq_confint_lim)
  colnames(matlim) <- c("quantile", "lower","upper")
  invisible(matlim)
}



# @param x A fitted object of class \code{mev_gev}.
# @param main title for the Q-Q plot
# @param xlab x-axis label
# @param ylab y-axis label
# @param ... additional argument passed to \code{matplot}.
#' @export
plot.mev_gev <- function(x, which = 1:2, main, xlab = "Theoretical quantiles", ylab = "Sample quantiles", ...) {
  if (!is.vector(x$xdat)) {
    stop("Object \"x\" does not contain \"exceedances\", or else the latter is not a vector")
  }
  if (!is.numeric(which) || any(which < 1) || any(which > 2)){
    stop("\"which\" must be in 1:2")
  }
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  old.par <- par(no.readonly = TRUE)
  if(sum(show) == 2){
    par(mfrow = c(1,2), mar = c(5,5,4,1))
  }
  on.exit(par(old.par))
  if (missing(main)) {
    main <- c("Probability-probability plot", "Quantile-quantile plot")
  } else{
    if(length(main) != sum(show)){
      stop("Invalid input: \"main\" must be of the same length as \"which\"")
    }
  }
  dat <- sort(x$xdat)
  pars <- x$param
  n <- length(dat)
  pp_confint_lim <- t(sapply(1:n, function(i) {qbeta(c(0.025, 0.975), i, n - i + 1) }))
  qq_confint_lim <- apply(pp_confint_lim, 2, function(y){
    evd::qgev(y, loc = pars[1], scale = pars[2], shape = pars[3])})
  pobs <- (1:n)/(n + 1)
  quant <- evd::qgev(pobs, loc = pars[1], scale = pars[2], shape = pars[3])
  if(show[1]){
    matplot(pobs, cbind(pp_confint_lim,
                        evd::pgev(dat, loc = pars[1], scale = pars[2], shape = pars[3])
                        ), main = main[1], xlab = xlab, ylab = ylab, type = "llp",
            pch = 20, col = c("grey", "grey", 1), ylim = c(0,1), xlim = c(0,1),
            lty = c(2, 2, 1), bty = "l", pty = "s", first={abline(0, 1, col="grey")}, ...)
  }

  if(show[2]){
    limqq <- c(min(c(quant[1], dat[1])), max(c(quant[n], dat[n])))
    matplot(quant, cbind(qq_confint_lim, dat), main = main[2], xlab = xlab, ylab = ylab, type = "llp",
            pch = 20, col = c("grey", "grey", 1), ylim = limqq, xlim = limqq,
            lty = c(2, 2, 1), bty = "l", pty = "s", first={abline(0, 1, col="grey")}, ...)
  }

  matlim <- cbind(quant, qq_confint_lim)
  colnames(matlim) <- c("quantile", "lower","upper")
  invisible(matlim)
}

# @param x A fitted object of class \code{mev_rlarg}.
# @param main title for the Q-Q plot
# @param xlab x-axis label
# @param ylab y-axis label
# @param ... additional argument passed to \code{matplot}.
#' @export
plot.mev_rlarg <- function(x, which = 1:2, main, xlab = "Theoretical quantiles", ylab = "Sample quantiles", ...) {
  if(!isTRUE(all.equal(x$param[3], 0, check.attributes = FALSE))){
  ppdat <- (1+x$param[3]*(as.matrix(x$xdat)-x$param[1])/x$param[2])^(-1/x$param[3])
  } else{
    ppdat <- exp(-(as.matrix(x$xdat)-x$param[1])/x$param[2])
  }
  if(ncol(ppdat) == 1){
    dat <- as.vector(ppdat)
  } else{
    dat <- sort(as.vector(c(ppdat[,1], apply(ppdat, 1, diff))))
  }
  n <- length(dat)
  if (!is.numeric(which) || any(which < 1) || any(which > 2)){
    stop("\"which\" must be in 1:2")
  }
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  old.par <- par(no.readonly = TRUE)
  if(sum(show) == 2){
    par(mfrow = c(1,2), mar = c(5,5,4,1))
  }
  on.exit(par(old.par))
  if (missing(main)) {
    main <- c("Probability-probability plot", "Quantile-quantile plot")
  } else{
    if(length(main) != sum(show)){
      stop("Invalid input: \"main\" must be of the same length as \"which\"")
    }
  }
  pp_confint_lim <- t(sapply(1:n, function(i) {qbeta(c(0.025, 0.975), i, n - i + 1) }))
  qq_confint_lim <- apply(pp_confint_lim, 2, function(y){qexp(y)})
  pobs <- (1:n)/(n + 1)
  quant <- qexp(pobs)
  if(show[1]){
    matplot(pobs, cbind(pp_confint_lim, pexp(dat)), main = main[1], xlab = xlab, ylab = ylab, type = "llp",
            pch = 20, col = c("grey", "grey", 1), ylim = c(0,1), xlim = c(0,1),
            lty = c(2, 2, 1), bty = "l", pty = "s", first={abline(0, 1, col="grey")}, ...)
  }

  if(show[2]){
    limqq <- c(0, max(c(quant[n], dat[n])))
    matplot(quant, cbind(qq_confint_lim, dat), main = main[2], xlab = xlab, ylab = ylab, type = "llp",
            pch = 20, col = c("grey", "grey", 1), ylim = limqq, xlim = limqq,
            lty = c(2, 2, 1), bty = "l", pty = "s", first={abline(0, 1, col="grey")}, ...)
  }

  matlim <- cbind(quant, qq_confint_lim)
  colnames(matlim) <- c("quantile", "lower","upper")
  invisible(matlim)
}


# @param x A fitted object of class \code{mev_pp}.
# @param main title for the Q-Q plot
# @param xlab x-axis label
# @param ylab y-axis label
# @param ... additional argument passed to \code{matplot}.
#' @export
plot.mev_pp <- function(x, which = 1:2, main = "Quantile-quantile plot", xlab = "Theoretical quantiles", ylab = "Sample quantiles", ...) {
  if (!is.vector(x$exceedances)) {
    stop("Object \"x\" does not contain \"exceedances\", or else the latter is not a vector")
  }
  if (!is.numeric(which) || any(which < 1) || any(which > 2)){
    stop("\"which\" must be in 1:2")
  }
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  old.par <- par(no.readonly = TRUE)
  if(sum(show) == 2){
    par(mfrow = c(1,2), mar = c(5,5,4,1))
  }
  on.exit(par(old.par))
  if (missing(main)) {
    main <- c("Probability-probability plot", "Quantile-quantile plot")
  } else{
    if(length(main) != sum(show)){
      stop("Invalid input: \"main\" must be of the same length as \"which\"")
    }
  }

  dat <- sort(x$exceedances)
  scalet <- x$param['scale'] + x$param['shape']*(x$threshold - x$param['loc'])
  n <- length(dat)
  pp_confint_lim <- t(sapply(1:n, function(i) {qbeta(c(0.025, 0.975), i, n - i + 1) }))
  qq_confint_lim <- apply(pp_confint_lim, 2, function(y){
    evd::qgpd(y, loc = x$threshold, scale = scalet, shape = x$param['shape'])})
  pobs <- (1:n)/(n + 1)
  quant <- evd::qgpd(pobs, loc = x$threshold, scale = scalet, shape = x$param['shape'])
  if(show[1]){
    matplot(pobs, cbind(pp_confint_lim,
                        evd::pgpd(dat, loc = x$threshold, scale = scalet, shape = x$param['shape'])),
            main = main[1], xlab = xlab, ylab = ylab, type = "llp",
            pch = 20, col = c("grey", "grey", 1), ylim = c(0,1), xlim = c(0,1),
            lty = c(2, 2, 1), bty = "l", pty = "s", first={abline(0, 1, col="grey")}, ...)
  }

  if(show[2]){
    limqq <- c(x$threshold, max(c(quant[n], dat[n])))
    matplot(quant, cbind(qq_confint_lim, dat), main = main[2], xlab = xlab, ylab = ylab, type = "llp",
            pch = 20, col = c("grey", "grey", 1), ylim = limqq, xlim = limqq,
            lty = c(2, 2, 1), bty = "l", pty = "s", first={abline(0, 1, col="grey")}, ...)
  }

  matlim <- cbind(quant, qq_confint_lim)
  colnames(matlim) <- c("quantile", "lower","upper")
  invisible(matlim)
}


# @param x A fitted object of class \code{mev_gpd}.
# @param digits Number of digits to display in \code{print} call.
# @param ... Additional argument passed to \code{print}.
#' @export
print.mev_gpd <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Method:", x$method, "\n")
  cat("Log-likelihood:", round(-x$nllh, digits), "\n")

  cat("\nThreshold:", round(x$threshold, digits), "\n")
  cat("Number Above:", x$nat, "\n")
  cat("Proportion Above:", round(x$pat, digits), "\n")

  cat("\nEstimates\n")
  print.default(format(x$estimate, digits = digits), print.gap = 2, quote = FALSE, ...)
  if (!is.null(x$std.err) && x$estimate[1] > -0.5) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, quote = FALSE, ...)
  }
  if(length(x$estimate) != length(x$param)){
    cat("\nParameters\n")
    print.default(format(x$param, digits = digits), print.gap = 2, quote = FALSE, ...)
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  if (x$method != "Grimshaw") {
    cat("  Function Evaluations:", x$counts["function"], "\n")
    if ((!is.na(x$counts["gradient"])) && x$method != "nlm")
      cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
    cat("\n")
  }
  invisible(x)
}

# @param x A fitted object of class \code{mev_gpdbayes}.
# @param digits Number of digits to display in \code{print} call.
# @param ... Additional argument passed to \code{print}.
#' @export
print.mev_gpdbayes <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nMethod:", switch(x$method, zs = "Zhang and Stephens", zhang = "Zhang"), "\n")
  cat("\nThreshold:", round(x$threshold, digits), "\n")
  cat("Number Above:", x$nat, "\n")
  cat("Proportion Above:", round(x$pat, digits), "\n")

  cat("\nApproximate posterior mean estimates\n")
  print.default(format(x$approx.mean, digits = 3), print.gap = 2, quote = FALSE)
  if (!is.null(x$post.mean)) {
    cat("\nPosterior mean estimates\n")
    print.default(format(x$post.mean, digits = 3), print.gap = 2, quote = FALSE)
    cat("\nMonte Carlo standard errors\n")
    print.default(format(x$post.se, digits = 3), print.gap = 2, quote = FALSE)
    cat("\nEstimates based on an adaptive MCMC\n Runs:   ", x$niter, "\n Burnin: ", x$burnin, "\n Acceptance rate:", round(x$accept.rate,
                                                                                                                           digits = 2), "\n Thinning:", x$thin, "\n")

  }
}

# @param x A fitted object of class \code{mev_gev}.
# @param digits Number of digits to display in \code{print} call.
# @param ... Additional argument passed to \code{print}.
#' @export
print.mev_gev <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Log-likelihood:", -x$nllh, "\n")

  cat("\nEstimates\n")
  print.default(format(x$estimate, digits = digits), print.gap = 2, quote = FALSE, ...)
  if (!isTRUE(any(is.na(x$std.err))) && x$estimate[1] > -0.5) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, quote = FALSE, ...)
  }
  if(length(x$estimate) != length(x$param)){
    cat("\nParameters\n")
    print.default(format(x$param, digits = digits), print.gap = 2, quote = FALSE, ...)
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  cat("  Function Evaluations:", x$counts["function"], "\n")
  cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
  cat("\n")
  invisible(x)
}


# @param x A fitted object of class \code{mev_gev}.
# @param digits Number of digits to display in \code{print} call.
# @param ... Additional argument passed to \code{print}.
#' @export
print.mev_pp <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Log-likelihood:", -x$nllh, "\n")

  cat("\nThreshold:", round(x$threshold, digits), "\n")
  cat("Number Above:", x$nat, "\n")
  cat("Proportion Above:", round(x$pat, digits), "\n")
  cat("Number of periods:", round(x$np, digits), "\n")

  cat("\nEstimates\n")
  print.default(format(x$estimate, digits = digits), print.gap = 2, quote = FALSE, ...)
  if (!isTRUE(any(is.na(x$std.err))) && x$estimate[1] > -0.5) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, quote = FALSE, ...)
  }
  if(length(x$estimate) != length(x$param)){
    cat("\nParameters\n")
    print.default(format(x$param, digits = digits), print.gap = 2, quote = FALSE, ...)
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  cat("  Function Evaluations:", x$counts["function"], "\n")
  cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
  cat("\n")
  invisible(x)
}


  # Methods for class "mev_gev", returned by mev::fit.gev()

  #' @export
  nobs.mev_gev <- function(object, ...) {
    return(object$nobs)
  }

  #' @export
  coef.mev_gev <- function(object, ...) {
    return(object$estimate)
  }

  #' @export
  vcov.mev_gev <- function(object, ...) {
    return(object$vcov)
  }

  #' @export
  logLik.mev_gev <- function(object, ...) {
    val <- -object$nllh
    attr(val, "nobs") <- nobs(object)
    attr(val, "df") <- length(coef(object))
    class(val) <- "logLik"
    return(val)
  }



  # Methods for class "mev_gpd", returned by mev::fit.gpd()

  #' @export
  nobs.mev_gpd <- function(object, ...) {
    return(object$nat)
  }

  #' @export
  coef.mev_gpd <- function(object, ...) {
    return(object$estimate)
  }

  #' @export
  vcov.mev_gpd <- function(object, ...) {
    return(object$vcov)
  }

  #' @export
  logLik.mev_gpd <- function(object, ...) {
    val <- -object$nllh
    attr(val, "nobs") <- nobs(object)
    attr(val, "df") <- length(coef(object))
    class(val) <- "logLik"
    return(val)
  }

  #' @export
  nobs.mev_rlarg <- function(object, ...) {
    #total number of observations is prod(xdat), so length(xdat)
    #n by r for a matrix
    return(object$nobs)
  }

  #' @export
  coef.mev_rlarg <- function(object, ...) {
    return(object$estimate)
  }

  #' @export
  vcov.mev_rlarg <- function(object, ...) {
    return(object$vcov)
  }

  #' @export
  logLik.mev_rlarg <- function(object, ...) {
    val <- -object$nllh
    attr(val, "nobs") <- nobs(object)
    attr(val, "df") <- length(coef(object))
    class(val) <- "logLik"
    return(val)
  }


  #' @export
  nobs.mev_pp <- function(object, ...) {
    return(length(object$xdat))
  }

  #' @export
  coef.mev_pp <- function(object, ...) {
    return(object$estimate)
  }

  #' @export
  vcov.mev_pp <- function(object, ...) {
    return(object$vcov)
  }

  #' @export
  logLik.mev_pp <- function(object, ...) {
    val <- -object$nllh
    attr(val, "nobs") <- nobs(object)
    attr(val, "df") <- length(coef(object))
    class(val) <- "logLik"
    return(val)
  }

  #' @export
  nobs.mev_egp <- function(object, ...) {
    return(object$nat)
  }

  #' @export
  coef.mev_egp <- function(object, ...) {
    return(object$estimate)
  }

  #' @export
  vcov.mev_egp <- function(object, ...) {
    return(object$vcov)
  }

  #' @export
  logLik.mev_egp <- function(object, ...) {
    val <- -object$nllh
    attr(val, "nobs") <- nobs(object)
    attr(val, "df") <- length(coef(object))
    class(val) <- "logLik"
    return(val)
  }
  #' @export
  anova.mev_gpd <- function(object, object2, ...){
    if (any(missing(object), missing(object2))){
      stop("At least two models must be specified.")
    }

    model1 <- deparse(substitute(object))
    model2 <- deparse(substitute(object2))
    models <- c(model1, model2)
    narg <- length(models)
    for (i in 1:narg) {
      if (!inherits(get(models[i], envir = parent.frame()), "mev_gpd")){
        stop("Invalid input: use only with objects of class 'mev_gpd'.")
      }
    }

    npar <- rep(0, length(models))
    dev <- rep(0, length(models))
    for (i in 1:narg) {
      evmod <- get(models[i], envir = parent.frame())
      if(evmod$method %in% c("obre","zs","zhang")){
        stop("Model comparison not supported for the chosen method.")
      }
      dev[i] <- 2*evmod$nllh
      npar[i] <- length(evmod$estimate)
    }
    if(sum(object2$wfixed) <= sum(object$wfixed)){
      stop("Invalid order: \"object2\" must contain a model with restrictions")
    }
    if(!isTRUE(all.equal(object$xdat, object2$xdat))){
      stop("Invalid arguments: the samples should be the same")
    }
    nested <- FALSE
    if(isTRUE(all(which(object$wfixed) %in% which(object2$wfixed)))){
      if(isTRUE(all.equal(object$param[which(object$wfixed)],
                          object2$param[which(object$wfixed)]))){
        nested <- TRUE
      }
    }
    if(!nested){
      stop("Invalid input: models are not nested.")
    }
    df <- -diff(npar)
    dvdiff <- diff(dev)
    pval <- pchisq(dvdiff, df = df, lower.tail = FALSE)
    table <- data.frame(npar, dev, c(NA, df), c(NA, dvdiff), c(NA,
                                                               pval))
    dimnames(table) <- list(models, c("npar", "Deviance", "Df",
                                      "Chisq", "Pr(>Chisq)"))
    structure(table, heading = c("Analysis of Deviance Table\n"),
              class = c("anova", "data.frame"))
  }
  #' @export
  anova.mev_gev <- function(object, object2, ...){
    if (any(missing(object), missing(object2))){
      stop("At least two models must be specified.")
    }

    model1 <- deparse(substitute(object))
    model2 <- deparse(substitute(object2))
    models <- c(model1, model2)
    narg <- length(models)
    for (i in 1:narg) {
      if (!inherits(get(models[i], envir = parent.frame()), "mev_gev"))
        stop("Invalid input: use only with objects of class 'mev_gev'.")
    }

    npar <- rep(0, length(models))
    dev <- rep(0, length(models))
    for (i in 1:narg) {
      evmod <- get(models[i], envir = parent.frame())
      dev[i] <- 2*evmod$nllh
      npar[i] <- length(evmod$estimate)
    }
    if(sum(object2$wfixed) <= sum(object$wfixed)){
      stop("Invalid order: \"object2\" must contain a model with restrictions")
    }
    if(!isTRUE(all.equal(object$xdat, object2$xdat))){
      stop("Invalid arguments: the samples should be the same")
    }
    nested <- FALSE
    if(isTRUE(all(which(object$wfixed) %in% which(object2$wfixed)))){
      if(isTRUE(all.equal(object$param[which(object$wfixed)],
                          object2$param[which(object$wfixed)]))){
        nested <- TRUE
      }
    }
    if(!nested){
      stop("Invalid input: models are not nested.")
    }
    df <- -diff(npar)
    dvdiff <- diff(dev)
    pval <- pchisq(dvdiff, df = df, lower.tail = FALSE)
    table <- data.frame(npar, dev, c(NA, df), c(NA, dvdiff), c(NA,
                                                               pval))
    dimnames(table) <- list(models, c("npar", "Deviance", "Df",
                                      "Chisq", "Pr(>Chisq)"))
    structure(table, heading = c("Analysis of Deviance Table\n"),
              class = c("anova", "data.frame"))
  }
  #' @export
  anova.mev_pp <- function(object, object2, ...){
    if (any(missing(object), missing(object2))){
      stop("At least two models must be specified.")
    }

    model1 <- deparse(substitute(object))
    model2 <- deparse(substitute(object2))
    models <- c(model1, model2)
    narg <- length(models)
    for (i in 1:narg) {
      if (!inherits(get(models[i], envir = parent.frame()), "mev_pp"))
        stop("Invalid input: use only with objects of class 'mev_pp'.")
    }

    npar <- rep(0, length(models))
    dev <- rep(0, length(models))
    for (i in 1:narg) {
      evmod <- get(models[i], envir = parent.frame())
      dev[i] <- 2*evmod$nllh
      npar[i] <- length(evmod$estimate)
    }
    if(sum(object2$wfixed) <= sum(object$wfixed)){
      stop("Invalid order: \"object2\" must contain a model with restrictions")
    }
    if(!isTRUE(all.equal(object$xdat, object2$xdat))){
      stop("Invalid arguments: the samples should be the same")
    }
    nested <- FALSE
    if(isTRUE(all(which(object$wfixed) %in% which(object2$wfixed)))){
      if(isTRUE(all.equal(object$param[which(object$wfixed)],
                          object2$param[which(object$wfixed)]))){
        nested <- TRUE
      }
    }
    if(!nested){
      stop("Invalid input: models are not nested.")
    }
    df <- -diff(npar)
    dvdiff <- diff(dev)
    pval <- pchisq(dvdiff, df = df, lower.tail = FALSE)
    table <- data.frame(npar, dev, c(NA, df), c(NA, dvdiff), c(NA,
                                                               pval))
    dimnames(table) <- list(models, c("npar", "Deviance", "Df",
                                      "Chisq", "Pr(>Chisq)"))
    structure(table, heading = c("Analysis of Deviance Table\n"),
              class = c("anova", "data.frame"))
  }
