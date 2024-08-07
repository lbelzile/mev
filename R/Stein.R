#' Stein's weighted generalized Pareto likelihood
#'
#' @param pars vector of length with scale and shape parameters
#' @param xdat vector of exceedances in decreasing order if \code{sort=FALSE}
#' @param weights vector of weights
#' @param ... additional arguments, currently ignored
#' @keywords internal
#' @export
#' @references Stein, M.L. (2023). A weighted composite log-likelihood approach to parametric estimation of the extreme quantiles of a distribution. \emph{Extremes} \bold{26}, 469-507 <doi:10.1007/s10687-023-00466-w>
stein_gp_lik <- function(pars, xdat, weights = rep(1, length(xdat)), ...){
  xdat <- as.numeric(xdat)
  pars <- as.numeric(pars) # strip names
  stopifnot(length(xdat) == length(weights))
  # Extract parameters
  scale <- pars[1]
  shape <- pars[2]
  # Get order statistics
  n <- length(xdat)
  args <- list(...)
  sorted <- isTRUE(args$sorted)
  if(!sorted){
   xdat <- sort(xdat, decreasing = TRUE)
  }
  # ydat <- rev(c(0, diff(xdat)))
  if(xdat[n] <= 0){
   stop("Invalid vector of threshold exceedances \"xdat\".")
  }
  if((scale < 0) | ((shape < 0) & (1+shape*xdat[1]/scale < 0)) | shape < -1){
    return(-1e10)
  }
  lp <- log1p(shape*xdat/scale)
  if(!isTRUE(all.equal(shape, 0, tol = 1e-5))){
  -sum(weights)*log(scale) -
    sum(weights * (seq_len(n)/shape + 1) * lp) +
    sum(weights[-n] * seq_len(n-1) * lp[-1])/shape
  }
}

#' Stein's vector of weights
#'
#' Computes a vector of decreasing weights for exceedances.
#'
#' @param n integer, sample size
#' @param gamma positive scalar for power
#' @return a vector of positive weights of length \code{n}
#' @export
#' @keywords internal
#' @references Stein, M.L. (2023). A weighted composite log-likelihood approach to parametric estimation of the extreme quantiles of a distribution. \emph{Extremes} \bold{26}, 469-507 <doi:10.1007/s10687-023-00466-w>
Stein_weights <- function(n, gamma = 1){
  stopifnot(length(gamma) == 1, gamma > 0)
  n <- as.integer(n)
  stopifnot(length(n) == 1L, n > 0)
  (gamma + 1)/gamma * (1-((0:(n-1))/n)^gamma)
}

#' Maximum likelihood estimation for weighted generalized Pareto distribution
#'
#' Weighted maximum likelihood estimation, with user-specified vector of weights.
#'
#' @export
#' @param xdat vector of observations
#' @param threshold numeric, value of the threshold
#' @param weightfun function whose first argument is the length of the weight vector
#' @return a list with components
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
#' \item \code{convergence} logical indicator of convergence.
#' \item \code{weights} vector of weights for exceedances.
#' \item \code{exceedances} excess over the threshold, sorted in decreasing order.
#' }
fit.wgpd <- function(xdat, threshold = 0, weightfun = Stein_weights, ...){
  xdat <- as.numeric(xdat[is.finite(xdat)])
  ntot <- length(xdat)
  # Extract exceedances
  xdat <- xdat[xdat > threshold] - threshold
  stopifnot(class(weightfun) == "function")
  weights <- try(weightfun(length(xdat), ...), silent = TRUE)
  if(inherits(weights, "try-error")){
    weights <- rep(1, length(xdat))
  }
  stopifnot(length(xdat) == length(weights))
  # Scale data for optimization
  sc <- sd(xdat)
  # Sort data to avoid doing this repeatedly in the optimization
  xdat <- sort(xdat, decreasing = TRUE)
  exc <- xdat
  xdat <- exc/sc
 # Define inequality and support constraints for the parameters
  hin <- function(par, maxdat = NULL, thresh = 0, ...) {
    stopifnot(`Argument "maxdat" is missing, with no default value.` = !is.null(maxdat))
    c(par[1], par[2], ifelse(par[2] < 0, thresh - par[1]/par[2] -
                               maxdat, 1e-05))
  }
  ineqLB <- c(0, -1, 0)
  ineqUB <- c(Inf, 10, Inf)
  LB <- c(0, -1)
  UB <- c(Inf, 2)
  start <- c(1, 0.1)
  opt <- Rsolnp::solnp(pars = start,
                fun = function(pars, ...){
                  -stein_gp_lik(pars = pars, ...)},
                ineqfun = hin,
                ineqLB = ineqLB,
                ineqUB = ineqUB,
                LB = LB,
                UB = UB,
                xdat = xdat,
                maxdat = xdat[1],
                weights = weights,
                sorted = TRUE,
                control = list(trace = 0))
  # Scale back parameters
  mle <- opt$pars * c(sc, 1)
  se <- "fail"
  cov <- try(solve(opt$hessian[-(1:length(ineqLB)),-(1:length(ineqLB))]),
             silent = TRUE)
  if(!inherits(cov, "try-error")){
    se <- try(sqrt(diag(cov)),
              silent = TRUE)
  }
  if(!is.character(se)){
    if(isTRUE(all(is.finite(se))) & isTRUE(all(se > 0))){
     se <- se * c(sc, 1)
    }
  }
  if(is.character(se)){
    se <- rep(NA, 2)
    cov <- NULL
  }
  return(list(
    estimate = mle,
    param = mle,
    std.err = se,
    vcov = cov,
    threshold = threshold,
    nllh = -(-opt$values[length(opt$values)] - sum(weights)*log(sc)),
    #stein_gp_lik(pars = mle, xdat = exc, weights = weights),
    convergence = isTRUE(opt$convergence == 0),
    exceedances = xdat * sc,
    nat = length(xdat),
    pat = length(xdat)/ntot,
    weights = weights
  ))
}

