.auglag <- function (x0, fn, gr = NULL, lower = NULL, upper = NULL, hin = NULL,
            hinjac = NULL, heq = NULL, heqjac = NULL, localsolver = c("COBYLA"),
            localtol = 1e-06, ineq2local = FALSE, nl.info = FALSE, control = list(),
            ...)
  {
    if (ineq2local) {
      stop("Inequalities to local solver: feature not yet implemented.")
    }
    localsolver <- toupper(localsolver)
    if (localsolver %in% c("COBYLA", "BOBYQA")) {
      dfree <- TRUE
      gsolver <- "NLOPT_LN_AUGLAG"
      lsolver <- paste("NLOPT_LN_", localsolver, sep = "")
    }
    else if (localsolver %in% c("LBFGS", "MMA", "SLSQP")) {
      dfree <- FALSE
      gsolver <- "NLOPT_LD_AUGLAG"
      lsolver <- paste("NLOPT_LD_", localsolver, sep = "")
    }
    else {
      stop("Only local solvers allowed: BOBYQA, COBYLA, LBFGS, MMA, SLSQP.")
    }
    .fn <- match.fun(fn)
    fn <- function(x) .fn(x, ...)
    if (!dfree && is.null(gr)) {
      gr <- function(x) nloptr::nl.grad(x, fn)
    }
    opts <- nloptr::nl.opts(control)
    opts$algorithm <- gsolver
    local_opts <- list(algorithm = lsolver, xtol_rel = localtol,
                       eval_grad_f = if (!dfree) gr else NULL)
    opts$local_opts <- local_opts
    if (!is.null(hin)) {
      .hin <- match.fun(hin)
      hin <- function(x) (-1) * .hin(x)
    }
    if (!dfree) {
      if (is.null(hinjac)) {
        hinjac <- function(x) nloptr::nl.jacobian(x, hin)
      }
      else {
        .hinjac <- match.fun(hinjac)
        hinjac <- function(x) (-1) * .hinjac(x)
      }
    }
    if (!is.null(heq)) {
      .heq <- match.fun(heq)
      heq <- function(x) .heq(x)
    }
    if (!dfree) {
      if (is.null(heqjac)) {
        heqjac <- function(x) nloptr::nl.jacobian(x, heq)
      }
      else {
        .heqjac <- match.fun(heqjac)
        heqjac <- function(x) .heqjac(x)
      }
    }
    S0 <- nloptr::nloptr(x0, eval_f = fn, eval_grad_f = gr, lb = lower,
                 ub = upper, eval_g_ineq = hin, eval_jac_g_ineq = hinjac,
                 eval_g_eq = heq, eval_jac_g_eq = heqjac, opts = opts)
    if (nl.info)
      print(S0)
    S1 <- list(par = S0$solution, value = S0$objective, iter = S0$iterations,
               global_solver = gsolver, local_solver = lsolver, convergence = S0$status,
               message = S0$message)
    return(S1)
  }

###################################################################
#####             Additional routines, July 2017             ######
###################################################################
#######################################
###  GEV in terms of return levels  ###
#######################################


#' Generalized Pareto maximum likelihood estimates for various quantities of interest
#'
#' This function calls the \code{gp.fit} routine on the sample of excesses and returns maximum likelihood
#' estimates for all quantities of interest, including scale and shape parameters, quantiles and value-at-risk,
#' expected shortfall and mean and quantiles of maxima of \code{N} threshold exceedances
#'
#' @param dat sample vector of excesses
#' @param args vector of strings indicating which arguments to return the maximum likelihood values for
#' @param m number of observations of interest for return levels. Required only for \code{args} values \code{"VaR"} or \code{"ES"}
#' @param N size of block over which to take maxima. Required only for \code{args} \code{Nmean} and \code{Nquant}.
#' @param p tail probability, equivalent to \eqn{1/m}. Required only for \code{args} \code{quant}.
#' @param q level of quantile for N-block maxima. Required only for \code{args} \code{Nquant}.
#' @return named vector with maximum likelihood values for arguments \code{args}
#' @export
#' @examples
#' dat <- evd::rgpd(n = 30, shape = 0.2)
#' gpd.mle(dat = dat, N = 100, p = 0.01, q = 0.5, m = 100)
gpd.mle <- function(dat, args = c("scale", "shape", "quant", "VaR", "ES", "Nmean", "Nquant"), m, N, p, q){
  args <- match.arg(args, c("scale", "shape", "quant", "VaR", "ES", "Nmean", "Nquant"), several.ok = TRUE)
  fitted <- try(gp.fit(xdat = dat, threshold = 0, method = "Grimshaw"))
  sigma <- fitted$estimate[1]; xi <- fitted$estimate[2]
  #Does not handle the case xi=0 because the optimizer does not return this value!
  a <- sapply(args, switch,
              scale = sigma,
              shape = xi,
              quant = sigma/xi*(p^(-xi)-1),
              Nquant = sigma/xi*((1-q^(1/N))^(-xi)-1),
              Nmean = (exp(lgamma(N + 1) + lgamma(1-xi) - lgamma(N + 1-xi))-1)*sigma/xi,
              VaR = sigma/xi*(m^xi-1),
              ES = ifelse(xi < 1, (sigma/xi*(m^xi-1) + sigma)/(1-xi), Inf)
  )
  a <- as.vector(unlist(a))
  names(a) = args
  return(a)
}




#'  Generalized extreme value maximum likelihood estimates for various quantities of interest
#'
#' #' This function calls the \code{fgev} routine on the sample of excesses and returns maximum likelihood
#' estimates for all quantities of interest, including scale and shape parameters, quantiles and value-at-risk,
#' expected shortfall and mean and quantiles of maxima of \code{N} threshold exceedances
#' @export
#' @param dat sample vector of excesses
#' @param args vector of strings indicating which arguments to return the maximum likelihood values for.
#' @param N size of block over which to take maxima. Required only for \code{args} \code{Nmean} and \code{Nquant}.
#' @param p tail probability. Required only for \code{arg} \code{quant}.
#' @param q level of quantile for maxima of \code{N} exceedances. Required only for \code{args} \code{Nquant}.
#' @return named vector with maximum likelihood estimated parameter values for arguments \code{args}
#' @examples
#' dat <- evd::rgev(n = 100, shape = 0.2)
#' gev.mle(dat = dat, N = 100, p = 0.01, q = 0.5)
gev.mle <- function(dat, args = c("loc", "scale", "shape", "quant", "Nmean", "Nquant"), N, p, q){
  args <- match.arg(args, c("loc", "scale", "shape", "quant", "Nmean", "Nquant"), several.ok = TRUE)
  fitted <- ismev::gev.fit(dat, method = "Nelder-Mead", show = FALSE)$mle
  xmax <- max(dat); xmin = min(dat)
  fitted <- nloptr::slsqp(x0 = fitted, fn = function(par){-gev.ll(par, dat = dat)},
                          gr = function(par){-gev.score(par, dat = dat)},
                          hin = function(par){ c(par[2] + par[3]*(xmax-par[1]), par[2] + par[3]*(xmin-par[1]))})$par
  mu <- fitted[1]; sigma <- fitted[2]; xi <- fitted[3]
  #Does not handle the case xi=0 because the optimizer does not return this value!
  a <- sapply(args, switch,
              loc = mu,
              scale = sigma,
              shape = xi,
              quant = evd::qgev(p = 1-p, loc = mu, scale = sigma, shape = xi),
              Nquant = ifelse(xi !=  0, mu - sigma / xi *(1 - (N / log(1/q))^xi),
                              mu + sigma * (log(N) - log(log(1/q)))),
              Nmean = ifelse(xi  !=  0, mu - sigma / xi *(1 - N^xi*gamma(1-xi)),
                             mu + sigma * (log(N) -psigamma(1)))
  )

  a <- as.vector(unlist(a))
  names(a) = args
  a
}






#' Confidence intervals for profile likelihood objects
#'
#' This function uses spline interpolation to derive \code{level} confidence intervals.
#'
#' @param object an object of class \code{extprof}, normally the output of \link{gpd.pll} or \link{gev.pll}.
#' @param parm a specification of which parameters are to be given confidence intervals,
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level, with default 0.95
#' @param ... additional arguments passed to functions. Providing a logical \code{warn=FALSE} turns off warning messages when the lower or upper confidence interval for \code{psi} are extrapolated beyond the provided calculations.
#' @return a 2 by 3 matrix containing point estimates, lower and upper confidence intervals based on the likelihood root and modified version thereof
#' @export
confint.extprof <- function(object, parm, level = 0.95, ...){
  args <- list(...)
  if("warn" %in% names(args) && is.logical(args$warn)){
    warn <- args$warn
  } else{
    warn <- TRUE
  }
  if(missing(parm)){
    parm <- NULL
    ind <- args$ind
    if(!is.null(object$pll) || !is.null(object$r)){
      parm <- c(parm, "profile")
      ind <- c(ind, 1)
    }
    if(!is.null(object$rstar)){
      parm <- c(parm, "tem")
      ind <- c(ind, 2)
    }
    if(!is.null(object$tem.pll)){
      parm <- c(parm, "modif.tem")
      ind <- c(ind, 3)
    }
    if(!is.null(object$empcov.pll)){
      parm <- c(parm, "modif.empcov")
      ind <- c(ind, 4)
    }

  } else {
   if(is.numeric(parm)){
     ind <- parm
     parm <-  c("profile","tem","modif.tem","modif.empcov")[ind]
   } else{
    parm <- match.arg(arg = parm, choices = c("profile","tem","modif.tem","modif.empcov","r","rstar"), several.ok = TRUE)
    parm[parm %in% "r"] <- "profile"
    parm[parm %in% "rstar"] <- "tem"
    ind <- which(c("profile","tem","modif.tem","modif.empcov") %in% parm)
   }
   parm <- unique(parm)
   ind <- unique(ind[ind %in% 1:4])
   }
  if(length(ind) == 0){
    stop("Invalid `parm` argument.")
  }
  conf <- NULL
  for(i in ind){
  if(i == 1){
    if(is.null(object$pll) && is.null(object$r)){ break;}
    if(is.null(object$r)){ #no r object, but must have pll + maxpll
      object$r <- sign(object$psi.max-object$psi)*sqrt(2*(object$maxpll-object$pll))
    }
    if(is.null(object$normal)){
      object$normal <- c(object$psi.max, object$std.error)
    }
    if(requireNamespace("cobs", quietly = TRUE)){
      fit.r <- cobs::cobs(x=object$r, y=object$psi, constraint = "decrease", lambda = 0, ic="SIC", pointwise = cbind(0,0,object$normal[1]),
                          knots.add = TRUE, repeat.delete.add = TRUE, print.mesg = FALSE, print.warn = FALSE)
      pr <- predict(fit.r,c(0,sqrt(qchisq(level,1)),-sqrt(qchisq(level,1))))[,2]
    } else{
      fit.r <- stats::smooth.spline(x=na.omit(cbind(object$r,object$psi)), cv=FALSE)
      pr <- predict(fit.r,c(0,sqrt(qchisq(level,1)),-sqrt(qchisq(level,1))))$y
      pr[1] <- object$normal[1]
    }
    conf <- cbind(conf, pr)
    if(warn){
      if(!any(object$r >  sqrt(qchisq(level,1)))){warning("Extrapolating the lower confidence interval for the profile likelihood ratio test")}
      if(!any(object$r < -sqrt(qchisq(level,1)))){warning("Extrapolating the upper confidence interval for the profile likelihood ratio test")}
    }
  } else if (i == 2){
    if(is.null(object$rstar)){ break;}
    if(requireNamespace("cobs", quietly = TRUE)){
      fit.rst <- cobs::cobs(x = object$rstar, y = object$psi, constraint = "decrease",
                            lambda = 0, ic = "SIC", knots.add = TRUE, repeat.delete.add = TRUE,
                            print.mesg = FALSE, print.warn = FALSE)
      prst <- predict(fit.rst,c(0, sqrt(qchisq(level,1)), -sqrt(qchisq(level,1))))[,2]
    } else{
      fit.rst <- stats::smooth.spline(x=na.omit(cbind(object$rstar,object$psi)), cv=FALSE)
      prst <- predict(fit.rst,c(0,sqrt(qchisq(level,1)),-sqrt(qchisq(level,1))))$y
    }
    if(!is.null(object$tem.psimax)){
      prst[1] <- object$tem.psimax
    }
    conf <- cbind(conf, prst)
    #lines(x=object$rstar,fit.rst$fitted,col=2,pch=19)
    if(warn){
      if(!any(object$rstar >  sqrt(qchisq(level,1)))){warning("Extrapolating the adjusted lower confidence interval for rstar.")}
      if(!any(object$rstar < -sqrt(qchisq(level,1)))){warning("Extrapolating the adjusted upper confidence interval for rstar")}
    }
  } else if (i == 3){
    if(is.null(object$tem.pll)){ break;}
    if(requireNamespace("cobs", quietly = TRUE)){
      fit.mtem <- cobs::cobs(x = sign(object$tem.mle-object$psi)*sqrt(-2*(object$tem.pll-object$tem.maxpll)),
                            y = object$psi, constraint = "decrease",
                            lambda = 0, ic = "SIC", knots.add = TRUE, repeat.delete.add = TRUE,
                            print.mesg = FALSE, print.warn = FALSE)
      ptem <- predict(fit.mtem,c(0, sqrt(qchisq(level,1)), -sqrt(qchisq(level,1))))[,2]
    } else{
      fit.mtem <- stats::smooth.spline(x=na.omit(cbind(
        sign(object$tem.mle-object$psi)*sqrt(-2*(object$tem.pll-object$tem.maxpll)),
        object$psi)), cv=FALSE)
      ptem <- predict(fit.mtem,c(0,sqrt(qchisq(level,1)),-sqrt(qchisq(level,1))))$y
    }
    ptem[1] <- object$tem.mle
    conf <- cbind(conf, ptem)
  } else if (i == 4){
    if(is.null(object$empcov.pll)){ break;}
    if(requireNamespace("cobs", quietly = TRUE)){
      fit.mempcov <- cobs::cobs(x = sign(object$empcov.mle-object$psi)*sqrt(-2*(object$empcov.pll-object$empcov.maxpll)),
                             y = object$psi, constraint = "decrease",
                             lambda = 0, ic = "SIC", knots.add = TRUE, repeat.delete.add = TRUE,
                             print.mesg = FALSE, print.warn = FALSE)
      pempcov <- predict(fit.mempcov,c(0, sqrt(qchisq(level,1)), -sqrt(qchisq(level,1))))[,2]
    } else{
      fit.mempcov <- stats::smooth.spline(x=na.omit(cbind(
        sign(object$empcov.mle-object$psi)*sqrt(-2*(object$empcov.pll-object$empcov.maxpll)),
        object$psi)), cv=FALSE)
      pempcov <- predict(fit.mempcov,c(0,sqrt(qchisq(level,1)),-sqrt(qchisq(level,1))))$y
    }
    pempcov[1] <- object$empcov.mle
    conf <- cbind(conf, pempcov)
  }
  }
  if(!is.null(conf)){
    colnames(conf) <- c("Profile LRT","TEM", "Modif. (TEM)","Modif. (emp. cov.)")[ind]
    rownames(conf) <- c("Estimate","Lower CI","Upper CI")
    return(conf)
  }
 }


#' Plot of (modified) profile likelihood
#'
#' The function plots the (modified) profile likelihood and the tangent exponential profile likelihood
#'
#' @param x an object of class \code{extprof} returned by \code{\link{gpd.pll}} or \code{\link{gev.pll}}.
#' @param ... further arguments to \code{plot}.
#' @return a graph of the (modified) profile likelihoods
#' @references Brazzale, A. R., Davison, A. C. and Reid, N. (2007). \emph{Applied Asymptotics: Case Studies in Small-Sample Statistics}. Cambridge University Press, Cambridge.
#' @references Severini, T. A. (2000). \emph{Likelihood Methods in Statistics}. Oxford University Press, Oxford.
#' @export
plot.extprof <- function(x, ...)
{ # plot the profile log-likelihoods
  old.pars <- par(no.readonly = TRUE)
    lik <- list()
    if(is.null(x$pll) && !is.null(x$r)){
      lik$npll <- -x$r^2/2
    } else if(!is.null(x$pll)){
      lik$npll <- x$pll-x$maxpll
    } else{
     stop("Invalid object provided")
    }
    if(!is.null(x$tem.pll)){
     lik$tem.npll <-  x$tem.pll-x$tem.maxpll
    }
    if(!is.null(x$empcov.pll)){
      lik$empcov.npll <-  x$empcov.pll-x$empcov.maxpll
    }
    ylim <- c(min(unlist(lapply(lik, min, na.rm = TRUE))), 0)
    args <- list(...)
    tikz <- FALSE
    level <- c(0.95,0.99)
    if(!is.null(args$level)){
     level <- args$level[1]
    }
    if(!is.null(args$tikz)){
     if(isTRUE(args$tikz)){
       tikz <- TRUE
     }
    }
    if(any(!is.null(args$which), !is.null(args$ind), !is.null(args$parm))){
      if(!is.null(args$parm)){
      parm <- match.arg(arg = args$parm, choices = c("profile","tem","modif.tem","modif.empcov","r","rstar"), several.ok = TRUE)
      parm[parm %in% "r"] <- "profile"
      parm[parm %in% "rstar"] <- "tem"
      ind <- which(c("profile","tem","modif.tem","modif.empcov") %in% parm)
    } else if(!is.null(args$ind)){
      ind <- args$ind[args$ind %in% 1:4]
      parm <- c("profile","tem","modif.tem","modif.empcov")[ind]
    } else if(!is.null(args$which)){
      ind <- args$which[args$which %in% 1:4]
      parm <- c("profile","tem","modif.tem","modif.empcov")[ind]
    }
    parm <- unique(parm)
    ind <- unique(ind[ind %in% 1:4])
    } else{
      ind <- 1:4
      parm <- c("profile","tem","modif.tem","modif.empcov")
    }
    plot(NULL,type = "n", bty = "l", xlim = c(min(x$psi), max(x$psi)), ylim = ylim,
         xlab = ifelse(!tikz, expression(psi), "$\\psi$"), ylab = "Profile log likelihood")
    abline(h = -qchisq(level,1)/2, col = 'gray')
    #Legend
    lcols <- NULL; llty <- NULL; llwd <- NULL; llegend <- NULL

    if(4 %in% ind && !is.null(x$empcov.mle)){
      abline(v = x$empcov.mle, lwd = 0.5, col = 4, lty = 4)
    }
    if(3 %in% ind && !is.null(x$tem.mle)){
      abline(v = x$tem.mle, lwd = 0.5, col = 2, lty = 2)
    }
    if(2 %in% ind && !is.null(x$tem.psimax)){
      abline(v = x$tem.psimax, lwd = 0.5, col = 3)
    }
    abline(v = x$psi.max, lwd = 0.5)

    if(4 %in% ind && !is.null(lik$empcov.npll)){
      lines(x$psi, lik$empcov.npll, lty = 4, col = 4, lwd=2)
      lcols <- c(lcols, 4); llty <- c(llty, 4);
      llwd <- c(llwd, 2); llegend <- c(llegend, "modif. emp. cov.")
    }
    if(3 %in% ind && !is.null(lik$tem.npll)){
      lines(x$psi, lik$tem.npll, lty = 2, col = 2, lwd = 2)
      lcols <- c(lcols, 2); llty <- c(llty, 2);
      llwd <- c(llwd, 2); llegend <- c(llegend, "modif. tem.")
    }
    if(2 %in% ind && !is.null(x$rstar)){
      lines(x$psi, -x$rstar^2/2, lwd = 2, col = 3, lty = 5)
      lcols <- c(lcols, 3); llty <- c(llty, 5);
      llwd <- c(llwd, 2); llegend <- c(llegend, "tem")
    }
    lcols <- c(lcols, 1); llty <- c(llty, 1);
    llwd <- c(llwd, 2); llegend <- c(llegend, "profile")
    lines(x$psi, lik$npll, lwd = 2)
    #add the legend in the top right corner
    legend(x = "topright", legend = rev(llegend),  lty = rev(llty), lwd = rev(llwd),
           col = rev(lcols), bty="n", x.intersp = 0.2, seg.len = 0.5, cex = 0.9)
  par(old.pars)

}



#' Modified profile likelihood for the generalized extreme value distribution
#'
#' This function calculates the profile likelihood along with two small-sample corrections
#' based on Severini's (1999) empirical covariance and the Fraser and Reid tangent exponential
#' model approximation.
#'
#' @details The two \code{mod} available are \code{tem}, the tangent exponential model (TEM) approximation and
#' \code{modif} for the penalized profile likelihood based on \eqn{p^*} approximation proposed by Severini.
#' For the latter, the penalization is based on the TEM or an empirical covariance adjustment term.
#'
#' @param psi parameter vector over which to profile (unidimensional)
#' @param param string indicating the parameter to profile over
#' @param mod string indicating the model. See \bold{Details}.
#' @param dat sample vector
#' @param N size of block over which to take maxima. Required only for \code{param} \code{Nmean} and \code{Nquant}.
#' @param p tail probability. Required only for \code{param} \code{quant}.
#' @param q probability level of quantile. Required only for \code{param} \code{Nquant}.
#' @param correction logical indicating whether to use \code{spline.corr} to smooth the tem approximation.
#' @param ... additional arguments such as output from call to \code{Vfun} if \code{mode="tem"}.
#'
#' @return a list with components
#' \itemize{
#' \item{\code{mle}:} maximum likelihood estimate
#' \item{\code{psi.max}:} maximum profile likelihood estimate
#' \item{\code{param}:} string indicating the parameter to profile over
#' \item{\code{std.error}:} standard error of \code{psi.max}
#' \item{\code{psi}:} vector of parameter \eqn{psi} given in \code{psi}
#' \item{\code{pll}:} values of the profile log likelihood at \code{psi}
#' \item{\code{maxpll}:} value of maximum profile log likelihood
#' }
#'
#'
#' In addition, if \code{mod} includes \code{tem}
#' \itemize{
#' \item{\code{normal}:}{maximum likelihood estimate and standard error of the interest parameter \eqn{psi}}
#' \item{\code{r}:}{values of likelihood root corresponding to \eqn{\psi}}
#' \item{\code{q}:}{vector of likelihood modifications}
#' \item{\code{rstar}:}{modified likelihood root vector}
#' \item{\code{rstar.old}:}{uncorrected modified likelihood root vector}
#' \item{\code{tem.psimax}:}{maximum of the tangent exponential model likelihood}
#' }
#' In addition, if \code{mod} includes \code{modif}
#' \itemize{
#' \item{\code{tem.mle}:} maximum of tangent exponential modified profile log likelihood
#' \item{\code{tem.profll}:} values of the modified profile log likelihood at \code{psi}
#' \item{\code{tem.maxpll}:} value of maximum modified profile log likelihood
#' \item{\code{empcov.mle}:} maximum of Severini's empirical covariance modified profile log likelihood
#' \item{\code{empcov.profll}:} values of the modified profile log likelihood at \code{psi}
#' \item{\code{empcov.maxpll}:} value of maximum modified profile log likelihood
#' }
#'
#' @references Fraser, D. A. S., Reid, N. and Wu, J. (1999), A simple general formula for tail probabilities for frequentist and Bayesian inference. \emph{Biometrika}, \bold{86}(2), 249--264.
#' @references Severini, T. (2000) Likelihood Methods in Statistics. Oxford University Press. ISBN 9780198506508.
#' @references Brazzale, A. R., Davison, A. C. and Reid, N. (2007) Applied asymptotics: case studies in small-sample statistics. Cambridge University Press, Cambridge. ISBN 978-0-521-84703-2
#'
#' @export
#' @examples
#' \dontrun{
#' dat <- evd::rgev(n = 100, loc = 0, scale = 2, shape = 0.3)
#' gev.pll(psi = seq(-0.5,1, by=0.01), param = "shape", dat = dat)
#' gev.pll(psi = seq(-3, 3, length = 50), param = "loc", dat = dat)
#' gev.pll(psi = seq(10, 30, by = 0.1), param = "quant", dat = dat, p = 0.01)
#' gev.pll(psi = seq(12, 100, by=1), param = "Nmean", N = 100, dat = dat)
#' gev.pll(psi = seq(12, 90, by=1), param = "Nquant", N = 100, dat = dat, q = 0.5)
#' }
gev.pll <- function(psi, param = c("loc", "scale", "shape", "quant", "Nmean", "Nquant"),
                    mod = c("tem","modif"), dat, N = NULL, p = NULL, q = NULL, correction = TRUE, ...){

  oldpar <- param <- match.arg(param, c("loc", "scale", "shape", "quant", "Nmean", "Nquant"))[1]
  mod <- match.arg(mod, c("tem","modif"), several.ok = TRUE)
  #Parametrization profiling over quant over scale is more numerically stable
  if(param == "quant"){
    stopifnot(!is.null(p))
    q <- 1 - p; N = 1; param <- "Nquant"
  }

  #Arguments for parametrization of the log likelihood
  if(param %in% c("loc","scale","shape")){
    args <- c("loc","scale","shape")
  } else if (param == "quant"){
    args <- c(param, "scale", "shape")
  } else{
    args <- c("loc", param, "shape")
  }
  #Sanity checks to ensure all arguments are provided
  if(is.null(N)){
    if(param %in% c("Nmean","Nquant")){
      stop("Argument `N` missing. Procedure aborted")
    } else { N <- NA
    }
  }
  if(is.null(q)){
    if(param %in% c("Nquant")){
      stop("Argument `q` missing. Procedure aborted")
    } else { q <- NA
    }
  }
  if(is.null(p)){
    if(param %in% c("quant")){
      stop("Argument `p` missing. Procedure aborted")
    } else { p <- NA
    }
  }
  xmin <- min(dat); xmax <- max(dat)
  #Find maximum likelihood estimates
  mle <- gev.mle(dat = dat, args = args, q = q, N = N, p = p)

  # if(missing(psi) || is.null(psi) || is.na(psi)){
  #   psi <- mle[param]
  # }
  # psi <- as.vector(psi)


  #Extract the components, notably V for model `tem`.
  #Keep other components for optimization
  Vprovided <- FALSE
  extra.args <- list(...)
  if("V" %in% names(extra.args)){
    V <- extra.args$V
    extra.args$V <- NULL
    if(isTRUE(all.equal(dim(V), c(length(dat), 2)))){
      Vprovided <- TRUE
    }
  }
  if(!Vprovided){
    V <- switch(param,
                loc = gev.Vfun(par = mle, dat = dat),
                scale = gev.Vfun(par = mle, dat = dat),
                shape = gev.Vfun(par = mle, dat = dat),
                quant = gevr.Vfun(par = mle, dat = dat,  p = p),
                Nmean = gevN.Vfun(par = mle, dat = dat, N = N, qty = "mean"),
                Nquant = gevN.Vfun(par = mle, dat = dat,  q = q, N = N, qty = "quantile"))
  }


  #Obtained constrained maximum likelihood estimates for given value of psi
  if(param %in% c("loc", "scale", "shape")){

    #Define observation-wise gradient
    gev.score.f <-  function(par, dat){
      dat <- as.vector(dat)
      mu = par[1]; sigma = par[2]; xi = as.vector(par[3])
      if(!isTRUE(all.equal(xi, 0))){
        cbind(-(-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)/sigma - xi*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)),
              -(dat - mu)*((dat - mu)*xi/sigma + 1)^(-1/xi - 1)/sigma^2 + (dat - mu)*xi*(1/xi + 1)/(sigma^2*((dat - mu)*xi/sigma + 1)) - 1/sigma,
              -(mu - dat)*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) -
                (log(-(mu - dat)*xi/sigma + 1)/xi^2 - (mu - dat)/(sigma*((mu - dat)*xi/sigma - 1)*xi))/(-(mu - dat)*xi/sigma + 1)^(1/xi)
              + log(-(mu - dat)*xi/sigma + 1)/xi^2
        )
      } else {
        cbind(-exp(mu/sigma - dat/sigma)/sigma + 1/sigma,
              mu*exp(mu/sigma - dat/sigma)/sigma^2 - dat*exp(mu/sigma - dat/sigma)/sigma^2 - mu/sigma^2 - 1/sigma + dat/sigma^2,
              rep(0,length(dat)))
      }
    }
    ind <- switch(param, loc=1, scale=2, shape=3)
    maxll <- gev.ll(mle, dat = dat)
    std.error <- sqrt(solve(gev.infomat(par=mle, dat = dat, method="exp"))[ind,ind])
    constr.mle.scale <- function(sigmat, dat = dat){
        x0 = c(median(dat/sigmat), 0.05)
      if(is.nan(gev.ll(c(x0[1], sigmat, x0[2]), dat = dat))){
        constr_fit <- try(evd::fgev(x = dat, std.err = FALSE, scale = sigmat, shape = x0[2]))
        if(!is.character(constr_fit)){
          if(constr_fit$convergence == "successful"){
            x0 <- as.vector(c(constr_fit$estimate['loc'], 0.05))
          } else{
            stop("Could not find starting values for optimization routine")
          }
        } else{
          stop("Could not find starting values for optimization routine")
        }
      }
      opt <- nloptr::sbplx(x0 = x0, fn = function(par){-gev.ll(c(par[1], sigmat, par[2]), dat = dat)})
      opt2 <- nloptr::slsqp(x0 = opt$par, fn = function(par){-gev.ll(c(par[1], sigmat, par[2]), dat = dat)},
                            gr = function(par){-gev.score(c(par[1], sigmat, par[2]), dat = dat)[-2]})
      if(opt2$convergence > 0){
        return(c(opt2$par, opt2$value))
      } else{
        return(rep(NA, 3))
      }
    }
    constr.mle.loc <- function(mut, dat = dat){
      opt <- suppressWarnings(nloptr::auglag(x0 = c(mad(dat, constant = 1),0.1), fn = function(par){
        val <- -gev.ll(c(mut, par[1:2]), dat = dat); ifelse(is.infinite(val), 1e10,val)},
        gr = function(par){-gev.score(c(mut, par[1:2]), dat = dat)[-ind]},
        hin = function(par){ifelse(par[2] <= 0, par[1] + par[2]*(xmax-mut), par[1] + par[2]*(xmin-mut))},
        localsolver = "SLSQP"))
      if(opt$convergence > 0){
        return(c(opt$par, opt$value))
      } else{
        return(rep(NA, 3))
      }
    }
    # constr.mle.shape <- function(xit, dat = dat){
    #   start.scale <- max(1e-2,mad(dat, constant = 1)/median(abs(evd::rgev(10000, shape=xit)-evd::qgev(0.5, shape=xit))))
    #   start.loc <- median(dat) - start.scale*ifelse(xit == 0, log(log(2)), (log(2)^(-xit)-1)/xit)
    #   if(any(c(start.scale + xit*(xmax-start.loc) <= 0, start.scale + xit*(xmin-start.loc) <= 0))){
    #     if(xit < 0){
    #       start.loc <-  start.scale/xit + xmax + 1e-3
    #     } else {
    #       start.loc <-  start.scale/xit + xmin + 1e-3
    #     }
    #   }
    #   #Subplex - simplex type algorithm, more robust than Nelder-Mead
    #   opt <- nloptr::sbplx(x0 = c(start.loc, start.scale), fn = function(par){-gev.ll(c(par,xit), dat=dat)}, lower= c(-Inf, 0),
    #                        control = list(xtol_rel = 1e-8, maxeval = 2000, ftol_rel = 1e-10))
    #   #If solution is not within the region
    #   if(ifelse(xit < 0, opt$par[2] + xit*(xmax-opt$par[1]) <= 0, opt$par[2] + xit*(xmin-opt$par[1]) <= 0)){
    #   opt <- nloptr::auglag(x0 = c(start.loc, start.scale),
    #                         fn = function(par){
    #                           val <- gev.ll.optim(c(par[1:2], xit), dat = dat); ifelse(is.infinite(val) || is.na(val), 1e10, val)},
    #                         gr = function(par){ val <- -gev.score(c(par[1], exp(par[2]), xit), dat = dat)[-ind]},
    #                         hin = function(par){c(ifelse(xit <= 0, exp(par[2]) + xit*(xmax-par[1]), exp(par[2]) + xit*(xmin-par[1])))}
    #     )
    #   }
    #   if(!is.character(opt)){
    #     if(all(c(opt$convergence > 0, abs(gev.score(c(opt$par, xit), dat = dat)[1:2]) < 5e-4))){
    #       return(c(opt$par, opt$value))
    #     } else {
    #       #evd::fgev(start = list(loc = opt$par[1], scale = opt$par[2]), shape = xit,
    #       # x = dat, method = "BFGS", control=list(reltol=1e-10, abstol = 1e-9))
    #     opt2 <- suppressWarnings(nloptr::slsqp(x0 = opt$par, fn = function(par){
    #       val <- gev.ll.optim(c(par[1:2], xit), dat = dat); ifelse(is.infinite(val) || is.na(val), 1e10, val)},
    #       gr = function(par){ val <- -gev.score(c(par[1], exp(par[2]), xit), dat = dat)[-ind]},
    #       hin = function(par){c(ifelse(xit <= 0, exp(par[2]) + xit*(xmax-par[1]), exp(par[2]) + xit*(xmin-par[1])))}
    #     ))
    #     opt2$par[2] <- exp(opt2$par[2])
    #     if(all(c(opt2$convergence > 0, !isTRUE(all.equal(opt2$value, 1e10)), abs(gev.score(c(opt2$par, xit), dat = dat)[1:2]) < 1e-4))){
    #        return(c(opt2$par, opt2$value))
    #     }
    #     }
    #   }
    #   return(rep(NA, 3))
    # }
    constr.mle.shape <- function(xit, dat = dat){
      if(abs(xit) < 1e-8){xit <- 0} #because rgev does not handle this case!
      start.scale <- mad(dat, constant = 1)/median(abs(evd::rgev(2000, shape=xit)-evd::qgev(0.5, shape=xit)))
      start.loc <- median(dat) - start.scale*ifelse(xit == 0, log(log(2)), (log(2)^(-xit)-1)/xit)
      if(start.scale + xit*(xmax-start.loc) <= 0){
        if(xit < 0){
          start.loc <-  start.scale/xit + xmax + 1e-3
        } else {
          start.loc <-  start.scale/xit + xmin + 1e-3
        }
      }
      opt <- try(suppressWarnings(nloptr::auglag(x0 = c(start.loc, start.scale), fn = function(par){
        val <- -gev.ll(c(par[1:2], xit), dat = dat); ifelse(is.infinite(val) || is.na(val), 1e10, val)},
        hin = function(par){ifelse(xit <= 0, par[2] + xit*(xmax-par[1]), par[2] + xit*(xmin-par[1]))},
        localsolver="SLSQP"
      )))
      if(!is.character(opt)){
        if(opt$convergence > 0 && !isTRUE(all.equal(opt$value, 1e10))){
          opt2 <- suppressWarnings(nloptr::slsqp(x0 = opt$par, fn = function(par){
            val <- -gev.ll(c(par[1:2], xit), dat = dat); ifelse(is.infinite(val) || is.na(val), 1e10, val)},
            gr = function(par){-gev.score(c(par[1:2], xit), dat = dat)[-ind]},
            hin = function(par){ifelse(xit <= 0, par[2] + xit*(xmax-par[1]), par[2] + xit*(xmin-par[1]))}
          ))
          return(c(opt2$par, opt2$value))
        }
      }
      return(rep(NA, 3))
    }
    #Missing psi vector
    if(missing(psi) || is.null(psi) || is.na(psi)){
      if(ind==2){
        psirangelow <- unique(min(1e-5,seq(-4, -1.5, length = 10)*std.error + mle[param]))
      } else{
        psirangelow <- seq(-4, -1.5, length = 10)*std.error + mle[param]
      }
      lowvals <-  -t(sapply(psirangelow, function(par){
        switch(param, loc = constr.mle.loc(par, dat = dat),
               scale = constr.mle.scale(par, dat = dat),
               shape = constr.mle.shape(par, dat = dat))}))[,3]-maxll
      psirangehigh <- seq(1.5, 4, length = 10)*std.error + mle[param]
      highvals <-  -t(sapply(psirangehigh, function(par){
        switch(param, loc = constr.mle.loc(par, dat = dat),
               scale = constr.mle.scale(par, dat = dat),
               shape = constr.mle.shape(par, dat = dat))}))[,3]-maxll
      lo <- approx(x=lowvals, y=psirangelow, xout = -4)$y
      hi <- approx(x=highvals, y=psirangehigh, xout = -4)$y
      lo <- ifelse(is.na(lo), predict(lm(psirangelow~lowvals), newdata = data.frame(lowvals = -4))[[1]], lo)
      hi <- ifelse(is.na(hi), predict(lm(psirangehigh~highvals), newdata = data.frame(highvals = -4))[[1]], hi)
      psi <- seq(lo, hi, length = 55)
    }
    if(param == "shape" && min(psi) < -1){
      warning("psi includes shape parameters below -1. These will be automatically removed.")
      psi <- psi[psi > -1]
    }
    #Calculate profile likelihood at psi
    lambda <- t(sapply(psi, function(par){
      switch(param, loc = constr.mle.loc(par, dat = dat),
             scale = constr.mle.scale(par, dat = dat),
             shape = constr.mle.shape(par, dat = dat))}))
    pars <- switch(param,
                   loc = cbind(psi, lambda[,1:2, drop=FALSE]),
                   scale = cbind(lambda[,1, drop=FALSE], psi, lambda[,2, drop=FALSE]),
                   shape = cbind(lambda[,1:2, drop=FALSE], psi))
    #Profile log likelihood values for psi
    profll <- -lambda[,3]
    if("tem" %in% mod){
      #Tangent exponential model approximation of Fraser and Reid to the profile likelihood
      r <- sign(mle[param]-psi)*sqrt(2*(maxll - profll))
      phi.mle <- gev.phi(par = mle, dat = dat, V = V)
      q2num <- ifelse(ind%%2==0, -1, 1)*
        apply(pars, 1, function(par){ det(rbind(c(c(phi.mle) - gev.phi(par = par, dat = dat, V = V)),
                                                  gev.dphi(par = par, dat = dat, V = V)[-ind,]))})
      if(isTRUE(any(sign(q2num)*sign(r) < 0, na.rm = TRUE))){
        warning("Correction factor and likelihood root are of opposite sign - check output")
      }
      logq <- apply(pars, 1, function(par){
        -0.5*log(det(gev.infomat(par = par, dat = dat, method = "obs")[-ind, -ind]))}) + log(abs(q2num))
      qmlecontrib <-  -log(det(gev.dphi(par = mle, dat = dat, V = V))) +
        0.5*log(det(gev.infomat(par = mle, dat = dat, method = "obs")))
      logq <- logq + qmlecontrib
      qcor <- sign(q2num)*exp(logq)
      rstar <- ifelse(r==0, 0, r + (logq-log(abs(r)))/r)
      ###

      tem.max.opt <- function(psi, dat = dat){
        lambda <- switch(param, loc = constr.mle.loc(psi, dat = dat),
                 scale = constr.mle.scale(psi, dat = dat),
                 shape = constr.mle.shape(psi, dat = dat))[1:2]
        para <- switch(ind, c(psi, lambda), c(lambda[1], psi, lambda[2]), c(lambda, psi))
        pll <- gev.ll(par = para, dat = dat)
        rs <- 2*(maxll - pll)
        logq <- -0.5*log(det(gev.infomat(par = para, dat = dat, method = "obs")[-ind, -ind])) + qmlecontrib +
          log(abs(det(rbind(c(c(phi.mle) - gev.phi(par = para, dat = dat, V = V)),
                            gev.dphi(par = para, dat = dat, V = V)[-ind,]))))
        abs(rs + logq - log(sqrt(abs(rs))))
      }
       tem.max <- optim(par = mle[ind], fn = tem.max.opt, method = "Brent", dat = dat,
                       lower = ifelse(ind == 2, max(1e-5, mle[ind] - std.error), mle[ind] - std.error),
                       upper = mle[ind] + std.error, control = list(abstol=1e-10))$par
      ###
    }
    #Modified profile likelihood based on p* approximation, two modifications due to Severini
    if("modif" %in% mod){
      tem.objfunc <- function(par, dat = dat){
        0.5*log(det(gev.infomat(par, dat = dat, method = "obs")[-ind,-ind])) -
          log(det(gev.dphi(par = par, dat = dat, V = V[, -ind])[-ind,]))}
      optim.tem.fn <- function(psi, dat = dat, param = param){
        theta.psi.opt <- switch(param,
                                loc = constr.mle.loc(psi, dat = dat),
                                scale = constr.mle.scale(psi, dat = dat),
                                shape = constr.mle.shape(psi, dat = dat))
        para <- switch(param,
                        loc = c(psi, theta.psi.opt[1:2]),
                        scale = c(theta.psi.opt[1], psi, theta.psi.opt[2]),
                        shape = c(theta.psi.opt[1:2], psi)
        )
        ll <- -theta.psi.opt[3]
        ll + tem.objfunc(para, dat = dat)
      }
      #TEM profile log likelihood values for psi
      proflltem <- profll + apply(pars, 1, tem.objfunc, dat = dat)
      #Maximum objective function for TEM (line search in neighborhood of the MLE)
      tem.mle.opt <- optim(par = mle[ind], fn = optim.tem.fn, method = "Brent", dat = dat, param = param,
                           lower = ifelse(param == "scale", max(1e-5, mle[ind] - std.error), mle[ind] - std.error),
                           upper = mle[ind] + std.error, control = list(fnscale=-1))


      #Severini empirical covariance function adjustment to profile likelihood
      #Score function - observation-wise

      #Score at MLE (sums to zero)
      score.scale.mle <- gev.score.f(mle, dat)[,-ind] #keep s_lambda
      empcov.objfunc <- function(par, dat){
        0.5*log(det(gev.infomat(par = par, dat = dat, method = "obs")[-ind,-ind])) -
          log(sum(score.scale.mle*gev.score.f(par, dat)[,-ind]))
      }
      profllempcov <- profll + apply(pars, 1, empcov.objfunc, dat = dat)
      optim.empcov.fn <- function(psi, param = param, dat = dat){
        theta.psi.opt <- switch(param,
                                loc = constr.mle.loc(psi, dat = dat),
                                scale = constr.mle.scale(psi, dat = dat),
                                shape = constr.mle.shape(psi, dat = dat))
        para <- switch(param,
                        loc = c(psi, theta.psi.opt[1:2]),
                        scale = c(theta.psi.opt[1], psi, theta.psi.opt[2]),
                        shape = c(theta.psi.opt[1:2], psi)
        )
        ll <- -theta.psi.opt[3]
        ll + empcov.objfunc(para, dat = dat)
      }

      empcov.mle.opt <- optim(par = mle[ind], fn=optim.empcov.fn, method = "Brent", dat = dat, param = param,
                              lower = ifelse(param == "scale", max(1e-5, mle[ind] - std.error), mle[ind]-std.error),
                              upper = mle[ind] + std.error, control=list(fnscale=-1))
    }

    #Return levels, quantiles or value-at-risk
  } else if(param == "quant") {
    maxll <- gevr.ll(mle, dat = dat, p = p)
    std.error <- sqrt(solve(gevr.infomat(par=mle, dat = dat, method="exp", p = p))[1])
    constr.mle.quant <- function(quant){
      fitted <- try(evd::fgev(x = dat, prob = p, quantile = quant, method = "Neld", std.err=FALSE))
      if(!is.character(fitted)){
          fitted <- optim(par = list(scale = log(fitted$estimate[1]), shape = fitted$estimate[2]),
                            dat = dat, p = p, fn = function(lambda, p, dat){gevr.ll.optim(c(quant, lambda), dat = dat, p = p)},
                          method = "BFGS", control = list(maxit = 500))
        return(c(exp(fitted$par[1]), fitted$par[2], -fitted$value))
      } else{ return(rep(NA,3))}
    }

    #Missing psi vector
    if(missing(psi) || is.null(psi) || is.na(psi)){
      psirangelow <- seq(-3, -1.5, length = 6)*std.error + mle[1]
      lowvals <-  sapply(psirangelow, constr.mle.quant)[3,]-maxll
      psirangehigh <- seq(1.5, 4, length = 10)*std.error + mle[1]
      highvals <-  sapply(psirangehigh, constr.mle.quant)[3,]-maxll
      lo <- approx(x=lowvals, y=psirangelow, xout = -4)$y
      hi <- approx(x=highvals, y=psirangehigh, xout = -4)$y
      lo <- ifelse(is.na(lo), lm(psirangelow~lowvals)$coef[2]*-4 + mle[2], lo)
      hi <- ifelse(is.na(hi), lm(psirangehigh~highvals)$coef[2]*-4 + mle[2], hi)
      psi <- seq(lo, hi, length = 55)
    }



    evals <- t(sapply(psi, constr.mle.quant))
    pars <- cbind(psi, evals[,1:2, drop=FALSE])
    #Profile log likelihood values for psi
    profll <- evals[,3] # apply(pars, 1, function(par){gevr.ll(par = par, dat = dat, p = p)})
    if("tem" %in% mod){
      r <- sign(mle[param]-psi)*sqrt(2*(maxll - profll))
      phi.mle <- gevr.phi(par = mle, dat = dat, V = V, p = p)
      q2num <- apply(pars, 1, function(par){ det(rbind(c(c(phi.mle) - gevr.phi(par = par, dat = dat, V = V, p = p)),
                                                       gevr.dphi(par = par, dat = dat, V = V, p = p)[-1,]))})
      if(isTRUE(any(sign(q2num)*sign(r) < 0, na.rm=TRUE))){
        warning("Correction factor and likelihood root are of opposite sign - check output")
      }
      logq <- apply(pars, 1, function(par){ -0.5*log(det(gevr.infomat(par = par, dat = dat, method = "obs", p = p)[-1, -1]))})  + log(abs(q2num))
      qmlecontrib <- -log(det(gevr.dphi(par = mle, dat = dat, V = V, p = p))) +
        0.5*log(det(gevr.infomat(par = mle, dat = dat, method = "obs", p = p)))
      logq <- logq + qmlecontrib
      qcor <- sign(q2num)*exp(logq)
      rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r))) / r)

      tem.max.opt <- function(psi, dat = dat){
        lam <- constr.mle.quant(psi)
        para <- c(psi, lam[1:2])
        pll <- lam[3]
        rs <- 2*(maxll - pll)
        logq <- -0.5*log(det(gevr.infomat(par = para, dat = dat, method = "obs", p = p)[-1, -1])) + qmlecontrib +
          log(abs(det(rbind(c(c(phi.mle) - gevr.phi(par = para, dat = dat, V = V, p = p)),
                            gevr.dphi(par = para, dat = dat, V = V, p = p)[-1,]))))
        abs(rs + logq - log(sqrt(abs(rs))))
      }
      tem.max <- optim(par = mle[1], fn = tem.max.opt, method = "Brent",
                       lower = max(1e-5, mle[1] - std.error), dat = dat,
                       upper = mle[1] + std.error, control = list(abstol=1e-10))$par
      }
    if("modif" %in% mod){
      #Tangent exponential model approximation of Fraser and Reid to the profile likelihood
      tem.objfunc.quant <- function(par){
        0.5*log(det(gevr.infomat(par = par, dat = dat, method = "obs", p = p)[-1,-1])) -
          log(det(gevr.dphi(par = par, dat = dat, p = p, V = V[, -1])[-1,]))
      }
      optim.tem.fn.quant <- function(psi){
        theta.psi.opt <- constr.mle.quant(psi)
        param <- c(psi, theta.psi.opt[1:2])
        ll <- theta.psi.opt[3] # gevr.ll(param, dat = dat, m = m)
        ll + tem.objfunc.quant(param)
      }
      #TEM profile log likelihood values for psi
      proflltem <- profll + suppressWarnings(apply(pars, 1, tem.objfunc.quant))
      #Maximum objective function for TEM
      tem.mle.opt <- optim(par=mle[1], fn=optim.tem.fn.quant, method = "Brent",
                           lower = max(1e-5,mle[1]-std.error), upper = mle[1] + std.error, control=list(fnscale=-1))

      #Severini empirical covariance function adjustment to profile likelihood
      #Score function - observation-wise for xi
      gevr.score.f <-  function(par, dat, p){
        sigma <- par[2]; xi <- par[3];  z <- par[1]
        cbind(((((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(-1/xi - 2)*xi*(1/xi + 1)*exp(-1/((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi))/sigma^2 - ((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(-2/xi - 2)*exp(-1/((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi))/sigma^2)*sigma*((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi + 1)*exp(1/(((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi)))),
              (-(dat*(-log(-p + 1))^xi - z*(-log(-p + 1))^xi - (dat*(-log(-p + 1))^xi - z*(-log(-p + 1))^xi - sigma)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))/((sigma*dat*xi*(-log(-p + 1))^xi - sigma*xi*z*(-log(-p + 1))^xi + sigma^2)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))),
              (-(xi*z*(-log(-p + 1))^xi - (dat*(-log(-p + 1))^xi - sigma*log(-log(-p + 1)))*xi + ((dat*(-log(-p + 1))^xi - sigma*log(-log(-p + 1)))*xi^2 + (dat*(-log(-p + 1))^xi - sigma*log(-log(-p + 1)))*xi - (xi^2*(-log(-p + 1))^xi + xi*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + (dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi - (dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + sigma)*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi)))/((dat*xi^3*(-log(-p + 1))^xi - xi^3*z*(-log(-p + 1))^xi + sigma*xi^2)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))
              ))

      }
      #Score at MLE (sums to zero)
      score.quant.mle <- gevr.score.f(mle, dat, p = p)
      empcov.objfunc.quant <- function(par){
        0.5*log(det(gevr.infomat(par = par, dat = dat, method = "obs", p = p)[-1,-1])) - log(sum(score.quant.mle*gevr.score.f(par, dat, p = p)))
      }
      profllempcov <- profll + suppressWarnings(apply(pars, 1, empcov.objfunc.quant))
      optim.empcov.fn.quant <- function(psi){
        theta.psi.opt <- constr.mle.quant(psi)
        param <- c(psi, theta.psi.opt[1:2])
        ll <- theta.psi.opt[3] #gevr.ll(param, dat = dat, p = p)
        ll +  empcov.objfunc.quant(param)
      }
      empcov.mle.opt <- optim(par=mle[1], fn=optim.empcov.fn.quant, method = "Brent",
                              lower = max(1e-5,mle[1]-std.error), upper = mle[1] + std.error, control=list(fnscale=-1))
    }
  }  else if(param %in% c("Nquant","Nmean")){

    qty <- switch(param, Nquant = "quantile", Nmean = "mean")
    maxll <- gevN.ll(mle, dat = dat, N = N, q = q, qty = qty)
    std.error <- sqrt(solve(gevN.infomat(par=mle, dat = dat, method="exp", N = N, q = q, qty = qty))[2,2])
    constr.mle.N <- function(zt, dat = dat){
      #browser();
      st_vals <- c(median(dat), 0.1)
      if(isTRUE(as.vector(mle["shape"]>0))){
        st_vals <- mle[c("loc","shape")]
      }
      opt <- try(suppressWarnings(
        .auglag(x0 = st_vals, fn = function(par){
        val <- -gevN.ll(par = c(par[1], zt, par[2]), dat = dat, q = q, qty = qty, N = N);
        ifelse(isTRUE(any(is.infinite(val), is.na(val), val <= - maxll)), 1e10, val)},
        hin = function(par){
          sigma <- switch(qty, quantile = (zt-par[1])*par[2]/(N^par[2]*(log(1/q))^(-par[2])-1),
                          mean = (zt-par[1])*par[2]/(N^par[2]*gamma(1-par[2])-1))
          c(sigma, sigma + par[2]*(xmax-par[1]), sigma + par[2]*(xmin-par[1]))},
        localsolver = "BOBYQA") #COBYLA unfortunately hangs from time to time.
        #so it cannot be used in optimization routine without risking hanging
        #despite the algorithm begin more robust and faster for this problem...
      ))
      # #With `alabama` package
      # opt <- try(suppressWarnings(
      #   alabama::auglag(par = c(median(dat), 0.1), fn = function(par){
      #     val <- -gevN.ll(par = c(par[1], zt, par[2]), dat = dat, q = q, qty = qty, N = N);
      #     ifelse(is.infinite(val) || is.na(val), 1e10, val)},
      #     gr = function(par){-gevN.score(par = c(par[1], zt, par[2]), dat = dat, q = q, qty = qty, N = N)[-2]},
      #     hin = function(par){
      #       sigma <- switch(qty, quantile = (zt-par[1])*par[2]/(N^par[2]*(log(1/q))^(-par[2])-1),
      #                       mean = (zt-par[1])*par[2]/(N^par[2]*gamma(1-par[2])-1))
      #       c(sigma, sigma + par[2]*(xmax-par[1]), sigma + par[2]*(xmin-par[1]))},
      #     control.outer = list (trace = FALSE)
      #   )))
      if(!is.character(opt)){
        # if(opt$convergence == 0  && !isTRUE(all.equal(opt$value, 1e10))){
        if(opt$convergence > 0 && !isTRUE(all.equal(opt$value, 1e10))){
        opt2 <- suppressWarnings(nloptr::slsqp(x0 = opt$par, fn = function(par){
          val <- -gevN.ll(par = c(par[1], zt, par[2]), dat = dat, q = q, qty = qty, N = N);
          ifelse(is.infinite(val) || is.na(val), 1e10, val)},
          gr = function(par){
            -gevN.score(par = c(par[1], zt, par[2]), dat = dat, q = q, qty = qty, N = N)[-2]},
          hin = function(par){
            sigma <- switch(qty, quantile = (zt-par[1])*par[2]/(N^par[2]*(log(1/q))^(-par[2])-1),
                            mean = (zt-par[1])*par[2]/(N^par[2]*gamma(1-par[2])-1))
            c(sigma, sigma + par[2]*(xmax-par[1]), sigma + par[2]*(xmin-par[1]))}
        ))
        if(!is.character(opt2)){
          if(opt2$convergence > 0 && !isTRUE(all.equal(opt2$value, 1e10))){
            return(c(opt2$par, opt2$value))
          }
        }
         return(c(opt$par, opt$value))
        }
      }
      return(rep(NA, 3))
    }

    #Missing psi vector
    if(missing(psi) || is.null(psi) || is.na(psi)){
      psirangelow <- pmax(1e-5, seq(ifelse(mle[3] < 0, -3.5,-2.5), -0.5, length = 6)*std.error + mle[2])
      lowvals <-  -sapply(psirangelow, constr.mle.N, dat = dat)[3,]-maxll
      psirangehigh <- seq(2.5, (1+mle[3])*8, length = 10)*std.error + mle[2]
      highvals <-  -sapply(psirangehigh, constr.mle.N, dat = dat)[3,]-maxll

      lo <- approx(x=lowvals, y=psirangelow, xout = -4)$y
      hi <- approx(x=highvals, y=psirangehigh, xout = -4)$y
      lo <- ifelse(is.na(lo), lm(psirangelow~lowvals)$coef[2]*-4 + mle[2], lo)
      hi <- ifelse(is.na(hi), lm(psirangehigh~highvals)$coef[2]*-4 + mle[2], hi)
      psi <- c(seq(lo, mle[2], length = 20)[-20], seq(mle[2], hi, length = 55)[-1])
    }

    lambda <- t(sapply(psi, constr.mle.N, dat = dat))
    pars <- cbind(lambda[,1, drop=FALSE], psi, lambda[,2, drop=FALSE])
    #Profile log likelihood values for psi
    profll <- -lambda[,3]
    if("tem" %in% mod){
      r <- sign(mle[param]-psi)*sqrt(2*(maxll - profll))
      phi.mle <- gevN.phi(par = mle, dat = dat, V = V, N = N, q = q, qty = qty)
      q2num <- apply(pars, 1, function(par){ -det(rbind(c(c(phi.mle) - gevN.phi(par = par, dat = dat, V = V, N = N, q = q, qty = qty)),
                                                       gevN.dphi(par = par, dat = dat, V = V, N = N, q = q, qty = qty)[-2,]))})
      if(isTRUE(any(sign(q2num)*sign(r) < 0, na.rm = TRUE))){
        warning("Correction factor and likelihood root are of opposite sign - check output")
      }

      logq <- apply(pars, 1, function(par){
        -0.5*log(det(gevN.infomat(par = par, dat = dat, method = "obs", N = N, q = q, qty = qty)[-2, -2]))}) +
        log(abs(q2num))
      qmlecontrib <-  -log(det(gevN.dphi(par = mle, dat = dat, V = V, N = N, q = q, qty = qty))) +
        0.5*log(det(gevN.infomat(par = mle, dat = dat, method = "obs", N = N, q = q, qty = qty)))
      logq <- logq + qmlecontrib
      qcor <- sign(q2num)*exp(logq)
      rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r))) / r)

      tem.max.opt <- function(psi, dat = dat){
        lam <- constr.mle.N(psi, dat = dat)
        para <- c(lam[1], psi, lam[2])
        pll <- -lam[3]
        rs <- 2*(maxll - pll)
        logq <- -0.5*log(det(gevN.infomat(par = para, dat = dat, method = "obs", N = N, q = q, qty = qty)[-2, -2])) + qmlecontrib +
          log(abs(det(rbind(c(c(phi.mle) - gevN.phi(par = para, dat = dat, V = V, q = q, N = N, qty = qty)),
                            gevN.dphi(par = para, dat = dat, V = V, q = q, N = N, qty = qty)[-2,]))))
        abs(rs + logq - log(sqrt(abs(rs))))
      }
      tem.max <- optim(par = mle[2], fn = tem.max.opt, method = "Brent", dat = dat,
                       lower = max(1e-5, mle[2] - std.error),
                       upper = mle[2] + std.error, control = list(abstol=1e-10))$par
       names(mle)[2] <- oldpar
    }
    if("modif" %in% mod){
      #Tangent exponential model approximation of Fraser and Reid to the profile likelihood
      tem.objfunc.N <- function(par, dat = dat, N = N, qty = qty, q = q){
        0.5*log(det(gevN.infomat(par = par, dat = dat, method = "obs", N = N, qty = qty, q = q)[-2,-2])) -
          log(det(gevN.dphi(par = par, dat = dat, N = N, V = V[, -2], qty = qty, q = q)[-2,]))
      }
      optim.tem.fn.N <- function(psi, dat = dat, q = q, qty = qty, N = N){
        theta.psi.opt <- constr.mle.N(psi, dat = dat)
        param <- c(theta.psi.opt[1], psi, theta.psi.opt[2])
        ll <- gevN.ll(param, dat = dat, N = N, q = q, qty = qty)
        ll + tem.objfunc.N(param, dat = dat, N = N, qty = qty, q = q)
      }
      #TEM profile log likelihood values for psi
      proflltem <- profll + suppressWarnings(apply(pars, 1, tem.objfunc.N, dat = dat, N = N, qty = qty, q = q))
      #Maximum objective function for TEM
      tem.mle.opt <- optim(par = mle[2], fn=optim.tem.fn.N, method = "Brent", dat = dat, q = q, qty = qty, N = N,
                           lower = max(min(dat),mle[2]-std.error), upper = mle[2] + std.error, control=list(fnscale=-1))
      #Severini empirical covariance function adjustment to profile likelihood
      #Score function - observation-wise for xi
      gevN.score.f <-  function(par, dat, N, q = q, qty = qty){
        qty <- match.arg(qty, c("mean", "quantile"))
        mu = par[1]; z = par[2]; xi = par[3];
        if(qty == "quantile"){ #quantiles at prob. q
          cbind(((-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*log(1/q)^(-xi) - 1)/(mu - z))/xi + ((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*log(1/q)^(-xi) - 1)/(mu - z))*(1/xi + 1)/((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1) - 1/(mu - z)),
                (-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)*(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)^2*xi) - (N^xi*log(1/q)^(-xi) - 1)*(dat - mu)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)^2) + 1/(mu - z)),
                ((-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi)*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(dat - mu)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)*xi) - log(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)/xi^2) - (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(dat - mu)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)) + (N^xi*log(1/q)^(-xi) - 1)*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - z)*xi/(N^xi*log(1/q)^(-xi) - 1)^2 - (mu - z)/(N^xi*log(1/q)^(-xi) - 1))/((mu - z)*xi) + log(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)/xi^2
                ))
        } else { #Mean
          cbind(((-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*gamma(-xi + 1) - 1)/(mu - z))/xi + ((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*gamma(-xi + 1) - 1)/(mu - z))*(1/xi + 1)/((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) - 1) - 1/(mu - z)),
                (-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)*(-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)^2*xi) - (N^xi*gamma(-xi + 1) - 1)*(dat - mu)*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)^2) + 1/(mu - z)),
                ((-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi)*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(dat - mu)/(((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)*xi) - log(-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)/xi^2) - (N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(dat - mu)*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)) + (N^xi*gamma(-xi + 1) - 1)*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - z)*xi/(N^xi*gamma(-xi + 1) - 1)^2 - (mu - z)/(N^xi*gamma(-xi + 1) - 1))/((mu - z)*xi) + log(-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)/xi^2))
        }
      }
    #Score at MLE (sums to zero)
    score.N.mle <- gevN.score.f(mle, dat, N, qty = qty, q = q)
    empcov.objfunc.N <- function(par, dat = dat, q = q, qty = qty, N = N){
      0.5*log(det(gevN.infomat(par = par, dat = dat, method = "obs", N = N, q = q, qty = qty)[-2,-2])) -
        log(sum(score.N.mle*gevN.score.f(par, dat, N = N, qty = qty, q = q)))
    }
    profllempcov <- profll + suppressWarnings(apply(pars, 1, empcov.objfunc.N, N = N, q = q, dat = dat, qty = qty))
    optim.empcov.fn.N <- function(psi, dat = dat, q = q, qty = qty, N = N){
      theta.psi.opt <- constr.mle.N(psi, dat = dat)
      param <- c(theta.psi.opt[1], psi, theta.psi.opt[2])
      ll <- gevN.ll(param, dat = dat, N = N, q = q, qty = qty)
      ll +  empcov.objfunc.N(param,dat = dat, q = q, qty = qty, N = N)
    }
    empcov.mle.opt <- optim(par=mle[2], fn=optim.empcov.fn.N, method = "Brent", dat = dat, qty = qty, q = q, N = N,
                            lower = max(min(dat), mle[2]-std.error), upper = mle[2] + std.error, control=list(fnscale=-1))
    }
  }
  #Return profile likelihood and quantities of interest (modified likelihoods)
  ans <- list(mle = mle, pars = pars, psi.max = as.vector(mle[oldpar]),
              param = oldpar, std.error = std.error, psi = psi, pll = profll, maxpll = maxll)
  if("tem" %in% mod){
    ans$r <- r
    ans$q <- qcor
    ans$rstar <- rstar
    ans$normal <- c(ans$psi.max, ans$std.error)
    if(correction && length(psi) > 10){
      ans <- spline.corr(ans)
    }
    ans$tem.psimax <- tem.max
  }
  if("modif" %in% mod){
    ans$tem.mle <- tem.mle.opt$par
    ans$tem.pll <- proflltem
    ans$tem.maxpll <- tem.mle.opt$value
    ans$empcov.mle <- empcov.mle.opt$par
    ans$empcov.pll <- profllempcov
    ans$empcov.maxpll <- empcov.mle.opt$value
  }

  if("tem" %in% mod){
    class(ans) <- c("extprof", "fr")
  } else{
    class(ans) <- "extprof"
  }
  ans
}

#' Modified profile likelihood for the generalized Pareto distribution
#'
#' This function calculates the (modified) profile likelihood based on the \eqn{p^*} formula.
#' There are two small-sample corrections that use a proxy for
#' \eqn{\ell_{\lambda; \hat{\lambda}}}{the sample space derivative of the nuisance},
#' which are based on Severini's (1999) empirical covariance
#' and the Fraser and Reid tangent exponential model approximation.
#' @details The two \code{mod} available are \code{tem}, the tangent exponential model (TEM) approximation and
#' \code{modif} for the penalized profile likelihood based on \eqn{p^*} approximation proposed by Severini.
#' For the latter, the penalization is based on the TEM or an empirical covariance adjustment term.
#'
#' @param psi parameter vector over which to profile (unidimensional)
#' @param param string indicating the parameter to profile over
#' @param mod string indicating the model. See \bold{Details}.
#' @param mle maximum likelihood estimate in \eqn{(\psi, \xi)} parametrization if \eqn{\psi \neq \xi} and \eqn{(\sigma, \xi)} otherwise (optional).
#' @param dat sample vector of excesses
#' @param m number of observations of interest for return levels. Required only for \code{args} values \code{"VaR"} or \code{"ES"}
#' @param N size of block over which to take maxima. Required only for \code{args} \code{Nmean} and \code{Nquant}.
#' @param p tail probability, equivalent to \eqn{1/m}. Required only for \code{args} \code{quant}.
#' @param q level of quantile for N-block maxima. Required only for \code{args} \code{Nquant}.
#' @param correction logical indicating whether to use \code{spline.corr} to smooth the tem approximation.
#' @param ... additional arguments such as output from call to \code{Vfun} if \code{mode="tem"}.
#'
#' @return a list with components
#' \itemize{
#' \item{\code{mle}:} maximum likelihood estimate
#' \item{\code{psi.max}:} maximum profile likelihood estimate
#' \item{\code{param}:} string indicating the parameter to profile over
#' \item{\code{std.error}:} standard error of \code{psi.max}
#' \item{\code{psi}:} vector of parameter \eqn{psi} given in \code{psi}
#' \item{\code{pll}:} values of the profile log likelihood at \code{psi}
#' \item{\code{maxpll}:} value of maximum profile log likelihood
#' }
#' In addition, if \code{mod} includes \code{tem}
#' \itemize{
#' \item{\code{normal}:}{maximum likelihood estimate and standard error of the interest parameter \eqn{psi}}
#' \item{\code{r}:}{values of likelihood root corresponding to \eqn{\psi}}
#' \item{\code{q}:}{vector of likelihood modifications}
#' \item{\code{rstar}:}{modified likelihood root vector}
#' \item{\code{rstar.old}:}{uncorrected modified likelihood root vector}
#' \item{\code{tem.psimax}:}{maximum of the tangent exponential model likelihood}
#' }
#' In addition, if \code{mod} includes \code{modif}
#' \itemize{
#' \item{\code{tem.mle}:} maximum of tangent exponential modified profile log likelihood
#' \item{\code{tem.profll}:} values of the modified profile log likelihood at \code{psi}
#' \item{\code{tem.maxpll}:} value of maximum modified profile log likelihood
#' \item{\code{empcov.mle}:} maximum of Severini's empirical covariance modified profile log likelihood
#' \item{\code{empcov.profll}:} values of the modified profile log likelihood at \code{psi}
#' \item{\code{empcov.maxpll}:} value of maximum modified profile log likelihood
#' }
#' @export
#' @examples
#' \dontrun{
#' dat <- evd::rgpd(n = 100, scale = 2, shape = 0.3)
#' gpd.pll(psi = seq(-0.5, 1, by=0.01), param = "shape", dat = dat)
#' gpd.pll(psi = seq(0.1, 5, by=0.1), param = "scale", dat = dat)
#' gpd.pll(psi = seq(20, 35, by=0.1), param = "quant", dat = dat, p = 0.01)
#' gpd.pll(psi = seq(20, 80, by=0.1), param = "ES", dat = dat, m = 100)
#' gpd.pll(psi = seq(15, 100, by=1), param = "Nmean", N = 100, dat = dat)
#' gpd.pll(psi = seq(15, 90, by=1), param = "Nquant", N = 100, dat = dat, q = 0.5)
#' }
gpd.pll <- function(psi, param = c("scale", "shape", "quant", "VaR", "ES", "Nmean", "Nquant"),
                    mod = c("tem","modif"), mle = NULL, dat, m = NULL, N = NULL, p = NULL, q = NULL, correction = TRUE, ...){
  param <- match.arg(param, c("scale", "shape", "quant", "VaR", "ES", "Nmean", "Nquant"))
  mod <- match.arg(mod, c("tem","modif"), several.ok = TRUE)
  #Arguments for parametrization of the log likelihood
  if(param == "shape"){
    args <- c("scale","shape")
  } else{
    args <- c(param, "shape")
  }
  #Sanity checks to ensure all arguments are provided
  if(missing(N)){
    if(param %in% c("Nmean","Nquant")){
      stop("Argument `N` missing. Procedure aborted")
    } else { N <- NA
    }
  }
  if(missing(m)){
    if(param %in% c("VaR","ES")){
      stop("Argument `m` missing. Procedure aborted")
    } else { m <- NA
    }
  }
  if(missing(p)){
    if(param == "quant"){
      stop("Argument `p` missing. Procedure aborted")
    } else { p <- NA
    }
  }
  if(missing(q)){
    if(param == "Nquant"){
      stop("Argument `q` missing. Procedure aborted")
    } else { q <- NA
    }
  }
  xmin <- min(dat); xmax <- max(dat)

  #If maximum likelihood estimates are not provided, find them
  if(is.null(mle)){
    mle <- gpd.mle(dat = dat, args = args, m = m, N = N, p = p, q = q)
  }

  #Extract the components, notably V for model `tem`.
  #Keep other components for optimization
  Vprovided <- FALSE
  extra.args <- list(...)
  if("V" %in% names(extra.args)){
    V <- extra.args$V
    extra.args$V <- NULL
    if(isTRUE(all.equal(dim(V), c(length(dat), 1)))){
      Vprovided <- TRUE
    }
  }
  if(!Vprovided){
    V <- switch(param,
                scale = gpd.Vfun(par = mle, dat = dat),
                shape = gpd.Vfun(par = mle, dat = dat),
                quant = gpdr.Vfun(par = mle, dat = dat,  m = 1/p),
                VaR = gpdr.Vfun(par = mle, dat = dat,  m = m),
                ES = gpde.Vfun(par = mle, dat = dat,  m = m),
                Nmean = gpdN.Vfun(par = mle, dat = dat,  N = N),
                Nquant = gpdr.Vfun(par = mle, dat = dat,  m = 1/(1-q^(1/N))))
  }
  #Obtained constrained maximum likelihood estimates for given value of psi
  if(param == "scale"){

    maxll <- gpd.ll(mle, dat = dat)
    std.error <- sqrt(solve(gpd.infomat(par = mle, dat = dat, method = "exp"))[1,1])
    constr.mle.scale <- function(sigmat){
      as.vector(evd::fpot(x = dat, threshold = 0,  model = "gpd", std.err = FALSE,
                          scale = sigmat, method = "Brent", lower = max(-1, -sigmat/xmax), upper = min(10, sigmat/xmin))$estimate)}

    #Missing psi vector
    if(missing(psi) || is.null(psi) || is.na(psi)){
      psirangelow <- pmax(1e-5, seq(-3, -1.5, length = 6)*std.error + mle[1])
      lowvals <-  sapply(psirangelow, function(par){gpd.ll(c(par,constr.mle.scale(par)), dat = dat)})-maxll
      psirangehigh <- seq(2.5, 4, length = 6)*std.error + mle[1]
      highvals <-  sapply(psirangehigh, function(par){gpd.ll(c(par,constr.mle.scale(par)), dat = dat)})-maxll

      lo <- approx(x=lowvals, y=psirangelow, xout = -4)$y
      hi <- approx(x=highvals, y=psirangehigh, xout = -4)$y
      lo <- ifelse(is.na(lo), lm(psirangelow~lowvals)$coef[2]*-4 + mle[1], lo)
      hi <- ifelse(is.na(hi), lm(psirangehigh~highvals)$coef[2]*-4 + mle[1], hi)
      psi <- seq(lo,hi, length = 55)
      }
    if(any(as.vector(psi) < 0)){
      warning("Negative scale values provided.");
      psi <- psi[psi>0]
      if(length(psi)==0){
       psi <- mle[1]
      }
    }

    pars <- cbind(psi, sapply(psi, constr.mle.scale))
    #Profile log likelihood values for psi
    profll <- apply(pars, 1, function(par){gpd.ll(par = par, dat = dat)})

      if("tem" %in% mod){
        r <- sign(mle[param]-psi)*sqrt(2*(maxll - profll))
        phi.mle <- gpd.phi(par = mle, dat = dat, V = V)
        q2num <- apply(pars, 1, function(par){ det(rbind(c(c(phi.mle) - gpd.phi(par = par, dat = dat, V = V)),
                                                         gpd.dphi(par = par, dat = dat, V = V)[-1,]))})
        if(isTRUE(any(sign(q2num)*sign(r) < 0, na.rm = TRUE))){
          warning("Correction factor and likelihood root are of opposite sign - check output")
        }
        logq <- apply(pars, 1, function(par){
          -0.5*log(gpd.infomat(par = par, dat = dat, method = "obs")[-1, -1])}) + log(abs(q2num))
        qmlecontrib <- -log(det(gpd.dphi(par = mle, dat = dat, V = V))) +
          0.5*log(det(gpd.infomat(par = mle, dat = dat, method = "obs")))
        logq <- logq + qmlecontrib
        qcor <- sign(q2num)*exp(logq)
        rstar <- ifelse(r==0, 0, r + (logq-log(abs(r)))/r)
        tem.max.opt <- function(psi, dat = dat){
          para <- c(psi, constr.mle.scale(psi))
          pll <- gpd.ll(par = para, dat = dat)
          rs <- 2*(maxll - pll)
          logq <- -0.5*log(gpd.infomat(par = para, dat = dat, method = "obs")[-1, -1]) + qmlecontrib +
        log(abs(det(rbind(c(c(phi.mle) - gpd.phi(par = para, dat = dat, V = V)),
                          gpd.dphi(par = para, dat = dat, V = V)[-1,]))))
          rs + logq - log(sqrt(abs(rs)))
        }
        tem.max <- optim(par=mle[1], fn=tem.max.opt, method = "Brent", dat = dat,
                             lower = max(1e-5,mle[1]-std.error), upper = mle[1] + std.error, control = list(abstol=1e-10))$par
      }

    if("modif" %in% mod){
      #Tangent exponential model approximation of Fraser and Reid to the profile likelihood
      tem.objfunc.scale <- function(par){
        0.5*log(gpd.infomat(par = par, dat = dat, method = "obs")[2,2]) -
          log(gpd.dphi(par = par, dat = dat, V = V[, 2, drop=FALSE])[2,1])}
      optim.tem.fn.scale <- function(psi){
        theta.psi.opt <- constr.mle.scale(psi)
        param <- c(psi, theta.psi.opt)          #ll <- -theta.psi.opt$deviance/2
        ll <- gpd.ll(param, dat = dat)
        ll + tem.objfunc.scale(param)
      }
      #TEM profile log likelihood values for psi
      proflltem <- profll + apply(pars, 1, tem.objfunc.scale)
      #Maximum objective function for TEM (line search in neighborhood of the MLE)

      tem.mle.opt <- optim(par=mle[1], fn=optim.tem.fn.scale, method = "Brent",
                           lower = max(1e-5,mle[1]-std.error), upper = mle[1] + std.error, control=list(fnscale=-1))
      tem.mle <- c(tem.mle.opt$par, constr.mle.scale(tem.mle.opt$par))

      #Severini empirical covariance function adjustment to profile likelihood
      #Score function - observation-wise
      gpd.score.f <-  function(par, dat){
        sigma = par[1]; xi = par[2]
        if(!isTRUE(all.equal(0,xi))){
          cbind(dat*(xi + 1)/(sigma^2*(dat*xi/sigma + 1)) - 1/sigma,
                -dat*(1/xi + 1)/(sigma*(dat*xi/sigma + 1)) + log(pmax(dat*xi/sigma + 1,0))/xi^2)
        } else {
          cbind((dat - sigma)/sigma^2, 1/2*(dat - 2*sigma)*dat/sigma^2)
        }
      }
      #Score at MLE (sums to zero)
      score.scale.mle <- gpd.score.f(mle, dat)[,2] #keep s_lambda
      empcov.objfunc.scale <- function(par){
        0.5*log(gpd.infomat(par = par, dat = dat, method = "obs")[2,2]) -
          log(sum(score.scale.mle*gpd.score.f(par, dat)[,2]))
      }
      profllempcov <- profll + apply(pars, 1, empcov.objfunc.scale)
      optim.empcov.fn.scale <- function(psi){
        theta.psi.opt <- constr.mle.scale(psi)
        param <- c(psi, theta.psi.opt)
        ll <- gpd.ll(param, dat = dat)
        ll +  empcov.objfunc.scale(param)
      }
      empcov.mle.opt <- optim(par=mle[1], fn=optim.empcov.fn.scale, method = "Brent",
                              lower = max(1e-5,mle[1]-std.error), upper = mle[1] + std.error, control=list(fnscale=-1))
      empcov.mle <- c(empcov.mle.opt$par, constr.mle.scale(empcov.mle.opt$par))
    }
    #Shape parameter
  } else if(param == "shape"){
    maxll <- gpd.ll(mle, dat = dat)
    std.error <- sqrt(solve(gpd.infomat(par = mle, dat = dat, method="exp"))[2,2])
    constr.mle.shape <- function(xit){
      as.vector(suppressWarnings(evd::fpot(x = dat, threshold = 0,  model = "gpd", std.err = FALSE,
                                           shape = xit, method = "Brent", lower = ifelse(xit<0, abs(xit)*xmax + 1e-5,1e-5), upper = 1e10)$estimate))
    }

    #Missing psi vector
    if(missing(psi) || is.null(psi) || is.na(psi)){
      psirangelow <- seq(ifelse(mle[2] < 0,-7,-5), -1.5, length = 10)*std.error + mle[2]
      lowvals <-  sapply(psirangelow, function(par){gpd.ll(c(constr.mle.shape(par), par), dat = dat)})-maxll
      psirangehigh <- seq(1.5, 10, length = 10)*std.error + mle[2]
      highvals <-  sapply(psirangehigh, function(par){gpd.ll(c(constr.mle.shape(par), par), dat = dat)})-maxll

      lo <- approx(x=lowvals, y=psirangelow, xout = -4)$y
      hi <- approx(x=highvals, y=psirangehigh, xout = -4)$y
      lo <- ifelse(is.na(lo), lm(psirangelow~lowvals)$coef[2]*-4 + mle[2], lo)
      hi <- ifelse(is.na(hi), lm(psirangehigh~highvals)$coef[2]*-4 + mle[2], hi)
      psi <- seq(lo, hi, length = 55)
    }

    pars <- cbind(sapply(psi, constr.mle.shape), psi)
    #Profile log likelihood values for psi
    profll <- apply(pars, 1, function(par){gpd.ll(par = par, dat = dat)})
    if("tem" %in% mod){
      r <- sign(mle[param]-psi)*sqrt(2*(maxll - profll))
      phi.mle <- gpd.phi(par = mle, dat = dat, V = V)
      q2num <- apply(pars, 1, function(par){ det(rbind(gpd.dphi(par = par, dat = dat, V = V)[-2,],
                                                       c(c(phi.mle) - gpd.phi(par = par, dat = dat, V = V))))})
      if(isTRUE(any(sign(q2num)*sign(r) < 0, na.rm = TRUE))){
        warning("Correction factor and likelihood root are of opposite sign - check output")
      }
      logq <- apply(pars, 1, function(par){
        -0.5*log(gpd.infomat(par = par, dat = dat, method = "obs")[-2, -2])}) + log(abs(q2num))
      qmlecontrib <-  - log(det(gpd.dphi(par = mle, dat = dat, V = V))) +
        0.5*log(det(gpd.infomat(par = mle, dat = dat, method = "obs")))
      logq <- logq + qmlecontrib
      qcor <- sign(q2num)*exp(logq)
      rstar <- ifelse(r == 0, 0, r + (logq-log(abs(r)))/r)
      #Maximum of TEM likelihood - indirect estimation via rstar
      tem.max.opt <- function(psi, dat = dat){
        para <- c(constr.mle.shape(psi), psi)
        pll <- gpd.ll(par = para, dat = dat)
        rs <- 2*(maxll - pll)
        logq <- -0.5*log(gpd.infomat(par = para, dat = dat, method = "obs")[-2, -2]) + qmlecontrib +
          log(abs(det(rbind(c(c(phi.mle) - gpd.phi(par = para, dat = dat, V = V)), gpd.dphi(par = para, dat = dat, V = V)[-2,]))))
        rs + logq - log(sqrt(abs(rs)))
      }
      tem.max <- optim(par=mle[2], fn=tem.max.opt, method = "Brent", dat = dat,
                       lower = mle[2]-std.error, upper = mle[2] + std.error, control = list(abstol=1e-10))$par

    }
    if("modif" %in% mod){
      #Tangent exponential model approximation of Fraser and Reid to the profile likelihood
      tem.objfunc.shape <- function(par){
        0.5*log(gpd.infomat(par = par, dat = dat, method = "obs")[1,1]) -
          log(gpd.dphi(par = par, dat = dat, V = V[, 1, drop=FALSE])[1,1])}
      optim.tem.fn.shape <- function(psi){
        theta.psi.opt <- constr.mle.shape(psi)
        param <- c(theta.psi.opt, psi)
        ll <- gpd.ll(param, dat = dat)
        ll + tem.objfunc.shape(param)
      }
      #TEM profile log likelihood values for psi
      proflltem <- profll + apply(pars, 1, tem.objfunc.shape)
      #Maximum objective function for TEM (line search in neighborhood of the MLE)
      tem.mle.opt <- optim(par = mle[2], fn = optim.tem.fn.shape, method = "Brent",
                           lower = mle[2]-std.error, upper = mle[2] + std.error, control = list(fnscale = -1))
      tem.mle <- c(constr.mle.shape(tem.mle.opt$par), tem.mle.opt$par)

      #Severini empirical covariance function adjustment to profile likelihood
      #Score function - observation-wise
      gpd.score.f <-  function(par, dat){
        sigma = par[1]; xi = par[2]
        if(!isTRUE(all.equal(0,xi))){
          cbind(dat*(xi + 1)/(sigma^2*(dat*xi/sigma + 1)) - 1/sigma,
                -dat*(1/xi + 1)/(sigma*(dat*xi/sigma + 1)) + log(pmax(dat*xi/sigma + 1,0))/xi^2)
        } else {
          cbind((dat - sigma)/sigma^2, 1/2*(dat - 2*sigma)*dat/sigma^2)
        }
      }
      #Score at MLE (sums to zero)
      score.shape.mle <- gpd.score.f(mle, dat)[,1] #keep s_lambda
      empcov.objfunc.shape <- function(par){
        0.5*log(gpd.infomat(par = par, dat = dat, method = "obs")[1,1]) -
          log(sum(score.shape.mle*gpd.score.f(par, dat)[,1]))
      }
      profllempcov <- profll + apply(pars, 1, empcov.objfunc.shape)
      optim.empcov.fn.shape <- function(psi){
        theta.psi.opt <- constr.mle.shape(psi)
        param <- c(theta.psi.opt, psi)
        ll <- gpd.ll(param, dat = dat)
        ll +  empcov.objfunc.shape(param)
      }
      empcov.mle.opt <- optim(par=mle[2], fn=optim.empcov.fn.shape, method = "Brent",
                              lower = mle[2]-std.error, upper = mle[2] + std.error, control=list(fnscale=-1))
      empcov.mle <- c(constr.mle.shape(empcov.mle.opt$par), empcov.mle.opt$par)
    }

    #Return levels, quantiles or value-at-risk
  } else if(param %in% c("quant", "VaR", "Nquant")) {
    if(param == "quant"){ m <- 1/p; }
    if(param == "Nquant"){ m <- 1/(1-q^(1/N)); }
    maxll <- gpdr.ll(mle, dat = dat, m = m)
    std.error <- sqrt(solve(gpdr.infomat(par = mle, dat = dat, method = "exp", m = m))[1,1])
    constr.mle.quant <- function(quant){
      nloptr::sbplx(x0=0.01, function(lambda, psi, m){
        -gpdr.ll(par=c(psi, lambda), dat = dat, m = m)}, psi = quant, m = m, lower = -1, upper=5)$par}


    #Missing psi vector
    if(missing(psi) || is.null(psi) || is.na(psi)){
      psirangelow <- unique(pmax(mean(dat), seq(-3, -0.5, length = 12)*std.error + mle[1]))
      lowvals <-  sapply(psirangelow, function(par){gpdr.ll(c(par, constr.mle.quant(par)), m = m, dat = dat)}) - maxll
      psirangehigh <- seq(2, 10, length = 12)*std.error + mle[1]
      highvals <-  sapply(psirangehigh, function(par){gpdr.ll(c(par, constr.mle.quant(par)), m = m, dat = dat)}) - maxll

      lo <- approx(x=lowvals, y=psirangelow, xout = -4)$y
      hi <- approx(x=highvals, y=psirangehigh, xout = -4)$y
      lo <- ifelse(is.na(lo), lm(psirangelow~lowvals)$coef[2]*-4 + mle[1], lo)
      hi <- ifelse(is.na(hi), lm(psirangehigh~highvals)$coef[2]*-4 + mle[1], hi)
      psi <- seq(lo,hi, length = 55)
    }


    pars <- cbind(psi, sapply(psi, constr.mle.quant))
    #Profile log likelihood values for psi
    profll <- apply(pars, 1, function(par){gpdr.ll(par = par, dat = dat, m = m)})
    if("tem" %in% mod){
      r <- sign(mle[param]-psi)*sqrt(2*(maxll - profll))
      phi.mle <- gpdr.phi(par = mle, dat = dat, m = m, V = V)
      q2num <- apply(pars, 1, function(par){ det(rbind(c(c(phi.mle) - gpdr.phi(par = par, dat = dat, V = V, m = m)),
                                                       gpdr.dphi(par = par, dat = dat, V = V, m = m)[-1,]))})
      if(isTRUE(any(sign(q2num)*sign(r) < 0, na.rm = TRUE))){
        warning("Correction factor and likelihood root are of opposite sign - check output")
      }
      logq <- apply(pars, 1, function(par){
        -0.5*log(gpdr.infomat(par = par, dat = dat, method = "obs", m = m)[-1, -1])}) + log(abs(q2num))
      qmlecontrib <- - log(det(gpdr.dphi(par = mle, dat = dat, V = V, m = m))) +
        0.5*log(det(gpdr.infomat(par = mle, dat = dat, method = "obs", m = m)))
      logq <- logq + qmlecontrib
      qcor <- sign(q2num)*exp(logq)
      rstar <- ifelse(r==0, 0, r + (logq-log(abs(r)))/r)

      tem.max.opt <- function(psi, dat = dat){
        para <- c(psi, constr.mle.quant(psi))
        pll <- gpdr.ll(par = para, dat = dat, m = m)
        rs <- 2*(maxll - pll)
        logq <- -0.5*log(gpdr.infomat(par = para, dat = dat, method = "obs", m = m)[-1, -1]) + qmlecontrib +
          log(abs(det(rbind(c(c(phi.mle) - gpdr.phi(par = para, dat = dat, V = V, m = m)), gpdr.dphi(par = para, dat = dat, V = V, m = m)[-1,]))))
        rs + logq - log(sqrt(abs(rs)))
      }
      tem.max <- optim(par=mle[1], fn=tem.max.opt, method = "Brent", dat = dat,
                       lower = max(1e-5,mle[1]-std.error), upper = mle[1] + std.error, control = list(abstol=1e-10))$par

    }
    if("modif" %in% mod){
      #Tangent exponential model approximation of Fraser and Reid to the profile likelihood
      tem.objfunc.quant <- function(par){
        0.5*log(gpdr.infomat(par = par, dat = dat, method = "obs", m = m)[2,2]) -
          log(gpdr.dphi(par = par, dat = dat, m = m, V = V[, 2, drop=FALSE])[2,1])
      }
      optim.tem.fn.quant <- function(psi){
        theta.psi.opt <- constr.mle.quant(psi)
        param <- c(psi, theta.psi.opt)
        ll <- gpdr.ll(param, dat = dat, m = m)
        ll + tem.objfunc.quant(param)
      }
      #TEM profile log likelihood values for psi
      proflltem <- profll + suppressWarnings(apply(pars, 1, tem.objfunc.quant))
      #Maximum objective function for TEM
      tem.mle.opt <- optim(par = mle[1], fn = optim.tem.fn.quant, method = "Brent",
                           lower = max(1e-5, mle[1] - std.error), upper = mle[1] + std.error, control = list(fnscale = -1))
      tem.mle <- c(tem.mle.opt$par, constr.mle.quant(tem.mle.opt$par))

      #Severini empirical covariance function adjustment to profile likelihood
      #Score function - observation-wise for xi
      gpdr.score.f <- function(par, dat, m){
        xi = par[2];  r= par[1]
        -dat*m^xi*(1/xi + 1)*log(m)/((dat*(m^xi - 1)/r + 1)*r) + (m^xi*r*xi*log(m)/(m^xi - 1)^2 - r/(m^xi - 1))*(m^xi - 1)/(r*xi) +
          log(dat*(m^xi - 1)/r + 1)/xi^2
      }
      #Score at MLE (sums to zero)
      score.quant.mle <- gpdr.score.f(mle, dat, m)
      empcov.objfunc.quant <- function(par){
        0.5*log(gpdr.infomat(par = par, dat = dat, method = "obs", m = m)[2,2]) - log(sum(score.quant.mle*gpdr.score.f(par, dat, m = m)))
      }
      profllempcov <- profll + suppressWarnings(apply(pars, 1, empcov.objfunc.quant))
      optim.empcov.fn.quant <- function(psi){
        theta.psi.opt <- constr.mle.quant(psi)
        param <- c(psi, theta.psi.opt)
        ll <- gpdr.ll(param, dat = dat, m = m)
        ll +  empcov.objfunc.quant(param)
      }
      empcov.mle.opt <- optim(par=mle[1], fn=optim.empcov.fn.quant, method = "Brent",
                              lower = max(1e-5,mle[1]-std.error), upper = mle[1] + std.error, control=list(fnscale=-1))
      empcov.mle <- c(empcov.mle.opt$par, constr.mle.quant(empcov.mle.opt$par))
    }
  } else if(param == "ES"){
    maxll <- gpde.ll(mle, dat = dat, m = m)
    std.error <- sqrt(solve(gpde.infomat(par=mle, dat = dat, method="exp", m = m))[1,1])
    constr.mle.es <- function(psif){
    # nloptr::auglag(x0 = 0.1, fn = function(x){-gpde.ll(par = c(psif, x), dat = dat, m = m)},
    #                            hin = function(x){ c(psif*(1-x)*((m^x-1)/x+1)^(-1),
    #                                                  1+psif*(1-x)*((m^x-1)/x+1)^(-1)/x*xmin,
    #                                                  1+psif*(1-x)*((m^x-1)/x+1)^(-1)/x*xmax)},
    #                lower = -1+1e-10, upper = 1-1e-5)$par}
      optim(par = 0.1, fn = function(x){-gpde.ll(par = c(psif, x), dat = dat, m = m)},
                       method="Brent", lower = -1+1e-10, upper = 1-1e-5)$par}

    #Missing psi vector
    if(missing(psi) || is.null(psi) || is.na(psi)){
      psirangelow <- unique(pmax(mean(dat), seq(-3, -0.1, length = 15)*std.error + mle[1]))
      lowvals <-  sapply(psirangelow, function(par){gpde.ll(c(par, constr.mle.es(par)), m = m, dat = dat)}) - maxll
      psirangehigh <- seq(2.5, 10, length = 12)*std.error + mle[1]
      highvals <-  sapply(psirangehigh, function(par){gpde.ll(c(par, constr.mle.es(par)), m = m, dat = dat)}) - maxll

      lo <- approx(x=lowvals, y=psirangelow, xout = -4)$y
      hi <- approx(x=highvals, y=psirangehigh, xout = -4)$y
      lo <- ifelse(is.na(lo), lm(psirangelow~lowvals)$coef[2]*-4 + mle[1], lo)
      hi <- ifelse(is.na(hi), lm(psirangehigh~highvals)$coef[2]*-4.5 + mle[1], hi)
      psi <- c(seq(lo,mle[1], length=15)[-15], seq(mle[1], min(mle[1]+2*std.error,hi), length=25)[-1])
      if(mle[1]+2*std.error < hi){
        psi <- c(psi, seq(mle[1]+2*std.error, hi, length = 20))
      }
    }

    pars <- cbind(psi, sapply(psi, constr.mle.es))
    #Profile log likelihood values for psi
    profll <- apply(pars, 1, function(par){gpde.ll(par = par, dat = dat, m = m)})
    profll[profll == 1e10] <- NA
    if("tem" %in% mod){
      r <- sign(mle[param]-psi)*sqrt(2*(maxll - profll))
      phi.mle <- gpde.phi(par = mle, dat = dat, m = m, V = V)
      q2num <- apply(pars, 1, function(par){ det(rbind(c(c(phi.mle) - gpde.phi(par = par, dat = dat, V = V, m = m)),
                                                       gpde.dphi(par = par, dat = dat, V = V, m = m)[-1,]))})
      if(isTRUE(any(sign(q2num)*sign(r) < 0, na.rm = TRUE))){
        warning("Correction factor and likelihood root are of opposite sign - check output")
      }
      logq <- apply(pars, 1, function(par){
        -0.5*log(gpde.infomat(par = par, dat = dat, method = "obs", m = m)[-1, -1])})  + log(abs(q2num))
      qmlecontrib <- - log(det(gpde.dphi(par = mle, dat = dat, V = V, m = m))) +
        0.5*log(det(gpde.infomat(par = mle, dat = dat, method = "obs", m = m)))
      logq <- logq + qmlecontrib
      qcor <- sign(q2num)*exp(logq)
      rstar <- ifelse(r == 0, 0, r + (logq-log(abs(r)))/r)

      tem.max.opt <- function(psi, dat = dat){
        para <- c(psi, constr.mle.es(psi))
        pll <- gpde.ll(par = para, dat = dat, m = m)
        rs <- 2*(maxll - pll)
        logq <- -0.5*log(gpde.infomat(par = para, dat = dat, method = "obs", m = m)[-1, -1]) + qmlecontrib +
          log(abs(det(rbind(c(c(phi.mle) - gpde.phi(par = para, dat = dat, V = V, m = m)), gpde.dphi(par = para, dat = dat, V = V, m = m)[-1,]))))
        rs + logq - log(sqrt(abs(rs)))
      }
      tem.max <- optim(par = mle[1], fn=tem.max.opt, method = "Brent", dat = dat,
                       lower = max(1e-5,mle[1]-std.error), upper = mle[1] + std.error, control = list(abstol=1e-10))$par

    }
    if("modif" %in% mod){
      #Tangent exponential model approximation of Fraser and Reid to the profile likelihood
      tem.objfunc.es <- function(par){
        0.5*log(gpde.infomat(par = par, dat = dat, method = "obs", m = m)[2,2]) -
          log(gpde.dphi(par = par, dat = dat, m = m, V = V[, 2, drop=FALSE])[2,1])
      }
      optim.tem.fn.es <- function(psi){
        theta.psi.opt <- constr.mle.es(psi)
        param <- c(psi, theta.psi.opt)
        ll <- gpde.ll(param, dat = dat, m = m)
        ll + tem.objfunc.es(param)
      }
      #TEM profile log likelihood values for psi
      proflltem <- profll + suppressWarnings(apply(pars, 1, tem.objfunc.es))
      #Maximum objective function for TEM

      tem.mle.opt <- optim(par=mle[1], fn=optim.tem.fn.es, method = "Brent",
                           lower = max(quantile(dat, 1-1/m), mle[1] - std.error), upper = mle[1] + std.error, control=list(fnscale=-1))
      tem.mle <- c(tem.mle.opt$par, constr.mle.es(tem.mle.opt$par))

      #Severini empirical covariance function adjustment to profile likelihood
      #Score function - observation-wise for xi
      gpde.score.f <-  function(par, dat, m){
        xi = par[2];  es = par[1]
        -((m^xi*log(m) + 1)*dat/(es*(xi - 1)) - dat*(m^xi + xi - 1)/(es*(xi - 1)^2))*(1/xi + 1)/(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1) +
          (m^xi + xi - 1)*((m^xi*log(m) + 1)*(xi - 1)*xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi + xi - 1))/((xi - 1)*xi) +
          log(pmax(0,-dat*(m^xi + xi - 1)/(es*(xi - 1)) + 1))/xi^2
      }
      #Score at MLE (sums to zero)
      score.es.mle <- gpde.score.f(mle, dat, m)
      empcov.objfunc.es <- function(par){
        0.5*log(gpde.infomat(par = par, dat = dat, method = "obs", m = m)[2,2]) - log(sum(score.es.mle*gpde.score.f(par, dat, m = m)))
      }
      profllempcov <- profll + suppressWarnings(apply(pars, 1, empcov.objfunc.es))
      optim.empcov.fn.es <- function(psi){
        theta.psi.opt <- constr.mle.es(psi)
        param <- c(psi, theta.psi.opt)
        ll <- gpde.ll(param, dat = dat, m = m)
        ll +  empcov.objfunc.es(param)
      }
      empcov.mle.opt <- optim(par=mle[1], fn=optim.empcov.fn.es, method = "Brent",
                              lower = max(quantile(dat, 1-1/m), mle[1]-std.error), upper = mle[1] + std.error, control=list(fnscale=-1))
      empcov.mle <- c(empcov.mle.opt$par, constr.mle.es(empcov.mle.opt$par))
    }
  } else if(param == "Nmean"){
    maxll <- gpdN.ll(mle, dat = dat, N = N)
    std.error <- sqrt(solve(gpdN.infomat(par=mle, dat = dat, method = "exp", N = N))[1])
    constr.mle.Nmean <- function(Nmeant){
      nloptr::sbplx(x0=0.01, function(lambda, psi, N){
        -gpdN.ll(par = c(psi, lambda), dat = dat, N = N)}, psi = Nmeant, N = N, lower = -1+1e-10, upper=1-1e-10)$par
    }

    #Missing psi vector
    if(missing(psi) || is.null(psi) || is.na(psi)){
      psirangelow <- unique(pmax(mean(dat), seq(-3, -0.25, length = 20)*std.error + mle[1]))
      lowvals <-  sapply(psirangelow, function(par){gpdN.ll(c(par,constr.mle.Nmean(par)), dat = dat, N=N)})-maxll
      psirangehigh <- seq(2.5, 8, length = 6)*std.error + mle[1]
      highvals <-  sapply(psirangehigh, function(par){gpdN.ll(c(par,constr.mle.Nmean(par)), dat = dat, N=N)})-maxll

      lo <- approx(x=lowvals, y=psirangelow, xout = -4)$y
      hi <- approx(x=highvals, y=psirangehigh, xout = -4)$y
      lo <- ifelse(is.na(lo), lm(psirangelow~lowvals)$coef[2]*-4 + mle[1], lo)
      hi <- ifelse(is.na(hi), lm(psirangehigh~highvals)$coef[2]*-4 + mle[1], hi)
      psi <- c(seq(lo,mle[1], length=15)[-15], seq(mle[1], min(mle[1]+2*std.error,hi), length=25)[-1])
      if(mle[1]+2*std.error < hi){
        psi <- c(psi, seq(mle[1]+2*std.error, hi, length = 20))
      }
    }
    if(any(as.vector(psi) < 0)){
      warning("Negative scale values provided.");
      psi <- psi[psi>0]
      if(length(psi)==0){
        psi <- mle[1]
      }
    }

    pars <- cbind(psi, sapply(psi, constr.mle.Nmean))
    #Profile log likelihood values for psi
    profll <- apply(pars, 1, function(par){gpdN.ll(par = par, dat = dat, N = N)})
    if("tem" %in% mod){
      r <- sign(mle[param]-psi)*sqrt(2*(maxll - profll))
      phi.mle <- gpdN.phi(par = mle, dat = dat, N = N, V = V)
      q2num <- apply(pars, 1, function(par){ det(rbind(c(c(phi.mle) - gpdN.phi(par = par, dat = dat, V = V, N = N)),
                                                       gpdN.dphi(par = par, dat = dat, V = V, N = N)[-1,]))})
      if(isTRUE(any(sign(q2num)*sign(r) < 0, na.rm = TRUE))){
        warning("Correction factor and likelihood root are of opposite sign - check output")
      }

      logq <- apply(pars, 1, function(par){
        -0.5*log(gpdN.infomat(par = par, dat = dat, method = "obs", N = N)[-1, -1])})  + log(abs(q2num))
      qmlecontrib <- - log(det(gpdN.dphi(par = mle, dat = dat, V = V, N = N))) +
        0.5*log(det(gpdN.infomat(par = mle, dat = dat, method = "obs", N = N)))
      logq <- logq + qmlecontrib
      qcor <- sign(q2num)*exp(logq)
      rstar <- ifelse(r == 0, 0, r + (logq-log(abs(r)))/r)

      tem.max.opt <- function(psi, dat = dat){
        para <- c(psi, constr.mle.Nmean(psi))
        pll <- gpdN.ll(par = para, dat = dat, N = N)
        rs <- 2*(maxll - pll)
        logq <- -0.5*log(gpdN.infomat(par = para, dat = dat, method = "obs", N = N)[-1, -1]) + qmlecontrib +
          log(abs(det(rbind(c(c(phi.mle) - gpdN.phi(par = para, dat = dat, V = V, N = N)),
                            gpdN.dphi(par = para, dat = dat, V = V, N = N)[-1,]))))
        rs + logq - log(sqrt(abs(rs)))
      }
      tem.max <- optim(par = mle[1], fn=tem.max.opt, method = "Brent", dat = dat,
                       lower = max(1e-5,mle[1]-std.error), upper = mle[1] + std.error,
                       control = list(abstol=1e-10))$par

    }
    if("modif" %in% mod){
      #Tangent exponential model approximation of Fraser and Reid to the profile likelihood
      tem.objfunc.Nmean <- function(par){
        0.5*log(gpdN.infomat(par = par, dat = dat, method = "obs", N = N)[2,2]) -
          log(gpdN.dphi(par = par, dat = dat, N = N, V = V[, 2, drop=FALSE])[2,1])
      }
      optim.tem.fn.Nmean <- function(psi){
        theta.psi.opt <- constr.mle.Nmean(psi)
        param <- c(psi, theta.psi.opt)
        ll <- gpdN.ll(param, dat = dat, N = N)
        ll + tem.objfunc.Nmean(param)
      }
      #TEM profile log likelihood values for psi
      proflltem <- profll + suppressWarnings(apply(pars, 1, tem.objfunc.Nmean))
      #Maximum objective function for TEM
      tem.mle.opt <- optim(par=mle[1], fn=optim.tem.fn.Nmean, method = "Brent",
                           lower = max(1e-5, mle[1]-std.error), upper = mle[1] + std.error, control=list(fnscale=-1))
      tem.mle <- c(tem.mle.opt$par, constr.mle.Nmean(tem.mle.opt$par))

      #Severini empirical covariance function adjustment to profile likelihood
      #Score function - observation-wise for xi
      gpdN.score.f <-  function(par, dat, N){
        z = par[1]; xi = par[2]
        cst <- exp(lgamma(N + 1) + lgamma(1-xi)-lgamma(N + 1-xi))
        -(psigamma(N - xi + 1)*cst - psigamma(-xi + 1)*cst)*dat*(1/xi + 1)/(z*(dat*(cst - 1)/z + 1)) +
          ((psigamma(N - xi + 1)*cst - psigamma(-xi + 1)*cst)*xi*z/(cst - 1)^2 - z/(cst - 1))*(cst - 1)/(xi*z) + log(dat*(cst - 1)/z + 1)/xi^2
      }
      #Score at MLE (sums to zero)
      score.Nmean.mle <- gpdN.score.f(mle, dat, N)
      empcov.objfunc.Nmean <- function(par){
        0.5*log(gpdN.infomat(par = par, dat = dat, method = "obs", N = N)[2,2]) - log(sum(score.Nmean.mle*gpdN.score.f(par, dat, N = N)))
      }
      profllempcov <- profll + suppressWarnings(apply(pars, 1, empcov.objfunc.Nmean))
      optim.empcov.fn.Nmean <- function(psi){
        theta.psi.opt <- constr.mle.Nmean(psi)
        param <- c(psi, theta.psi.opt)
        ll <- gpdN.ll(param, dat = dat, N = N)
        ll +  empcov.objfunc.Nmean(param)
      }
      empcov.mle.opt <- optim(par=mle[1], fn=optim.empcov.fn.Nmean, method = "Brent",
                              lower = max(1e-5,mle[1]-std.error), upper = mle[1] + std.error, control=list(fnscale=-1))
      empcov.mle <- c(empcov.mle.opt$par, constr.mle.Nmean(empcov.mle.opt$par))
    }
  }
  #Return profile likelihood and quantities of interest (modified likelihoods)

  ans <- list(mle = mle, pars = pars, psi.max = as.vector(mle[param]), param = param,
              std.error = std.error, psi = psi, pll = profll, maxpll = maxll)
  if("tem" %in% mod){
    ans$r <- r
    ans$q <- qcor
    ans$rstar <- rstar
    ans$normal <- c(ans$psi.max, ans$std.error)
    if(correction && length(psi) > 10){
      ans <- spline.corr(ans)
    }
    ans$tem.psimax <- as.vector(tem.max)
  }
  if("modif" %in% mod){
    ans$tem.mle <- ifelse(param == "shape", tem.mle[2], tem.mle[1])
    ans$tem.pll <- proflltem
    ans$tem.maxpll <- as.vector(tem.mle.opt$value)
    ans$empcov.mle <- ifelse(param == "shape", empcov.mle[2], empcov.mle[1])
    ans$empcov.pll <- as.vector(profllempcov)
    ans$empcov.maxpll <- as.vector(empcov.mle.opt$value)
  }
  if("tem" %in% mod){
    class(ans) <- c("extprof", "fr")
  } else{
    class(ans) <- "extprof"
  }
  ans
}



#' Plot of tangent exponential model profile likelihood
#'
#' This function is adapted from \code{\link[hoa]{plot.fr}}. It differs mostly in
#' the placement of legends.
#'
#' @param x an object of class \code{fr} returned by \code{\link{gpd.tem}} or \code{\link{gev.tem}}.
#' @param ... further arguments to \code{plot} currently ignored. Providing a numeric vector \code{which} allows for custom selection of the plots. A logical \code{all}. See \strong{Details}.
#' @return graphs depending on argument \code{which}
#' @details Plots produced depend on the integers provided in \code{which}. \code{1} displays the Wald pivot, the likelihood root \code{r}, the modified likelihood root \code{rstar} and the likelihood modification \code{q} as functions of the parameter \code{psi}. \code{2} gives the renormalized profile log likelihood and adjusted form, with the maximum likelihood having ordinate value of zero. \code{3} provides the significance function, a transformation of \code{1}. Lastly, \code{4} plots the correction factor as a function of the likelihood root; it is a diagnostic plot aimed for detecting failure of
#' the asymptotic approximation, often due to poor numerics in a neighborhood of \code{r=0}; the function should be smooth. The function \code{\link{spline.corr}} is designed to handle this by correcting numerically unstable estimates, replacing outliers and missing values with the fitted values from the fit.
#'
#'
#' @references Brazzale, A. R., Davison, A. C. and Reid, N. (2007). \emph{Applied Asymptotics: Case Studies in Small-Sample Statistics}. Cambridge University Press, Cambridge.
#' @export
plot.fr <- function(x, ...)
{ # plot a fraser-reid object
  whichPlot <- c(1:4) #default
  if(length(list(...))>0){
    if("which" %in% names(list(...))){
      whichPlot <- list(...)$which
      whichPlot <- (1:4)[c(1:4 %in% whichPlot)]
    } else if("all" %in% names(list(...))){
      if(!is.logical(all[1])){stop("Invalid `all' parameter")}
      if(list(...)$all){  whichPlot <- 1:4} else{ whichPlot <- 1:2}
    }
  }
  old.pars <- par(no.readonly = TRUE)
  if(sum(c(1,2,3,4) %in% whichPlot)>2){
    par(mfrow=c(2,2), mar=c(4.8,4.8,1,0.1))
  } else if(sum(c(1,2,3,4) %in% whichPlot)==2){
    par(mfrow=c(1,2))
  }

  fr <- x
  xl <- fr$param

  if(1 %in% whichPlot){
    plot(fr$psi,fr$r,type="l",xlab=xl,ylab="Value of pivot",ylim=c(-4, 4),
         panel.first=abline(h=qnorm(c(0.005,0.025,0.05,0.5,0.95,0.975,0.995)),col="grey",lwd=0.7), bty = "l")
    lines(fr$psi,(fr$normal[1]-fr$psi)/fr$normal[2],col="green", lwd=1.5)
    lines(fr$psi,fr$q,col="red", lwd=1.5)
    lines(fr$psi,fr$r, lwd=1.5)
    lines(fr$psi,fr$rstar,col="blue")
    legend(x="topright",c("Wald pivot","lik. root","modif. root", expression(q(psi))),
           lty=c(1,1,1,1), x.intersp = 0.2, lwd = 1.5, seg.len = 0.5,
           col=c("green","black","blue","red"), bty="n", cex=0.9, xjust = 1)
  }
  # top right: log likelihood (and adjusted version, I think?) as a function of psi
  if(2 %in% whichPlot){
    plot(fr$psi,-fr$r^2/2,type="l",xlab=xl,ylab="Profile log likelihood",ylim=c(-8, 0),
         panel.first=abline(h=-qchisq(c(0.95,0.99),df=1)/2,col="grey", lwd = 0.7),lwd=1.5, bty = "l")
    lines(fr$psi,-fr$rstar^2/2,col="blue",lwd=1.5)
    legend(x="bottomright",c(expression(paste("\u2113"["p"])),expression(paste("\u2113"["tem"]))),
           lty=c(1,1), x.intersp = 0.2, lwd = 1.5, seg.len = 0.5, col=c("black","blue"),bty="n")
    # optional: add diagnostic panels
  }
  if(3 %in% whichPlot){

    # lower left: plot of Phi(pivot) as a function of psi

    plot(fr$psi,pnorm(fr$r),type="l",xlab=xl,ylab="Significance function",ylim=c(0,1),
         panel.first=abline(h=c(0.025,0.05,0.5,0.95,0.975),col="grey",lwd=0.7), lwd = 1.5, bty = "l")
    lines(fr$psi,pnorm(fr$q),col="red", lwd = 1.5)
    lines(fr$psi,pnorm(fr$rstar),col="blue", lwd = 1.5)
    legend(x="topright",c("lik. root","modif. root",expression(q(psi))),
           lty=c(1,1,1),col=c("black","blue","red"),bty="n", x.intersp = 0.2, lwd = 1.5, seg.len = 0.5, cex = 0.9)
  }
  # lower right: log(q/r)/r as a function of r (should be smooth)
  if(4 %in% whichPlot){
    plot(fr$r,fr$rstar,type="l",xlab="Likelihood root r", ylab=expression(paste("Modified root r"^"*")),
         panel.first={ abline(h=0,col="grey", lwd = 0.7); abline(v=0,col="grey", lwd = 0.7)}, lwd=1.5, bty = "l")
  }


  par(old.pars)

}

#' Spline correction for Fraser-Reid approximations
#'
#' The tangent exponential model can be numerically unstable for values close to \eqn{r = 0}.
#' This function corrects these incorrect values, which are interpolated using splines.
#' The function takes as input an object of class \code{fr} and returns the same object with
#' different \code{rstar} values.
#' @section Warning:
#'
#' While penalized (robust) splines often do a good job at capturing and correcting for numerical outliers and \code{NA}, it
#' may also be driven by unusual features of the curve or fail at detecting outliers (or falsely identifying `correct' values as outliers). The user should always validate by comparing the plots of both the uncorrected (raw) output of the object with that of \code{spline.corr}.
#' @details If available, the function uses \code{cobs} from the eponym package. The latter handles constraints and smoothness penalties, and is more robust than the equivalent \code{\link[stats]{smooth.spline}}.
#'
#' @param fr an object of class \code{fr}, normally the output of \link{gpd.tem} or \link{gev.tem}.
#'
#' @return an object of class \code{fr}, containing as additional arguments \code{spline} and a modified \code{rstar} argument.
#' @export
#' @importFrom stats smooth.spline
spline.corr <- function(fr){
  #Step 1: fit a smoothing spline to rstar

  #If fit failed for some values (for example when shape forced to be < 1)
  #Remove those values
  fitfailed <- which(is.na(fr$r))
  if(length(fitfailed)>0){
    fr$r <- fr$r[-fitfailed]
    fr$rstar <- fr$rstar[-fitfailed]
    fr$q <- fr$q[-fitfailed]
    fr$psi <- fr$psi[-fitfailed]
  }
  w <- pchisq(fr$r^2, 0.5)
  #If any correction for q failed and returned NA
  corfailed <- which(is.na(fr$rstar))
  #If equispaced values for psi between MLE and other, then we have r = 0
  corfailed <- c(corfailed, which(fr$r == 0))
  if(length(corfailed)>0){
    resp <- (fr$rstar-fr$r)[-corfailed]
    regr <- fr$r[-corfailed]
    w <- w[-corfailed]
  } else{
    resp <-  (fr$rstar-fr$r)
    regr <- fr$r
  }
  if(requireNamespace("cobs", quietly = TRUE)){
    spline <- cobs::cobs(y = resp, x = regr, w = w, constraint = "none", lambda = 1, print.mesg = FALSE, print.warn = FALSE)$fitted
  } else{
    spline <- rev(stats::smooth.spline(y = resp, x = regr, w = w, spar = 0.9, cv = TRUE)$y)
  }
  #Compute difference between fitted values and rstar
  departure <- spline-resp
  #Ad-hoc fix of the values close to MLE where the numerical precision causes difficulty
  #Outlier detection via chi-square test
  #From package outliers, (c)Lukasz Komsta
  scores <- function (x, prob = NA){
    abs((x - mean(x))^2/var(x)) > qchisq(prob, 1)
  }
  bad <- which(scores(departure, prob = 0.95))

  if(length(bad)>0){
    #Exclude those values if they are in the end of the distribution
    bad <- bad[which(bad < 0.85*length(departure) && bad > 0.15*length(departure))]
    #Remove outliers and fit again (with less smoothness)
  }
  if(length(bad)>0){
    resp[bad] <- NA
    w <- w[-bad]
  }
  if(requireNamespace("cobs", quietly = TRUE)){
    spline <- cobs::cobs(y = resp, x = regr, constraint = "none", w = w, lambda = -1, ic = "SIC",
                         knots.add = TRUE, repeat.delete.add = TRUE, print.mesg = FALSE, print.warn = FALSE)
    fr$spline <- spline
    fr$rstar <- predict(spline, fr$r, interval = "none")[, 2]+fr$r
  } else{
    spline <- stats::smooth.spline(x = na.omit(cbind(regr, resp)), w = w, cv = FALSE, all.knots = TRUE)
    fr$spline <- spline
    fr$rstar <- predict(spline, fr$r)$y+fr$r
  }


  return(fr)
}

#' Confidence intervals for profile likelihood derived from TEM
#'
#' This function uses spline interpolation to derive \code{level} confidence intervals
#' using the output of either \link{gev.tem} or \link{gpd.tem}.
#'
#' @param object an object of class \code{fr}, normally the output of \link{gpd.tem} or \link{gev.tem}.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level, with default 0.95
#' @param ... additional arguments passed to functions. Providing a logical \code{warn = FALSE} turns off warning messages when the lower or upper confidence interval for \code{psi} are extrapolated beyond the provided calculations.
#' @return a 2 by 3 matrix containing point estimates, lower and upper confidence intervals based on \code{r} and \code{rstar}
#' @export
confint.fr <- function(object, parm, level = 0.95, ...){
  args <- list(...)
  if("warn" %in% names(args) && is.logical(args$warn)){
    warn <- args$warn
  } else{
    warn <- TRUE
  }
  #plot(object$psi~object$r)
  if(missing(parm)){
    ind  <- c(1, 2)
  } else if(is.numeric(parm)){
    ind <- c(1, 2)[c(1, 2) %in% parm]
  } else{
    ind <- c(1, 2)[c("r", "rstar") %in% parm]
  }
  if(length(ind) == 0){
    stop("Invalid `parm` argument.")
  }
  if(1 %in% ind){
    if(requireNamespace("cobs", quietly = TRUE)){
      fit.r <- cobs::cobs(x = object$r, y = object$psi, constraint = "none", lambda = 0, ic = "SIC", pointwise = cbind(0, 0, object$normal[1]),
                          knots.add = TRUE, repeat.delete.add = TRUE, print.mesg = FALSE, print.warn = FALSE)
      pr <- predict(fit.r, c(0, sqrt(qchisq(level, 1)), -sqrt(qchisq(level, 1))))[, 2]
    } else{
      fit.r <- stats::smooth.spline(x = na.omit(cbind(object$r, object$psi)), cv = FALSE)
      pr <- predict(fit.r, c(0, sqrt(qchisq(level, 1)), -sqrt(qchisq(level, 1))))$y
      pr[1] <- object$normal[1]
    }
    #lines(object$r, fit.r$fitted, col = 2)
    if(warn){
      if(!any(object$r >  sqrt(qchisq(level, 1)))){warning("Extrapolating the lower confidence interval for psi")}
      if(!any(object$r < -sqrt(qchisq(level, 1)))){warning("Extrapolating the upper confidence interval for psi")}
    }
  }
  if(2 %in% ind){
    #plot(object$psi~object$rstar)
    if(requireNamespace("cobs", quietly = TRUE)){
      fit.rst <- cobs::cobs(x = object$rstar, y = object$psi, constraint = "none", lambda = 0, ic = "SIC",
                            knots.add = TRUE, repeat.delete.add = TRUE, print.mesg = FALSE, print.warn = FALSE)
      prst <- predict(fit.rst, c(0, sqrt(qchisq(level, 1)), -sqrt(qchisq(level, 1))))[, 2]
    } else{
      fit.rst <- stats::smooth.spline(x = na.omit(cbind(object$rstar, object$psi)), cv = FALSE)
      prst <- predict(fit.rst, c(0, sqrt(qchisq(level, 1)), -sqrt(qchisq(level, 1))))$y
    }
    #lines(x = object$rstar, fit.rst$fitted, col = 2, pch = 19)
    if(warn){
      if(!any(object$r >  sqrt(qchisq(level, 1)))){warning("Extrapolating the adjusted lower confidence interval for psi.")}
      if(!any(object$r < -sqrt(qchisq(level, 1)))){warning("Extrapolating the adjusted upper confidence interval for psi")}
    }
  }

  if(all(c(1, 2) %in% ind)){
    conf <- cbind(pr, prst)
    colnames(conf) <- c("r", "rstar")
    rownames(conf) <- c("Estimate", "Lower CI", "Upper CI")
    return(conf)
  } else if(1 %in% ind){
    names(pr) <- c("Estimate", "Lower CI", "Upper CI")
    return(pr)
  } else if(2 %in% ind){
    names(prst) <- c("Estimate", "Lower CI", "Upper CI")
    return(prst)
  }

}



#' Tangent exponential model approximation for the GEV distribution
#'
#' The function \code{gev.tem} provides a tangent exponential model (TEM) approximation
#' for higher order likelihood inference for a scalar parameter for the generalized extreme value distribution.
#' Options include location scale and shape parameters as well as value-at-risk (or return levels).
#' The function attempts to find good values for \code{psi} that will
#' cover the range of options, but the fail may fit and return an error.
#'
#'
#' @param param parameter over which to profile
#' @param psi scalar or ordered vector of values for the interest parameter. If \code{NULL} (default), a grid of values centered at the MLE is selected
#' @param N size of block over which to take maxima. Required only for \code{param} \code{Nmean} and \code{Nquant}.
#' @param p tail probability for the (1-p)th quantile (return levels). Required only if \code{param = "retlev"}
#' @param q probability level of quantile. Required only for \code{param} \code{Nquant}.
#' @param dat sample vector for the GEV distribution
#' @param n.psi number of values of \code{psi} at which the likelihood is computed, if \code{psi} is not supplied (\code{NULL}). Odd values are more prone to give rise to numerical instabilities near the MLE. If \code{psi} is a vector of length 2 and \code{n.psi} is greater than 2, these are taken to be endpoints of the sequence.
#' @param plot logical indicating whether \code{plot.fr} should be called upon exit
#' @param correction logical indicating whether \link{spline.corr} should be called.
#' @author Leo Belzile, from code by A. Davison extracted from the \code{hoa} package bundle.
#' @importFrom ismev gev.fit
#' @return an invisible object of class \code{fr} (see \code{\link[hoa]{tem}}) with elements
#' \itemize{
#' \item{\code{normal}: }{maximum likelihood estimate and standard error of the interest parameter \eqn{psi}}
#' \item{\code{par.hat}: }{maximum likelihood estimates}
#' \item{\code{par.hat.se}: }{standard errors of maximum likelihood estimates}
#' \item{\code{th.rest}: }{estimated maximum profile likelihood at (\eqn{psi}, \eqn{\hat{\lambda}})}
#' \item{\code{r}: }{values of likelihood root corresponding to \eqn{\psi}}
#' \item{\code{psi}: }{vector of interest parameter}
#' \item{\code{q}: }{vector of likelihood modifications}
#' \item{\code{rstar}: }{modified likelihood root vector}
#' \item{\code{rstar.old}: }{uncorrected modified likelihood root vector}
#' \item{\code{param}: }{parameter}
#' }
#' @export
#' @examples
#' \dontrun{
#' dat <- evd::rgev(n = 40, loc = 0, scale = 2, shape = -0.1)
#' gev.tem("shape", dat = dat, plot = TRUE)
#' gev.tem("quant", dat = dat, p = 0.01, plot = TRUE)
#' gev.tem("scale", psi = seq(1, 4, by = 0.1), dat = dat, plot = TRUE)
#' dat <- evd::rgev(n = 40, loc = 0, scale = 2, shape = 0.2)
#' gev.tem("loc", dat = dat, plot = TRUE)
#' gev.tem("Nmean", dat = dat, p = 0.01, N=100, plot = TRUE)
#' gev.tem("Nquant", dat = dat, q = 0.5, N=100, plot = TRUE)
#' }
gev.tem <- function (param = c("loc", "scale", "shape", "quant", "Nmean", "Nquant"), dat, psi = NULL,
                     p = NULL, q = 0.5, N = NULL, n.psi = 50, plot = TRUE, correction = TRUE){
  if(param %in% c("VaR","retlev")){ param <- "quant"} #Compatibility following change of notation 13/07/2017
  if(param  ==  "scale" &&!is.null(psi)){
    if(isTRUE(any(psi<0))){
      stop("Invalid argument: scale parameter must be positive")
    }
  }
  tem_out <- gev.pll(psi = psi, param = param, mod = "tem", dat = dat, N = N, p = p, q = q, correction = correction)
  if(plot){
    plot.fr(tem_out)
  }
  return(invisible(tem_out))
}

#' Tangent exponential model approximation for the GP distribution
#'
#' The function \code{gpd.tem} provides a tangent exponential model (TEM) approximation
#' for higher order likelihood inference for a scalar parameter for the generalized Pareto distribution. Options include
#' scale and shape parameters as well as value-at-risk (also referred to as quantiles, or return levels)
#' and expected shortfall. The function attempts to find good values for \code{psi} that will
#' cover the range of options, but the fit may fail and return an error. In such cases, the user can try to find good
#' grid of starting values and provide them to the routine.
#'
#' As of version 1.11, this function is a wrapper around \code{gpd.pll}.
#'
#' @details The interpretation for \code{m} is as follows: if there are on average \eqn{m_y} observations per year above the threshold, then  \eqn{m = Tm_y} corresponds to \eqn{T}-year return level.
#'
#' @param param parameter over which to profile
#' @param threshold threshold value corresponding to the lower bound of the support or the location parameter of the generalized Pareto distribution.
#' @param psi scalar or ordered vector of values for the interest parameter. If \code{NULL} (default), a grid of values centered at the MLE is selected. If \code{psi} is of length 2 and \code{n.psi}>2, it is assumed to be the minimal and maximal values at which to evaluate the profile log likelihood.
#' @param m number of observations of interest for return levels. See \strong{Details}. Required only for \code{param = "VaR"} or \code{param = "ES"}.
#' @param N size of block over which to take maxima. Required only for \code{args} \code{Nmean} and \code{Nquant}.
#' @param p tail probability, equivalent to \eqn{1/m}. Required only for \code{args} \code{quant}.
#' @param q level of quantile for N-block maxima. Required only for \code{args} \code{Nquant}.
#' @param dat sample vector for the GP distribution
#' @param n.psi number of values of \code{psi} at which the likelihood is computed, if \code{psi} is not supplied (\code{NULL}). Odd values are more prone to give rise to numerical instabilities near the MLE
#' @param plot logical indiating whether \code{plot.fr} should be called upon exit
#' @param correction logical indicating whether \link{spline.corr} should be called.
#' @author Leo Belzile
#' @return an invisible object of class \code{fr} (see \code{\link[hoa]{tem}}) with elements
#' \itemize{
#' \item{\code{normal}: }{maximum likelihood estimate and standard error of the interest parameter \eqn{psi}}
#' \item{\code{par.hat}: }{maximum likelihood estimates}
#' \item{\code{par.hat.se}: }{standard errors of maximum likelihood estimates}
#' \item{\code{th.rest}: }{estimated maximum profile likelihood at (\eqn{psi}, \eqn{\hat{\lambda}})}
#' \item{\code{r}: }{values of likelihood root corresponding to \eqn{\psi}}
#' \item{\code{psi}: }{vector of interest parameter}
#' \item{\code{q}: }{vector of likelihood modifications}
#' \item{\code{rstar}: }{modified likelihood root vector}
#' \item{\code{rstar.old}: }{uncorrected modified likelihood root vector}
#' \item{\code{param}: }{parameter}
#' }
#' @export
#' @examples
#' set.seed(123)
#' dat <- evd::rgpd(n = 40, scale = 1, shape = -0.1)
#' #with plots
#' m1 <- gpd.tem(param = "shape", n.psi = 50, dat = dat, plot = TRUE)
#' m2 <- gpd.tem(param = "scale", n.psi = 50, dat = dat)
#' m3 <- gpd.tem(param = "VaR", n.psi = 50, dat = dat, m = 100)
#' #Providing psi
#' \dontrun{
#' psi <- c(seq(2, 5, length = 15), seq(5, 35, length = 45))
#' m4 <- gpd.tem(param = "ES", dat = dat, m = 100, psi = psi, correction = FALSE)
#' plot.fr(m4, which = c(2, 4))
#' plot(fr4 <- spline.corr(m4))
#' confint(m1)
#' confint(m4, parm = 2, warn = FALSE)
#' m5 <- gpd.tem(param = "Nmean", dat = dat, N = 100, psi = psi, correction = FALSE)
#' m6 <- gpd.tem(param = "Nquant", dat = dat, N = 100, q = 0.7, correction = FALSE)
#' }
gpd.tem <- function (dat, param = c("scale", "shape", "quant", "VaR", "ES", "Nmean", "Nquant"), psi = NULL,
                     m = NULL, threshold = 0, n.psi = 50, N = NULL, p = NULL,
                     q = NULL, plot = FALSE, correction = TRUE){
  if(param %in% c("VaR", "ES") && is.null(m)){
    stop("Parameter `m' missing")
  }
  dat <- dat - threshold
  tem <- gpd.pll(psi = psi, param = param, mod = "tem", dat = dat, N = N, m = m, mle = NULL,
                 q = q, p = p, correction = correction)
  if(plot){
    plot.fr(tem)
  }
  return(invisible(tem))
}
