#' Generalized Pareto distribution
#'
#' Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized Pareto distribution
#'
#' @author Leo Belzile
#' @name gpd
#' @param par vector of \code{scale} and \code{shape}
#' @param dat sample vector
#' @param tol numerical tolerance for the exponential model
#' @param method string indicating whether to use the expected  (\code{"exp"}) or the observed (\code{"obs"} - the default) information matrix.
#' @param V vector calculated by \code{gpd.Vfun}
#' @param n sample size
#' @section Usage: \preformatted{gpd.ll(par, dat, tol=1e-5)
#' gpd.ll.optim(par, dat, tol=1e-5)
#' gpd.score(par, dat)
#' gpd.infomat(par,dat,method = c("obs","exp"))
#' gpd.bias(par, n)
#' gpd.Fscore(par, dat, method=c("obs","exp"))
#' gpd.Vfun(par, dat)
#' gpd.phi(par, dat, V)
#' gpd.dphi(par, dat, V)}
#'
#' @section Functions:
#'
#' \itemize{
#' \item{\code{gpd.ll}:} {log-likelihood}
#' \item{\code{gpd.ll.optim}:} {negative log-likelihood parametrized in terms of \code{log(scale)} and shape
#' in order to perform unconstrained optimization}
#' \item{\code{gpd.score}:} {score vector}
#' \item{\code{gpd.infomat}:} {observed or expected information matrix}
#' \item{\code{gpd.bias}:} {Cox-Snell first order bias}
#' \item{\code{gpd.Fscore}:} {Firth's modified score equation}
#' \item{\code{gpd.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gpd.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gpd.dphi}:} {derivative matrix of the canonical parameter in the local
#' exponential family approximation}
#' }
#' @references Firth, D. (1993). Bias reduction of maximum likelihood estimates, \emph{Biometrika}, \strong{80}(1), 27--38.
#' @references Coles, S. (2001). \emph{An Introduction to Statistical Modeling of Extreme Values}, Springer, 209 p.
#' @references Cox, D. R. and E. J. Snell (1968). A general definition of residuals, \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \strong{30}, 248--275.
#' @references Cordeiro, G. M. and R. Klein (1994). Bias correction in ARMA models, \emph{Statistics and Probability Letters}, \strong{19}(3), 169--176.
#' @references Giles, D. E., Feng, H. and R. T. Godwin (2016).  Bias-corrected maximum likelihood estimation of the  parameters of the generalized Pareto distribution, \emph{Communications in Statistics - Theory and Methods}, \strong{45}(8), 2465--2483.
#'
NULL



#' @title Generalized extreme value distribution
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized extreme value distribution
#'
#' @name gev
#' @param par vector of \code{loc}, \code{scale} and \code{shape}
#' @param dat sample vector
#' @param method string indicating whether to use the expected  (\code{"exp"}) or the observed (\code{"obs"} - the default) information matrix.
#' @param V vector calculated by \code{gev.Vfun}
#' @param n sample size
#' @param p vector of probabilities
#' @section Usage: \preformatted{gev.ll(par, dat)
#' gev.ll.optim(par, dat)
#' gev.score(par, dat)
#' gev.infomat(par, dat, method = c("obs","exp"))
#' gev.retlev(par, p)
#' gev.bias(par, n)
#' gev.Fscore(par, dat, method=c("obs","exp"))
#' gev.Vfun(par, dat)
#' gev.phi(par, dat, V)
#' gev.dphi(par, dat, V)}
#'
#' @section Note:
#' The Gumbel case is not currently handled.
#'
#' @section Functions:
#' \itemize{
#' \item{\code{gev.ll}:} {log-likelihood}
#' \item{\code{gev.ll.optim}:} {negative log-likelihood parametrized in terms of location, \code{log(scale)} and shape
#' in order to perform unconstrained optimization}
#' \item{\code{gev.score}:} {score vector}
#' \item{\code{gev.infomat}:} {observed or expected information matrix}
#' \item{\code{gev.retlev}:} {return level, corresponding to the \eqn{(1-p)}th quantile}
#' \item{\code{gev.bias}:} {Cox-Snell first order bias}
#' \item{\code{gev.Fscore}:} {Firth's modified score equation}
#' \item{\code{gev.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gev.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gev.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
#' @references Firth, D. (1993). Bias reduction of maximum likelihood estimates, \emph{Biometrika}, \strong{80}(1), 27--38.
#' @references Coles, S. (2001). \emph{An Introduction to Statistical Modeling of Extreme Values}, Springer, 209 p.
#' @references Cox, D. R. and E. J. Snell (1968). A general definition of residuals, \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \strong{30}, 248--275.
#' @references Cordeiro, G. M. and R. Klein (1994). Bias correction in ARMA models, \emph{Statistics and Probability Letters}, \strong{19}(3), 169--176.
NULL


#' Log-likelihood for the generalized Pareto distribution
#'
#' Function returning the density of an \code{n} sample from the GP distribution.
#' \code{gpd.ll.optim} returnsthe negative log-likelihood parametrized in terms of \code{log(scale)} and shape
#' in order to perform unconstrained optimization
#' @seealso \code{\link{gpd}}
#' @inheritParams gpd
#' @export
#' @keywords internal
gpd.ll <- function(par, dat, tol=1e-5){
  sigma = par[1]; xi = par[2]
  if(abs(xi)>tol){
    -length(dat)*log(sigma)-(1+1/xi)*sum(log(pmax(1+xi/sigma*dat,0)))
  } else{
    -length(dat)*log(sigma)-sum(dat)/sigma
  }
}



#' @rdname gpd.ll
#' @inheritParams gpd
#' @export
gpd.ll.optim <- function(par, dat,tol=1e-5){
  sigma = exp(par[1]); xi = par[2]
  if(abs(xi)>tol){
    length(dat)*log(sigma)+(1+1/xi)*sum(log(pmax(1+xi/sigma*dat,0)))
  } else{
    length(dat)*log(sigma)+sum(dat)/sigma
  }
}
#' Score vector for the generalized Pareto distribution
#'
#' @seealso \code{\link{gpd}}
#' @inheritParams gpd
#' @export
#' @keywords internal
gpd.score <- function(par, dat){
  sigma = par[1]; xi = par[2]
  c(sum(dat*xi*(1/xi + 1)/(sigma^2*(dat*xi/sigma + 1)) - 1/sigma),
    sum(-dat*(1/xi + 1)/(sigma*(dat*xi/sigma + 1)) + log(pmax(dat*xi/sigma + 1,0))/xi^2))
}
#????
# Jf.gpd.grad <- function(par){
#   xi =par[2]; sigma = par[1];
#   c(-1/(sigma^2*sqrt(2*xi + 1)*(xi + 1)),
#   -1/(sigma*sqrt(2*xi + 1)*(xi + 1)^2) - 1/(sigma*(2*xi + 1)^(3/2)*(xi + 1))
# )}


#' Information matrix for the generalized Pareto distribution
#'
#' The function returns the expected or observed information matrix.
#' @seealso \code{\link{gpd}}
#' @inheritParams gpd
#' @export
#' @keywords internal
gpd.infomat <- function(par,dat,method = c("obs","exp")){
  if(missing(method)){method = "obs"}
  sigma = par[1]; xi = par[2]
  if(method=="obs"){
    -cbind(c( (length(dat)-(1+xi)*sum(dat*(2*sigma+xi*dat)/(sigma+xi*dat)^2))/(sigma^2),
              (1+xi)*sum(dat/(sigma+xi*dat)^2)/xi-sum(dat/(sigma+xi*dat))/(sigma*xi)),
           c((1+xi)*sum(dat/(sigma+xi*dat)^2)/xi-sum(dat/(sigma+xi*dat))/(sigma*xi),
             (2/xi^2)*sum(dat/(sigma+xi*dat))-2*sum(log(1+xi*dat/sigma))/(xi^3)+(1+1/xi)*sum(dat^2/(sigma+xi*dat)^2)
           ))
  } else if(method=="exp"){
    k11 = -2/((1+xi)*(1+2*xi))
    k22 = -1/(sigma^2*(1+2*xi))
    k12 = -1/(sigma*(1+xi)*(1+2*xi))
    -cbind(c(k11,k12),c(k12,k22))
  }
}
#' Tangent exponential model statistics for the generalized Pareto distribution
#'
#' Matrix of approximate ancillary statistics, sample space derivative of the
#' log-likelihood and mixed derivative for the generalized Pareto distribution.
#' @seealso \code{\link{gpd}}
#' @inheritParams gpd
#' @export
#' @name gpd.temstat
#' @keywords internal
gpd.Vfun <- function(par, dat){
  sigma = par[1]; xi = par[2]
  cbind(dat/sigma,
        sigma*(dat*xi/sigma + 1)^(-1/xi)*(log(dat*xi/sigma + 1)/xi^2 - dat/(sigma*(dat*xi/sigma + 1)*xi))/(dat*xi/sigma + 1)^(-1/xi - 1)
  )
}

#' @inheritParams gpd
#' @rdname gpd.temstat
#' @export
gpd.phi <- function(par, dat, V){
  sigma = par[1]; xi = par[2]
  rbind(-xi*(1/xi + 1)/(sigma*(dat*xi/sigma + 1))	)%*%V
}

#' @inheritParams gpd
#' @export
#' @rdname gpd.temstat
gpd.dphi <- function(par, dat, V){
  sigma = par[1]; xi = par[2]
  rbind(xi*(1/xi + 1)/(sigma^2*(dat*xi/sigma + 1)) - dat*xi^2*(1/xi + 1)/(sigma^3*(dat*xi/sigma + 1)^2),
        -(1/xi + 1)/(sigma*(dat*xi/sigma + 1)) + dat*xi*(1/xi + 1)/(sigma^2*(dat*xi/sigma + 1)^2) + 1/(sigma*(dat*xi/sigma + 1)*xi))%*%V
}


#' Log-likelihood for the generalized extreme value distribution
#'
#' Function returning the density of an \code{n} sample from the GEV distribution.
#'
#' \code{gev.ll.optim} returns the negative log-likelihood parametrized in terms of location, \code{log(scale)} and shape in order to perform unconstrained optimization
#'
#' @inheritParams gev
#' @export
#' @keywords internal
#' @seealso \code{\link{gev}}
gev.ll <- function(par, dat){
  dat <- as.vector(dat)

  tx <- pmax(1+par[3]/par[2]*(dat-par[1]),0)
  sum(-log(par[2])-(1/par[3]+1)*log(tx)-tx^(-1/par[3]))
}

#' @rdname gev.ll
#' @inheritParams gev
#' @export
#' @keywords internal
gev.ll.optim <- function(par, dat){
  tpar = par; tpar[2] = exp(par[2])
  #parameters with log-scale
  nll = -gev.ll(tpar,dat)
  return(nll)
}

#' Score vector for the generalized extreme value distribution
#'
#' @inheritParams gev
#' @export
#' @keywords internal
gev.score <- function(par, dat){
  mu = par[1]; sigma = par[2]; xi = par[3]
  c(sum(-(-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)/sigma - xi*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1))),
    sum(-(dat - mu)*((dat - mu)*xi/sigma + 1)^(-1/xi - 1)/sigma^2 + (dat - mu)*xi*(1/xi + 1)/(sigma^2*((dat - mu)*xi/sigma + 1)) - 1/sigma),
    sum(-(mu - dat)*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) - (log(-(mu - dat)*xi/sigma + 1)/xi^2 - (mu - dat)/(sigma*((mu - dat)*xi/sigma - 1)*xi))/(-(mu - dat)*xi/sigma + 1)^(1/xi) + log(-(mu - dat)*xi/sigma + 1)/xi^2)
  )
}

#' Information matrix for the generalized extreme value distribution
#'
#' The function returns the expected or observed information matrix.
#' @inheritParams gev
#' @export
#' @keywords internal
gev.infomat <- function(par, dat, method=c("obs","exp")){
  method <- match.arg(method,c("obs","exp"))
  if(missing(method)){
    method="obs"
  }
  n <- length(dat)
  if(method=="exp"){
    #(Expected) Fisher information does not depend on location parameter
    if(length(par)==3){sigma = par[2]; xi = par[3];
    } else{	sigma=par[1]; xi=par[2];
    }
    p = (1+xi)^2*gamma(1+2*xi)
    q = gamma(2+xi)*(digamma(1+xi)+(1+xi)/xi)
    return(n*cbind(
      c(p/sigma^2, -(p-gamma(2+xi))/(sigma^2*xi), (p/xi-q)/(sigma*xi)),
      c(-(p-gamma(2+xi))/(sigma^2*xi), (1-2*gamma(2+xi)+p)/(sigma^2*xi^2), -(1+digamma(1)+(1-gamma(2+xi))/xi-q+p/xi)/(sigma*xi^2)),
      c((p/xi-q)/(sigma*xi), -(1+digamma(1)+(1-gamma(2+xi))/xi-q+p/xi)/(sigma*xi^2), (pi^2/6+(1+digamma(1)+1/xi)^2-2*q/xi+p/xi^2)/xi^2)
    ))
  } else if(method=="obs"){
    mu = par[1]; sigma = par[2]; xi = par[3];
    infomat <- matrix(0, ncol=3,nrow=3)
    infomat[1,1] <- sum(-((dat - mu)*xi/sigma + 1)^(-1/xi - 2)*xi*(1/xi + 1)/sigma^2 + xi^2*(1/xi + 1)/(sigma^2*((dat - mu)*xi/sigma + 1)^2))
    infomat[1,2] <- infomat[2,1] <- sum(-(dat - mu)*((dat - mu)*xi/sigma + 1)^(-1/xi - 2)*xi*(1/xi + 1)/sigma^3 + ((dat - mu)*xi/sigma + 1)^(-1/xi - 1)/sigma^2 - xi*(1/xi + 1)/(sigma^2*((dat - mu)*xi/sigma + 1)) + (dat - mu)*xi^2*(1/xi + 1)/(sigma^3*((dat - mu)*xi/sigma + 1)^2))
    infomat[1,3] <- infomat[3,1] <- sum((-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)*((mu - dat)*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) - log(-(mu - dat)*xi/sigma + 1)/xi^2)/sigma - (1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) + (mu - dat)*xi*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)^2) + 1/(sigma*((mu - dat)*xi/sigma - 1)*xi))
    infomat[2,3] <- infomat[3,2] <- sum((mu - dat)*(-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)*(log(-(mu - dat)*xi/sigma + 1)/xi^2 - (mu - dat)/(sigma*((mu - dat)*xi/sigma - 1)*xi))/sigma^2 + (mu - dat)*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)) - (mu - dat)^2*xi*(1/xi + 1)/(sigma^3*((mu - dat)*xi/sigma - 1)^2) - (mu - dat)/(sigma^2*((mu - dat)*xi/sigma - 1)*xi) + (mu - dat)^2/(sigma^3*(-(mu - dat)*xi/sigma + 1)^(1/xi)*((mu - dat)*xi/sigma - 1)^2))
    infomat[3,3] <- sum(-(log(-(mu - dat)*xi/sigma + 1)/xi^2 - (mu - dat)/(sigma*((mu - dat)*xi/sigma - 1)*xi))^2/(-(mu - dat)*xi/sigma + 1)^(1/xi) + (2*log(-(mu - dat)*xi/sigma + 1)/xi^3 - 2*(mu - dat)/(sigma*((mu - dat)*xi/sigma - 1)*xi^2) - (mu - dat)^2/(sigma^2*((mu - dat)*xi/sigma - 1)^2*xi))/(-(mu - dat)*xi/sigma + 1)^(1/xi) + (mu - dat)^2*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)^2) - 2*log(-(mu - dat)*xi/sigma + 1)/xi^3 + 2*(mu - dat)/(sigma*((mu - dat)*xi/sigma - 1)*xi^2))
    infomat[2,2] <- sum(-(mu - dat)^2*(-(mu - dat)*xi/sigma + 1)^(-1/xi - 2)*xi*(1/xi + 1)/sigma^4 - 2*(mu - dat)*(-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)/sigma^3 - 2*(mu - dat)*xi*(1/xi + 1)/(sigma^3*((mu - dat)*xi/sigma - 1)) + (mu - dat)^2*xi^2*(1/xi + 1)/(sigma^4*((mu - dat)*xi/sigma - 1)^2) + 1/sigma^2)

    return(-infomat)
  }
}


#' Tangent exponential model statistics for the generalized extreme value distribution
#'
#' Matrix of approximate ancillary statistics, sample space derivative of the
#' log-likelihood and mixed derivative for the generalized extreme value distribution.
#' @seealso \code{\link{gev}}
#' @inheritParams gev
#' @export
#' @name gev.temstat
#' @keywords internal
gev.Vfun <- function(par, dat){
  cbind(1,
        (dat-par[1])/par[2],
        par[2]*(-(par[1] - dat)*par[3]/par[2] + 1)^(-1/par[3])*(log(-(par[1] - dat)*par[3]/par[2] + 1)/par[3]^2 - (par[1] - dat)/(par[2]*((par[1] - dat)*par[3]/par[2] - 1)*par[3]))/(-(par[1] - dat)*par[3]/par[2] + 1)^(-1/par[3] - 1))
}

#' @rdname gev.temstat
#' @inheritParams gev
#' @export
#' @keywords internal
gev.phi <- function(par, dat, V){
  mu = par[1]; sigma = par[2]; xi = par[3]
  t(((dat - mu)*xi/sigma + 1)^(-1/xi - 1)/sigma + xi*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)))%*%V
}


#' @rdname gev.temstat
#' @inheritParams gev
#' @export
#' @keywords internal
gev.dphi <- function(par, dat, V){
  mu = par[1]; sigma = par[2]; xi = par[3]
  rbind((-(mu - dat)*xi/sigma + 1)^(-1/xi - 2)*xi*(1/xi + 1)/sigma^2 - xi^2*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)^2),
        -(mu - dat)*(-(mu - dat)*xi/sigma + 1)^(-1/xi - 2)*xi*(1/xi + 1)/sigma^3 - (-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)/sigma^2 - xi*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)) + (mu - dat)*xi^2*(1/xi + 1)/(sigma^3*((mu - dat)*xi/sigma - 1)^2),
        -(-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)*((mu - dat)*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) - log(-(mu - dat)*xi/sigma + 1)/xi^2)/sigma + (1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) - (mu - dat)*xi*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)^2) - 1/(sigma*((mu - dat)*xi/sigma - 1)*xi))%*%V
}






#' Cox-Snell first order bias expression for the GEV distribution
#'
#' Bias vector for the GEV distribution based on an \code{n} sample.
#' @inheritParams gev
#' @export
#' @keywords internal
#' @seealso \code{\link{gev}}
gev.bias <- function(par, n){
	if(length(n)>1){stop("Invalid argument for sample size")}
	sigma <- par[2] ; xi <- par[3]
	if(xi < -1/3){stop("Cox-Snell correction only valid if the shape is greater than -1/3")}
	zeta3 = 1.20205690315959428539973816151144999076498629234049888179227155
	k111 <- ((1 + xi)^2*(1 + 4*xi)*gamma(1 + 3*xi))/sigma^3
	k112 <- (xi+1)*(gamma(2*xi+2)-(xi+1)*(4*xi+1)*gamma(3*xi+1))/(sigma^3*xi)
	k113 <- (1 + xi)*((1 + xi)*(1 + 4*xi)*gamma(1+3*xi)- gamma(1+2*xi)*(1 + 2*xi*(2 + xi) + xi*(1 + 2*xi)*psigamma(2+2*xi,0)))/(sigma^2*xi^2)
	k122 <- ((1 - xi)*gamma(2 + xi) - gamma(3 + 2*xi) + (1 + xi)^2*(1 + 4*xi)*gamma(1 + 3*xi))/(sigma^3*xi^2)
	k123 <- ((-gamma(2 + xi))*(1 + 2*xi + xi*psigamma(2+xi,0)) + (1 + xi)*((-(1 + 5*xi + 4*xi^2))*gamma(1 + 3*xi) +  gamma(1 + 2*xi)*(2 + 7*xi + 3*xi^2 + xi*(1 + 2*xi)* psigamma(2+2*xi,0))))/(sigma^2*xi^3)
	k133 <- ((1 + xi)^2*(1+6*xi+8*xi^2)*gamma(1+3*xi)- gamma(3 + 2*xi)*(1 + 5*xi + 3*xi^2 +  xi*(1 + 2*xi)*psigamma(2+2*xi,0)) + (1 + 2*xi)*gamma(1+xi)*(1 + 6*xi + 5*xi^2 + 2*xi^3 + 2*xi*(1 + 3*xi + 2* xi^2)*psigamma(2+xi,0) + xi^2*(1 + xi)*psigamma(2+xi,0)^2+ (xi^2)*(1 + xi)*psigamma(2+xi,1)))/(sigma*xi^4*(1 + 2*xi))
	k222 <- (1 - 3*xi + 3*(xi-1)*gamma(2+xi)+ 1.5*gamma(3 + 2*xi)-(1 + xi)^2*(1 + 4*xi)*gamma(1+3*xi))/(sigma^3*xi^3)
	k223 <- -(1 + 2*xi +digamma(1)*xi - xi^2 - digamma(1)*xi^2 - 27*xi^3*gamma(3*xi)-12*xi^4*gamma(3*xi)-3*gamma(2+xi) - 5*xi*gamma(2+xi) + 3*(1 + xi)*gamma(1+2*xi)+ 10*xi*(1 + xi)*gamma(1+2*xi)+4*xi^2*(1 + xi)*gamma(1+2*xi)-gamma(1+3*xi)-6*xi*gamma(1+3*xi)- 2*xi*gamma(2+xi)*psigamma(2+xi,0)+xi*(1 + xi)*(1 + 2*xi)*gamma(1+2*xi)*psigamma(2+2*xi,0))/(sigma^2*xi^4)
	k233 <- (1 + 7*xi + 2*digamma(1)*xi + 4*xi^2 + 6*digamma(1)*xi^2 + digamma(1)^2*xi^2 + (pi^2*xi^2)/6 -3*xi^3*(9 + 4*xi)*gamma(3*xi) + 3*gamma(1+2*xi) + 17*xi*gamma(1+2*xi)+22*xi^2*gamma(1+2*xi)+8*xi^3*gamma(1 + 2*xi)-(1+6*xi)*gamma(1+3*xi)+(1+3*xi+2*xi^2)*2*xi*gamma(1+2*xi)*psigamma(2+2*xi,0)-gamma(1+xi)*(3 + 16*xi + 13*xi^2 + 4*xi^3 + 2*xi*(2 + 5*xi + 3*xi^2)*psigamma(2+xi,0)+xi^2*(1 + xi)*psigamma(2+xi,0)^2 + xi^2*(1 + xi)*psigamma(2+xi,1)))/(sigma*xi^5)
	k333 <- (-3*gamma(3+2*xi)*(1+6*xi+4*xi^2+xi*(1+2*xi)*psigamma(2+2*xi,0))+6*(1+2*xi)*gamma(1+xi)*(1+8*xi+7*xi^2+4*xi^3+2*xi*(1+4*xi+3* xi^2)*psigamma(2+xi,0)+xi^2*(1+xi)*psigamma(2+xi,0)^2+xi^2*(1+xi)*psigamma(2+xi,1))+(1+2*xi)*(-2-6*(4+digamma(1))*xi-(30+42*digamma(1)+6*digamma(1)^2+pi^2)*xi^2+2*(1+xi)^2*(1+4*xi)*gamma(1+3*xi)+xi^3*(-8-18*digamma(1)^2-2*digamma(1)^3-3*pi^2-digamma(1)*(24+pi^2)+4*zeta3)))/(xi^6*(2+4*xi))

	k11.2 <- -2*(xi + 1)^2*gamma(2*xi + 1)/sigma^3
	k11.3 <- 2*(xi + 1)^2*psigamma(2*xi + 1)*gamma(2*xi + 1)/sigma^2 + 2*(xi + 1)*gamma(2*xi + 1)/sigma^2
	k12.2 <- 2*((xi + 1)^2*gamma(2*xi + 1) - gamma(xi + 2))/(sigma^3*xi)
	k12.3 <- -(2*(xi + 1)^2*psigamma(2*xi + 1)*gamma(2*xi + 1) + 2*(xi + 1)*gamma(2*xi + 1) - psigamma(xi + 2)*gamma(xi + 2))/(sigma^2*xi) + ((xi + 1)^2*gamma(2*xi + 1) - gamma(xi + 2))/(sigma^2*xi^2)
	k13.2 <- -(((xi + 1)^2*gamma(2*xi + 1) - xi*((xi + 1)/xi + psigamma(xi + 1))*gamma(xi + 2))/(sigma^2*xi^2))
	k13.3 <- -(-(2*(xi + 1)^2*psigamma(2*xi + 1)*gamma(2*xi + 1) - xi*((xi + 1)/xi + psigamma(xi + 1))*psigamma(xi + 2)*gamma(xi + 2) + xi*((xi + 1)/xi^2 - 1/xi - psigamma(xi + 1,1))*gamma(xi + 2) + 2*(xi + 1)*gamma(2*xi + 1) - ((xi + 1)/xi + psigamma(xi + 1))*gamma(xi + 2))/(sigma*xi^2) + 2*((xi + 1)^2*gamma(2*xi + 1) - xi*((xi + 1)/xi + psigamma(xi + 1))*gamma(xi + 2))/(sigma*xi^3))
	k22.2 <- -2*((xi + 1)^2*gamma(2*xi + 1) - 2*gamma(xi + 2) + 1)/(sigma^3*xi^2)
	k22.3 <- 2*((xi + 1)^2*psigamma(2*xi + 1)*gamma(2*xi + 1) + (xi + 1)*gamma(2*xi + 1) - psigamma(xi + 2)*gamma(xi + 2))/(sigma^2*xi^2) - 2*((xi + 1)^2*gamma(2*xi + 1) - 2*gamma(xi + 2) + 1)/(sigma^2*xi^3)
	k23.2 <- -(-((digamma(1)) + (xi + 1)^2*gamma(2*xi + 1)/xi - ((xi + 1)/xi + psigamma(xi + 1))*gamma(xi + 2) - (gamma(xi + 2) - 1)/xi + 1)/(sigma^2*xi^2))
	k23.3 <- -((2*(xi + 1)^2*psigamma(2*xi + 1)*gamma(2*xi + 1)/xi - ((xi + 1)/xi + psigamma(xi + 1))*psigamma(xi + 2)*gamma(xi + 2) + ((xi + 1)/xi^2 - 1/xi - psigamma(xi + 1,1))*gamma(xi + 2) - (xi + 1)^2*gamma(2*xi + 1)/xi^2 + 2*(xi + 1)*gamma(2*xi + 1)/xi - psigamma(xi + 2)*gamma(xi + 2)/xi + (gamma(xi + 2) - 1)/xi^2)/(sigma*xi^2) - 2*((digamma(1)) + (xi + 1)^2*gamma(2*xi + 1)/xi - ((xi + 1)/xi + psigamma(xi + 1))*gamma(xi + 2) - (gamma(xi + 2) - 1)/xi + 1)/(sigma*xi^3))
	k33.2 <- 0
	k33.3 <- 2*((xi + 1)^2*psigamma(2*xi + 1)*gamma(2*xi + 1)/xi^2 - ((xi + 1)/xi + psigamma(xi + 1))*psigamma(xi + 2)*gamma(xi + 2)/xi + ((xi + 1)/xi^2 - 1/xi - psigamma(xi + 1,1))*gamma(xi + 2)/xi - (xi + 1)^2*gamma(2*xi + 1)/xi^3 + (xi + 1)*gamma(2*xi + 1)/xi^2 + ((xi + 1)/xi + psigamma(xi + 1))*gamma(xi + 2)/xi^2 - (digamma(1) + 1/xi + 1)/xi^2)/xi^2 - 1/3*(pi^2 + 6*(digamma(1) + 1/xi + 1)^2 + 6*(xi + 1)^2*gamma(2*xi + 1)/xi^2 - 12*((xi + 1)/xi + psigamma(xi + 1))*gamma(xi + 2)/xi)/xi^3

	#Derivatives of information matrix
	A1 <- 0.5*cbind(c(k111,k112,k113),c(k112,k122,k123),c(k113,k123,k133))
	A2 <- -cbind(c(k11.2,k12.2,k13.2),c(k12.2,k22.2,k23.2),c(k13.2,k23.2,k33.2))+0.5*cbind(c(k112,k122,k123),c(k122,k222,k223),c(k123,k223,k233))
	A3 <- -cbind(c(k11.3,k12.3,k13.3),c(k12.3,k22.3,k23.3),c(k13.3,k23.3,k33.3))+0.5*cbind(c(k113,k123,k133),c(k123,k223,k233),c(k133,k233,k333))

	#Information matrix
	p = (1+xi)^2*gamma(1+2*xi)
	q = gamma(2+xi)*(digamma(1+xi)+(1+xi)/xi)
	infomat <- cbind(
		c(p/sigma^2, -(p-gamma(2+xi))/(sigma^2*xi), (p/xi-q)/(sigma*xi)),
		c(-(p-gamma(2+xi))/(sigma^2*xi), (1-2*gamma(2+xi)+p)/(sigma^2*xi^2),
			-(1+digamma(1)+(1-gamma(2+xi))/xi-q+p/xi)/(sigma*xi^2)),
		c((p/xi-q)/(sigma*xi), -(1+digamma(1)+(1-gamma(2+xi))/xi-q+p/xi)/(sigma*xi^2),
			(pi^2/6+(1+digamma(1)+1/xi)^2-2*q/xi+p/xi^2)/xi^2)
	)

	infoinv <- solve(infomat)

	return(infoinv%*%cbind(A1,A2,A3)%*%c(infoinv)/n)
}

#
# GPDbias <- function(par=c(2,0.1), n=1){
#   sigma=par[1]; xi=par[2]
#
#   k111 = 24/((1+xi)*(1+2*xi)*(1+3*xi))
#   k112 = 8/(sigma*(1+xi)*(1+2*xi)*(1+3*xi))
#   k122 = 4/(sigma^2*(1+2*xi)*(1+3*xi))
#   k222 = 4/(sigma^3*(1+3*xi))
#   k11d1 = 2*(3+4*xi)/((1+xi)^2*(1+2*xi)^2)
#   k11d2 = 0
#   k22d1 = 2/(sigma^2*(1+2*xi)^2)
#   k22d2 = 2/(sigma^3*(1+2*xi))
#   k12d1 = (3+4*xi)/(sigma*(1+xi)^2*(1+2*xi)^2)
#   k12d2 = 1/(sigma^2*(1+xi)*(1+2*xi))
#   k11 = -2/((1+xi)*(1+2*xi))
#   k22 = -1/(sigma^2*(1+2*xi))
#   k12 = -1/(sigma*(1+xi)*(1+2*xi))
#   infomat.gpd <- -cbind(c(k11,k12),c(k12,k22))
#   A1g <- cbind(c(k11d1,k12d1),c(k12d1,k22d1))-0.5*cbind(c(k111,k112),c(k112,k122))
#   A2g <- cbind(c(k11d2,k12d2),c(k12d2,k22d2))-0.5*cbind(c(k112,k122),c(k122,k222))
#   invinfomat.gpd <- function(par,n){
#     sigma = par[1]; xi = par[2]
#     cbind(c((1+xi)^2, - sigma*(1+xi)),c(-sigma*(1+xi),2*(1+xi)*sigma^2))
#   }
#   c(invinfomat.gpd(par,1)%*%cbind(A1g, A2g)%*%c(invinfomat.gpd(par,1))/n)
# }
#   invinfomat.gpd <- function(par,n){
#     sigma = par[1]; xi = par[2]
#     cbind(c((1+xi)^2, - sigma*(1+xi)),c(-sigma*(1+xi),2*(1+xi)*sigma^2))
#   }


#' Cox-Snell first order bias expression for the generalized Pareto distribution
#'
#' Bias vector for the GP distribution based on an \code{n} sample.
#' @inheritParams gpd
#' @references Coles, S. (2001). \emph{An Introduction to Statistical Modeling of Extreme Values}, Springer, 209 p.
#'@references Cox, D. R. and E. J. Snell (1968). A general definition of residuals, \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \strong{30}, 248--275.
#' @references Cordeiro, G. M. and R. Klein (1994). Bias correction in ARMA models, \emph{Statistics and Probability Letters}, \strong{19}(3), 169--176.
#' @references Giles, D. E., Feng, H. and R. T. Godwin (2016).  Bias-corrected maximum likelihood estimation of the  parameters of the generalized Pareto distribution, \emph{Communications in Statistics - Theory and Methods}, \strong{45}(8), 2465--2483.
#' @export
#' @keywords internal
#' @seealso \code{\link{gpd}}, \code{\link{gpd.bcor}}
gpd.bias <- function(par, n){ #scale, shape
	if(length(par)!=2){stop("Invalid input for correction")}
	if(length(n)>1){stop("Invalid argument for sample size")}
	if(3*par[2]< -1){stop("Invalid bias correction for GPD; need shape > -1/3");}
	c(par[1]*(3+5*par[2]+4*par[2]^2), -(1+par[2])*(3+par[2]))/(n+3*n*par[2])
}

#'  Firth's modified score equation for the generalized Pareto distribution
#'
#' @inheritParams gpd
#' @references Firth, D. (1993). Bias reduction of maximum likelihood estimates, \emph{Biometrika}, \strong{80}(1), 27--38.
#' @seealso \code{\link{gpd}}, \code{\link{gpd.bcor}}
#' @export
#' @keywords internal
gpd.Fscore <- function(par, dat, method=c("obs","exp")){
	if(missing(method) || method!="exp"){	method="obs"}
	gpd.score(par, dat) - gpd.infomat(par, dat, method)%*%gpd.bias(par, length(dat))
}

#'  Firth's modified score equation for the generalized extreme value distribution
#'
#' @seealso \code{\link{gev}}
#' @inheritParams gev
#' @param method string indicating whether to use the expected ("exp") or the observed ("obs" - the default) information matrix.
#' @references Firth, D. (1993). Bias reduction of maximum likelihood estimates, \emph{Biometrika}, \strong{80}(1), 27--38.
#' @export
#' @keywords internal
gev.Fscore <- function(par, dat, method="obs"){
	if(missing(method) || method!="exp"){	method="obs"}
	gev.score(par, dat) - gev.infomat(par, dat, method)%*%gev.bias(par, length(dat))
}

#' Bias correction using Firth's modified score function or bias substraction
#'
#' The routine uses the MLE (bias-corrected) as starting values and proceeds
#' to find the solution using a root finding algorithm.
#' Since the bias-correction is not valid for \eqn{xi < -1/3}, any solution that is unbounded
#' will return a vector of \code{NA} - additionally, passing a \code{par} argument with shape less than -1/3
#' will return an error if \code{method="subtract"} is selected, as the bias correction does not exist then.
#' For small samples, expected and observed information can return very different estimates.
#' @importFrom nleqslv nleqslv
#' @importFrom rootSolve multiroot
#' @param par parameter vector (\code{scale}, \code{shape})
#' @param dat sample of observations
#' @param corr string indicating which correction to employ either \code{subtract} or \code{firth}
#' @param method string indicating whether to use the expected  (\code{"exp"}) or the observed (\code{"obs"} --- the default) information matrix. Used only if \code{corr="firth"}
#' @return vector of bias-corrected parameters
#' @export
#' @examples
#' set.seed(1)
#' dat <- evd::rgpd(n=40, scale=1, shape=-0.2)
#' par <- mev::gp.fit(dat, threshold=0, show=FALSE)$estimate
#' gpd.bcor(par,dat, "subtract")
#' gpd.bcor(par,dat, "firth") #observed information
#' gpd.bcor(par,dat, "firth","exp")
gpd.bcor <- function(par, dat, corr=c("subtract","firth"), method=c("obs","exp")){
	corr <- match.arg(corr, c("subtract","firth"))
#Basic bias correction - substract bias at MLE parbc=par-bias(par)
#bcor1 <- function(par, dat){ par-gpd.bias(par,length(dat))}
#Other bias correction - find bias corrected that solves implicit eqn parbc=par-bias(parbc)
bcor <-  function(par, dat){
		if(par[2]< -1/3){
			warning("Invalid bias correction for GPD; need shape > -1/3")
			return(rep(NA,2))
		} else{
		bcor.rootfind <- nleqslv::nleqslv(x=par, fn=function(parbc, par, dat){
			parbc-par+gpd.bias(parbc, length(dat))}, par=par, dat=dat)
		if(bcor.rootfind$termcd == 1 || (bcor.rootfind$termcd==2 && isTRUE(all.equal(bcor.rootfind$fvec,c(0,0),tolerance=1.5e-8)))){
			return(bcor.rootfind$x)
		} else{
			return(rep(NA,2))
		}
	}
}
#Firth correction, shifted shape so as to take advantage of the positive argument of multiroot
#It is easier to check the output for failed convergence as it is known not be valid there.
bcorF <-  function(par, dat, method=c("obs","exp")){
	# if(! "rootSolve" %in% installed.packages()[,1]){ #TODO check this for package
	# 	stop("Please install package `rootSolve' to use Firth correction")
	# } else{
		method <- match.arg(method, c("obs","exp"))
		#Score function, location transformed so that failed fit lies on boundary without errors
		#should in principle return the same numerical estimate
		gpd.Fscore.plus <- function(par, dat,method=c("obs","exp")){
			parcopy <- c(par[1], par[2]-0.3)
			if(missing(method)){method="obs"}
			gpd.score(parcopy, dat) - gpd.infomat(parcopy, dat, method)%*%gpd.bias(parcopy, length(dat))
		}

		firthplus = try(multiroot(gpd.Fscore.plus,start=par+c(0,0.3),
															dat=dat, method=method, positive=TRUE), silent=TRUE)
		if(is.character(firthplus) ||
			 any(c(isTRUE(all.equal(firthplus$root[2],target=0, check.names=FALSE)),
					isTRUE(all.equal(firthplus$root[1],target=0, check.names=FALSE)),
			 			is.nan(c(firthplus$f.root)),
					firthplus$root[2]>2,
					!isTRUE(all.equal(c(firthplus$f.root),target=rep(0,2), check.names=FALSE))
					))){
			#Can fail if sigma+xi*x_max < 0 - error message
			#Or can reach the boundary and not be able to evaluate the root
			firthplus <- rep(NA,2)
		} else{
			firthplus <- firthplus$root-c(0,1/3)
		}
		return(firthplus)

		# st <- get(ifelse(par[2] > -1/3, "parbc", "par0"))
		# firth = try(rootSolve::multiroot(gpd.Fscore, start=st,dat=dat, method=method)$root, silent=TRUE)
		# if(is.character(firth)){
		#   firth <- rep(NA,2)
		# }
		# return(firth)
	# }
}
#Return values
if(corr=="subtract"){
	return(bcor(par=par, dat=dat))
	}
if(corr=="firth"){
	return(bcorF(par=par, dat=dat, method=method))
	}
}
