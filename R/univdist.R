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
#' gpd.infomat(par, dat, method = c("obs","exp"))
#' gpd.bias(par, n)
#' gpd.Fscore(par, dat, method = c("obs","exp"))
#' gpd.Vfun(par, dat)
#' gpd.phi(par, dat, V)
#' gpd.dphi(par, dat, V)}
#'
#' @section Functions:
#'
#' \itemize{
#' \item{\code{gpd.ll}:} {log likelihood}
#' \item{\code{gpd.ll.optim}:} {negative log likelihood parametrized in terms of \code{log(scale)} and shape
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
#'
#' @section Functions:
#' \itemize{
#' \item{\code{gev.ll}:} {log likelihood}
#' \item{\code{gev.ll.optim}:} {negative log likelihood parametrized in terms of location, \code{log(scale)} and shape
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


#' log likelihood for the generalized Pareto distribution
#'
#' Function returning the density of an \code{n} sample from the GP distribution.
#' \code{gpd.ll.optim} returns the negative log likelihood parametrized in terms of \code{log(scale)} and shape
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
  sigma = par[1]; xi = as.vector(par[2])
  if(!isTRUE(all.equal(0,xi))){
    c(sum(dat*(xi + 1)/(sigma^2*(dat*xi/sigma + 1)) - 1/sigma),
      sum(-dat*(1/xi + 1)/(sigma*(dat*xi/sigma + 1)) + log(pmax(dat*xi/sigma + 1,0))/xi^2))
  } else {
    c(sum((dat - sigma)/sigma^2), sum(1/2*(dat - 2*sigma)*dat/sigma^2))
  }
}

#' Information matrix for the generalized Pareto distribution
#'
#' The function returns the expected or observed information matrix.
#' @seealso \code{\link{gpd}}
#' @inheritParams gpd
#' @param nobs number of observations
#' @export
#' @keywords internal
gpd.infomat <- function(par, dat, method = c("obs", "exp"), nobs = length(dat)){
  if(missing(method)){method <- "obs"}
  sigma <- as.vector(par[1]); xi <- as.vector(par[2])
  if(method=="obs"){
    if(any((1 + xi*dat/sigma)<0)){
      stop("Data outside of range specified by parameter, yielding a zero likelihood")
    }
    c11 <- (length(dat)-(1+xi)*sum(dat*(2*sigma+xi*dat)/(sigma+xi*dat)^2))/(sigma^2)
    if(!isTRUE(all.equal(xi,0))){
      c12 <- (1+xi)*sum(dat/(sigma+xi*dat)^2)/xi-sum(dat/(sigma+xi*dat))/(sigma*xi)
      c22 <- (2/xi^2)*sum(dat/(sigma+xi*dat))-2*sum(log(pmax(dat*xi/sigma + 1,0)))/(xi^3)+(1+1/xi)*sum(dat^2/(sigma+xi*dat)^2)
    } else{
      c12 <- sum(-(dat^2 - dat*sigma)/sigma^3)
      c22 <- sum(-1/3*(2*dat - 3*sigma)*dat^2/sigma^3)
    }
    -matrix(c(c11,c12,c12,c22), nrow=2,ncol=2,byrow=TRUE)
  } else if(method=="exp"){
    k22 <- -2/((1+xi)*(1+2*xi))
    k11 <- -1/(sigma^2*(1+2*xi))
    k12 <- -1/(sigma*(1+xi)*(1+2*xi))
    -nobs*cbind(c(k11,k12),c(k12,k22)) #fixed 12-10-2016
  }
}
#' Tangent exponential model statistics for the generalized Pareto distribution
#'
#' Matrix of approximate ancillary statistics, sample space derivative of the
#' log likelihood and mixed derivative for the generalized Pareto distribution.
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
  rbind(-(xi + 1)/(sigma*(dat*xi/sigma + 1))	)%*%V
 }

#' @inheritParams gpd
#' @export
#' @rdname gpd.temstat
gpd.dphi <- function(par, dat, V){
  sigma = par[1]; xi = as.vector(par[2])
  if(!isTRUE(all.equal(xi,0))){
    rbind(xi*(1/xi + 1)/(sigma^2*(dat*xi/sigma + 1)) - dat*xi^2*(1/xi + 1)/(sigma^3*(dat*xi/sigma + 1)^2),
          -(1/xi + 1)/(sigma*(dat*xi/sigma + 1)) + dat*(xi + 1)/(sigma^2*(dat*xi/sigma + 1)^2) + 1/(sigma*(dat*xi/sigma + 1)*xi))%*%V
  } else{
    rbind(rep(sigma^(-2), length(dat)), (dat - sigma)/sigma^2) %*% V
  }
}


#' log likelihood for the generalized extreme value distribution
#'
#' Function returning the density of an \code{n} sample from the GEV distribution.
#'
#' \code{gev.ll.optim} returns the negative log likelihood parametrized in terms of location, \code{log(scale)} and shape in order to perform unconstrained optimization
#'
#' @inheritParams gev
#' @export
#' @keywords internal
#' @seealso \code{\link{gev}}
gev.ll <- function(par, dat){
  dat <- as.vector(dat)
  par <- as.vector(par)
  tx <- pmax(1+par[3]/par[2]*(dat-par[1]),0)
  if(!isTRUE(all.equal(as.vector(par[3]), 0, tolerance = 1e-7))){
    sum(-log(par[2])-(1/par[3]+1)*log(tx)-tx^(-1/par[3]))
  } else{
    sum(-log(par[2])-(dat-par[1])/par[2]-exp(-(dat-par[1])/par[2]))
  }
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
  mu = par[1]; sigma = par[2]; xi = as.vector(par[3])
  if(!isTRUE(all.equal(xi, 0))){
    c(sum(-(-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)/sigma - xi*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1))),
      sum(-(dat - mu)*((dat - mu)*xi/sigma + 1)^(-1/xi - 1)/sigma^2 + (dat - mu)*(xi + 1)/(sigma^2*((dat - mu)*xi/sigma + 1)) - 1/sigma),
      sum(-(mu - dat)*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) -
            (log(-(mu - dat)*xi/sigma + 1)/xi^2 - (mu - dat)/(sigma*((mu - dat)*xi/sigma - 1)*xi))/(-(mu - dat)*xi/sigma + 1)^(1/xi) +
            log(-(mu - dat)*xi/sigma + 1)/xi^2))
  } else{
    c(sum(-exp(mu/sigma - dat/sigma)/sigma + 1/sigma),
      sum(mu*exp(mu/sigma - dat/sigma)/sigma^2 - dat*exp(mu/sigma - dat/sigma)/sigma^2 - mu/sigma^2 - 1/sigma + dat/sigma^2),
      0)
  }

}

#' Information matrix for the generalized extreme value distribution
#'
#' The function returns the expected or observed information matrix.
#' @inheritParams gev
#' @param nobs number of observations
#' @export
#' @keywords internal
gev.infomat <- function(par, dat, method=c("obs","exp"), nobs=length(dat)){
  method <- match.arg(method,c("obs","exp"))
  if(missing(method)){
    method="obs"
  }
  if(method=="exp"){
    #(Expected) Fisher information does not depend on location parameter
    if(length(par)==3){sigma = par[2]; xi = as.vector(par[3]);
    } else{	sigma=par[1]; xi=as.vector(par[2]);
    }
    #Limiting case when xi=0
    if(isTRUE(all.equal(xi,0))){
      return(nobs*cbind(c(1/sigma^2,
                          -0.42278433509846713939348790991759756895784066/sigma^2,
                          0.41184033042643969478888356141823227689702419/sigma),
                        c(-0.42278433509846713939348790991759756895784066/sigma^2,
                          1.8236806608528793895777671228364645537940484/sigma^2,
                          0.33248490716027406147000563764932365201044312/sigma),
                        c(0.41184033042643969478888356141823227689702419/sigma,
                          0.33248490716027406147000563764932365201044312/sigma,
                          2.4236060551770285007097120629947238284701404)))
    }
    p = (1+xi)^2*gamma(1+2*xi)
    q = gamma(2+xi)*(digamma(1+xi)+(1+xi)/xi)
    infomat <- nobs*cbind(
      c(p/sigma^2, -(p-gamma(2+xi))/(sigma^2*xi), (p/xi-q)/(sigma*xi)),
      c(-(p-gamma(2+xi))/(sigma^2*xi), (1-2*gamma(2+xi)+p)/(sigma^2*xi^2), -(1+digamma(1)+(1-gamma(2+xi))/xi-q+p/xi)/(sigma*xi^2)),
      c((p/xi-q)/(sigma*xi), -(1+digamma(1)+(1-gamma(2+xi))/xi-q+p/xi)/(sigma*xi^2), (pi^2/6+(1+digamma(1)+1/xi)^2-2*q/xi+p/xi^2)/xi^2)
    )
    if(xi<0.003 & xi>0){
      #Linearization because numerically unstable in this region
      y2 <- (pi^2/6+(1+digamma(1)+1/0.003)^2-2*(gamma(2+0.003)*(digamma(1+0.003)+(1+0.003)/0.003))/0.003+(1+0.003)^2*gamma(1+2*0.003)/0.003^2)/0.003^2
      y1 <- 2.4236060551770285007097120629947238284701404
      infomat[3,3] <- nobs*((y2-y1)/0.003*xi+y1)
    } else if(xi > -0.003 & xi < 0){
      y2 <- (pi^2/6+(1+digamma(1)-1/0.003)^2+2*(gamma(2-0.003)*(digamma(1-0.003)-(1-0.003)/0.003))/0.003+(1-0.003)^2*gamma(1-2*0.003)/0.003^2)/0.003^2
      y1 <- 2.4236060551770285007097120629947238284701404
      infomat[3,3] <- nobs*(-(y2-y1)/0.003*xi+y1)
    }
    return(infomat)
  } else if(method=="obs"){
    if(length(par)!=3){
      stop("Invalid parameter vector")
    }
    mu = par[1]; sigma = par[2]; xi = as.vector(par[3]);
    #Bug fixed 21-10-2016 (parameter were defined after they were used).
    if(any((1 + xi*(dat-mu)/sigma)<0)){
      stop("Data outside of range specified by parameter, yielding a zero likelihood")
    }
    infomat <- matrix(0, ncol=3, nrow=3)
    if(!isTRUE(all.equal(xi,0))){
      infomat[1,1] <- sum(-((dat - mu)*xi/sigma + 1)^(-1/xi - 2)*(xi + 1)/sigma^2 + xi^2*(1/xi + 1)/(sigma^2*((dat - mu)*xi/sigma + 1)^2))
      infomat[1,2] <- infomat[2,1] <- sum(-(dat - mu)*((dat - mu)*xi/sigma + 1)^(-1/xi - 2)*(xi + 1)/sigma^3 + ((dat - mu)*xi/sigma + 1)^(-1/xi - 1)/sigma^2 - xi*(1/xi + 1)/(sigma^2*((dat - mu)*xi/sigma + 1)) + (dat - mu)*xi^2*(1/xi + 1)/(sigma^3*((dat - mu)*xi/sigma + 1)^2))
      infomat[1,3] <- infomat[3,1] <- sum((-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)*((mu - dat)*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) - log(-(mu - dat)*xi/sigma + 1)/xi^2)/sigma - (1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) + (mu - dat)*(xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)^2) + 1/(sigma*((mu - dat)*xi/sigma - 1)*xi))
      infomat[2,3] <- infomat[3,2] <- sum((mu - dat)*(-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)*(log(-(mu - dat)*xi/sigma + 1)/xi^2 - (mu - dat)/(sigma*((mu - dat)*xi/sigma - 1)*xi))/sigma^2 + (mu - dat)*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)) - (mu - dat)^2*(xi + 1)/(sigma^3*((mu - dat)*xi/sigma - 1)^2) - (mu - dat)/(sigma^2*((mu - dat)*xi/sigma - 1)*xi) + (mu - dat)^2/(sigma^3*(-(mu - dat)*xi/sigma + 1)^(1/xi)*((mu - dat)*xi/sigma - 1)^2))
      infomat[3,3] <- sum(-(log(-(mu - dat)*xi/sigma + 1)/xi^2 - (mu - dat)/(sigma*((mu - dat)*xi/sigma - 1)*xi))^2/(-(mu - dat)*xi/sigma + 1)^(1/xi) + (2*log(-(mu - dat)*xi/sigma + 1)/xi^3 - 2*(mu - dat)/(sigma*((mu - dat)*xi/sigma - 1)*xi^2) - (mu - dat)^2/(sigma^2*((mu - dat)*xi/sigma - 1)^2*xi))/(-(mu - dat)*xi/sigma + 1)^(1/xi) + (mu - dat)^2*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)^2) - 2*log(-(mu - dat)*xi/sigma + 1)/xi^3 + 2*(mu - dat)/(sigma*((mu - dat)*xi/sigma - 1)*xi^2))
      infomat[2,2] <- sum(-(mu - dat)^2*(-(mu - dat)*xi/sigma + 1)^(-1/xi - 2)*(xi + 1)/sigma^4 - 2*(mu - dat)*(-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)/sigma^3 - 2*(mu - dat)*(xi + 1)/(sigma^3*((mu - dat)*xi/sigma - 1)) + (mu - dat)^2*xi^2*(1/xi + 1)/(sigma^4*((mu - dat)*xi/sigma - 1)^2) + 1/sigma^2)
    } else {
      infomat[1,1] <- sum(-exp(mu/sigma - dat/sigma)/sigma^2)
      infomat[1,2] <- infomat[2,1] <- sum(((mu + sigma)*exp(mu/sigma) - dat*exp(mu/sigma) - sigma*exp(dat/sigma))*exp(-dat/sigma)/sigma^3)
      infomat[2,2] <- sum(-2*mu*exp(mu/sigma - dat/sigma)/sigma^3 + 2*dat*exp(mu/sigma - dat/sigma)/sigma^3 + 2*mu/sigma^3 + 1/sigma^2 - 2*dat/sigma^3 - (mu^2*exp(mu/sigma) - 2*mu*dat*exp(mu/sigma) + dat^2*exp(mu/sigma))*exp(-dat/sigma)/sigma^4)
      infomat[1,3] <- infomat[3,1] <- sum(-1/2*(mu^2*exp(mu/sigma) + 2*mu*sigma*exp(mu/sigma) - 2*mu*dat*exp(mu/sigma) - 2*sigma*dat*exp(mu/sigma) + dat^2*exp(mu/sigma) - 2*mu*sigma*exp(dat/sigma) - 2*sigma^2*exp(dat/sigma) + 2*sigma*dat*exp(dat/sigma))*exp(-dat/sigma)/sigma^3)
      infomat[2,3] <- infomat[3,2] <- sum(1/2*(mu^2*exp(mu/sigma) + 2*mu*sigma*exp(mu/sigma) - 2*mu*dat*exp(mu/sigma) - 2*sigma*dat*exp(mu/sigma) + dat^2*exp(mu/sigma) - 2*mu*sigma*exp(dat/sigma) - 2*sigma^2*exp(dat/sigma) + 2*sigma*dat*exp(dat/sigma))*(mu - dat)*exp(-dat/sigma)/sigma^4)
      infomat[3,3] <- sum(-1/12*(3*mu^2*exp(mu/sigma) + 8*mu*sigma*exp(mu/sigma) - 6*mu*dat*exp(mu/sigma) - 8*sigma*dat*exp(mu/sigma) + 3*dat^2*exp(mu/sigma) - 8*mu*sigma*exp(dat/sigma) - 12*sigma^2*exp(dat/sigma) + 8*sigma*dat*exp(dat/sigma))*(mu - dat)^2*exp(-dat/sigma)/sigma^4)
    }

    return(-infomat)
  }
}


#' Tangent exponential model statistics for the generalized extreme value distribution
#'
#' Matrix of approximate ancillary statistics, sample space derivative of the
#' log likelihood and mixed derivative for the generalized extreme value distribution.
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
  mu = par[1]; sigma = par[2]; xi = as.vector(par[3])
  if(!isTRUE(all.equal(xi, 0))){
    t(((dat - mu)*xi/sigma + 1)^(-1/xi - 1)/sigma + xi*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)))%*%V
  } else{
    t(exp((mu-dat)/sigma)/sigma - 1/sigma)%*%V

  }
}


#' @rdname gev.temstat
#' @inheritParams gev
#' @export
#' @keywords internal
gev.dphi <- function(par, dat, V){
  mu = par[1]; sigma = par[2]; xi = as.vector(par[3])
  if(!isTRUE(all.equal(xi,0))){
    rbind((-(mu - dat)*xi/sigma + 1)^(-1/xi - 2)*(xi + 1)/sigma^2 - xi^2*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)^2),
          -(mu - dat)*(-(mu - dat)*xi/sigma + 1)^(-1/xi - 2)*(xi + 1)/sigma^3 - (-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)/sigma^2 - xi*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)) + (mu - dat)*xi^2*(1/xi + 1)/(sigma^3*((mu - dat)*xi/sigma - 1)^2),
          -(-(mu - dat)*xi/sigma + 1)^(-1/xi - 1)*((mu - dat)*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) - log(-(mu - dat)*xi/sigma + 1)/xi^2)/sigma + (1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) - (mu - dat)*(xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1)^2) - 1/(sigma*((mu - dat)*xi/sigma - 1)*xi))%*%V
  } else{
    rbind(exp(-dat/sigma + mu/sigma)/sigma^2,
          (sigma*exp(dat/sigma) + (dat - mu - sigma)*exp(mu/sigma))*exp(-dat/sigma)/sigma^3,
          1/2*(2*dat*sigma*exp(dat/sigma) - 2*mu*sigma*exp(dat/sigma) - 2*sigma^2*exp(dat/sigma) + dat^2*exp(mu/sigma) - 2*dat*mu*exp(mu/sigma) + mu^2*exp(mu/sigma) - 2*dat*sigma*exp(mu/sigma) + 2*mu*sigma*exp(mu/sigma))*exp(-dat/sigma)/sigma^3)%*%V
  }
}

#' @title Generalized Pareto distribution (expected shortfall parametrization)
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized Pareto distribution parametrized in terms of expected shortfall.
#'
#' The parameter \code{m} corresponds to \eqn{\zeta_u}/(1-\eqn{\alpha}), where \eqn{\zeta_u} is the rate of exceedance over the threshold
#' \code{u} and \eqn{\alpha} is the percentile of the expected shortfall.
#' Note that the actual parametrization is in terms of excess expected shortfall, meaning expected shortfall minus threshold.
#'
#' @details The observed information matrix was calculated from the Hessian using symbolic calculus in Sage.
#'
#' @author Leo Belzile
#' @name gpde
#' @param par vector of length 2 containing \eqn{e_m} and \eqn{\xi}, respectively the expected shortfall at probability 1/(1-\eqn{\alpha}) and the shape parameter.
#' @param dat sample vector
#' @param m number of observations of interest for return levels. See \strong{Details}
#' @param tol numerical tolerance for the exponential model
#' @param method string indicating whether to use the expected  (\code{"exp"}) or the observed (\code{"obs"} - the default) information matrix.
#' @param nobs number of observations
#' @param V vector calculated by \code{gpde.Vfun}
#'
#' @section Usage: \preformatted{gpde.ll(par, dat, m, tol=1e-5)
#' gpde.ll.optim(par, dat, m, tol=1e-5)
#' gpde.score(par, dat, m)
#' gpde.infomat(par, dat, m, method = c("obs", "exp"), nobs = length(dat))
#' gpde.Vfun(par, dat, m)
#' gpde.phi(par, dat, V, m)
#' gpde.dphi(par, dat, V, m)}
#'
#' @section Functions:
#'
#' \itemize{
#' \item{\code{gpde.ll}:} {log likelihood}
#' \item{\code{gpde.ll.optim}:} {negative log likelihood parametrized in terms of log expected
#' shortfall and shape in order to perform unconstrained optimization}
#' \item{\code{gpde.score}:} {score vector}
#' \item{\code{gpde.infomat}:} {observed information matrix for GPD parametrized in terms of rate of expected shortfall and shape}
#' \item{\code{gpde.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gpde.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gpde.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
NULL

#' @title Generalized Pareto distribution (return level parametrization)
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized Pareto distribution parametrized in terms of return levels.
#'
#' @details The observed information matrix was calculated from the Hessian using symbolic calculus in Sage.
#' @details The interpretation for \code{m} is as follows: if there are on average \eqn{m_y} observations per year above the threshold, then  \eqn{m=Tm_y} corresponds to \eqn{T}-year return level.
#'
#' @author Leo Belzile
#' @name gpdr
#' @param par vector of length 2 containing \eqn{y_m} and \eqn{\xi}, respectively the \eqn{m}-year return level and the shape parameter.
#' @param dat sample vector
#' @param m number of observations of interest for return levels. See \strong{Details}
#' @param tol numerical tolerance for the exponential model
#' @param method string indicating whether to use the expected  (\code{"exp"}) or the observed (\code{"obs"} - the default) information matrix.
#' @param nobs number of observations
#' @param V vector calculated by \code{gpdr.Vfun}
#'
#' @section Usage: \preformatted{gpdr.ll(par, dat, m, tol=1e-5)
#' gpdr.ll.optim(par, dat, m, tol=1e-5)
#' gpdr.score(par, dat, m)
#' gpdr.infomat(par, dat, m, method = c("obs", "exp"), nobs = length(dat))
#' gpdr.Vfun(par, dat, m)
#' gpdr.phi(par, V, dat, m)
#' gpdr.dphi(par, V, dat, m)}
#' @section Functions:
#'
#' \itemize{
#' \item{\code{gpdr.ll}:} {log likelihood}
#' \item{\code{gpdr.ll.optim}:} {negative log likelihood parametrized in terms of \code{log(scale)} and shape
#' in order to perform unconstrained optimization}
#' \item{\code{gpdr.score}:} {score vector}
#' \item{\code{gpdr.infomat}:} {observed information matrix for GPD parametrized in terms of rate of \eqn{m}-year return level and shape}
#' \item{\code{gpdr.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gpdr.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gpdr.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
NULL

#' @title Generalized extreme value distribution (return level parametrization)
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized extreme value distribution  parametrized in terms of the return level \eqn{z}, scale and shape.
#'
#' @author Leo Belzile
#' @name gevr
#' @param par vector of \code{retlev}, \code{scale} and \code{shape}
#' @param dat sample vector
#' @param p tail probability, corresponding to \eqn{(1-p)}th quantile for \eqn{z}
#' @param method string indicating whether to use the expected  (\code{"exp"}) or the observed (\code{"obs"} - the default) information matrix.
#' @param nobs number of observations
#' @param V vector calculated by \code{gevr.Vfun}
#'
#' @section Usage: \preformatted{gevr.ll(par, dat, p)
#' gevr.ll.optim(par, dat, p)
#' gevr.score(par, dat, p)
#' gevr.infomat(par, dat, p, method = c("obs", "exp"), nobs = length(dat))
#' gevr.Vfun(par, dat, p)
#' gevr.phi(par, dat, p, V)
#' gevr.dphi(par, dat, p, V)}
#'
#' @section Functions:
#' \itemize{
#' \item{\code{gevr.ll}:} {log likelihood}
#' \item{\code{gevr.ll.optim}:} {negative log likelihood parametrized in terms of return levels, \code{log(scale)} and shape in order to perform unconstrained optimization}
#' \item{\code{gevr.score}:} {score vector}
#' \item{\code{gevr.infomat}:} {observed information matrix}
#' \item{\code{gevr.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gevr.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gevr.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
NULL

#' Negative log likelihood of the generalized Pareto distribution (expected shortfall)
#'
#' @seealso \code{\link{gpde}}
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.ll <- function(par, dat, m){
  es = par[1]; xi = par[2]
  if(any(xi > 1, es < 0, min(1 + (m^xi - 1 + xi)/((1 - xi) * es) * dat) < 0)){ return(-1e10)}
  ifelse(!isTRUE(all.equal(xi, 0, tolerance=1e-7)),
  sum(-log(xi*(1-xi)/(m^xi-1+xi))-log(es)-(1+1/xi)*log(1+(m^xi-1+xi)/((1-xi)*es)*dat)),
   sum(-log(1/(log(m) + 1))-log(es)-exp(dat*log(m)/es + dat/es)))
}

#' Negative log likelihood of the generalized Pareto distribution (expected shortfall) - optimization
#' The negative log likelihood is parametrized in terms of log expected shortfall and shape in order to perform unconstrained optimization
#' @rdname gpde.ll
#' @seealso \code{\link{gpde}}
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.ll.optim <- function(par, dat, m){
  es = exp(par[1]); xi = par[2]
  if(xi>1){ return(1e10)}
  -sum(-log(xi*(1-xi)/(m^xi-1+xi))-log(es)-(1+1/xi)*log(1+(m^xi-1+xi)/((1-xi)*es)*dat))
}
#' Score vector for the GP distribution (expected shortfall)
#' @seealso \code{\link{gpde}}
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.score <- function(par, dat, m){
  es = par[1]; xi = par[2]
  if(xi > 1){ return(rep(NA,2))}
  c(sum(-1/es + dat*(m^xi + xi - 1)*(1/xi + 1)/(es^2*(xi - 1)*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1))),
    sum(-((m^xi*log(m) + 1)*dat/(es*(xi - 1)) - dat*(m^xi + xi - 1)/(es*(xi - 1)^2))*(1/xi + 1)/(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1) + (m^xi + xi - 1)*((m^xi*log(m) + 1)*(xi - 1)*xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi + xi - 1))/((xi - 1)*xi) + log(-dat*(m^xi + xi - 1)/(es*(xi - 1)) + 1)/xi^2))
}

#' Observed information matrix for the GP distribution (expected shortfall)
#'
#' The information matrix is parametrized in terms of excess expected shortfall and shape
#' @seealso \code{\link{gpde}}
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.infomat <- function(par, dat, m, method = c("obs", "exp"), nobs = length(dat)){
  method <- match.arg(method, c("obs", "exp")) #default to observed information
  es = as.vector(par[1]); xi = as.vector(par[2])
  if(xi >= 1 || xi <= -0.5){ return(matrix(NA, 2, 2))}
  if(method == "obs"){
    k11 = sum(-1/es^2 + 2*dat*(m^xi + xi - 1)*(1/xi + 1)/(es^3*(xi - 1)*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)) - dat^2*(m^xi + xi - 1)^2*(1/xi + 1)/(es^4*(xi - 1)^2*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)^2))
    k22 = sum(-((m^xi*log(m) + 1)*dat/(es*(xi - 1)) - dat*(m^xi + xi - 1)/(es*(xi - 1)^2))^2*(1/xi + 1)/(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)^2 + (dat*m^xi*log(m)^2/(es*(xi - 1)) - 2*(m^xi*log(m) + 1)*dat/(es*(xi - 1)^2) + 2*dat*(m^xi + xi - 1)/(es*(xi - 1)^3))*(1/xi + 1)/(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1) - (m^xi*(xi - 1)*xi*log(m)^2/(m^xi + xi - 1)^2 - 2*(m^xi*log(m) + 1)^2*(xi - 1)*xi/(m^xi + xi - 1)^3 + 2*(m^xi*log(m) + 1)*(xi - 1)/(m^xi + xi - 1)^2 + 2*(m^xi*log(m) + 1)*xi/(m^xi + xi - 1)^2 - 2/(m^xi + xi - 1))*(m^xi + xi - 1)/((xi - 1)*xi) - (m^xi*log(m) + 1)*((m^xi*log(m) + 1)*(xi - 1)*xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi + xi - 1))/((xi - 1)*xi) + (m^xi + xi - 1)*((m^xi*log(m) + 1)*(xi - 1)*xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi + xi - 1))/((xi - 1)*xi^2) + (m^xi + xi - 1)*((m^xi*log(m) + 1)*(xi - 1)*xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi + xi - 1))/((xi - 1)^2*xi) - 2*((m^xi*log(m) + 1)*dat/(es*(xi - 1)) - dat*(m^xi + xi - 1)/(es*(xi - 1)^2))/(xi^2*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)) + 2*log(-dat*(m^xi + xi - 1)/(es*(xi - 1)) + 1)/xi^3)
    k12 = sum(-((m^xi*log(m) + 1)*dat/(es^2*(xi - 1)) - dat*(m^xi + xi - 1)/(es^2*(xi - 1)^2))*(1/xi + 1)/(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1) + dat*(m^xi + xi - 1)*((m^xi*log(m) + 1)*dat/(es*(xi - 1)) - dat*(m^xi + xi - 1)/(es*(xi - 1)^2))*(1/xi + 1)/(es^2*(xi - 1)*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)^2) + dat*(m^xi + xi - 1)/(es^2*(xi - 1)*xi^2*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)))
    return(cbind(c(k11,k12),c(k12,k22)))
  } else if (method == "exp"){
    sigmae = es*(1-xi)*xi/(m^xi-1+xi)
    Jac <- rbind(c((1-xi)*xi/(m^xi-1+xi),
                   (m^xi*log(m) + 1)*es*(xi - 1)*xi/(m^xi + xi - 1)^2 - es*(xi - 1)/(m^xi + xi - 1) - es*xi/(m^xi + xi - 1)),
                 c(0,1))
    return(t(Jac) %*% gpd.infomat(par = c(sigmae, xi), dat = dat, method = "exp", nobs = nobs) %*% Jac)
  }
}
#' Tangent exponential model statistics for the generalized Pareto distribution (expected shortfall)
#'
#' Vector implementing conditioning on approximate ancillary statistics for the TEM
#' @seealso \code{\link{gpde}}
#' @name gpde.temstat
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.Vfun <- function(par, dat, m){
  es = par[1]; xi = par[2]
  cbind(dat/es,
        es*(xi - 1)*xi*(-dat*(m^xi + xi - 1)/(es*(xi - 1)) + 1)^(-1/xi)*(((m^xi*log(m) + 1)*dat/(es*(xi - 1)) - dat*(m^xi + xi - 1)/(es*(xi - 1)^2))/(xi*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)) - log(-dat*(m^xi + xi - 1)/(es*(xi - 1)) + 1)/xi^2)/((m^xi + xi - 1)*(-dat*(m^xi + xi - 1)/(es*(xi - 1)) + 1)^(-1/xi - 1))
  )
}

#' Canonical parameter in the local exponential family approximation
#'
#' @inheritParams gpde
#' @rdname gpde.temstat
#' @keywords internal
#' @export
gpde.phi <- function(par, dat, V, m){
  es = par[1]; xi = par[2]
  rbind(-(m^xi + xi - 1)*(1/xi + 1)/(es*(xi - 1)*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)))%*%V
}

## Derivative matrix of the canonical parameter in the local exponential family approximation
#' @inheritParams gpde
#' @rdname gpde.temstat
#' @keywords internal
#' @export
gpde.dphi <- function(par, dat, V, m){
  es = par[1]; xi = par[2]
  rbind((m^xi + xi - 1)*(1/xi + 1)/(es^2*(xi - 1)*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)) - dat*(m^xi + xi - 1)^2*(1/xi + 1)/(es^3*(xi - 1)^2*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)^2)
        ,(m^xi + xi - 1)*((m^xi*log(m) + 1)*dat/(es*(xi - 1)) - dat*(m^xi + xi - 1)/(es*(xi - 1)^2))*(1/xi + 1)/(es*(xi - 1)*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)^2) - (m^xi*log(m) + 1)*(1/xi + 1)/(es*(xi - 1)*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)) + (m^xi + xi - 1)*(1/xi + 1)/(es*(xi - 1)^2*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)) + (m^xi + xi - 1)/(es*(xi - 1)*xi^2*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)))%*%V
}

#' Negative log likelihood of the generalized Pareto distribution (return levels)
#'
#' @seealso \code{\link{gpdr}}
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.ll <- function(par, dat, m, tol=1e-5){
  ym = par[1]; xi = par[2]
  if(par[1]<0){ return(-1e15)}
  nn = length(dat)
  pr <- m^xi-1
  if(abs(xi)>tol){
    -nn*log(xi/pr)-nn*log(ym)-(1+1/xi)*sum(log(pmax(1+exp(log(dat)-log(ym))*pr,0)))
  } else{
    nn*log(log(m))-nn*log(ym)+sum(-dat*log(m)/ym)
  }
}

#' Negative log likelihood parametrized in terms of log return level and shape in order to perform unconstrained optimization
#' @seealso \code{\link{gpdr}}
#' @inheritParams gpdr
#' @keywords internal
#' @rdname gpdr.ll
#' @export
gpdr.ll.optim <- function(par, dat, m){
  tol <- 1e-5
  nn = length(dat)
  ym = exp(par[1]); xi = par[2]
  p <- m^xi-1
  if(abs(xi)>tol){
    -(-nn*log(xi/p)-nn*log(ym)-(1+1/xi)*sum(log(pmax(1+exp(log(dat)-log(ym))*p,0))))
  } else{
    -( nn*log(log(m))-nn*log(ym)+sum(-dat*log(m)/ym))
  }
}

#' Score of the profile log likelihood for the GP distribution (return levels parametrization)
#' @seealso \code{\link{gpdr}}
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.score <- function(par, dat, m){
  nn = length(dat); xi = par[2];  ym = par[1]
  p <- (m)^xi-1
  c(	-nn/ym+(1+1/xi)*sum(dat*p/(ym+dat*p))/ym,
     -nn/xi+exp(log(nn)+xi*log(m)+log(log(m)))/p+sum(log(1+dat/ym*p))/xi^2-(1+1/xi)*
       sum(exp(log(dat)+log(log(m))+xi*log(m)-log(ym+dat*p)))
  )
}

#' Observed information matrix for GP distribution (return levels)
#'
#'The information matrix is parametrized in terms of rate of \code{m}-year return level and shape
#' @seealso \code{\link{gpdr}}
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.infomat <- function(par, dat, m, method = c("obs", "exp"), nobs = length(dat)){
  xi = as.vector(par[2]); r = as.vector(par[1])
  method <- method[1]
  if(xi >= 1 || xi <= -0.5){ return(matrix(NA, 2, 2))}
  if(method == "obs"){
    info <- matrix(ncol=2,nrow=2)
    info[1,1] <- sum(dat^2*(m^xi - 1)^2*(1/xi + 1)/((dat*(m^xi - 1)/r + 1)^2*r^4) - 2*dat*(m^xi - 1)*(1/xi + 1)/((dat*(m^xi - 1)/r + 1)*r^3) + 1/r^2)
    info[2,2] <- sum(dat^2*(m^xi)^2*(1/xi + 1)*log(m)^2/((dat*(m^xi - 1)/r + 1)^2*r^2) - dat*m^xi*(1/xi + 1)*log(m)^2/((dat*(m^xi - 1)/r + 1)*r) + m^xi*(m^xi*xi*log(m)/(m^xi - 1)^2 - 1/(m^xi - 1))*log(m)/xi + (m^xi*xi*log(m)^2/(m^xi - 1)^2 - 2*(m^xi)^2*xi*log(m)^2/(m^xi - 1)^3 + 2*m^xi*log(m)/(m^xi - 1)^2)*(m^xi - 1)/xi - (m^xi - 1)*(m^xi*xi*log(m)/(m^xi - 1)^2 - 1/(m^xi - 1))/xi^2 + 2*dat*m^xi*log(m)/((dat*(m^xi - 1)/r + 1)*r*xi^2) - 2*log(dat*(m^xi - 1)/r + 1)/xi^3)
    info[2,1] <- info[1,2] <- sum(-dat^2*(m^xi - 1)*m^xi*(1/xi + 1)*log(m)/((dat*(m^xi - 1)/r + 1)^2*r^3) + dat*m^xi*(1/xi + 1)*log(m)/((dat*(m^xi - 1)/r + 1)*r^2) - dat*(m^xi - 1)/((dat*(m^xi - 1)/r + 1)*r^2*xi^2))
    return(-info)
  } else if(method == "exp"){
    sigmar = r*xi/(m^xi-1)
    Jac <- rbind(c(xi/(m^xi - 1), -m^xi*r*xi*log(m)/(m^xi - 1)^2 + r/(m^xi - 1)), c(0, 1))
    return(t(Jac) %*% gpd.infomat(par = c(sigmar, xi), method = "exp", dat = dat, nobs = nobs) %*% Jac)
  }

}

#' Tangent exponential model statistics for the generalized Pareto distribution (return level)
#'
#' Vector implementing conditioning on approximate ancillary statistics for the TEM
#' @seealso \code{\link{gpdr}}
#' @name gpdr.temstat
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.Vfun <- function(par, dat, m){
  xi = par[2]; ym = par[1]
  p <- m^xi-1
  cbind(dat/ym,
        -(ym+dat*p)/p*((m)^xi*log(m)/(p+ym/dat)-log(1+dat/ym*p)/xi)
  )
}

#' Canonical parameter in the local exponential family approximation
#' @seealso \code{\link{gpdr}}
#' @rdname gpdr.temstat
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.phi <- function(par, dat, V, m){
  xi = par[2]; r = par[1]
  matrix(-(1+1/xi)*(m^xi-1)/(r+dat*(m^xi-1)),nrow=1)%*%V
}

#' Derivative of the canonical parameter \eqn{\phi(\theta)} in the local exponential family approximation
#' @rdname gpdr.temstat
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.dphi <- function(par, dat, V, m){
  xi = par[2];  r = par[1];  p = m^xi-1
  rbind((1+1/xi)*p/(r+dat*p)^2,
        p/(xi^2*(r+dat*p))-(1+1/xi)*log(m)*m^xi*r/(r+dat*p)^2
  )%*%V

}
#######################################
###  GEV in terms of return levels  ###
#######################################


#' Return level for the generalized extreme value distribution
#'
#' This function returns the \eqn{1-p}th quantile of the GEV.
#' @inheritParams gev
#' @seealso \code{\link{gev}}
#' @keywords internal
#' @export
gev.retlev <- function(par, p){
  mu = par[1]; sigma = par[2]; xi = as.vector(par[3])
  if(!isTRUE(all.equal(xi, 0))){
    mu - sigma/xi*(1-(-log(1-p))^(-xi))
  } else{
    mu - sigma*log(-log(-p + 1))
  }
}

#' This function returns the mean of N observations from the GEV.
#' @inheritParams gevN
#' @seealso \code{\link{gevN}}
#' @keywords internal
#' @export
gevN.mean <- function(par, N){
  mu = par[1]; sigma = par[2]; xi = as.vector(par[3])
  if(!isTRUE(all.equal(xi, 0))){
    mu - sigma / xi *(1 - N^xi*gamma(1-xi))
  } else{
    mu + sigma * (log(N) -psigamma(1))
  }
}

gevN.quant <- function(par, N, q){
  mu = par[1]; sigma = par[2]; xi = as.vector(par[3])
  if(!isTRUE(all.equal(xi, 0))){
    mu - sigma / xi *(1 - (N / log(1/q))^xi)
  } else {
    mu + sigma * (log(N) - log(log(1/q)))
  }
}

#' This function returns the mean of N maxima from the GP.
#' @inheritParams gpdN
#' @seealso \code{\link{gpdN}}
#' @keywords internal
#' @export
gpdN.mean <- function(par, N){
  sigma = par[1]; xi = as.vector(par[2])
  if(!isTRUE(all.equal(xi, 0))){
    (exp(lgamma(N + 1) + lgamma(1-xi) - lgamma(N + 1-xi))-1)*sigma/xi
  } else{
    (-psigamma(1) + psigamma(N + 1))*sigma
  }
}

#' This function returns the qth percentile of N maxima from the GP.
#' @inheritParams gpdN
#' @seealso \code{\link{gpdN}}
#' @keywords internal
#' @export
gpdN.quant <- function(par, q, N){
  sigma = par[1]; xi = as.vector(par[2])
  if(!isTRUE(all.equal(xi, 0))){
    sigma/xi*((1-q^(1/N))^(-xi)-1)
  } else{
    sigma*(-log(-q^(1/N) + 1))
  }
}


#' Negative log likelihood of the generalized extreme value distribution (return levels)
#'
#' @seealso \code{\link{gevr}}
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.ll <- function(par, dat, p){
  z = par[1]; sigma = par[2]; xi = as.vector(par[3])
  muf = ifelse(!isTRUE(all.equal(xi,0)), z + sigma/xi*(1-(-log(1-p))^(-xi)), sigma*log(-log(-p + 1)) + z)
  if(!isTRUE(all.equal(xi, 0, tolerance = 1e-8))){
    sum(-log(sigma)-(1/xi+1)*log(pmax(1+xi*(dat-muf)/sigma,0))-pmax(1+xi*(dat-muf)/sigma,0)^(-1/xi))
  } else{
    sum((muf - dat)/sigma - exp((muf-dat)/sigma) - log(sigma))
  }

}

#' Negative log likelihood parametrized in terms of location, log return level and shape in order to perform unconstrained optimization
#' @rdname gevr.ll
#' @seealso \code{\link{gevr}}
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.ll.optim <- function(par, dat, p){
  tpar = par; tpar[2] = exp(par[2])
  #if(tpar[3] <= -1){return(1e20)}
  -gevr.ll(tpar,dat, p)
}

#' Score of the log likelihood for the GEV distribution (return levels)
#' @seealso \code{\link{gevr}}
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.score <- function(par, dat, p){
  z = par[1]; sigma = par[2]; xi = par[3];
  c(sum((((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(-1/xi - 2)*(xi + 1)*exp(-1/((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi))/sigma^2 - ((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(-2/xi - 2)*exp(-1/((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi))/sigma^2)*sigma*((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi + 1)*exp(1/(((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi)))),
    sum(-(dat*(-log(-p + 1))^xi - z*(-log(-p + 1))^xi - (dat*(-log(-p + 1))^xi - z*(-log(-p + 1))^xi - sigma)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))/((sigma*dat*xi*(-log(-p + 1))^xi - sigma*xi*z*(-log(-p + 1))^xi + sigma^2)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))),
    sum(-(xi*z*(-log(-p + 1))^xi - (dat*(-log(-p + 1))^xi - sigma*log(-log(-p + 1)))*xi + ((dat*(-log(-p + 1))^xi - sigma*log(-log(-p + 1)))*xi^2 + (dat*(-log(-p + 1))^xi - sigma*log(-log(-p + 1)))*xi - (xi^2*(-log(-p + 1))^xi + xi*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + (dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi - (dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + sigma)*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi)))/((dat*xi^3*(-log(-p + 1))^xi - xi^3*z*(-log(-p + 1))^xi + sigma*xi^2)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))
    ))
}

#' Observed information matrix for GEV distribution (return levels)
#'
#'The information matrix is parametrized in terms of return level (\eqn{(1-p)}th quantile), scale and shape.
#' @seealso \code{\link{gevr}}
#' @inheritParams gevr
#' @param nobs integer number of observations
#' @keywords internal
#' @export
gevr.infomat <- function(par, dat, method = c("obs", "exp"), p, nobs = length(dat)){
  method <- match.arg(method, c("obs", "exp"))
  z = par[1]; sigma = par[2]; xi = par[3];
  if(method == "obs"){
    infomat <- matrix(0, ncol=3,nrow=3)

    infomat[1,1] <- sum(-(xi*(-log(-p + 1))^(2*xi) - (xi^2*(-log(-p + 1))^(2*xi) + xi*(-log(-p + 1))^(2*xi))*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + (-log(-p + 1))^(2*xi))/((dat^2*xi^2*(-log(-p + 1))^(2*xi) + xi^2*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*xi*(-log(-p + 1))^xi + sigma^2 - 2*(dat*xi^2*(-log(-p + 1))^(2*xi) + sigma*xi*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))

    infomat[1,2] <- infomat[2,1] <-  sum(-(dat*(-log(-p + 1))^(2*xi) - z*(-log(-p + 1))^(2*xi) - sigma*(-log(-p + 1))^xi + (sigma*xi*(-log(-p + 1))^xi + sigma*(-log(-p + 1))^xi)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))/((sigma*dat^2*xi^2*(-log(-p + 1))^(2*xi) + sigma*xi^2*z^2*(-log(-p + 1))^(2*xi) + 2*sigma^2*dat*xi*(-log(-p + 1))^xi + sigma^3 - 2*(sigma*dat*xi^2*(-log(-p + 1))^(2*xi) + sigma^2*xi*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))
    infomat[1,3] <- infomat[3,1] <- sum(-((sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - dat*(-log(-p + 1))^(2*xi))*xi^2 + (sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - dat*(-log(-p + 1))^(2*xi))*xi + (xi^2*(-log(-p + 1))^(2*xi) + xi*(-log(-p + 1))^(2*xi))*z - (sigma*xi^3*(-log(-p + 1))^xi*log(-log(-p + 1)) + xi^2*z*(-log(-p + 1))^(2*xi) + (sigma*(-log(-p + 1))^xi*(log(-log(-p + 1)) + 1) - dat*(-log(-p + 1))^(2*xi))*xi^2)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + (dat*xi*(-log(-p + 1))^(2*xi) - xi*z*(-log(-p + 1))^(2*xi) + sigma*(-log(-p + 1))^xi)*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi)))/((dat^2*xi^4*(-log(-p + 1))^(2*xi) + xi^4*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*xi^3*(-log(-p + 1))^xi + sigma^2*xi^2 - 2*(dat*xi^4*(-log(-p + 1))^(2*xi) + sigma*xi^3*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))
    infomat[2,2] <- sum((dat^2*xi*(-log(-p + 1))^(2*xi) + (xi*(-log(-p + 1))^(2*xi) - (-log(-p + 1))^(2*xi))*z^2 - dat^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*(-log(-p + 1))^xi - 2*(dat*xi*(-log(-p + 1))^(2*xi) - dat*(-log(-p + 1))^(2*xi) + sigma*(-log(-p + 1))^xi)*z - (dat^2*xi*(-log(-p + 1))^(2*xi) + xi*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*(-log(-p + 1))^xi - sigma^2 - 2*(dat*xi*(-log(-p + 1))^(2*xi) + sigma*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))/((sigma^2*dat^2*xi^2*(-log(-p + 1))^(2*xi) + sigma^2*xi^2*z^2*(-log(-p + 1))^(2*xi) + 2*sigma^3*dat*xi*(-log(-p + 1))^xi + sigma^4 - 2*(sigma^2*dat*xi^2*(-log(-p + 1))^(2*xi) + sigma^3*xi*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))
    infomat[3,2] <- infomat[2,3] <- sum(-((sigma*dat*(-log(-p + 1))^xi*log(-log(-p + 1)) - dat^2*(-log(-p + 1))^(2*xi))*xi^2 - (xi^2*(-log(-p + 1))^(2*xi) + xi*(-log(-p + 1))^(2*xi))*z^2 + (sigma*dat*(-log(-p + 1))^xi*log(-log(-p + 1)) - dat^2*(-log(-p + 1))^(2*xi))*xi - ((sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - 2*dat*(-log(-p + 1))^(2*xi))*xi^2 + (sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - 2*dat*(-log(-p + 1))^(2*xi))*xi)*z - (sigma*dat*xi^3*(-log(-p + 1))^xi*log(-log(-p + 1)) - xi^2*z^2*(-log(-p + 1))^(2*xi) + (sigma*dat*(-log(-p + 1))^xi*(log(-log(-p + 1)) + 1) - dat^2*(-log(-p + 1))^(2*xi))*xi^2 - (sigma*xi^3*(-log(-p + 1))^xi*log(-log(-p + 1)) + (sigma*(-log(-p + 1))^xi*(log(-log(-p + 1)) + 1) - 2*dat*(-log(-p + 1))^(2*xi))*xi^2)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + (dat^2*xi*(-log(-p + 1))^(2*xi) + xi*z^2*(-log(-p + 1))^(2*xi) + sigma*dat*(-log(-p + 1))^xi - (2*dat*xi*(-log(-p + 1))^(2*xi) + sigma*(-log(-p + 1))^xi)*z)*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi)))/((sigma*dat^2*xi^4*(-log(-p + 1))^(2*xi) + sigma*xi^4*z^2*(-log(-p + 1))^(2*xi) + 2*sigma^2*dat*xi^3*(-log(-p + 1))^xi + sigma^3*xi^2 - 2*(sigma*dat*xi^4*(-log(-p + 1))^(2*xi) + sigma^2*xi^3*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))
    infomat[3,3] <- sum((sigma*dat*xi^4*(-log(-p + 1))^xi*log(-log(-p + 1))^2 + (4*sigma*dat*(-log(-p + 1))^xi*log(-log(-p + 1)) - 3*dat^2*(-log(-p + 1))^(2*xi))*xi^3 + (2*sigma*dat*(-log(-p + 1))^xi*(log(-log(-p + 1)) - 1) - (log(-log(-p + 1))^2 - 2*log(-log(-p + 1)))*sigma^2 - dat^2*(-log(-p + 1))^(2*xi))*xi^2 - (3*xi^3*(-log(-p + 1))^(2*xi) + xi^2*(-log(-p + 1))^(2*xi))*z^2 - (dat^2*xi^2*(-log(-p + 1))^(2*xi) + xi^2*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*xi*(-log(-p + 1))^xi + sigma^2 - 2*(dat*xi^2*(-log(-p + 1))^(2*xi) + sigma*xi*(-log(-p + 1))^xi)*z)*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^2 - (sigma*xi^4*(-log(-p + 1))^xi*log(-log(-p + 1))^2 + 2*(2*sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - 3*dat*(-log(-p + 1))^(2*xi))*xi^3 + 2*(sigma*(-log(-p + 1))^xi*(log(-log(-p + 1)) - 1) - dat*(-log(-p + 1))^(2*xi))*xi^2)*z - (sigma*dat*xi^5*(-log(-p + 1))^xi*log(-log(-p + 1))^2 + ((log(-log(-p + 1))^2 + 2*log(-log(-p + 1)))*sigma*dat*(-log(-p + 1))^xi - dat^2*(-log(-p + 1))^(2*xi))*xi^4 + (4*sigma*dat*(-log(-p + 1))^xi*log(-log(-p + 1)) - 3*dat^2*(-log(-p + 1))^(2*xi))*xi^3 - 2*(sigma*dat*(-log(-p + 1))^xi - sigma^2*log(-log(-p + 1)))*xi^2 - (xi^4*(-log(-p + 1))^(2*xi) + 3*xi^3*(-log(-p + 1))^(2*xi))*z^2 - (sigma*xi^5*(-log(-p + 1))^xi*log(-log(-p + 1))^2 + ((log(-log(-p + 1))^2 + 2*log(-log(-p + 1)))*sigma*(-log(-p + 1))^xi - 2*dat*(-log(-p + 1))^(2*xi))*xi^4 + 2*(2*sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - 3*dat*(-log(-p + 1))^(2*xi))*xi^3 - 2*sigma*xi^2*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + 2*(dat^2*xi^3*(-log(-p + 1))^(2*xi) - (sigma*dat*(-log(-p + 1))^xi*(log(-log(-p + 1)) - 2) - dat^2*(-log(-p + 1))^(2*xi))*xi^2 + (xi^3*(-log(-p + 1))^(2*xi) + xi^2*(-log(-p + 1))^(2*xi))*z^2 + (sigma*dat*(-log(-p + 1))^xi - sigma^2*(log(-log(-p + 1)) - 1))*xi - (2*dat*xi^3*(-log(-p + 1))^(2*xi) - (sigma*(-log(-p + 1))^xi*(log(-log(-p + 1)) - 2) - 2*dat*(-log(-p + 1))^(2*xi))*xi^2 + sigma*xi*(-log(-p + 1))^xi)*z - (dat^2*xi^3*(-log(-p + 1))^(2*xi) + xi^3*z^2*(-log(-p + 1))^(2*xi) +2*sigma*dat*xi^2*(-log(-p + 1))^xi + sigma^2*xi - 2*(dat*xi^3*(-log(-p + 1))^(2*xi) + sigma*xi^2*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi)))/((dat^2*xi^6*(-log(-p + 1))^(2*xi) + xi^6*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*xi^5*(-log(-p + 1))^xi + sigma^2*xi^4 - 2*(dat*xi^6*(-log(-p + 1))^(2*xi) + sigma*xi^5*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))
return(-infomat)
      #else expected information matrix
  } else{
    muf = z+sigma/xi*(1-(-log(1-p))^(-xi))
    if(!isTRUE(all.equal(xi,0, tolerance = 1e-7))){
      Jac <- rbind(c(1,-((-log(-p + 1))^(-xi) - 1)/xi, sigma*(-log(-p + 1))^(-xi)*log(-log(-p + 1))/xi + sigma*((-log(-p + 1))^(-xi) - 1)/xi^2),
                   c(0,1,0),c(0,0,1))
    } else{
      Jac <- rbind(c(1, log(-log(-p + 1)), -sigma*log(-log(-p + 1))^2 + 1/2*sigma*(-log(-p + 1))^(-xi)*log(-log(-p + 1))^2),
                   c(0,1,0),c(0,0,1))
    }
    return(t(Jac) %*% gev.infomat(dat = dat, method = "exp", par = c(muf, sigma, xi)) %*% Jac)

  }
}


#' Tangent exponential model statistics for the GEV distribution (return level)
#'
#' Vector implementing conditioning on approximate ancillary statistics for the TEM
#' @seealso \code{\link{gevr}}
#' @name gevr.temstat
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.Vfun <- function(par, dat, p) {
  z = par[1]
  sigma = par[2]
  xi = par[3]
  cbind(1, (dat - z) / sigma,
        sigma * ((dat - z) * xi / sigma + (-log(-p + 1)) ^ (-xi)) ^ (-1 / xi) *
          (((-log(-p + 1)) ^ (-xi) * log(-log(-p + 1)) - (dat - z) / sigma
          ) / (((dat - z) * xi / sigma + (-log(-p + 1)) ^ (-xi)) * xi) +
            log((dat - z) * xi / sigma + (-log(-p + 1)) ^ (-xi)) / xi ^ 2) / ((dat - z) * xi / sigma + (-log(-p + 1)) ^ (-xi))^(-1 / xi - 1)
  )
}

#' Canonical parameter in the local exponential family approximation
#' @rdname gevr.temstat
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.phi <- function(par, dat, p, V){
  z = par[1]; sigma = par[2]; xi = par[3];
  t(-((xi*(-log(-p + 1))^xi + (-log(-p + 1))^xi)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) - (-log(-p + 1))^xi)/((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))%*%V
}

#' Derivative of the canonical parameter \eqn{\phi(\theta)} in the local exponential family approximation
#' @rdname gevr.temstat
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.dphi <- function(par, dat, p, V){
  z = par[1]; sigma = par[2]; xi = par[3]
  rbind(
    (xi*(-log(-p + 1))^(2*xi) - (xi^2*(-log(-p + 1))^(2*xi) + xi*(-log(-p + 1))^(2*xi))*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + (-log(-p + 1))^(2*xi))/((dat^2*xi^2*(-log(-p + 1))^(2*xi) + xi^2*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*xi*(-log(-p + 1))^xi + sigma^2 - 2*(dat*xi^2*(-log(-p + 1))^(2*xi) + sigma*xi*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)),

    #dphi sigma
    (dat*(-log(-p + 1))^(2*xi) - z*(-log(-p + 1))^(2*xi) - sigma*(-log(-p + 1))^xi + (sigma*xi*(-log(-p + 1))^xi + sigma*(-log(-p + 1))^xi)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))/((sigma*dat^2*xi^2*(-log(-p + 1))^(2*xi) + sigma*xi^2*z^2*(-log(-p + 1))^(2*xi) + 2*sigma^2*dat*xi*(-log(-p + 1))^xi + sigma^3 - 2*(sigma*dat*xi^2*(-log(-p + 1))^(2*xi) + sigma^2*xi*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)),
    #dphi xi
    ((sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - dat*(-log(-p + 1))^(2*xi))*xi^2 + (sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - dat*(-log(-p + 1))^(2*xi))*xi + (xi^2*(-log(-p + 1))^(2*xi) + xi*(-log(-p + 1))^(2*xi))*z - (sigma*xi^3*(-log(-p + 1))^xi*log(-log(-p + 1)) + xi^2*z*(-log(-p + 1))^(2*xi) + (sigma*(-log(-p + 1))^xi*(log(-log(-p + 1)) + 1) - dat*(-log(-p + 1))^(2*xi))*xi^2)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + (dat*xi*(-log(-p + 1))^(2*xi) - xi*z*(-log(-p + 1))^(2*xi) + sigma*(-log(-p + 1))^xi)*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi)))/((dat^2*xi^4*(-log(-p + 1))^(2*xi) + xi^4*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*xi^3*(-log(-p + 1))^xi + sigma^2*xi^2 - 2*(dat*xi^4*(-log(-p + 1))^(2*xi) + sigma*xi^3*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))
  )%*%V
}


#' @title Generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized Pareto distribution parametrized in terms of average maximum of \code{N} exceedances.
#'
#' The parameter \code{N} corresponds to the number of threshold exceedances of interest over which the maxima is taken.
#' \eqn{z} is the corresponding expected value of this block maxima.
#' Note that the actual parametrization is in terms of excess expected mean, meaning expected mean minus threshold.
#'
#' @details The observed information matrix was calculated from the Hessian using symbolic calculus in Sage.
#'
#' @author Leo Belzile
#' @name gpdN
#' @param par vector of length 2 containing \eqn{z} and \eqn{\xi}, respectively the mean excess of the maxima of N exceedances above the threshold and the shape parameter.
#' @param dat sample vector
#' @param N block size for threshold exceedances.
#' @param tol numerical tolerance for the exponential model
#' @param V vector calculated by \code{gpdN.Vfun}
#' @section Usage: \preformatted{gpdN.ll(par, dat, N, tol=1e-5)
#' gpdN.score(par, dat, N)
#' gpdN.infomat(par, dat, N, method = c("obs", "exp"), nobs = length(dat))
#' gpdN.Vfun(par, dat, N)
#' gpdN.phi(par, dat, N, V)
#' gpdN.dphi(par, dat, N, V)}
#'
#' @section Functions:
#'
#' \itemize{
#' \item{\code{gpdN.ll}:} {log likelihood}
#' \item{\code{gpdN.score}:} {score vector}
#' \item{\code{gpdN.infomat}:} {observed information matrix for GP parametrized in terms of mean of the maximum of \code{N} exceedances and shape}
#' \item{\code{gpdN.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gpdN.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gpdN.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
NULL

#' Negative log likelihood of the generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#'
#' @seealso \code{\link{gpdN}}
#' @inheritParams gpdN
#' @keywords internal
#' @export
gpdN.ll <-  function(par, dat, N){
  xi = par[2]; z = par[1]
  sigma = z*xi/(exp(lgamma(N + 1) + lgamma(-xi + 1)-lgamma(N - xi + 1))-1)
  gpd.ll(par = c(sigma, xi), dat = dat)
}

#' Score vector of the generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#' @seealso \code{\link{gpdN}}
#' @inheritParams gpdN
#' @keywords internal
#' @export
gpdN.score <- function(par, dat, N){
  z = par[1]; xi = par[2]
  cst <- exp(lgamma(N + 1) + lgamma(1-xi)-lgamma(N + 1-xi))
  c(sum(dat*(cst - 1)*(1/xi + 1)/(z^2*(dat*(cst - 1)/z + 1)) - 1/z),
    sum(-(psigamma(N - xi + 1)*cst - psigamma(-xi + 1)*cst)*dat*(1/xi + 1)/(z*(dat*(cst - 1)/z + 1)) +
          ((psigamma(N - xi + 1)*cst - psigamma(-xi + 1)*cst)*xi*z/(cst - 1)^2 - z/(cst - 1))*(cst - 1)/(xi*z) + log(dat*(cst - 1)/z + 1)/xi^2))
}

#' Information matrix of the generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#' @seealso \code{\link{gpdN}}
#' @inheritParams gpdN
#' @keywords internal
#' @export
gpdN.infomat <- function(par, dat, N, method = c("obs", "exp"), nobs = length(dat)){
  z = par[1]; xi = par[2]
  if(xi >= 1 || xi <= -0.5){ return(matrix(NA,2,2))}
  method <- method[1] #default corresponds to observed information rather than Fisher information
  #Fisher information matrix
  cst <- exp(lgamma(N + 1) + lgamma(1-xi)-lgamma(N + 1-xi))
  if(method == "exp"){
    sigmaf <- z*xi/(cst-1)
    Jac <- rbind(c(xi/(cst - 1),
                   -(digamma(N - xi + 1)*cst - digamma(-xi + 1)*cst)*xi*z/(cst - 1)^2 + z/(cst - 1)),
                 c(0, 1))
    return(t(Jac) %*% gpd.infomat(par = c(sigmaf, xi), dat = dat, method = "exp") %*% Jac)

  } else if (method == "obs"){
    #Observed information
    k11 <- sum(2*dat*(cst - 1)*(1/xi + 1)/(z^3*(dat*(cst - 1)/z + 1)) - dat^2*(cst - 1)^2*(1/xi + 1)/(z^4*(dat*(cst - 1)/z + 1)^2) - 1/z^2)
    k12 <- sum(-(digamma(N - xi + 1)*cst - digamma(-xi + 1)*cst)*dat*(1/xi + 1)/(z^2*(dat*(cst - 1)/z + 1)) + (digamma(N - xi + 1)*cst - digamma(-xi + 1)*cst)*dat^2*(cst - 1)*(1/xi + 1)/(z^3*(dat*(cst - 1)/z + 1)^2) + dat*(cst - 1)/(xi^2*z^2*(dat*(cst - 1)/z + 1)))
    k22 <- sum(-(digamma(N - xi + 1)*cst - digamma(-xi + 1)*cst)^2*dat^2*(1/xi + 1)/(z^2*(dat*(cst - 1)/z + 1)^2) + (digamma(N - xi + 1)^2*cst - 2*digamma(N - xi + 1)*digamma(-xi + 1)*cst + digamma(-xi + 1)^2*cst - psigamma(deriv=1, N - xi + 1)*cst + psigamma(deriv=1, -xi + 1)*cst)*dat*(1/xi + 1)/(z*(dat*(cst - 1)/z + 1)) - (digamma(N - xi + 1)*cst - digamma(-xi + 1)*cst)*((digamma(N - xi + 1)*cst - digamma(-xi + 1)*cst)*xi*z/(cst - 1)^2 - z/(cst - 1))/(xi*z) + (2*(digamma(N - xi + 1)*cst - digamma(-xi + 1)*cst)^2*xi*z/(cst - 1)^3 - (digamma(N - xi + 1)^2*cst - 2*digamma(N - xi + 1)*digamma(-xi + 1)*cst + digamma(-xi + 1)^2*cst - psigamma(deriv=1, N - xi + 1)*cst + psigamma(deriv=1, -xi + 1)*cst)*xi*z/(cst - 1)^2 - 2*(digamma(N - xi + 1)*cst - digamma(-xi + 1)*cst)*z/(cst - 1)^2)*(cst - 1)/(xi*z) + ((digamma(N - xi + 1)*cst - digamma(-xi + 1)*cst)*xi*z/(cst - 1)^2 - z/(cst - 1))*(cst - 1)/(xi^2*z) - 2*(digamma(N - xi + 1)*cst - digamma(-xi + 1)*cst)*dat/(xi^2*z*(dat*(cst - 1)/z + 1)) + 2*log(dat*(cst - 1)/z + 1)/xi^3)
    return(cbind(c(k11,k12),c(k12,k22)))
  }
}

#' Tangent exponential model statistics for the generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#'
#' Vector implementing conditioning on approximate ancillary statistics for the TEM
#' @seealso \code{\link{gpdN}}
#' @name gpdN.temstat
#' @inheritParams gpdN
#' @keywords internal
#' @export
gpdN.Vfun <- function(par, dat, N){
  xi = par[2]; z = par[1]
  cst <- exp(lgamma(N + 1) + lgamma(-xi + 1)-lgamma(N - xi + 1))
  cbind(dat/z,-xi*z*(dat*(cst - 1)/z + 1)^(-1/xi)*((psigamma(N - xi + 1)*cst - psigamma(-xi + 1)*cst)*dat/(xi*z*(dat*(cst - 1)/z + 1)) - log(dat*(cst - 1)/z + 1)/xi^2)/((dat*(cst - 1)/z + 1)^(-1/xi - 1)*(cst - 1)))
}

#' Canonical parameter in the local exponential family approximation
#'
#' @inheritParams gpdN
#' @rdname gpdN.temstat
#' @keywords internal
#' @export
gpdN.phi <- function(par, dat, N, V){
  xi = par[2]; z = par[1]
  cst <- exp(lgamma(N + 1) + lgamma(-xi + 1)-lgamma(N - xi + 1))
  matrix(-(cst - 1)*(1/xi + 1)/(z*(dat*(cst - 1)/z + 1)), nrow=1) %*% V
}

## Derivative matrix of the canonical parameter in the local exponential family approximation
#' @inheritParams gpdN
#' @rdname gpdN.temstat
#' @keywords internal
#' @export
gpdN.dphi <- function(par, dat, N, V){
  xi = par[2]; z = par[1]
  cst <- exp(lgamma(N + 1) + lgamma(-xi + 1)-lgamma(N - xi + 1))
  rbind((cst - 1)*(1/xi + 1)/(z^2*(dat*(cst - 1)/z + 1)) - dat*(cst - 1)^2*(1/xi + 1)/(z^3*(dat*(cst - 1)/z + 1)^2),
        -(psigamma(N - xi + 1)*cst - psigamma(-xi + 1)*cst)*(1/xi + 1)/(z*(dat*(cst - 1)/z + 1)) + (psigamma(N - xi + 1)*cst - psigamma(-xi + 1)*cst)*dat*(cst - 1)*(1/xi + 1)/(z^2*(dat*(cst - 1)/z + 1)^2) + (cst - 1)/(xi^2*z*(dat*(cst - 1)/z + 1))
  ) %*% V
}



############################################################################################

############################################################################################

#' @title Generalized extreme value distribution (quantile/mean of N-block maxima parametrization)
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized extreme value distribution  parametrized in terms of the
#' quantiles/mean of N-block maxima parametrization \eqn{z}, scale and shape.
#'
#' @author Leo Belzile
#' @name gevN
#' @param par vector of \code{loc}, quantile/mean of N-block maximum and \code{shape}
#' @param dat sample vector
#' @param V vector calculated by \code{gevN.Vfun}
#' @param q probability, corresponding to \eqn{q}th quantile of the \code{N}-block maximum
#' @param qty string indicating whether to calculate the \code{q} quantile or the mean
#' @section Usage: \preformatted{gevN.ll(par, dat, N, q, qty = c("mean", "quantile"))
#' gevN.ll.optim(par, dat, N, q = 0.5, qty = c("mean", "quantile"))
#' gevN.score(par, dat, N, q = 0.5, qty = c("mean", "quantile"))
#' gevN.infomat(par, dat, qty = c("mean", "quantile"), method = c("obs", "exp"), N, q = 0.5, nobs = length(dat))
#' gevN.Vfun(par, dat, N, q = 0.5, qty = c("mean", "quantile"))
#' gevN.phi(par, dat, N, q = 0.5, qty = c("mean", "quantile"), V)
#' gevN.dphi(par, dat, N, q = 0.5, qty = c("mean", "quantile"), V)}
#'
#' @section Functions:
#' \itemize{
#' \item{\code{gevN.ll}:} {log likelihood}
#' \item{\code{gevN.score}:} {score vector}
#' \item{\code{gevN.infomat}:} {expected and observed information matrix}
#' \item{\code{gevN.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gevN.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gevN.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
NULL

#' Negative log likelihood of the generalized extreme value distribution (quantile/mean of N-block maxima parametrization)
#'
#' @seealso \code{\link{gevN}}
#' @inheritParams gevN
#' @keywords internal
#' @export
gevN.ll <- function(par, dat, N, q = 0.5, qty = c("mean", "quantile")){
  qty <- match.arg(qty, c("mean", "quantile"))
  mu = par[1]; z = par[2]; xi = as.vector(par[3])
  if(!isTRUE(all.equal(xi, 0))){
    sigma <- switch(qty,
                  quantile = (z-mu)*xi/(N^xi*(log(1/q))^(-xi)-1),
                  mean = (z-mu)*xi/(N^xi*gamma(1-xi)-1))
  } else{
    sigma <- switch(qty,
                    quantile = (z-mu)/(log(N) - log(-log(q))),
                    mean = (z-mu)/(-psigamma(1)+log(N)))
  }
  gev.ll(par = c(mu, sigma, xi), dat = dat)
}

#' Score of the generalized extreme value distribution (quantile/mean of N-block maxima parametrization)
#'
#' @seealso \code{\link{gevN}}
#' @inheritParams gevN
#' @keywords internal
#' @export
gevN.score <- function(par, dat, N, q = 0.5, qty = c("mean", "quantile")){
  qty <- match.arg(qty, c("mean", "quantile"))
  mu = par[1]; z = par[2]; xi = par[3];
  if(qty == "quantile"){ #quantiles at prob. q
    c(sum((-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*log(1/q)^(-xi) - 1)/(mu - z))/xi + ((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*log(1/q)^(-xi) - 1)/(mu - z))*(1/xi + 1)/((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1) - 1/(mu - z)),
      sum(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)*(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)^2*xi) - (N^xi*log(1/q)^(-xi) - 1)*(dat - mu)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)^2) + 1/(mu - z)),
      sum((-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi)*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(dat - mu)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)*xi) - log(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)/xi^2) - (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(dat - mu)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)) + (N^xi*log(1/q)^(-xi) - 1)*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - z)*xi/(N^xi*log(1/q)^(-xi) - 1)^2 - (mu - z)/(N^xi*log(1/q)^(-xi) - 1))/((mu - z)*xi) + log(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)/xi^2
      ))
  } else { #Mean
    c(sum((-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*gamma(-xi + 1) - 1)/(mu - z))/xi + ((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*gamma(-xi + 1) - 1)/(mu - z))*(1/xi + 1)/((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) - 1) - 1/(mu - z)),
      sum(-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)*(-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)^2*xi) - (N^xi*gamma(-xi + 1) - 1)*(dat - mu)*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)^2) + 1/(mu - z)),
      sum((-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi)*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(dat - mu)/(((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)*xi) - log(-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)/xi^2) - (N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(dat - mu)*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)) + (N^xi*gamma(-xi + 1) - 1)*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - z)*xi/(N^xi*gamma(-xi + 1) - 1)^2 - (mu - z)/(N^xi*gamma(-xi + 1) - 1))/((mu - z)*xi) + log(-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)/xi^2))

  }
}

#' Information matrix of the generalized extreme value distribution (quantile/mean of N-block maxima parametrization)
#' @seealso \code{\link{gevN}}
#' @inheritParams gevN
#' @keywords internal
#' @export
gevN.infomat <- function(par, dat, method = c("obs", "exp"), qty = c("mean", "quantile"), N, q = 0.5, nobs = length(dat)){
  mu = par[1]; z = par[2]; xi = par[3]
  qty <- match.arg(qty, c("mean", "quantile"))[1]
  xizero <- isTRUE(all.equal(xi,0))
  if(xizero && method == "obs"){
    stop("Unimplemented for the case xi=0 for observed information")
  }
  eulergamma <- 0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917

  if(qty == "mean"){
    if(!xizero){
      #z = mu + sigma/xi*(N^xi*gamma(1-xi)-1)
      sigmaq <- (z-mu)*xi/(N^xi*gamma(1-xi)-1)
    } else {
      #z = mu + sigma*(log(N) + eulergamma)
      sigmaq <- (z-mu)/(log(N) + eulergamma)
    }

    if(method == "obs"){
      k11 <- sum(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 2)*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*gamma(-xi + 1) - 1)/(mu - z))^2*(1/xi + 1)/xi - 2*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^3 - (N^xi*gamma(-xi + 1) - 1)/(mu - z)^2)/xi - ((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*gamma(-xi + 1) - 1)/(mu - z))^2*(1/xi + 1)/((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^2 + 2*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^3 - (N^xi*gamma(-xi + 1) - 1)/(mu - z)^2)*(1/xi + 1)/((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1) - 1/(mu - z)^2)
      k12 <- sum(-(N^xi*gamma(-xi + 1) - 1)*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 2)*(mu - dat)*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*gamma(-xi + 1) - 1)/(mu - z))*(1/xi + 1)/((mu - z)^2*xi) + ((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*(2*(N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^3 - (N^xi*gamma(-xi + 1) - 1)/(mu - z)^2)/xi - (2*(N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^3 - (N^xi*gamma(-xi + 1) - 1)/(mu - z)^2)*(1/xi + 1)/((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1) + (N^xi*gamma(-xi + 1) - 1)*(mu - dat)*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*gamma(-xi + 1) - 1)/(mu - z))*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)^2) + 1/(mu - z)^2)
      k13 <- sum(-((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - dat)*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)) - log((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)/xi^2)*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*gamma(-xi + 1) - 1)/(mu - z))/xi + ((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - dat)/(mu - z)^2 - (N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))/(mu - z))/xi - ((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - dat)/(mu - z)^2 - (N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))/(mu - z))*(1/xi + 1)/((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1) + (N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - dat)*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*gamma(-xi + 1) - 1)/(mu - z))*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)) - ((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*gamma(-xi + 1) - 1)/(mu - z))/xi^2 + ((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*gamma(-xi + 1) - 1)/(mu - z))/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)*xi^2))
      k22 <- sum((N^xi*gamma(-xi + 1) - 1)^2*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 2)*(mu - dat)^2*(1/xi + 1)/((mu - z)^4*xi) - 2*(N^xi*gamma(-xi + 1) - 1)*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*(mu - dat)/((mu - z)^3*xi) - (N^xi*gamma(-xi + 1) - 1)^2*(mu - dat)^2*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)^4) + 2*(N^xi*gamma(-xi + 1) - 1)*(mu - dat)*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)^3) - 1/(mu - z)^2)
      k23 <- sum((N^xi*gamma(-xi + 1) - 1)*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*(mu - dat)*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - dat)*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)) - log((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)/xi^2)/((mu - z)^2*xi) - (N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*(mu - dat)/((mu - z)^2*xi) - (N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(N^xi*gamma(-xi + 1) - 1)*(mu - dat)^2*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)^3) + (N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - dat)*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)^2) + (N^xi*gamma(-xi + 1) - 1)*((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*(mu - dat)/((mu - z)^2*xi^2) - (N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)^2*xi^2))
      k33 <- sum(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi)*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - dat)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)*xi) - log((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)/xi^2)^2 + ((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi)*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))^2*(mu - dat)^2/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)^2*xi) - (N^xi*log(N)^2*gamma(-xi + 1) - 2*N^xi*log(N)*psigamma(-xi + 1)*gamma(-xi + 1) + N^xi*psigamma(-xi + 1)^2*gamma(-xi + 1) + N^xi*psigamma(1, -xi + 1)*gamma(-xi + 1))*(mu - dat)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)*xi) + 2*(N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - dat)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)*xi^2) - 2*log((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)/xi^3) - (N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))^2*(mu - dat)^2*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)^2) + (N^xi*log(N)^2*gamma(-xi + 1) - 2*N^xi*log(N)*psigamma(-xi + 1)*gamma(-xi + 1) + N^xi*psigamma(-xi + 1)^2*gamma(-xi + 1) + N^xi*psigamma(1, -xi + 1)*gamma(-xi + 1))*(mu - dat)*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)) + (N^xi*gamma(-xi + 1) - 1)*(2*(N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))^2*(mu - z)*xi/(N^xi*gamma(-xi + 1) - 1)^3 - (N^xi*log(N)^2*gamma(-xi + 1) - 2*N^xi*log(N)*psigamma(-xi + 1)*gamma(-xi + 1) + N^xi*psigamma(-xi + 1)^2*gamma(-xi + 1) + N^xi*psigamma(1, -xi + 1)*gamma(-xi + 1))*(mu - z)*xi/(N^xi*gamma(-xi + 1) - 1)^2 - 2*(N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - z)/(N^xi*gamma(-xi + 1) - 1)^2)/((mu - z)*xi) - (N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - z)*xi/(N^xi*gamma(-xi + 1) - 1)^2 - (mu - z)/(N^xi*gamma(-xi + 1) - 1))/((mu - z)*xi) + (N^xi*gamma(-xi + 1) - 1)*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - z)*xi/(N^xi*gamma(-xi + 1) - 1)^2 - (mu - z)/(N^xi*gamma(-xi + 1) - 1))/((mu - z)*xi^2) - 2*(N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - dat)/(((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)*xi^2) + 2*log((N^xi*gamma(-xi + 1) - 1)*(mu - dat)/(mu - z) + 1)/xi^3)
      return( cbind(c(k11,k12,k13),c(k12,k22,k23),c(k13,k23,k33)))
    } else if (method == "exp"){
      if(!xizero){
        Jac <- rbind(c(1,0,0),
                     c(-xi/(N^xi*gamma(-xi + 1) - 1),
                       xi/(N^xi*gamma(-xi + 1) - 1),
                       (N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(mu - z)*xi/(N^xi*gamma(-xi + 1) - 1)^2 - (mu - z)/(N^xi*gamma(-xi + 1) - 1)),
                     c(0, 0, 1))
        return(t(Jac) %*% gev.infomat(par = c(mu, sigmaq, xi), method = "exp", nobs = nobs) %*% Jac)
      } else{
        Jac <- rbind(c(1,0,0),
                     c(-1/(eulergamma + log(N)),
                       1/(eulergamma + log(N)),
                       1/12*(6*eulergamma^2 + pi^2 + 12*eulergamma*log(N) + 6*log(N)^2)*(mu - z)/(eulergamma + log(N))^2),
                     c(0, 0, 1))

        return(t(Jac) %*% gev.infomat(par = c(mu, sigmaq, 0), method = "exp", nobs = nobs) %*% Jac)

      }
    }
  } else if (qty == "quantile"){
    if(!xizero){
      #z = mu + sigma/xi*(N^xi*log(1/q)^(-xi)-1)
      sigmaq <- (z - mu)*xi/(N^xi*(log(1/q))^(-xi)-1)
    } else{
      #z = mu + sigma*log(N/log(1/q))
      sigmaq <- (z - mu)/log(N/log(1/q))
    }
    #Quantiles, observed information
    if(method == "obs"){
      k11 <- sum(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 2)*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z))^2*(1/xi + 1)/xi - 2*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^3 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z)^2)/xi - ((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z))^2*(1/xi + 1)/((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^2 + 2*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^3 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z)^2)*(1/xi + 1)/((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1) - 1/(mu - z)^2)
      k12 <- sum(-(N^xi*log(1/q)^(-xi) - 1)*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 2)*(mu - dat)*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z))*(1/xi + 1)/((mu - z)^2*xi) + ((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*(2*(N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^3 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z)^2)/xi - (2*(N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^3 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z)^2)*(1/xi + 1)/((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1) + (N^xi*log(1/q)^(-xi) - 1)*(mu - dat)*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z))*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)^2) + 1/(mu - z)^2)
      k13 <- sum(-((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - dat)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)) - log((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)/xi^2)*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z))/xi + ((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - dat)/(mu - z)^2 - (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))/(mu - z))/xi - ((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - dat)/(mu - z)^2 - (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))/(mu - z))*(1/xi + 1)/((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1) + (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - dat)*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z))*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)) - ((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z))/xi^2 + ((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z)^2 - (N^xi*log(1/q)^(-xi) - 1)/(mu - z))/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)*xi^2))
      k22 <- sum((N^xi*log(1/q)^(-xi) - 1)^2*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 2)*(mu - dat)^2*(1/xi + 1)/((mu - z)^4*xi) - 2*(N^xi*log(1/q)^(-xi) - 1)*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*(mu - dat)/((mu - z)^3*xi) - (N^xi*log(1/q)^(-xi) - 1)^2*(mu - dat)^2*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)^4) + 2*(N^xi*log(1/q)^(-xi) - 1)*(mu - dat)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)^3) - 1/(mu - z)^2)
      k23 <- sum((N^xi*log(1/q)^(-xi) - 1)*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*(mu - dat)*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - dat)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)) - log((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)/xi^2)/((mu - z)^2*xi) - (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*(mu - dat)/((mu - z)^2*xi) - (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(N^xi*log(1/q)^(-xi) - 1)*(mu - dat)^2*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)^3) + (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - dat)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)^2) + (N^xi*log(1/q)^(-xi) - 1)*((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi - 1)*(mu - dat)/((mu - z)^2*xi^2) - (N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)^2*xi^2))
      k33 <- sum(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi)*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - dat)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)*xi) - log((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)/xi^2)^2 + ((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^(-1/xi)*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))^2*(mu - dat)^2/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)^2*xi) - (N^xi*log(1/q)^(-xi)*log(N)^2 - 2*N^xi*log(1/q)^(-xi)*log(N)*log(log(1/q)) + N^xi*log(1/q)^(-xi)*log(log(1/q))^2)*(mu - dat)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)*xi) + 2*(N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - dat)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)*xi^2) - 2*log((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)/xi^3) - (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))^2*(mu - dat)^2*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)^2*(mu - z)^2) + (N^xi*log(1/q)^(-xi)*log(N)^2 - 2*N^xi*log(1/q)^(-xi)*log(N)*log(log(1/q)) + N^xi*log(1/q)^(-xi)*log(log(1/q))^2)*(mu - dat)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)) + (N^xi*log(1/q)^(-xi) - 1)*(2*(N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))^2*(mu - z)*xi/(N^xi*log(1/q)^(-xi) - 1)^3 - (N^xi*log(1/q)^(-xi)*log(N)^2 - 2*N^xi*log(1/q)^(-xi)*log(N)*log(log(1/q)) + N^xi*log(1/q)^(-xi)*log(log(1/q))^2)*(mu - z)*xi/(N^xi*log(1/q)^(-xi) - 1)^2 - 2*(N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - z)/(N^xi*log(1/q)^(-xi) - 1)^2)/((mu - z)*xi) - (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - z)*xi/(N^xi*log(1/q)^(-xi) - 1)^2 - (mu - z)/(N^xi*log(1/q)^(-xi) - 1))/((mu - z)*xi) + (N^xi*log(1/q)^(-xi) - 1)*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - z)*xi/(N^xi*log(1/q)^(-xi) - 1)^2 - (mu - z)/(N^xi*log(1/q)^(-xi) - 1))/((mu - z)*xi^2) - 2*(N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - dat)/(((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)*(mu - z)*xi^2) + 2*log((N^xi*log(1/q)^(-xi) - 1)*(mu - dat)/(mu - z) + 1)/xi^3)
      return( cbind(c(k11,k12,k13),c(k12,k22,k23),c(k13,k23,k33)))

    } else if(method == "exp"){
      if(!xizero){
        Jac <- rbind(c(1, 0, 0),
                     c(-xi/(N^xi*log(1/q)^(-xi) - 1),
                       xi/(N^xi*log(1/q)^(-xi) - 1),
                       (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(mu - z)*xi/(N^xi*log(1/q)^(-xi) - 1)^2 - (mu - z)/(N^xi*log(1/q)^(-xi) - 1)),
                     c(0, 0, 1))
        return(t(Jac) %*% gev.infomat(par = c(mu, sigmaq, xi), method = "exp", nobs = nobs) %*% Jac )
      } else{
        Jac <- rbind(c(1, 0, 0), c(-1/(log(N) - log(-log(q))), 1/(log(N) - log(-log(q))), 1/2*(mu - z)), c(0, 0, 1))
        return(t(Jac) %*% gev.infomat(par = c(mu, sigmaq, 0), method = "exp", nobs = nobs) %*% Jac)
      }
    }
  }

}



#' Tangent exponential model statistics for the generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#'
#' Vector implementing conditioning on approximate ancillary statistics for the TEM
#' @seealso \code{\link{gevN}}
#' @name gevN.temstat
#' @inheritParams gevN
#' @keywords internal
#' @export
gevN.Vfun <- function(par, dat, N, q = 0.5, qty = c("mean", "quantile")){
  #Quantiles, profiling by substituting sigma by zq
  qty <- match.arg(qty, c("mean", "quantile"))
  mu = par[1]; z = par[2]; xi = par[3];
  if(qty == "quantile"){
    cbind((mu - z)*((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*log(1/q)^(-xi) - 1)/(mu - z))/(N^xi*log(1/q)^(-xi) - 1),
          -(dat - mu)/(mu - z),
          (-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi)*(mu - z)*xi*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(dat - mu)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)*xi) - log(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)/xi^2)/((N^xi*log(1/q)^(-xi) - 1)*(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1))
    )
  } else{
    cbind((mu - z)*((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*gamma(-xi + 1) - 1)/(mu - z))/(N^xi*gamma(-xi + 1) - 1),
          -(dat - mu)/(mu - z),
          (-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi)*(mu - z)*xi*((N^xi*log(N)*gamma(-xi + 1) - N^xi*psigamma(-xi + 1)*gamma(-xi + 1))*(dat - mu)/(((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)*xi) - log(-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)/xi^2)/((N^xi*gamma(-xi + 1) - 1)*(-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1))
    )
  }
}

#' Canonical parameter in the local exponential family approximation
#'
#' @inheritParams gevN
#' @rdname gevN.temstat
#' @keywords internal
#' @export
gevN.phi <- function(par, dat, N, q = 0.5, qty = c("mean", "quantile"), V){
  qty <- match.arg(qty, c("mean", "quantile"))
  mu = par[1]; z = par[2]; xi = par[3];
  if(qty == "mean"){
    matrix(-(N^xi*gamma(-xi + 1) - 1)*(-(N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)*xi) - (N^xi*gamma(-xi + 1) - 1)*(1/xi + 1)/(((N^xi*gamma(-xi + 1) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)), nrow = 1) %*% V
  } else{
    matrix(-(N^xi*log(1/q)^(-xi) - 1)*(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)*xi) - (N^xi*log(1/q)^(-xi) - 1)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)), nrow = 1) %*% V
  }
}

  ## Derivative matrix of the canonical parameter in the local exponential family approximation
  #' @inheritParams gevN
  #' @rdname gevN.temstat
  #' @keywords internal
  #' @export
  gevN.dphi <- function(par, dat, N, q = 0.5, qty = c("mean", "quantile"), V){
    qty <- match.arg(qty, c("mean", "quantile"))
    mu = par[1]; z = par[2]; xi = par[3];
    if(qty == "mean"){
      rbind(
        -(N^xi*mu*xi^2*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*gamma(-xi + 1) - N^xi*xi^2*z*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*gamma(-xi + 1) + N^xi*mu*xi*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*gamma(-xi + 1) - N^xi*xi*z*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*gamma(-xi + 1) - N^xi*mu*xi*gamma(-xi + 1) + N^xi*xi*z*gamma(-xi + 1) - N^xi*dat*gamma(-xi + 1) + N^xi*z*gamma(-xi + 1) + dat - z)*(N^xi*gamma(-xi + 1) - 1)/((N^xi*dat*gamma(-xi + 1) - N^xi*mu*gamma(-xi + 1) - dat + z)^2*(mu - z)*xi^2*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)),
        (mu*xi^2*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi) - xi^2*z*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi) + mu*xi*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi) - xi*z*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi) - N^xi*dat*gamma(-xi + 1) + N^xi*mu*gamma(-xi + 1) - mu*xi + xi*z + dat - mu)*(N^xi*gamma(-xi + 1) - 1)/((N^xi*dat*gamma(-xi + 1) - N^xi*mu*gamma(-xi + 1) - dat + z)^2*(mu - z)*xi^2*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)),
        (N^xi*mu*xi^3*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*log(N)*gamma(-xi + 1) - N^xi*xi^3*z*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*log(N)*gamma(-xi + 1) - N^xi*mu*xi^3*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*psigamma(-xi + 1)*gamma(-xi + 1) + N^xi*xi^3*z*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*psigamma(-xi + 1)*gamma(-xi + 1) + N^xi*mu*xi^2*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*log(N)*gamma(-xi + 1) - N^xi*xi^2*z*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*log(N)*gamma(-xi + 1) - N^xi*mu*xi^2*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*psigamma(-xi + 1)*gamma(-xi + 1) + N^xi*xi^2*z*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*psigamma(-xi + 1)*gamma(-xi + 1) - N^xi*mu*xi^2*log(N)*gamma(-xi + 1) + N^xi*xi^2*z*log(N)*gamma(-xi + 1) + N^xi*mu*xi^2*psigamma(-xi + 1)*gamma(-xi + 1) - N^xi*xi^2*z*psigamma(-xi + 1)*gamma(-xi + 1) + N^(2*xi)*dat*xi*((N^xi*mu*gamma(-xi + 1) -
(N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*gamma(-xi + 1)^2 - N^(2*xi)*mu*xi*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*gamma(-xi + 1)^2 - N^(2*xi)*dat*xi*log(N)*gamma(-xi + 1)^2 + N^(2*xi)*mu*xi*log(N)*gamma(-xi + 1)^2 + N^(2*xi)*dat*xi*psigamma(-xi + 1)*gamma(-xi + 1)^2 - N^(2*xi)*mu*xi*psigamma(-xi + 1)*gamma(-xi + 1)^2 - 2*N^xi*dat*xi*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*gamma(-xi + 1) + N^xi*mu*xi*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*gamma(-xi + 1) + N^xi*xi*z*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi)*gamma(-xi + 1) + N^xi*dat*xi*log(N)*gamma(-xi + 1) - N^xi*mu*xi*log(N)*gamma(-xi + 1) - N^xi*dat*xi*psigamma(-xi + 1)*gamma(-xi + 1) + N^xi*mu*xi*psigamma(-xi + 1)*gamma(-xi + 1) - N^(2*xi)*dat*xi*gamma(-xi + 1)^2 + N^(2*xi)*mu*xi*gamma(-xi + 1)^2 + N^(2*xi)*dat*log((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))*gamma(-xi + 1)^2 - N^(2*xi)*mu*log((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))*gamma(-xi + 1)^2 + 2*N^xi*dat*xi*gamma(-xi + 1) - N^xi*mu*xi*gamma(-xi + 1) - N^xi*xi*z*gamma(-xi + 1) - 2*N^xi*dat*log((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))*gamma(-xi + 1) + N^xi*mu*log((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))*gamma(-xi + 1) + N^xi*z*log((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))*gamma(-xi + 1) + dat*xi*((N^xi*mu*gamma(-xi + 1) -
(N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi) - xi*z*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi) - dat*xi + xi*z + dat*log((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z)) - z*log((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z)))/((N^xi*dat*gamma(-xi + 1) - N^xi*mu*gamma(-xi + 1) - dat + z)^2*xi^3*((N^xi*mu*gamma(-xi + 1) - (N^xi*gamma(-xi + 1) - 1)*dat - z)/(mu - z))^(1/xi))) %*% V
    } else{
      rbind((N^xi*log(1/q)^(-xi) - 1)*(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 2)*((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*log(1/q)^(-xi) - 1)/(mu - z))*(1/xi + 1)/((mu - z)*xi) - (N^xi*log(1/q)^(-xi) - 1)*((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z)^2 + (N^xi*log(1/q)^(-xi) - 1)/(mu - z))*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)^2*(mu - z)) + (N^xi*log(1/q)^(-xi) - 1)*(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)^2*xi) + (N^xi*log(1/q)^(-xi) - 1)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)^2),
            -(N^xi*log(1/q)^(-xi) - 1)^2*(dat - mu)*(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 2)*(1/xi + 1)/((mu - z)^3*xi) - (N^xi*log(1/q)^(-xi) - 1)*(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)^2*xi) + (N^xi*log(1/q)^(-xi) - 1)^2*(dat - mu)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)^2*(mu - z)^3) - (N^xi*log(1/q)^(-xi) - 1)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)^2),
            (N^xi*log(1/q)^(-xi) - 1)*(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)*((N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(dat - mu)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)) - log(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)/xi^2)/((mu - z)*xi) - (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)*xi) + (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)^2*(mu - z)^2) - (N^xi*log(1/q)^(-xi)*log(N) - N^xi*log(1/q)^(-xi)*log(log(1/q)))*(1/xi + 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)) + (N^xi*log(1/q)^(-xi) - 1)*(-(N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)*xi^2) + (N^xi*log(1/q)^(-xi) - 1)/(((N^xi*log(1/q)^(-xi) - 1)*(dat - mu)/(mu - z) - 1)*(mu - z)*xi^2)) %*% V
    }
  }
