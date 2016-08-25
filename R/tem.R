#' @title Generalized Pareto distribution (expected shortfall parametrization)
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized Pareto distribution parametrized in terms of expected shortfall.
#'
#' @details The observed information matrix was calculated using the second Bartlett identity as the negative of the second derivative of the log-likelihood function, using symbolic calculus in Sage.
#' @details The interpretation for \code{m} is as follows: if there are on average \eqn{m_y} observations per year above the threshold, then  \eqn{m=Tm_y} corresponds to \eqn{T}-year return level.
#'
#' @author Leo Belzile
#' @name gpde
#' @param par vector of length 2 containing \eqn{y_m} and \eqn{\xi}, respectively the \eqn{m}-year return level and the shape parameter.
#' @param dat sample vector
#' @param m number of observations of interest for return levels. See \strong{Details}
#' @param tol numerical tolerance for the exponential model
#' @param V vector calculated by \code{gpde.Vfun}
#' @section Usage: \preformatted{gpde.ll(par, dat, m, tol=1e-5)
#' gpde.ll.optim(par, dat, m, tol=1e-5)
#' gpde.score(par, dat, m)
#' gpde.infomat(par, dat, m
#' gpde.Vfun(par, dat, m)
#' gpde.phi(par, dat, V, m)
#' gpde.dphi(par, dat, V, m)}
#'
#' @section Functions:
#'
#' \itemize{
#' \item{\code{gpde.ll}:} {log-likelihood}
#' \item{\code{gpde.ll.optim}:} {negative log-likelihood parametrized in terms of log expected
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
#' @details The observed information matrix was calculated using the second Bartlett identity as the negative of the second derivative of the log-likelihood function, using symbolic calculus in Sage.
#' @details The interpretation for \code{m} is as follows: if there are on average \eqn{m_y} observations per year above the threshold, then  \eqn{m=Tm_y} corresponds to \eqn{T}-year return level.
#'
#' @author Leo Belzile
#' @name gpdr
#' @param par vector of length 2 containing \eqn{y_m} and \eqn{\xi}, respectively the \eqn{m}-year return level and the shape parameter.
#' @param dat sample vector
#' @param m number of observations of interest for return levels. See \strong{Details}
#' @param tol numerical tolerance for the exponential model
#' @param V vector calculated by \code{gpdr.Vfun}
#'
#' @section Usage: \preformatted{gpdr.ll(par, dat, m, tol=1e-5)
#' gpdr.ll.optim(par, dat, m, tol=1e-5)
#' gpdr.score(par, dat, m)
#' gpdr.infomat(par, dat, m
#' gpdr.Vfun(par, dat, m)
#' gpdr.phi(par, V, dat, m)
#' gpdr.dphi(par, V, dat, m)}
#' @section Functions:
#'
#' \itemize{
#' \item{\code{gpdr.ll}:} {log-likelihood}
#' \item{\code{gpdr.ll.optim}:} {negative log-likelihood parametrized in terms of \code{log(scale)} and shape
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
#' @param V vector calculated by \code{gevr.Vfun}
#' @param p probability, corresponding to \eqn{(1-p)}th quantile for \eqn{z}
#' @section Usage: \preformatted{gevr.ll(par, dat, p)
#' gevr.ll.optim(par, dat, p)
#' gevr.score(par, dat, p)
#' gevr.infomat(par, dat, p)
#' gevr.Vfun(par, dat, p)
#' gevr.phi(par, dat, p, V)
#' gevr.dphi(par, dat, p, V)}
#'
#' @section Functions:
#' \itemize{
#' \item{\code{gevr.ll}:} {log-likelihood}
#' \item{\code{gevr.ll.optim}:} {negative log-likelihood parametrized in terms of return levels, \code{log(scale)} and shape in order to perform unconstrained optimization}
#' \item{\code{gevr.score}:} {score vector}
#' \item{\code{gevr.infomat}:} {observed information matrix}
#' \item{\code{gevr.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gevr.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gevr.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
NULL

#' Negative log-likelihood of the generalized Pareto distribution (expected shortfall)
#'
#' @seealso \code{\link{gpde}}
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.ll <- function(par, dat, m){
	es = par[1]; xi = par[2]
	if(xi>1){ return(1e10)}
	sum(-log(xi*(1-xi)/(m^xi-1+xi))-log(es)-(1+1/xi)*log(1+(m^xi-1+xi)/((1-xi)*es)*dat))
}


#' The negative log-likelihood is parametrized in terms of log expected shortfall and shape in order to perform unconstrained optimization
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
	if(xi>1){ return(rep(NA,2))}
	c(sum(-1/es + dat*(m^xi + xi - 1)*(1/xi + 1)/(es^2*(xi - 1)*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1))),
		sum(-((m^xi*log(m) + 1)*dat/(es*(xi - 1)) - dat*(m^xi + xi - 1)/(es*(xi - 1)^2))*(1/xi + 1)/(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1) + (m^xi + xi - 1)*((m^xi*log(m) + 1)*(xi - 1)*xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi + xi - 1))/((xi - 1)*xi) + log(-dat*(m^xi + xi - 1)/(es*(xi - 1)) + 1)/xi^2))
}

#' Observed information matrix for the GP distribution (expected shortfall)
#'
#' The information matrix is parametrized in terms of rate of expected shortfall and shape
#' @seealso \code{\link{gpde}}
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.infomat <- function(par, dat, m){
	es = par[1]; xi = par[2]
	if(xi>1){ return(matrix(NA,2,2))}
	k11 = sum(-1/es^2 + 2*dat*(m^xi + xi - 1)*(1/xi + 1)/(es^3*(xi - 1)*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)) - dat^2*(m^xi + xi - 1)^2*(1/xi + 1)/(es^4*(xi - 1)^2*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)^2))
	k22 = sum(-((m^xi*log(m) + 1)*dat/(es*(xi - 1)) - dat*(m^xi + xi - 1)/(es*(xi - 1)^2))^2*(1/xi + 1)/(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)^2 + (dat*m^xi*log(m)^2/(es*(xi - 1)) - 2*(m^xi*log(m) + 1)*dat/(es*(xi - 1)^2) + 2*dat*(m^xi + xi - 1)/(es*(xi - 1)^3))*(1/xi + 1)/(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1) - (m^xi*(xi - 1)*xi*log(m)^2/(m^xi + xi - 1)^2 - 2*(m^xi*log(m) + 1)^2*(xi - 1)*xi/(m^xi + xi - 1)^3 + 2*(m^xi*log(m) + 1)*(xi - 1)/(m^xi + xi - 1)^2 + 2*(m^xi*log(m) + 1)*xi/(m^xi + xi - 1)^2 - 2/(m^xi + xi - 1))*(m^xi + xi - 1)/((xi - 1)*xi) - (m^xi*log(m) + 1)*((m^xi*log(m) + 1)*(xi - 1)*xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi + xi - 1))/((xi - 1)*xi) + (m^xi + xi - 1)*((m^xi*log(m) + 1)*(xi - 1)*xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi + xi - 1))/((xi - 1)*xi^2) + (m^xi + xi - 1)*((m^xi*log(m) + 1)*(xi - 1)*xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi + xi - 1))/((xi - 1)^2*xi) - 2*((m^xi*log(m) + 1)*dat/(es*(xi - 1)) - dat*(m^xi + xi - 1)/(es*(xi - 1)^2))/(xi^2*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)) + 2*log(-dat*(m^xi + xi - 1)/(es*(xi - 1)) + 1)/xi^3)
	k12 = sum(-((m^xi*log(m) + 1)*dat/(es^2*(xi - 1)) - dat*(m^xi + xi - 1)/(es^2*(xi - 1)^2))*(1/xi + 1)/(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1) + dat*(m^xi + xi - 1)*((m^xi*log(m) + 1)*dat/(es*(xi - 1)) - dat*(m^xi + xi - 1)/(es*(xi - 1)^2))*(1/xi + 1)/(es^2*(xi - 1)*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)^2) + dat*(m^xi + xi - 1)/(es^2*(xi - 1)*xi^2*(dat*(m^xi + xi - 1)/(es*(xi - 1)) - 1)))
	cbind(c(k11,k12),c(k12,k22))
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

#' Negative log-likelihood of the generalized Pareto distribution (return levels)
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

#' Negative log-likelihood parametrized in terms of log return level and shape in order to perform unconstrained optimization
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

#' Score of the profile log-likelihood for the GP distribution (return levels parametrization)
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
gpdr.infomat <- function(par, dat, m){
	xi = par[2]; r  = par[1]
	info <- matrix(ncol=2,nrow=2)
	info[1,1] <- sum(dat^2*(m^xi - 1)^2*(1/xi + 1)/((dat*(m^xi - 1)/r + 1)^2*r^4) - 2*dat*(m^xi - 1)*(1/xi + 1)/((dat*(m^xi - 1)/r + 1)*r^3) + 1/r^2)
	info[2,2] <- sum(dat^2*(m^xi)^2*(1/xi + 1)*log(m)^2/((dat*(m^xi - 1)/r + 1)^2*r^2) - dat*m^xi*(1/xi + 1)*log(m)^2/((dat*(m^xi - 1)/r + 1)*r) + m^xi*(m^xi*xi*log(m)/(m^xi - 1)^2 - 1/(m^xi - 1))*log(m)/xi + (m^xi*xi*log(m)^2/(m^xi - 1)^2 - 2*(m^xi)^2*xi*log(m)^2/(m^xi - 1)^3 + 2*m^xi*log(m)/(m^xi - 1)^2)*(m^xi - 1)/xi - (m^xi - 1)*(m^xi*xi*log(m)/(m^xi - 1)^2 - 1/(m^xi - 1))/xi^2 + 2*dat*m^xi*log(m)/((dat*(m^xi - 1)/r + 1)*r*xi^2) - 2*log(dat*(m^xi - 1)/r + 1)/xi^3)
	info[2,1] <- info[1,2] <- sum(-dat^2*(m^xi - 1)*m^xi*(1/xi + 1)*log(m)/((dat*(m^xi - 1)/r + 1)^2*r^3) + dat*m^xi*(1/xi + 1)*log(m)/((dat*(m^xi - 1)/r + 1)*r^2) - dat*(m^xi - 1)/((dat*(m^xi - 1)/r + 1)*r^2*xi^2))
	-info
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


#' Derivative of the canonical parameter \eqn{\phi(\theta)} in the local exponential family approximation
#' @rdname gpdr.temstat
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.dphi <- function(par, dat, V, m){ #TODO fix
	xi = par[2];  r = par[1];  p = m^xi-1
	rbind((1+1/xi)*p/(r+dat*p)^2,
				p/(xi^2*(r+dat*p))-(1+1/xi)*log(m)*m^xi*r/(r+dat*p)^2
	)%*%V

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


###########################################################################################################################

#' Tangent exponential model approximation for the GP distribution
#'
#' The function \code{gpd.tem} provides a tangent exponential model (TEM) approximation
#' for higher order likelihood inference for a scalar parameter for the generalized Pareto distribution. Options include
#' scale and shape parameters as well as value-at-risk (also referred to as quantiles, or return levels)
#' and expected shortfall. The function attempts to find good values for \code{psi} that will
#' cover the range of options, but the fail may fit and return an error.
#'
#'
#' @param param parameter over which to profile
#' @param psi scalar or ordered vector of values for the interest parameter. If \code{NULL} (default), a grid of values centered at the MLE is selected
#' @param m number of observations (year) of interest. Required only for \code{param="VaR"} or \code{param="ES"}.
#' @param dat sample vector for the GP distribution
#' @param n.psi number of values of \code{psi} at which the likelihood is computed, if \code{psi} is not supplied (\code{NULL}). Odd values are more prone to give rise to numerical instabilities near the MLE
#' @param plot logical indiating whether \code{plot.fr} should be called upon exit
#' @author Leo Belzile, from code by A. Davison from the \code{hoa} package
#' @return an object of class \code{fr} (see \code{\link[hoa]{tem}}) with elements
#' \itemize{
#' \item{\code{normal}: }{maximum likelihood estimate and standard error of the interest parameter \eqn{psi}}
#' \item{\code{par.hat}: }{maximum likelihood estimates}
#' \item{\code{par.hat.se}: }{standard errors of maximum likelihood estimates}
#' \item{\code{th.rest}: }{estimated maximum profile likelihood at (\eqn{psi},\eqn{\hat{\lambda}})}
#' \item{\code{r}: }{values of likelihood root corresponding to \eqn{\psi}}
#' \item{\code{psi}: }{vector of interest parameter}
#' \item{\code{q}: }{vector of likelihood modifications}
#' \item{\code{rstar}: }{modified likelihood root vector}
#' \item{\code{param}: }{parameter}
#' }
#' @export
#' @examples
#' set.seed(123)
#' dat <- evd::rgpd(n=40, scale=1, shape=-0.1)
#' #with plots
#' m1 = gpd.tem(param="shape", n.psi=50, dat=dat,plot=TRUE)
#' m2 = gpd.tem(param="scale", n.psi=50, dat=dat)
#' m3 = gpd.tem(param="VaR", n.psi=50, dat=dat, m=100)
#' #Providing psi
#' m4 = gpd.tem(param="ES", dat=dat, m=100, psi = c(seq(2,5,length=15),seq(5, 35, length=45)))
#' plot(m4, c(2,4)) #displays numerical instability for values of r around zero
#' plot(fr4 <- spline.corr(m4), which=c(2,4))
#' confint(m1)
#' confint(m4, parm=2, warn=FALSE)
gpd.tem <- function (param=c("scale","shape","VaR","ES"), psi = NULL,
										 m=NULL, dat, n.psi = 50, plot=FALSE){
	if(param %in% c("VaR","ES") && is.null(m)){
		stop("Parameter `m' missing")
	}
	gpd.startvals <- function(dat, m=1000){
		#dat <- rgpd(n=n, scale=scale, shape=shape)
		start <- suppressWarnings(try(as.vector(mev::gp.fit(dat, 0, method="Grimshaw",show=FALSE)$estimate)))
		# if(is.character(start) || start[2] < -1){ #actually, always returns warnings rather than errors
		# 	start <- as.vector(extRemes::fevd(dat, threshold=0, type="GP",method="Lmoments")$results)
		# }
		scale0 <- start[1]
		shape0 <- start[2]
		ym0 <- scale0/shape0*(m^shape0-1)
		ES0 <- ifelse(shape0 < 1, (ym0+scale0)/(1-shape0), Inf)
		return(c("scale"=scale0, "shape"=shape0, "VaR"=ym0, "ES"=ES0))
	}

	th.init <- gpd.startvals(dat=dat,ifelse(is.null(m),10,m))
	tr.par <- function(par){ c(exp(par[1]),par[2])}
	itr.par <- function(par){ c(log(par[1]),par[2])}
	if(param %in% c("scale","shape")){
		xmax <- max(dat)
		th.init <- th.init[1:2]
		pin <- switch(param,scale=1,shape=2)
		nlogL.full <- function(par, dat){ gpd.ll.optim(par, dat)};
		nlogL.rest <- function(lam, psi, dat){ gpd.ll.optim(switch(pin,c(psi,lam),c(lam,psi)), dat)}
		gr.full <- function(par, dat){ -gpd.score(tr.par(par), dat)}
		gr.rest <- function(lam, psi, dat){
			parsc <- rep(0,2); parsc[pin]=psi; parsc[-pin] <- lam; parsc <- tr.par(parsc);
			-gpd.score(parsc, dat)[-pin]
		}
		# suppressWarnings(full <- optim(par=itr.par(th.init), method="BFGS",
		#                                fn=nlogL.full,gr=gr.full,dat=dat, control=list(maxit=500)))
		# if(full$convergence!=0){
		#   stop("Could not find maximum likelihood estimates. Change starting values")
		# }
		# suppressWarnings(full <- optim(par=full$par, method="Nelder-Mead",
		#                                fn=nlogL.full,dat=dat,control=list(parscale=full$par, abstol=1e-10)))
		# full$estimate <- tr.par(full$par);
		#if(full$par[2] <= -0.5){stop("Fisher information matrix requires shape > -0.5")}
		#L.full <- -full$value
		## May 13th, alternative method based on fit using Grimshaw
		full <- mev::gp.fit(dat,0,show = FALSE)
		if(full$estimate[2] <= -0.5){stop("Fisher information matrix requires shape > -0.5")}
		full$hessian <- gpd.infomat(full$estimate, dat=dat)
		th.full <- full$estimate
		th.se <- sqrt(diag(base::solve(full$hessian)))
		psi.se <- th.se[pin]
		L.full <- -full$deviance/2
		out <- NULL
		out$normal <- c(th.full[pin], psi.se)
		out$par.hat <- th.full
		out$par.hat.se <- th.se
		V.f <- gpd.Vfun(full$estimate, dat=dat)
		if (is.null(psi)) {
			if(pin==1){
				psi <- seq(from = th.full[pin] - 2*psi.se, to = th.full[pin] + 3*psi.se, length = n.psi) #TODO
			} else if(pin==2){
				psi <- seq(from = max(-0.49,th.full[pin] - 2.3*psi.se), to = th.full[pin] + 3.5*psi.se, length = n.psi) #TODO
				#Cannot have shape less than -0.5, because the correction is based on the information matrix
			}
		}
		#psi <- seq(from = max(max(dat),th.full[1] - 2*psi.se), to = th.full[1] + 6*psi.se, length = n.psi) #TODO
		#Adjusted to return the correct value for psi, meaning strictly positive
		n.psi <- length(psi)
		K <- ceiling(n.psi/2)
		if (n.psi%%2==0){   K <- K + 1}
		L.rest <- J.rest <- psi
		out$th.rest <- matrix(th.full, n.psi, length(th.full), byrow = TRUE)
		if(pin==1){psi <- log(psi)}
		for (j in (K - 1):1) {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,
											par=itr.par(out$th.rest[(j+1),]-diff(out$th.rest[(j+1):(j+2),]))[-pin],
											psi=psi[j], dat=dat, lower=switch(pin, -exp(psi[j])/xmax,
																												log(ifelse(psi[j]<0, -psi[j]*xmax, 1e-8)) ),
											upper=switch(pin,6*out$par.hat.se[1]+out$par.hat[1],3),
											method="Brent", hessian=TRUE)
			)
			rest <- optim(fn=nlogL.rest,gr=gr.rest,
										par=itr.par(out$th.rest[(j+1),]-diff(out$th.rest[(j+1):(j+2),]))[-pin],
										psi=psi[j], dat=dat, lower=switch(pin, -exp(psi[j])/xmax,
																											log(ifelse(psi[j]<0, -psi[j]*xmax, 1e-8))),
										upper=switch(pin,6*out$par.hat.se[1]+out$par.hat[1],3),
										method="Brent", hessian=TRUE,
										control=list(parscale=rest$par, abstol=1e-10))
			if(rest$convergence==0){
				out$th.rest[j, -pin] <- rest$par
				out$th.rest[j, pin] <- psi[j]
				out$th.rest[j,] <- tr.par(out$th.rest[j,])
				L.rest[j] <- -rest$value
				J.rest[j] <- gpd.infomat(out$th.rest[j,],dat=dat)[-pin,-pin]

			}
		}
		# if (n.psi%%2==0) {
		suppressWarnings(
			rest <- optim(fn=nlogL.rest,gr=gr.rest,
										par=itr.par(out$th.rest[K-1,]-diff(out$th.rest[(K-1):(K-2),]))[-pin],
										psi=psi[K], dat=dat, lower=switch(pin, -exp(psi[K])/xmax,
																											log(ifelse(psi[K]<0, -psi[K]*xmax, 1e-8))),
										upper=switch(pin,6*out$par.hat.se[1]+out$par.hat[1],3),
										method="Brent", hessian=TRUE,
										control=list(parscale=rest$par, abstol=1e-10))
		)
		# rest <- optim(fn=nlogL.rest,gr=gr.rest,par=itr.par(th.full)[-pin], psi=psi[K], dat=dat, method="BFGS"))
		#  }  else {
		# suppressWarnings(
		#   rest <- optim(fn=nlogL.rest,gr=gr.rest,par=itr.par(out$th.rest)[K - 1, -pin], psi=psi[K], dat=dat, method="BFGS"))
		#  }
		# rest <- optim(fn=nlogL.rest,gr=gr.rest,par=rest$par, psi=psi[K], dat=dat, method="BFGS",
		#               control=list(parscale=rest$par, abstol=1e-10))
		out$th.rest[K, pin] <- psi[K]
		out$th.rest[K, -pin] <- rest$par
		out$th.rest[K,] <- tr.par(out$th.rest[K,])
		L.rest[K] <- -rest$value
		J.rest[K] <- gpd.infomat(out$th.rest[K,],dat=dat)[-pin,-pin]
		for (j in (K + 1):n.psi) {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,
											par=itr.par(out$th.rest[(j-1),]-diff(out$th.rest[(j-1):(j-2),]))[-pin],
											psi=psi[j], dat=dat, lower=switch(pin, -exp(psi[j])/xmax,
																												log(ifelse(psi[j]<0, -psi[j]*xmax, 1e-8))),
											upper=switch(pin,6*out$par.hat.se[1]+out$par.hat[1],3),
											method="Brent", hessian=TRUE)
			)
			rest <- optim(fn=nlogL.rest,gr=gr.rest,
										par=itr.par(out$th.rest[(j-1),]-diff(out$th.rest[(j-1):(j-2),]))[-pin],
										psi=psi[j], dat=dat, lower=switch(pin, -exp(psi[j])/xmax,
																											log(ifelse(psi[j]<0, -psi[j]*xmax, 1e-8))),
										upper=switch(pin,6*out$par.hat.se[1]+out$par.hat[1],3),
										method="Brent", hessian=TRUE,
										control=list(parscale=rest$par, abstol=1e-10))
			# suppressWarnings(
			#   rest <- optim(fn=nlogL.rest,gr=gr.rest,par=sapply(out$th.rest[j - 1, -pin],jitter),
			#                 psi=psi[j], dat=dat, method="BFGS"))
			# rest <- optim(fn=nlogL.rest,gr=gr.rest,par=rest$par, psi=psi[j], dat=dat, method="BFGS",
			#               control=list(parscale=rest$par, abstol=1e-10))
			out$th.rest[j, -pin] <- rest$par
			out$th.rest[j, pin] <- psi[j]
			out$th.rest[j,] <- tr.par(out$th.rest[j,])
			L.rest[j] <- -rest$value
			J.rest[j] <- gpd.infomat(out$th.rest[j,],dat=dat)[-pin,-pin]
		}
		if(pin==1){psi <- exp(psi)}
		out$r <- sign(out$normal[1] - out$th.rest[, pin]) * sqrt(2 *(L.full - L.rest))
		dphi.dth.full <- gpd.dphi(par=th.full,V=V.f, dat=dat)
		D.bot <-  det(dphi.dth.full)
		j.th.th <- det(full$hessian)
		out$q <- out$psi <- psi
		for (j in 1:n.psi) {
			dphi.dth.rest <- gpd.dphi(par=out$th.rest[j,], V=V.f,dat=dat)
			dphi.dth.rest[pin, ] <- gpd.phi(out$par.hat, V=V.f,dat=dat) -
				gpd.phi(out$th.rest[j,], V=V.f,dat=dat)
			D.top <- det(dphi.dth.rest)
			j.lam.lam <- J.rest[j]
			out$q[j] <- sign(D.top)*sign(D.bot)*exp(log(abs(D.top))-log(abs(D.bot))) * sqrt(j.th.th/j.lam.lam)
		}
		####################################################################################################
	} else if(param=="VaR"){
		if(th.init[2]<0){
			ub <- -th.init[1]/th.init[2]
		} else{
			ub <- Inf
		}
		th.init <- th.init[c(3,2)]
		pin <- 1
		nlogL.full <- function(par, dat){ gpdr.ll.optim(par, dat, m)};
		nlogL.rest <- function(lam, psi, dat){ gpdr.ll.optim(c(psi,lam), dat, m)}
		gr.full <- function(par, dat){ -gpdr.score(tr.par(par), dat, m)}
		gr.rest <- function(lam, psi, dat){
			parsc <- rep(0,2); parsc[pin]=psi; parsc[-pin] <- lam; parsc <- tr.par(parsc);
			-gpdr.score(parsc, dat, m)[-pin]
		}
		suppressWarnings(full <- optim(par=itr.par(th.init), method="BFGS",
																	 fn=nlogL.full,gr=gr.full,dat=dat, control=list(maxit=500)))
		if(full$convergence!=0){
			stop("Could not find maximum likelihood estimates. Change starting values")
		}

		suppressWarnings(full <- optim(par=full$par, method="Nelder-Mead",
																	 fn=nlogL.full,dat=dat,control=list(parscale=full$par, abstol=1e-10)))
		full$estimate <- tr.par(full$par);
		if(full$estimate[2] <= -0.5){stop("Fisher information matrix requires shape > -0.5")}
		full$hessian <- gpdr.infomat(full$estimate, dat=dat, m)
		th.full <- full$estimate
		th.se <- sqrt(diag(base::solve(full$hessian)))
		psi.se <- th.se[pin]
		L.full <- -full$value
		out <- NULL
		out$normal <- c(th.full[pin], psi.se)
		out$par.hat <- th.full
		out$par.hat.se <- th.se
		V.f <- gpdr.Vfun(full$estimate, dat=dat,m=m)

		if (is.null(psi)){

			psi <- seq(from = max(quantile(dat,0.5),th.full[pin] - 2*psi.se),
								 to = min(ub,th.full[pin] + 5*psi.se), length = n.psi)
		}

		#TODO
		#psi <- seq(from = max(max(dat),th.full[1] - 2*psi.se), to = th.full[1] + 6*psi.se, length = n.psi) #TODO
		#Adjusted to return the correct value for psi, meaning strictly positive
		n.psi <- length(psi)
		K <- ceiling(n.psi/2)
		if (n.psi%%2==0){   K <- K + 1}
		L.rest <- J.rest <- psi
		out$th.rest <- matrix(th.full, n.psi, length(th.full), byrow = TRUE)
		if(pin==1){psi <- log(psi)}
		for (j in (K - 1):1) {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,par=itr.par(out$th.rest[j+1,])[-pin],
											psi=psi[j], dat=dat, method="BFGS", hessian=TRUE)
			)
			rest <- optim(fn=nlogL.rest,gr=gr.rest,par=rest$par, psi=psi[j], dat=dat, method="BFGS",
										control=list(parscale=rest$par, abstol=1e-10))
			if(rest$convergence==0){
				out$th.rest[j, -pin] <- rest$par
				out$th.rest[j, pin] <- psi[j]
				out$th.rest[j,] <- tr.par(out$th.rest[j,])
				L.rest[j] <- -rest$value
				J.rest[j] <- gpdr.infomat(out$th.rest[j,],dat=dat, m=m)[-pin,-pin]

			}
		}
		if (n.psi%%2==0) {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,par=itr.par(th.full)[-pin], psi=psi[K], dat=dat, method="BFGS"))
		}  else {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,par=itr.par(out$th.rest[K - 1,])[-pin], psi=psi[K], dat=dat, method="BFGS"))
		}
		rest <- optim(fn=nlogL.rest,gr=gr.rest,par=rest$par, psi=psi[K], dat=dat, method="BFGS",
									control=list(parscale=rest$par, abstol=1e-10))
		out$th.rest[K, pin] <- psi[K]
		out$th.rest[K, -pin] <- rest$par
		out$th.rest[K,] <- tr.par(out$th.rest[K,])
		L.rest[K] <- -rest$value
		J.rest[K] <- gpdr.infomat(out$th.rest[K,],dat=dat, m=m)[-pin,-pin]
		for (j in (K + 1):n.psi) {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,par=itr.par(out$th.rest[j - 1,])[-pin], psi=psi[j], dat=dat, method="BFGS"))
			rest <- optim(fn=nlogL.rest,gr=gr.rest,par=rest$par, psi=psi[j], dat=dat, method="BFGS",
										control=list(parscale=rest$par, abstol=1e-10))
			out$th.rest[j, -pin] <- rest$par
			out$th.rest[j, pin] <- psi[j]
			out$th.rest[j,] <- tr.par(out$th.rest[j,])
			L.rest[j] <- -rest$value
			J.rest[j] <- gpdr.infomat(out$th.rest[j,],dat=dat, m=m)[-pin,-pin]
		}
		if(pin==1){psi <- exp(psi)}
		out$r <- sign(out$normal[1] - out$th.rest[, pin]) * sqrt(2 *(L.full - L.rest))
		dphi.dth.full <- gpdr.dphi(par=th.full,V=V.f, dat=dat, m=m)
		D.bot <-  det(dphi.dth.full)
		j.th.th <- det(full$hessian)
		out$q <- out$psi <- psi
		for (j in 1:n.psi) {
			dphi.dth.rest <- gpdr.dphi(par=out$th.rest[j,], V=V.f,dat=dat, m=m)
			dphi.dth.rest[pin, ] <- gpdr.phi(out$par.hat, V=V.f,dat=dat, m=m) -
				gpdr.phi(out$th.rest[j,], V=V.f,dat=dat, m=m)
			D.top <- det(dphi.dth.rest)
			j.lam.lam <- J.rest[j]
			out$q[j] <- sign(D.top)*sign(D.bot)*exp(log(abs(D.top))-log(abs(D.bot))) * sqrt(j.th.th/j.lam.lam)
		}
		####################################################################################################
	} else if(param=="ES"){
		th.init <- th.init[c(4,2)]
		if(th.init[2]>1){stop("Infinite expected shortfall")}
		pin <- 1
		nlogL.full <- function(par, dat){ gpde.ll.optim(par, dat, m)};
		nlogL.rest <- function(lam, psi, dat){ gpde.ll.optim(switch(pin,c(psi,lam),c(lam,psi)), dat, m)}
		gr.full <- function(par, dat){ -gpde.score(tr.par(par), dat, m)}
		gr.rest <- function(lam, psi, dat){
			parsc <- rep(0,2); parsc[pin]=psi; parsc[-pin] <- lam; parsc <- tr.par(parsc);
			-gpde.score(parsc, dat, m)[-pin]
		}
		suppressWarnings(full <- optim(par=sapply(itr.par(th.init), jitter), method="BFGS",
																	 fn=nlogL.full,gr=gr.full,dat=dat, control=list(maxit=500)))
		if(full$convergence!=0){
			stop("Could not find maximum likelihood estimates. Change starting values")
		}
		suppressWarnings(full <- optim(par=full$par, method="Nelder-Mead",
																	 fn=nlogL.full,dat=dat,control=list(parscale=full$par, abstol=1e-10)))
		full$estimate <- tr.par(full$par);
		if(full$estimate[2] <= -0.5){stop("Fisher information matrix requires shape > -0.5")}
		full$hessian <- gpde.infomat(full$estimate, dat=dat, m)
		th.full <- full$estimate
		th.se <- sqrt(diag(base::solve(full$hessian)))
		psi.se <- th.se[pin]
		L.full <- -full$value
		out <- NULL
		out$normal <- c(th.full[pin], psi.se)
		out$par.hat <- th.full
		out$par.hat.se <- th.se
		V.f <- gpde.Vfun(full$estimate, dat=dat,m=m)
		if (is.null(psi))
			psi <- c(seq(from = th.full[pin] - 1.25*psi.se, to=th.full[pin], length=floor(n.psi))[-floor(n.psi)],
							 seq(from=th.full[pin], to = th.full[pin] + 10*psi.se, length = n.psi)[-1]) #TODO
		#psi <- seq(from = max(max(dat),th.full[1] - 2*psi.se), to = th.full[1] + 6*psi.se, length = n.psi) #TODO
		#Adjusted to return the correct value for psi, meaning strictly positive
		n.psi <- length(psi)
		K <- ceiling(n.psi/2)
		if (n.psi%%2==0){   K <- K + 1}
		L.rest <- J.rest <- psi
		out$th.rest <- matrix(th.full, n.psi, length(th.full), byrow = TRUE)
		if(pin==1){psi <- log(psi)}
		for (j in (K - 1):1) {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,par=itr.par(out$th.rest[j+1,])[-pin],
											psi=psi[j], dat=dat, method="BFGS", hessian=TRUE)
			)
			rest <- optim(fn=nlogL.rest,gr=gr.rest,par=rest$par, psi=psi[j], dat=dat, method="BFGS",
										control=list(parscale=rest$par, abstol=1e-10))
			if(rest$convergence==0){
				out$th.rest[j, -pin] <- rest$par
				out$th.rest[j, pin] <- psi[j]
				out$th.rest[j,] <- tr.par(out$th.rest[j,])
				L.rest[j] <- -rest$value
				J.rest[j] <- gpde.infomat(out$th.rest[j,],dat=dat, m=m)[-pin,-pin]

			}
		}
		if (n.psi%%2==0) {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,par=itr.par(th.full)[-pin], psi=psi[K], dat=dat, method="BFGS"))
		}  else {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,par=itr.par(out$th.rest)[K - 1, -pin], psi=psi[K], dat=dat, method="BFGS"))
		}
		rest <- optim(fn=nlogL.rest,gr=gr.rest,par=rest$par, psi=psi[K], dat=dat, method="BFGS",
									control=list(parscale=rest$par, abstol=1e-10))
		out$th.rest[K, pin] <- psi[K]
		out$th.rest[K, -pin] <- rest$par
		out$th.rest[K,] <- tr.par(out$th.rest[K,])
		L.rest[K] <- -rest$value
		J.rest[K] <- gpde.infomat(out$th.rest[K,],dat=dat, m=m)[-pin,-pin]
		for (j in (K + 1):n.psi) {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,par=sapply(out$th.rest[j - 1, -pin],jitter),
											psi=psi[j], dat=dat, method="BFGS"))
			rest <- optim(fn=nlogL.rest,gr=gr.rest,par=rest$par, psi=psi[j], dat=dat, method="BFGS",
										control=list(parscale=rest$par, abstol=1e-10))
			out$th.rest[j, -pin] <- rest$par
			out$th.rest[j, pin] <- psi[j]
			out$th.rest[j,] <- tr.par(out$th.rest[j,])
			L.rest[j] <- -rest$value
			J.rest[j] <- gpde.infomat(out$th.rest[j,],dat=dat, m=m)[-pin,-pin]
		}
		if(pin==1){psi <- exp(psi)}
		out$r <- sign(out$normal[1] - out$th.rest[, pin]) * sqrt(2 *(L.full - L.rest))
		dphi.dth.full <- gpde.dphi(par=th.full,V=V.f, dat=dat, m=m)
		D.bot <-  det(dphi.dth.full)
		j.th.th <- det(full$hessian)
		out$q <- out$psi <- psi
		for (j in 1:n.psi) {
			dphi.dth.rest <- gpde.dphi(par=out$th.rest[j,], V=V.f,dat=dat, m=m)
			dphi.dth.rest[pin, ] <- gpde.phi(out$par.hat, V=V.f,dat=dat, m=m) -
				gpde.phi(out$th.rest[j,], V=V.f,dat=dat, m=m)
			D.top <- det(dphi.dth.rest)
			j.lam.lam <- J.rest[j]
			out$q[j] <- sign(D.top)*sign(D.bot)*exp(log(abs(D.top))-log(abs(D.bot))) * sqrt(j.th.th/j.lam.lam)
		}
	}
	out$rstar <- out$r + log(out$q/out$r)/out$r
	out$param <- param
	class(out) <- "fr"
	if(plot){
		plot(out)
	}
	return(out)
}



#' Plot of tangent exponential model profile likelihood
#'
#' This function is adapted from \code{\link[hoa]{plot.fr}}. It differs mostly in
#' the placement of legends.
#'
#' @param x an object of class \code{fr} returned by \code{\link{gpd.tem}} or \code{\link{gev.tem}}.
#' @param ... further arguments to \code{plot} currently ignored. Providing a numeric vector \code{which} allows for custom selection of the plots. A logical \code{all}. See \strong{Details}.
#' @return graphs depending on argument \code{which}
#' @details Plots produced depend on the integers provided in \code{which}. \code{1} displays the Wald pivot, the likelihood root \code{r}, the modified likelihood root \code{rstar} and the likelihood modification \code{q} as functions of the parameter \code{psi}. \code{2} gives the renormalized profile log likelihood and adjusted form, with the maximum likelihood having ordinate value of zero. \code{3} provides the significance function, a transformation of \code{1}. Lastly, \code{4} plots the correction factor as a function of the likelihood root; it is a diagnostic plot aimed for detecting failure of the asymptotic approximation, often
#' due to poor numerics in a neighborhood of \code{r=0}; the function should be smooth. The function \code{\link{spline.corr}} is designed to handle this by removing the incorrect corrections and providing smoothed estimates, replacing outliers and missing values with the fitted values from the fit.
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
		par(mfrow=c(2,2), mar=c(4.5,4.5,1,0.1))
	} else if(sum(c(1,2,3,4) %in% whichPlot)==2){
		par(mfrow=c(1,2))
	}

	# top left: plot of pivot as a function of psi
	#ad hoc
	fr <- x
	xl <- fr$param
	#   fr$q[whichPlot(abs(fr$q)<0.1)] <- NA
	#   fr$rstar <- predict(loess(c(fr$r+log(fr$q/fr$r)/fr$r,rep(0,100))~c(fr$r,rep(0,100))),fr$r)
	#   fr$q <- predict(loess(fr$q~fr$r),fr$r)

	if(1 %in% whichPlot){
		plot(fr$psi,fr$r,type="l",xlab=xl,ylab="Value of pivot",ylim=c(-4, 4),
				 panel.first=abline(h=qnorm(c(0.005,0.025,0.05,0.5,0.95,0.975,0.995)),col="grey",lwd=0.7))
		lines(fr$psi,(fr$normal[1]-fr$psi)/fr$normal[2],col="green")
		lines(fr$psi,fr$q,col="red")
		lines(fr$psi,fr$r)
		lines(fr$psi,fr$rstar,col="blue")
		legend(x="topright",c("Wald pivot","Lik root","Modif root",expression(q(psi))),
					 lty=c(1,1,1,1),col=c("green","black","blue","red"),bty="n",cex=0.9)
	}
	# top right: log likelihood (and adjusted version, I think?) as a function of psi
	if(2 %in% whichPlot){
		plot(fr$psi,-fr$r^2/2,type="l",xlab=xl,ylab="Profile log likelihood",ylim=c(-8, 0),
				 panel.first=abline(h=-qchisq(c(0.95,0.99),df=1)/2,col="grey"),lwd=0.7)
		lines(fr$psi,-fr$rstar^2/2,col="blue")
		legend(x="bottomright",c(expression(l[p]),expression(l[a])),
					 lty=c(1,1),col=c("black","blue"),bty="n",cex=0.9)
		# optional: add diagnostic panels
	}
	if(3 %in% whichPlot){

		# lower left: plot of Phi(pivot) as a function of psi

		plot(fr$psi,pnorm(fr$r),type="l",xlab=xl,ylab="Significance function",ylim=c(0,1),
				 panel.first=abline(h=c(0.025,0.05,0.5,0.95,0.975),col="grey",lwd=0.7))
		lines(fr$psi,pnorm(fr$q),col="red")
		lines(fr$psi,pnorm(fr$rstar),col="blue")
		legend(x="topright",c("Lik. root","Modif. root","q(psi)"),
					 lty=c(1,1,1),col=c("black","blue","red"),bty="n")
	}
		# lower right: log(q/r)/r as a function of r (should be smooth)
	if(4 %in% whichPlot){
		plot(fr$r,fr$rstar,type="l",xlab="Likelihood root r",ylab="Correction log(q/r)/r",
				 panel.first={ abline(h=0,col="grey"); abline(v=0,col="grey")})
	}


	par(old.pars)

}

#' Spline correction for Fraser-Reid approximations
#'
#' The tangent exponential model can be numerically unstable for values close to \eqn{r=0}.
#' This function corrects these incorrect values, which are interpolated using splines.
#' The function takes as input an object of class \code{fr} and returns the same object with
#' different \eqn{r^*}{\code{rstar}} values.
#' @section Warning:
#'
#' While penalized (robust) splines often do a good job at capturing and correcting for numerical outliers and \code{NA}, it
#' may also be driven by unusual curves or fail at detecting outliers (or falsely identifying `correct' values as outliers). The user should always validate by comparing the plots of both the uncorrected (raw) output of the \code{\link{gpd.tem}} or \code{\link{gev.tem}} with that of \code{spline.corr}.
#' @details If available, the function uses \code{cobs} from the eponym package. The latter handles constraints and smoothness penalties, and is more robust than the equivalent \code{\link[stats]{smooth.spline}}.
#' @param fr an object of class \code{fr}, normally the output of \link{gpd.tem} or \link{gev.tem}.
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
	w <- pchisq(fr$r^2,0.5)
	#If any correction for q failed and returned NA
	corfailed <- which(is.na(fr$rstar))
	#If equispaced values for psi between MLE and other, than have r=0
	corfailed <- c(corfailed,which(fr$r==0))
	if(length(corfailed)>0){
		resp <- (fr$rstar-fr$r)[-corfailed]
		regr <- fr$r[-corfailed]
		w <- w[-corfailed]
	} else{
		resp <-  (fr$rstar-fr$r)
		regr <- fr$r
	}
	if(requireNamespace("cobs", quietly = TRUE)){
		spline <- cobs::cobs(y=resp, x=regr, w=w, constraint="none", lambda=1, print.mesg=FALSE, print.warn=FALSE)$fitted
	} else{
		spline <- rev(stats::smooth.spline(y=resp, x=regr, w=w, spar=0.9)$y)
	}
	#Compute difference between fitted values and rstar
	departure <- spline-resp
	#Ad-hoc fix of the values close to MLE where the numerical precision causes difficulty
	#Outlier detection via chi-square test
	#From package outliers, (c)Lukasz Komsta
	scores <- function (x, prob = NA){
		abs((x - mean(x))^2/var(x)) > qchisq(prob, 1)
	}
	bad <- which(scores(departure,prob=0.95))

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
		spline <- cobs::cobs(y=resp, x=regr, constraint="none", w=w, lambda=-1, ic="SIC",
												 knots.add=TRUE,repeat.delete.add=TRUE, print.mesg=FALSE, print.warn=FALSE)
		fr$spline <- spline
		fr$rstar <- predict(spline, fr$r, interval="none")[,2]+fr$r
	} else{
		spline <- stats::smooth.spline(x=na.omit(cbind(regr,resp)),w=w, cv=FALSE, all.knots=TRUE)
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
#' @param ... additional arguments passed to functions. Providing a logical \code{warn=FALSE} turns off warning messages when the lower or upper confidence interval for \code{psi} are extrapolated beyond the provided calculations.
#' @return a 2 by 3 matrix containing point estimates, lower and upper confidence intervals based on \code{r} and \code{rstar}
#' @export
confint.fr <- function(object, parm, level=0.95, ...){
	args <- list(...)
	if("warn" %in% names(args) && is.logical(args$warn)){
		warn <- args$warn
	} else{
		warn <- TRUE
	}
	#plot(object$psi~object$r)
	if(missing(parm)){
		ind  <- c(1,2)
	} else if(is.numeric(parm)){
		ind <- c(1,2)[c(1,2) %in% parm]
	} else{
	  ind <- c(1,2)[c("r","rstar") %in% parm]
	}
	if(length(ind)==0){
		stop("Invalid `parm` argument.")
	}
	if(1 %in% ind){
		if(requireNamespace("cobs", quietly = TRUE)){
			fit.r <- cobs::cobs(x=object$r,y=object$psi,constraint="none",lambda = 0, ic="SIC",pointwise=cbind(0,0,object$normal[1]),
										knots.add=TRUE,repeat.delete.add=TRUE, print.mesg=FALSE, print.warn=FALSE)
			pr <- predict(fit.r,c(0,sqrt(qchisq(level,1)),-sqrt(qchisq(level,1))))[,2]
		} else{
			fit.r <- stats::smooth.spline(x=na.omit(cbind(object$r,object$psi)), cv=FALSE)
			pr <- predict(fit.r,c(0,sqrt(qchisq(level,1)),-sqrt(qchisq(level,1))))$y
			pr[1] <- object$normal[1]
		}
		#lines(object$r,fit.r$fitted,col=2)
		if(warn){
			if(!any(object$r >  sqrt(qchisq(level,1)))){warning("Extrapolating the lower confidence interval for psi")}
			if(!any(object$r < -sqrt(qchisq(level,1)))){warning("Extrapolating the upper confidence interval for psi")}
		}
	}
	if(2 %in% ind){
	#plot(object$psi~object$rstar)
		if(requireNamespace("cobs", quietly = TRUE)){
			fit.rst <- cobs::cobs(x=object$rstar,y=object$psi,constraint="none",lambda = 0, ic="SIC",
											knots.add=TRUE,repeat.delete.add=TRUE, print.mesg=FALSE, print.warn=FALSE)
			prst <- predict(fit.rst,c(0,sqrt(qchisq(level,1)),-sqrt(qchisq(level,1))))[,2]
		} else{
			fit.rst <- stats::smooth.spline(x=na.omit(cbind(object$rstar,object$psi)), cv=FALSE)
			prst <- predict(fit.rst,c(0,sqrt(qchisq(level,1)),-sqrt(qchisq(level,1))))$y
		}
		#lines(x=object$rstar,fit.rst$fitted,col=2,pch=19)
		if(warn){
			if(!any(object$r >  sqrt(qchisq(level,1)))){warning("Extrapolating the adjusted lower confidence interval for psi.")}
			if(!any(object$r < -sqrt(qchisq(level,1)))){warning("Extrapolating the adjusted upper confidence interval for psi")}
		}
	}

	if(all(c(1,2) %in% ind)){
		conf <- cbind(pr,prst)
		colnames(conf) <- c("r","rstar")
		rownames(conf) <- c("Estimate","Lower CI","Upper CI")
		return(conf)
	} else if(1 %in% ind){
		names(pr) <- c("Estimate","Lower CI","Upper CI")
		return(pr)
	} else if(2 %in% ind){
		names(prst) <- c("Estimate","Lower CI","Upper CI")
		return(prst)
	}

}

###########################################################################################

#Tangent exponential model approximation for the GEV distribution
#Script: 31-03-2016


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
  mu = par[1]; sigma = par[2]; xi = par[3]
  mu - sigma/xi*(1-(-log(1-p))^(-xi))
}


#' Negative log-likelihood of the generalized extreme value distribution (return levels)
#'
#' @seealso \code{\link{gevr}}
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.ll <- function(par, dat, p){
	z = par[1]; sigma = par[2]; xi = par[3]
	muf = z+sigma/xi*(1-(-log(1-p))^(-xi))
	sum(-log(sigma)-(1/xi+1)*log(pmax(1+xi*(dat-muf)/sigma,0))-pmax(1+xi*(dat-muf)/sigma,0)^(-1/xi))
	#sum(log(-((dat - z)*xi/sigma + (-log(-p + 1))^(-xi))^(-1/xi - 1)*exp(-((dat - z)*xi/sigma + (-log(-p + 1))^(-xi))^(-1/xi))/sigma))
}

#' Negative log-likelihood parametrized in terms of location, log return level and shape in order to perform unconstrained optimization
#' @rdname gevr.ll
#' @seealso \code{\link{gevr}}
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.ll.optim <- function(par, dat, p){
	tpar = par; tpar[2] = exp(par[2])
	-gevr.ll(tpar,dat, p)
}

#' Score of the log-likelihood for the GEV distribution (return levels)
#' @seealso \code{\link{gevr}}
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.score <- function(par, dat, p){
	z = par[1]; sigma = par[2]; xi = par[3];
	c(sum((((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(-1/xi - 2)*xi*(1/xi + 1)*exp(-1/((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi))/sigma^2 - ((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(-2/xi - 2)*exp(-1/((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi))/sigma^2)*sigma*((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi + 1)*exp(1/(((dat - z)*xi/sigma + 1/(-log(-p + 1))^xi)^(1/xi)))),
		sum(-(dat*(-log(-p + 1))^xi - z*(-log(-p + 1))^xi - (dat*(-log(-p + 1))^xi - z*(-log(-p + 1))^xi - sigma)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))/((sigma*dat*xi*(-log(-p + 1))^xi - sigma*xi*z*(-log(-p + 1))^xi + sigma^2)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))),
		sum(-(xi*z*(-log(-p + 1))^xi - (dat*(-log(-p + 1))^xi - sigma*log(-log(-p + 1)))*xi + ((dat*(-log(-p + 1))^xi - sigma*log(-log(-p + 1)))*xi^2 + (dat*(-log(-p + 1))^xi - sigma*log(-log(-p + 1)))*xi - (xi^2*(-log(-p + 1))^xi + xi*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + (dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi - (dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + sigma)*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi)))/((dat*xi^3*(-log(-p + 1))^xi - xi^3*z*(-log(-p + 1))^xi + sigma*xi^2)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))
		))
}

#' Observed information matrix for GEV distribution (return levels)
#'
#'The information matrix is parametrized in terms of location, \eqn{(1-p)}th quantile and shape.
#' @seealso \code{\link{gevr}}
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.infomat <- function(par, dat, p){
	z = par[1]; sigma = par[2]; xi = par[3];
	infomat <- matrix(0, ncol=3,nrow=3)

	infomat[1,1] <- sum(-(xi*(-log(-p + 1))^(2*xi) - (xi^2*(-log(-p + 1))^(2*xi) + xi*(-log(-p + 1))^(2*xi))*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + (-log(-p + 1))^(2*xi))/((dat^2*xi^2*(-log(-p + 1))^(2*xi) + xi^2*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*xi*(-log(-p + 1))^xi + sigma^2 - 2*(dat*xi^2*(-log(-p + 1))^(2*xi) + sigma*xi*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))

	infomat[1,2] <- infomat[2,1] <-  sum(-(dat*(-log(-p + 1))^(2*xi) - z*(-log(-p + 1))^(2*xi) - sigma*(-log(-p + 1))^xi + (sigma*xi*(-log(-p + 1))^xi + sigma*(-log(-p + 1))^xi)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))/((sigma*dat^2*xi^2*(-log(-p + 1))^(2*xi) + sigma*xi^2*z^2*(-log(-p + 1))^(2*xi) + 2*sigma^2*dat*xi*(-log(-p + 1))^xi + sigma^3 - 2*(sigma*dat*xi^2*(-log(-p + 1))^(2*xi) + sigma^2*xi*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))
	infomat[1,3] <- infomat[3,1] <- sum(-((sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - dat*(-log(-p + 1))^(2*xi))*xi^2 + (sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - dat*(-log(-p + 1))^(2*xi))*xi + (xi^2*(-log(-p + 1))^(2*xi) + xi*(-log(-p + 1))^(2*xi))*z - (sigma*xi^3*(-log(-p + 1))^xi*log(-log(-p + 1)) + xi^2*z*(-log(-p + 1))^(2*xi) + (sigma*(-log(-p + 1))^xi*(log(-log(-p + 1)) + 1) - dat*(-log(-p + 1))^(2*xi))*xi^2)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + (dat*xi*(-log(-p + 1))^(2*xi) - xi*z*(-log(-p + 1))^(2*xi) + sigma*(-log(-p + 1))^xi)*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi)))/((dat^2*xi^4*(-log(-p + 1))^(2*xi) + xi^4*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*xi^3*(-log(-p + 1))^xi + sigma^2*xi^2 - 2*(dat*xi^4*(-log(-p + 1))^(2*xi) + sigma*xi^3*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))
	infomat[2,2] <- sum((dat^2*xi*(-log(-p + 1))^(2*xi) + (xi*(-log(-p + 1))^(2*xi) - (-log(-p + 1))^(2*xi))*z^2 - dat^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*(-log(-p + 1))^xi - 2*(dat*xi*(-log(-p + 1))^(2*xi) - dat*(-log(-p + 1))^(2*xi) + sigma*(-log(-p + 1))^xi)*z - (dat^2*xi*(-log(-p + 1))^(2*xi) + xi*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*(-log(-p + 1))^xi - sigma^2 - 2*(dat*xi*(-log(-p + 1))^(2*xi) + sigma*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))/((sigma^2*dat^2*xi^2*(-log(-p + 1))^(2*xi) + sigma^2*xi^2*z^2*(-log(-p + 1))^(2*xi) + 2*sigma^3*dat*xi*(-log(-p + 1))^xi + sigma^4 - 2*(sigma^2*dat*xi^2*(-log(-p + 1))^(2*xi) + sigma^3*xi*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))
	infomat[3,2] <- infomat[2,3] <- sum(-((sigma*dat*(-log(-p + 1))^xi*log(-log(-p + 1)) - dat^2*(-log(-p + 1))^(2*xi))*xi^2 - (xi^2*(-log(-p + 1))^(2*xi) + xi*(-log(-p + 1))^(2*xi))*z^2 + (sigma*dat*(-log(-p + 1))^xi*log(-log(-p + 1)) - dat^2*(-log(-p + 1))^(2*xi))*xi - ((sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - 2*dat*(-log(-p + 1))^(2*xi))*xi^2 + (sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - 2*dat*(-log(-p + 1))^(2*xi))*xi)*z - (sigma*dat*xi^3*(-log(-p + 1))^xi*log(-log(-p + 1)) - xi^2*z^2*(-log(-p + 1))^(2*xi) + (sigma*dat*(-log(-p + 1))^xi*(log(-log(-p + 1)) + 1) - dat^2*(-log(-p + 1))^(2*xi))*xi^2 - (sigma*xi^3*(-log(-p + 1))^xi*log(-log(-p + 1)) + (sigma*(-log(-p + 1))^xi*(log(-log(-p + 1)) + 1) - 2*dat*(-log(-p + 1))^(2*xi))*xi^2)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + (dat^2*xi*(-log(-p + 1))^(2*xi) + xi*z^2*(-log(-p + 1))^(2*xi) + sigma*dat*(-log(-p + 1))^xi - (2*dat*xi*(-log(-p + 1))^(2*xi) + sigma*(-log(-p + 1))^xi)*z)*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi)))/((sigma*dat^2*xi^4*(-log(-p + 1))^(2*xi) + sigma*xi^4*z^2*(-log(-p + 1))^(2*xi) + 2*sigma^2*dat*xi^3*(-log(-p + 1))^xi + sigma^3*xi^2 - 2*(sigma*dat*xi^4*(-log(-p + 1))^(2*xi) + sigma^2*xi^3*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))
	infomat[3,3] <- sum((sigma*dat*xi^4*(-log(-p + 1))^xi*log(-log(-p + 1))^2 + (4*sigma*dat*(-log(-p + 1))^xi*log(-log(-p + 1)) - 3*dat^2*(-log(-p + 1))^(2*xi))*xi^3 + (2*sigma*dat*(-log(-p + 1))^xi*(log(-log(-p + 1)) - 1) - (log(-log(-p + 1))^2 - 2*log(-log(-p + 1)))*sigma^2 - dat^2*(-log(-p + 1))^(2*xi))*xi^2 - (3*xi^3*(-log(-p + 1))^(2*xi) + xi^2*(-log(-p + 1))^(2*xi))*z^2 - (dat^2*xi^2*(-log(-p + 1))^(2*xi) + xi^2*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*xi*(-log(-p + 1))^xi + sigma^2 - 2*(dat*xi^2*(-log(-p + 1))^(2*xi) + sigma*xi*(-log(-p + 1))^xi)*z)*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^2 - (sigma*xi^4*(-log(-p + 1))^xi*log(-log(-p + 1))^2 + 2*(2*sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - 3*dat*(-log(-p + 1))^(2*xi))*xi^3 + 2*(sigma*(-log(-p + 1))^xi*(log(-log(-p + 1)) - 1) - dat*(-log(-p + 1))^(2*xi))*xi^2)*z - (sigma*dat*xi^5*(-log(-p + 1))^xi*log(-log(-p + 1))^2 + ((log(-log(-p + 1))^2 + 2*log(-log(-p + 1)))*sigma*dat*(-log(-p + 1))^xi - dat^2*(-log(-p + 1))^(2*xi))*xi^4 + (4*sigma*dat*(-log(-p + 1))^xi*log(-log(-p + 1)) - 3*dat^2*(-log(-p + 1))^(2*xi))*xi^3 - 2*(sigma*dat*(-log(-p + 1))^xi - sigma^2*log(-log(-p + 1)))*xi^2 - (xi^4*(-log(-p + 1))^(2*xi) + 3*xi^3*(-log(-p + 1))^(2*xi))*z^2 - (sigma*xi^5*(-log(-p + 1))^xi*log(-log(-p + 1))^2 + ((log(-log(-p + 1))^2 + 2*log(-log(-p + 1)))*sigma*(-log(-p + 1))^xi - 2*dat*(-log(-p + 1))^(2*xi))*xi^4 + 2*(2*sigma*(-log(-p + 1))^xi*log(-log(-p + 1)) - 3*dat*(-log(-p + 1))^(2*xi))*xi^3 - 2*sigma*xi^2*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi) + 2*(dat^2*xi^3*(-log(-p + 1))^(2*xi) - (sigma*dat*(-log(-p + 1))^xi*(log(-log(-p + 1)) - 2) - dat^2*(-log(-p + 1))^(2*xi))*xi^2 + (xi^3*(-log(-p + 1))^(2*xi) + xi^2*(-log(-p + 1))^(2*xi))*z^2 + (sigma*dat*(-log(-p + 1))^xi - sigma^2*(log(-log(-p + 1)) - 1))*xi - (2*dat*xi^3*(-log(-p + 1))^(2*xi) - (sigma*(-log(-p + 1))^xi*(log(-log(-p + 1)) - 2) - 2*dat*(-log(-p + 1))^(2*xi))*xi^2 + sigma*xi*(-log(-p + 1))^xi)*z - (dat^2*xi^3*(-log(-p + 1))^(2*xi) + xi^3*z^2*(-log(-p + 1))^(2*xi) +2*sigma*dat*xi^2*(-log(-p + 1))^xi + sigma^2*xi - 2*(dat*xi^3*(-log(-p + 1))^(2*xi) + sigma*xi^2*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi))*log((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi)))/((dat^2*xi^6*(-log(-p + 1))^(2*xi) + xi^6*z^2*(-log(-p + 1))^(2*xi) + 2*sigma*dat*xi^5*(-log(-p + 1))^xi + sigma^2*xi^4 - 2*(dat*xi^6*(-log(-p + 1))^(2*xi) + sigma*xi^5*(-log(-p + 1))^xi)*z)*((dat*xi*(-log(-p + 1))^xi - xi*z*(-log(-p + 1))^xi + sigma)/(sigma*(-log(-p + 1))^xi))^(1/xi)))
	return(-infomat)
}


#' Tangent exponential model statistics for the GEV distribution (return level)
#'
#' Vector implementing conditioning on approximate ancillary statistics for the TEM
#' @seealso \code{\link{gevr}}
#' @name gevr.temstat
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.Vfun <- function(par, dat, p){
	z = par[1]; sigma = par[2]; xi = par[3];
	cbind(1,
				(dat-z)/sigma,
				sigma*((dat - z)*xi/sigma + (-log(-p + 1))^(-xi))^(-1/xi)*(((-log(-p + 1))^(-xi)*log(-log(-p + 1)) - (dat - z)/sigma)/(((dat - z)*xi/sigma + (-log(-p + 1))^(-xi))*xi) + log((dat - z)*xi/sigma + (-log(-p + 1))^(-xi))/xi^2)/((dat - z)*xi/sigma + (-log(-p + 1))^(-xi))^(-1/xi - 1))
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
#' @param p probability associated with the (1-p)th quantile for return levels if \code{param="VaR"}.
#' @param dat sample vector for the GEV distribution
#' @param n.psi number of values of \code{psi} at which the likelihood is computed, if \code{psi} is not supplied (\code{NULL}). Odd values are more prone to give rise to numerical instabilities near the MLE
#' @param plot logical indiating whether \code{plot.fr} should be called upon exit
#' @author Leo Belzile, from code by A. Davison from the \code{hoa} package
#' @importFrom ismev gev.fit
#' @return an object of class \code{fr} (see \code{\link[hoa]{tem}}) with elements
#' \itemize{
#' \item{\code{normal}: }{maximum likelihood estimate and standard error of the interest parameter \eqn{psi}}
#' \item{\code{par.hat}: }{maximum likelihood estimates}
#' \item{\code{par.hat.se}: }{standard errors of maximum likelihood estimates}
#' \item{\code{th.rest}: }{estimated maximum profile likelihood at (\eqn{psi},\eqn{\hat{\lambda}})}
#' \item{\code{r}: }{values of likelihood root corresponding to \eqn{\psi}}
#' \item{\code{psi}: }{vector of interest parameter}
#' \item{\code{q}: }{vector of likelihood modifications}
#' \item{\code{rstar}: }{modified likelihood root vector}
#' \item{\code{param}: }{parameter}
#' }
#' @export
gev.tem <- function (param=c("loc", "scale", "shape", "VaR"), dat, psi = NULL,
										 p=NULL, n.psi = 50, plot=TRUE){
	#n.psi = 50; psi = NULL; make.V=gevr.Vfun; th.init <- par0r; plot=TRUE
	#n.psi = 50; psi = NULL; make.V=gev.Vfun; th.init <- parstart; plot=TRUE
	tr.par <- function(par){ c(par[1],exp(par[2]),par[3])}
	itr.par <- function(par){ c(par[1],log(par[2]),par[3])}
	gev.startvals <- function(dat, p=0.01){
		start <- ismev::gev.fit(dat,show=FALSE)$mle
		return(c(start, gev.retlev(start,p)))
	}
	th.init <- gev.startvals(dat, ifelse(is.null(p),0.01,p))

	# if(!missing(th.init) && !length(th.init)==3){
	# 	stop("Invalid starting parameter. Either leave unspecified or provide correct input")
	# }
	if(param %in% c("loc","scale","shape")){
		th.init <- th.init[1:3]
		make.V = gev.Vfun
		pin <- switch(param,loc=1,scale=2,shape=3)
		nlogL.full <- function(th, dat) gev.ll.optim(th, dat);
		nlogL.rest <- function(lam, psi, dat) gev.ll.optim(switch(pin,c(psi,lam),c(lam[1],psi,lam[2]),c(lam,psi)), dat)
		gr.full <- function(th, dat){ -gev.score(c(th[1],exp(th[2]),th[3]), dat)}
		gr.rest <- function(lam, psi, dat){
			parsc <- rep(0,3); parsc[pin]=psi; parsc[-pin] <- lam; parsc[2] <- exp(parsc[2]);
			-gev.score(parsc, dat)[-pin]
		}
		suppressWarnings(full <- optim(par=sapply(itr.par(th.init), jitter), method="BFGS",
																	 fn=nlogL.full,gr=gr.full,dat=dat))
		if(full$convergence!=0){
			stop("Could not find maximum likelihood estimates. Change starting values")
		}
		suppressWarnings(full <- optim(par=full$par, method="Nelder-Mead",
																	 fn=nlogL.full,dat=dat,control=list(parscale=full$par, abstol=1e-10)))
		full$estimate <- tr.par(full$par);
		full$hessian <- gev.infomat(full$estimate, dat=dat, "obs")
		th.full <- full$estimate
		th.se <- sqrt(diag(base::solve(full$hessian)))
		psi.se <- th.se[pin]
		L.full <- -full$value
		out <- NULL
		out$normal <- c(th.full[pin], psi.se)
		out$th.hat <- th.full
		out$th.hat.se <- th.se
		V.f <- make.V(full$estimate, dat=dat)
		if (is.null(psi))
			psi <- seq(from = th.full[pin] - 2.5*psi.se, to = th.full[pin] + 3.5*psi.se, length = n.psi) #TODO
		#psi <- seq(from = max(max(dat),th.full[1] - 2*psi.se), to = th.full[1] + 6*psi.se, length = n.psi) #TODO
		#Adjusted to return the correct value for psi, meaning strictly positive
		n.psi <- length(psi)
		K <- ceiling(n.psi/2)
		if (n.psi%%2==0){   K <- K + 1}
		L.rest <- J.rest <- psi
		out$th.rest <- matrix(th.full, n.psi, length(th.full), byrow = TRUE)
		if(pin==2){psi <- log(psi)}
		for (j in (K - 1):1) {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,
											par=itr.par(out$th.rest[(j+1),]-diff(out$th.rest[(j+1):(j+2),]))[-pin],
											psi=psi[j], dat=dat, method="BFGS"))
			rest <- optim(fn=nlogL.rest,par=rest$par, psi=psi[j], dat=dat, method="Nelder-Mead",
										control=list(parscale=rest$par, abstol=1e-10))
			if(rest$convergence==0){
				out$th.rest[j, -pin] <- rest$par
				out$th.rest[j, pin] <- psi[j]
				out$th.rest[j,] <- tr.par(out$th.rest[j,])
				L.rest[j] <- -rest$value
				J.rest[j] <- det(gev.infomat(out$th.rest[j,],dat=dat)[-pin,-pin])

			}
		}
		if (n.psi%%2==0) {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,par=itr.par(th.full)[-pin], psi=psi[K],
											dat=dat, method="BFGS"))
		}  else {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,par=itr.par(out$th.rest)[K - 1, -pin],
											psi=psi[K], dat=dat, method="BFGS"))
		}
		rest <- optim(fn=nlogL.rest,par=rest$par, psi=psi[K], dat=dat, method="Nelder-Mead",
									control=list(parscale=rest$par, abstol=1e-10))
		out$th.rest[K, pin] <- psi[K]
		out$th.rest[K, -pin] <- rest$par
		out$th.rest[K,] <- tr.par(out$th.rest[K,])
		L.rest[K] <- -rest$value
		J.rest[K] <- det(gev.infomat(out$th.rest[K,],dat=dat)[-pin,-pin])
		for (j in (K + 1):n.psi) {
			suppressWarnings(
				rest <- optim(fn=nlogL.rest,gr=gr.rest,par=sapply(out$th.rest[j - 1, -pin],jitter),
											psi=psi[j], dat=dat, method="BFGS"))
			rest <- optim(fn=nlogL.rest,par=rest$par, psi=psi[j], dat=dat, method="Nelder-Mead",
										control=list(parscale=rest$par, abstol=1e-10))
			out$th.rest[j, -pin] <- rest$par
			out$th.rest[j, pin] <- psi[j]
			out$th.rest[j,] <- tr.par(out$th.rest[j,])
			L.rest[j] <- -rest$value
			J.rest[j] <- det(gev.infomat(out$th.rest[j,],dat=dat)[-pin,-pin])
		}
		if(pin==2){psi <- exp(psi)}
		out$r <- sign(out$normal[1] - out$th.rest[, pin]) * sqrt(2 *(L.full - L.rest))
		dphi.dth.full <- gev.dphi(par=th.full,V=V.f, dat=dat)
		D.bot <-  det(dphi.dth.full)
		j.th.th <- det(full$hessian)
		out$q <- out$psi <- psi
		for (j in 1:n.psi) {
			dphi.dth.rest <- gev.dphi(par=out$th.rest[j,], V=V.f,dat=dat)
			dphi.dth.rest[pin, ] <- gev.phi(out$th.hat, V=V.f,dat=dat) -
				gev.phi(out$th.rest[j,], V=V.f,dat=dat)
			D.top <- det(dphi.dth.rest)
			j.lam.lam <- J.rest[j]
			out$q[j] <- sign(D.top)*sign(D.bot)*exp(log(abs(D.top))-log(abs(D.bot))) * sqrt(j.th.th/j.lam.lam)
		}
		###################################################################################################
		#Quantiles, or value-at-risk
	} else if(param=="VaR") {
		if(is.null(p)){stop("Invalid period for return levels")}
		if(length(th.init)==4){ #th.init was not provided by user
			th.init <- th.init[c(4,2,3)]
		}
		make.V = gevr.Vfun
		nlogLr.full <- function(th, dat, p) gevr.ll.optim(th, dat, p);
		nlogLr.rest <- function(lam, psi, dat, p) gevr.ll.optim(c(psi,lam), dat, p)
		grr.full <- function(th, dat, p){ -gevr.score(c(th[1],exp(th[2]),th[3]), dat, p=p)}
		grr.rest <- function(lam, psi, dat,p){ -gevr.score(tr.par(c(psi,lam)), dat=dat,p=p)[-1]}
		suppressWarnings(full <- optim(par=sapply(itr.par(th.init), jitter), method="BFGS",
																	 fn=nlogLr.full,gr=grr.full,dat=dat,p=p))
		suppressWarnings(full <- optim(par=full$par, method="Nelder-Mead",
																	 fn=nlogLr.full,gr=grr.full,dat=dat,p=p, control=list(parscale=full$par,reltol=1e-10)))
		if(full$convergence!=0){
			stop("Could not find maximum likelihood estimates. Change starting values")
		}
		full$estimate <- tr.par(full$par);
		full$hessian <- gevr.infomat(full$estimate, dat=dat, p=p)
		th.full <- full$estimate
		th.se <- sqrt(diag(base::solve(full$hessian)))
		psi.se <- th.se[1]
		L.full <- -full$value
		out <- NULL
		out$normal <- c(th.full[1], psi.se)
		out$th.hat <- th.full
		out$th.hat.se <- th.se
		V.f <- make.V(full$estimate, dat=dat, p=p)
		if (is.null(psi))
			psi <- c(seq(from = max(th.full[1] - 2.5*psi.se,quantile(dat,0.9)),th.full[1],length=floor(n.psi/2))[-floor(n.psi/2)],
							 seq(from=th.full[1], to = th.full[1] + ifelse(th.full[3]<0,7,4)*psi.se, length = ceiling(n.psi/2))[-1]) #TODO
		#psi <- seq(from = max(max(dat),th.full[1] - 2*psi.se), to = th.full[1] + 6*psi.se, length = n.psi) #TODO
		#Adjusted to return the correct value for psi, meaning strictly positive
		n.psi <- length(psi)
		K <- ceiling(n.psi/2)
		if (n.psi%%2==0){   K <- K + 1 }
		L.rest <- J.rest <- psi
		out$th.rest <- matrix(th.full, n.psi, length(th.full), byrow = TRUE)
		for (j in (K - 1):1) {
			suppressWarnings(
				rest <- optim(fn=nlogLr.rest,gr=grr.rest,p=p, par=out$th.rest[j+1,][-1], psi=psi[j], dat=dat,
											method="BFGS",control=list(parscale=th.init[-1],maxit=250)))
			suppressWarnings(
				rest <- optim(fn=nlogLr.rest,par=rest$par, psi=psi[j], dat=dat, p=p,
											control=list(maxit=250), method="Nelder-Mead"))
			#if(is.character(rest)) next;

			if(rest$convergence==0){
				out$th.rest[j, -1] <- rest$par
				out$th.rest[j, 1] <- psi[j]
				out$th.rest[j,] <- tr.par(out$th.rest[j,])
				L.rest[j] <- -rest$value
				J.rest[j] <- det(gevr.infomat(out$th.rest[j,],dat=dat, p=p)[-1,-1])
			}
		}
		if (n.psi%%2==0) {
			suppressWarnings(
				rest <- optim(fn=nlogLr.rest,gr=grr.rest,par=th.full[-1], psi=psi[K], dat=dat, p=p, method="BFGS"))
		}  else {
			suppressWarnings(
				rest <- optim(fn=nlogLr.rest,gr=grr.rest,par=out$th.rest[K - 1, -1], psi=psi[K], dat=dat, p=p, method="BFGS"))
		}
		rest <- optim(fn=nlogLr.rest,gr=grr.rest,par=rest$par, psi=psi[K],
									dat=dat, p=p, method="Nelder-Mead",control=list(parscale=rest$par))
		out$th.rest[K, 1] <- psi[K]
		out$th.rest[K, -1] <- rest$par
		out$th.rest[K,] <- tr.par(out$th.rest[K,])
		L.rest[K] <- -rest$value
		J.rest[K] <- det(gevr.infomat(out$th.rest[K,],dat=dat, p=p)[-1,-1])
		for (j in (K + 1):n.psi) {
			suppressWarnings(
				rest <- optim(fn=nlogLr.rest,gr=grr.rest,par=out$th.rest[j - 1, -1], psi=psi[j], dat=dat, p=p,
											control=list(parscale=th.full[-1], maxit=250), method="BFGS"))
			#Poor man's fix for cases were convergence is dubious
			suppressWarnings(
				rest <- optim(fn=nlogLr.rest,par=rest$par, psi=psi[j], dat=dat, p=p,
											control=list(maxit=250), method="Nelder-Mead"))
			out$th.rest[j, -1] <- rest$par
			out$th.rest[j, 1] <- psi[j]
			out$th.rest[j,] <- tr.par(out$th.rest[j,])
			L.rest[j] <- -rest$value
			J.rest[j] <- det(gevr.infomat(out$th.rest[j,],dat=dat, p=p)[-1,-1])
		}
		out$Lrest <- L.rest
		out$Lfull <- L.full
		out$r <- sign(out$normal[1] - out$th.rest[, 1]) * sqrt(2 *(L.full - L.rest))
		dphi.dth.full <- gevr.dphi(par=th.full,V=V.f, dat=dat, p=p)
		D.bot <-  det(dphi.dth.full)
		j.th.th <- det(full$hessian)
		out$q <- out$psi <- psi
		for (j in 1:n.psi) {
			dphi.dth.rest <- gevr.dphi(par=out$th.rest[j,], V=V.f,dat=dat, p=p)
			dphi.dth.rest[1, ] <- gevr.phi(out$th.hat, V=V.f,dat=dat, p=p) -
				gevr.phi(out$th.rest[j,], V=V.f,dat=dat, p=p)
			D.top <- det(dphi.dth.rest)
			j.lam.lam <- J.rest[j]
			out$q[j] <- sign(D.top)*sign(D.bot)*exp(log(abs(D.top))-log(abs(D.bot))) * sqrt(j.th.th/j.lam.lam)
		}
		#End of for loop for "VaR"
	}

	out$rstar <- out$r + log(out$q/out$r)/out$r
	class(out) <- "fr"
	if(plot){
		plot(out)
	}
	return(out)
}
