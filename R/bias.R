#' Cox-Snell first order bias expression for the GEV distribution
#'
#' Bias vector for the GEV distribution based on an \code{n} sample.
#' @inheritParams gev
#' @export
#' @keywords internal
#' @seealso \code{\link{gev}}
gev.bias <- function(par, n){
	if(length(n)>1){
	  stop("Invalid argument for sample size")
	 }
  if(length(par) != 3){
    stop("Invalid argument for parameter vector, must be a vector of length 3.")
  }
	sigma <- par[2] ; xi <- par[3]
	if(xi < -1/3){
	  stop("Cox-Snell correction only valid if the shape is greater than -1/3")
	  }
	zeta3 <- 1.20205690315959428539973816151144999076498629234049888179227155
	zeta5 <- 1.0369277551433699263313654864570341680570809 #gsl::zeta(5)
	k111 <- ((1 + xi)^2*(1 + 4*xi)*gamma(1 + 3*xi))/sigma^3

	if(abs(xi) < 1e-3){ #Limiting case when xi=0, some of the calculations break down
	euler_gamma <- -psigamma(1)
	k112 <- (euler_gamma-3)/sigma^2
	k113 <- -1/12*(36*euler_gamma - 6*euler_gamma^2 - pi^2 - 24)/sigma^2
	k122 <- -1/6*(36*euler_gamma - 6*euler_gamma^2 - pi^2 - 24)/sigma^3
  k123 <- 1/12*(60*euler_gamma + 6*euler_gamma^3 - euler_gamma*pi^2 + 4*pi^2*(euler_gamma - 1) - 48*euler_gamma^2 - 4*pi^2 + 12*zeta3 - 12)/sigma^2
  k133 <- 0.10683192718888033249425142127224548061317544/sigma
  k222 <- 1/4*(48*euler_gamma + 4*euler_gamma^3 + 9*euler_gamma*pi^2 - 4*pi^2*(2*euler_gamma - 3) + pi^2*(euler_gamma - 1) - 36*euler_gamma^2 - 17*pi^2 + 8*zeta3 - 16)/sigma^2
  k223 <- 1/40*(20*euler_gamma^4 + 3*pi^4 - 200*euler_gamma^3 + 20*euler_gamma^2*(pi^2 + 18) + 60*pi^2 - 20*euler_gamma*(5*pi^2 - 8*zeta3 + 8) - 400*zeta3)/sigma^2
  k233 <- 1/48*(12*euler_gamma^5 - 140*euler_gamma^4 - 21*pi^4 + 20*euler_gamma^3*(pi^2 + 16) - 4*euler_gamma^2*(35*pi^2 - 60*zeta3 + 48) + 8*pi^2*(5*zeta3 - 4) + euler_gamma*(9*pi^4 + 160*pi^2 - 1120*zeta3) + 288*zeta5 + 640*zeta3)/sigma
  k333 <- -20.807671559558883514171052830537917750303231

  k11.2 <- -2*(xi + 1)^2*gamma(2*xi + 1)/sigma^3
  k11.3 <- 2*(xi + 1)^2*psigamma(2*xi + 1)*gamma(2*xi + 1)/sigma^2 + 2*(xi + 1)*gamma(2*xi + 1)/sigma^2
  k12.2 <- -2*(euler_gamma - 1)/sigma^3
  k12.3 <- (6*euler_gamma - 3*euler_gamma^2 - 1/2*pi^2 - 2)/(2*sigma^2)
  k13.2 <- 1/12*(12*euler_gamma - 6*euler_gamma^2 - pi^2)/sigma^2
  k13.3 <- (-6*euler_gamma - 4*euler_gamma^3 - 7/2*euler_gamma*pi^2 + 3/2*pi^2*(euler_gamma - 1) + 12*euler_gamma^2 + 7/2*pi^2 - 8*zeta3)/(6*sigma)
  k22.2 <- 1/3*(12*euler_gamma - 6*euler_gamma^2 - pi^2 - 6)/sigma^3
  k22.3 <- -1/6*(12*euler_gamma + 6*euler_gamma^3 + 4*euler_gamma*pi^2 - pi^2*(euler_gamma - 1) - 18*euler_gamma^2 - 4*pi^2 + 12*zeta3)/sigma^2
  k23.2 <- -1/12*(12*euler_gamma + 6*euler_gamma^3 - euler_gamma*pi^2 + 4*pi^2*(euler_gamma - 1) - 18*euler_gamma^2 + pi^2 + 12*zeta3)/sigma^2
  k23.3 <- -3.7096580935190566493843882211576614781371729/sigma
  k33.2 <- 0
  k33.3 <- -5.4502140978602180294657833995281271927087253

	} else{
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
	}

	#Derivatives of information matrix
	A1 <- 0.5*cbind(c(k111,k112,k113),c(k112,k122,k123),c(k113,k123,k133))
	A2 <- -cbind(c(k11.2,k12.2,k13.2),c(k12.2,k22.2,k23.2),c(k13.2,k23.2,k33.2))+0.5*cbind(c(k112,k122,k123),c(k122,k222,k223),c(k123,k223,k233))
	A3 <- -cbind(c(k11.3,k12.3,k13.3),c(k12.3,k22.3,k23.3),c(k13.3,k23.3,k33.3))+0.5*cbind(c(k113,k123,k133),c(k123,k223,k233),c(k133,k233,k333))

	#Information matrix
	infomat <- gev.infomat(par=c(0, sigma, xi), dat = 1, method = "exp", nobs = 1)
	infoinv <- solve(infomat)

	return(infoinv%*%cbind(A1,A2,A3)%*%c(infoinv)/n)
}


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
	if(missing(method) || method!="exp"){	method <- "obs"}
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
	if(missing(method) || method!="exp"){	method <- "obs"}
	gev.score(par, dat) - gev.infomat(par, dat, method)%*%gev.bias(par, length(dat))
}

#' Bias correction for GP distribution using Firth's modified score function or bias substraction
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
#' par <- gp.fit(dat, threshold=0, show=FALSE)$estimate
#' gpd.bcor(par,dat, "subtract")
#' gpd.bcor(par,dat, "firth") #observed information
#' gpd.bcor(par,dat, "firth","exp")
gpd.bcor <- function(par, dat, corr=c("subtract","firth"), method=c("obs","exp")){
	corr <- match.arg(corr, c("subtract","firth"))
#Basic bias correction - substract bias at MLE parbc=par-bias(par)
#bcor1 <- function(par, dat){ par-gpd.bias(par,length(dat))}
#Other bias correction - find bias corrected that solves implicit eqn parbc=par-bias(parbc)
	if(length(par)!=2){
	  stop("Invalid `par` argument.")
	}
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
			gpd.score(parcopy, dat) - gpd.infomat(parcopy, dat, method)%*%gpd.bias(parcopy, length(dat))
		}

		firthplus = try(rootSolve::multiroot(gpd.Fscore.plus,start=par+c(0,0.3),
															dat=dat, method=method, positive=TRUE,
															atol=1e-10, rtol=1e-8, ctol=1e-10), silent=TRUE)
		#Changed tolerance on 12-10-2016 to ensure that the root passes the all.equal test
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
			firthplus <- firthplus$root-c(0,0.3)
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



#' Bias correction for GEV distribution using Firth's modified score function or bias substraction
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
#' dat <- evd::rgev(n=40, loc = 1, scale=1, shape=-0.2)
#' par <- evd::fgev(dat)$estimate
#' gev.bcor(par,dat, "subtract")
#' gev.bcor(par,dat, "firth") #observed information
#' gev.bcor(par,dat, "firth","exp")
gev.bcor <- function(par, dat, corr=c("subtract","firth"), method=c("obs","exp")){
  corr <- match.arg(corr, c("subtract","firth"))
  #Basic bias correction - substract bias at MLE parbc=par-bias(par)
  #bcor1 <- function(par, dat){ par-gpd.bias(par,length(dat))}
  #Other bias correction - find bias corrected that solves implicit eqn parbc=par-bias(parbc)
  bcor <-  function(par, dat){
    if(par[3]< -1/3){
      warning("Invalid bias correction for GEV; need shape > -1/3")
      return(rep(NA,3))
    } else{
      bcor.rootfind <- nleqslv::nleqslv(x=par, fn=function(parbc, par, dat){
        parbc-par+gev.bias(parbc, length(dat))}, par=par, dat=dat) #termcd is termination code, 1 or 2 implies successful convergence
      if(bcor.rootfind$termcd == 1 || (bcor.rootfind$termcd==2 && isTRUE(all.equal(bcor.rootfind$fvec,rep(0,3),tolerance=1.5e-8)))){
        return(bcor.rootfind$x)
      } else{
        return(rep(NA,3))
      }
    }
  }
  bcorF <-  function(par, dat, method=c("obs","exp")){
    method <- match.arg(method, c("obs","exp"))
    firth = try(rootSolve::multiroot(gev.Fscore, start=par,
                                         dat=dat, method=method, positive=FALSE,
                                         atol=1e-10, rtol=1e-8, ctol=1e-10), silent=TRUE)
    #Changed tolerance on 12-10-2016 to ensure that the root passes the all.equal test
    if(is.character(firth) ||
       any(c(isTRUE(all.equal(firth$root[2],target=0, check.names=FALSE)),
             is.nan(c(firth$f.root)),
             firth$root[3]>3)
       )){
      #Can fail if combination rends observation to have a zero likelihood - error message
      #Or can reach the boundary and not be able to evaluate the root
      firth <- rep(NA,3)
    } else{
      firth <- firth$root
    }
    return(firth)
  }
  #Return values
  if(corr=="subtract"){
    return(bcor(par=par, dat=dat))
  }
  if(corr=="firth"){
    return(bcorF(par=par, dat=dat, method=method))
  }
}

#' Posterior predictive distribution and density for the GEV distribution
#'
#' This function calculates the posterior predictive density at points x
#' based on a matrix of posterior density parameters
.gev.postpred <- function(x, posterior, Nyr = 100, type = c("density","quantile")){
  rowMeans(cbind(apply(rbind(posterior), 1, function(par){
    switch(type,
    density = evd::dgev(x=x,loc=par[1]-par[2]*(1-Nyr^par[3])/par[3], scale=par[2]*Nyr^par[3], shape=par[3]),
    quantile = evd::qgev(x=x,loc=par[1]-par[2]*(1-Nyr^par[3])/par[3], scale=par[2]*Nyr^par[3], shape=par[3]))
  })))
}


#' N-year return levels, median and mean estimate
#'
#' @param par vector of location, scale and shape parameters for the GEV distribution
#' @param nobs integer number of observation on which the fit is based
#' @param N integer number of observations for return level. See \strong{Details}
#' @param type string indicating the statistic to be calculated (can be abbreviated).
#' @param p probability indicating the return level, corresponding to the quantile at 1-1/p
#'
#' @details If there are \eqn{n_y} observations per year, the \code{L}-year return level is obtained by taking
#' \code{N} equal to \eqn{n_yL}.
#' @export
#' @return a list with components
#' \itemize{
#' \item{est} point estimate
#' \item{var} variance estimate based on delta-method
#' \item{type} statistic
#' }
gev.Nyr <- function(par, nobs, N, type = c("retlev", "median", "mean"), p = 1/N){
  #Create copy of parameters
  mu <- par[1]; sigma <- par[2]; xi <- par[3]
  type <- match.arg(type, c("retlev", "median", "mean"))[1]
  #Check whether arguments are well defined
  if(type == "retlev"){
    stopifnot(p >= 0, p < 1)
    yp <- -log(1-p)
  } else{
    stopifnot(N >= 1)
  }

  #Euler-Masc. constant :
  emcst <- -psigamma(1)
  #Return levels, N-year median and mean for GEV
  estimate <- switch(type,
    retlev = ifelse(xi==0, mu - sigma * log(yp), mu - sigma / xi * (1 - yp^(-xi))),
    median = ifelse(xi==0, mu + sigma * (log(N) - log(log(2))), mu - sigma / xi *(1 - (N / log(2))^xi)),
    mean = ifelse(xi==0, mu + sigma * (log(N) + emcst), mu - sigma / xi *(1 - N^xi*gamma(1-xi)))
  )
  if(type == "retlev"){
    if(p > 0){
      if(xi == 0){
      grad_retlev <- c(1, -log(yp), 0.5*sigma*log(yp)^2)
      } else{
      grad_retlev <- c(1, -(1-yp^(-xi))/xi, sigma*(1-yp^(-xi))/xi^2-sigma/xi*yp^(-xi)*log(yp))
      }
    }
    if(p == 0){
      if(xi < 0){
        grad_retlev <- c(1, -1/xi, sigma/xi^2)
      } else{
        stop("Invalid argument; maximum likelihood estimate of the uptyper endpoint is Inf");
      }
    }
    #Variance estimate based on delta-method
    var_retlev <- t(grad_retlev) %*% solve(gev.infomat(par = par, dat = 1, method = "exp", nobs = nobs)) %*% grad_retlev
  } else if(type == "median"){
    #Gradient of N-years maxima median
    if(xi == 0){
      grad_Nmed <- c(1, log(N/log(2)), 0.5*sigma*log(N/log(2))^2)
    } else{
      grad_Nmed <- c(1, ((N/log(2))^xi - 1)/xi, sigma*(N/log(2))^xi*log(N/log(2))/xi - sigma*((N/log(2))^xi - 1)/xi^2)
    }
    #Delta-method covariance matrix
    var_Nmed <- t(grad_Nmed) %*% solve(gev.infomat(par = par, dat = 1, method = "exp", nobs = nobs)) %*% grad_Nmed
  } else if(type == "mean"){
    if(xi == 0){
      grad_Nmean <- c(1, log(N) + emcst, 0.5*sigma*(emcst^2 + pi^2 / 6 + 2 * emcst*log(N) + log(N)^2))
    } else{
      grad_Nmean <- c(1,(N^xi*gamma(-xi + 1) - 1)/xi, (N^xi*log(N)*gamma(-xi + 1) - N^xi*digamma(-xi + 1)*gamma(-xi + 1))*sigma/xi - (N^xi*gamma(-xi + 1) - 1)*sigma/xi^2)
    }
    var_Nmean <- t(grad_Nmean) %*% solve(gev.infomat(par = par, dat = 1, method = "exp", nobs = nobs)) %*% grad_Nmean
  }
  var_est <- switch(type, retlev = var_retlev, median = var_Nmed, mean = var_Nmean)
  return(list(est = estimate, var = var_est[1,1], type = type))
}


#' Asymptotic bias of block maxima for fixed sample sizes
#'
#' @param shape shape parameter
#' @param rho second-order parameter, non-positive
#' @references Dombry, C. and A. Ferreira (2017). Maximum likelihood estimators based on the block maxima method. \code{https://arxiv.org/abs/1705.00465}
#' @export
#' @return a vector of length three containing the bias for location, scale and shape (in this order)
gev.abias <- function(shape, rho){
  stopifnot(rho <= 0, shape > -0.5)
  if(shape != 0 && rho < 0){
    bmu <- (1+shape)/(shape*rho*(shape+rho))*(-(shape+rho)*gamma(1+shape)+(1+shape)*rho*gamma(1+2*shape)+shape*(1-rho)*gamma(1+shape-rho))
    bsigma <- (-shape-rho+(1+shape)*(shape+2*rho)*gamma(1+shape)-(1+shape)^2*rho*gamma(1+2*shape)+shape*gamma(2-rho)-shape*(1+shape)*(1-rho)*gamma(1+shape-rho))/(shape^2*rho*(shape+rho))
    bshape <- ((shape+rho)*(1+shape+psigamma(1)*shape)-(shape+shape^2*(1+rho)+2*rho*(1+shape))*gamma(1+shape)+(1+shape)^2*rho*gamma(1+2*shape)+shape^2*gamma(1-rho)-shape*(1+shape)*gamma(2-rho)+shape*(1+shape)*(1-rho)*gamma(1+shape-rho)-shape*rho*psigamma(2+shape)*gamma(2+shape)-shape^2*psigamma(2-rho)*gamma(2-rho))/(shape^3*rho*(shape+rho))
  } else if(rho == 0 && shape != 0){
    bmu <- (1+shape)/shape^2*((1+shape)*gamma(1+2*shape)-gamma(2+shape)-shape*psigamma(1+shape)*gamma(1+shape))
    bsigma <- (shape^2*gamma(1+shape)*psigamma(1+shape) - (shape + 1)^2*gamma(1+2*shape) + shape^2*gamma(1+shape) + shape*gamma(1+shape)*psigamma(1+shape) - (psigamma(1)+1)*shape + 3*shape*gamma(1+shape) + 2*gamma(1+shape) - 1)/shape^3
    bshape <- ((1+shape+shape*psigamma(1))^2+shape^2*pi^2/6+(1+shape)^2*gamma(1+2*shape)-2*(1+shape)*((1+shape)*gamma(1+shape)+shape*psigamma(1+shape)*gamma(1+shape)))/shape^4
  } else if(rho < 0 && shape == 0){
    bmu <- (-1 + rho + psigamma(1)* rho + (1-rho)*gamma(1 - rho))/rho^2
    bsigma <- (6 - 6*rho - 6*psigamma(1)^2*rho - pi^2*rho + 6*psigamma(1)*(1-2*rho)-6*(1-rho)*gamma(1-rho)*(1+psigamma(1-rho)))/(6*rho^2)
    bshape <- -1/(12*rho^3)*(-12*gamma(2-rho)*psigamma(2-rho) + 6*gamma(1-rho)*(2 + 2*(-1 + rho)^2*psigamma(1-rho)+ (-1 + rho) * rho * psigamma(1 - rho)^2 + (-1 + rho)*rho*psigamma(1 - rho, deriv=1)) +  rho*(psigamma(1)^2*(6 - 18*rho) + pi^2*(1 - 3*rho) + 6*(-psigamma(1)^3)*rho + 3*(-psigamma(1))*(-4 + (4 + pi^2)*rho) + 6*rho*(-2 - 2*psigamma(1, deriv = 2) + psigamma(2,2))))
  } else if(rho == 0 && shape == 0){
    #Last row of information matrix with (0, 1, shape)
    bmu <- 0.41184033042643969478888356141823227689702419
    bsigma <- 0.3324849071602740614700056376493236520104431
    bshape <- 2.4236060551770285007097120629947238284701404
  }
  c(solve(gev.infomat(c(0, 1, shape), dat = 1, method = "exp", nobs = 1)) %*% c(bmu, bsigma, bshape))
}

#' Asymptotic bias of threshold exceedances for k order statistics
#'
#' The formula given in de Haan and Ferreira, 2007 (Springer). Note that the latter differs from that found in Drees, Ferreira and de Haan.
#' @references Dombry, C. and A. Ferreira (2017). Maximum likelihood estimators based on the block maxima method. \code{https://arxiv.org/abs/1705.00465}
#' @param shape shape parameter
#' @param rho second-order parameter, non-positive
#' @export
#' @return a vector of length containing the bias for scale and shape (in this order)
gpd.abias <- function(shape, rho){
  stopifnot(rho <= 0, shape > -0.5)
  c(-rho/((1-rho)*(1+shape-rho)), (1+shape)/((1-rho)*(1+shape-rho)))
}
