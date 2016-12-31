##  Copyright (C) 2015-2017 Leo Belzile
##
##  This file is part of the "mev" package for R.  This program
##  is free software; you can redistribute it and/or modify it under the
##  terms of the GNU General Public License 3.0 as published by the Free
##  Software Foundation.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
##  MA 02111-1307 USA or look up the web page
##  http://www.gnu.org/copyleft/gpl.html.


#' Exact simulations of multivariate extreme value distributions
#'
#' Implementation of the random number generators for multivariate extreme-value distributions
#' and max-stable processes based on the two algorithms described in
#' Dombry, Engelke and Oesting (2016).
#'
#' @param n number of observations
#' @param d dimension of sample
#' @param param parameter vector for the logistic, bilogistic, negative bilogistic and extremal Dirichlet (Coles and Tawn) model.
#' Parameter matrix for the Dirichlet mixture. Degree of freedoms for extremal student model. See \bold{Details}.
#' @param sigma covariance matrix for Brown-Resnick and extremal Student-t distributions. Symmetric matrix of squared  coefficients \eqn{\lambda^2} for the Husler-Reiss model, with zero diagonal elements.
#' @param asy list of asymmetry parameters, as in \code{\link[evd]{rmvevd}}, of \eqn{2^d-1} vectors of size corresponding to the power set of \code{d}, with sum to one constraints.
#' @param alg algorithm, either simulation via extremal function (\code{'ef'}) or via the spectral measure (\code{'sm'}). Default to \code{ef}.
#' @param model for multivariate extreme value distributions, users can choose between 1-parameter logistic and negative logistic, asymmetric logistic and negative logistic, bilogistic, Husler-Reiss, extremal Dirichlet model (Coles and Tawn) or the Dirichlet mixture. Spatial models include
#' the Brown-Resnick, Smith, Schlather and extremal Student max-stable processes.
#' @param vario function specifying the variogram. Used only if provided in conjonction with \code{loc} and if \code{sigma} is missing
#' @param loc \code{d} by \code{k} matrix of location, used as input in the variogram \code{vario} or as parameter for the Smith model. If \code{grid} is \code{TRUE}, unique entries should be supplied.
#' @param weights vector of length \code{m} for the \code{m} mixture components. Must sum to one
#' @param grid Logical. \code{TRUE} if the coordinates are two-dimensional grid points (spatial models).
#' @param ... additional arguments for the \code{vario} function
#' @author Leo Belzile
#' @details The vector param differs depending on the model
#' \itemize{
#'  \item \code{log}: one dimensional parameter greater than 1
#'  \item \code{alog}: \eqn{2^d-d-1} dimensional parameter for \code{dep}. Values are recycled if needed.
#'  \item \code{neglog}: one dimensional positive parameter
#'  \item \code{aneglog}: \eqn{2^d-d-1} dimensional parameter for \code{dep}. Values are recycled if needed.
#'  \item \code{bilog}: \code{d}-dimensional vector of parameters in \eqn{[0,1]}
#'  \item \code{negbilog}: \code{d}-dimensional vector of negative parameters
#'  \item \code{ct, dir, negdir}: \code{d}-dimensional vector of positive (a)symmetry parameters. For \code{dir} and \code{negdir}, a \eqn{d+1}
#'  vector consisting of the \code{d} Dirichlet parameters and the last entry is an index of regular variation in \eqn{(0, 1]} treated as shape parameter
#'  \item \code{xstud}: one dimensional parameter corresponding to degrees of freedom \code{alpha}
#'  \item \code{dirmix}: \code{d} by \code{m}-dimensional matrix of positive (a)symmetry parameters
#' }
#'
#' Stephenson points out that the multivariate asymmetric negative logistic model is not a valid distribution function.
#' The implementation in \code{mev} uses the same construction as the asymmetric logistic distribution,
#' and as such it does not match the bivariate implementation of \link[evd]{rbvevd}.
#'
#' The dependence parameter of the \code{evd} package for the Husler-Reiss distribution can be recovered taking
#' for the Brown--Resnick model  \eqn{2/r=\sqrt(2\gamma(h))} where \eqn{h} is the lag vector between sites and \eqn{r=1/\lambda} for the Husler--Reiss.
#'
#' @section Warning:
#'As of version 1.8 (August 16, 2016), there is a distinction between models \code{hr} and \code{br}. The latter is meant to be used in conjonction with variograms. The parametrization differs between the two models.
#'
#' @return an \code{n} by \code{d} exact sample from the corresponding multivariate extreme value model
#'
#' @export
#' @references
#'Dombry, Engelke and Oesting (2016). Exact simulation of max-stable processes, \emph{Biometrika}, \bold{103}(2), 303--317.
#' @seealso \link{rmevspec}, \link[evd]{rmvevd}, \link[evd]{rbvevd}
#' @examples
#'set.seed(1)
#'rmev(n=100, d=3, param=2.5, model="log", alg="ef")
#'rmev(n=100, d=4, param=c(0.2,0.1,0.9,0.5), model="bilog", alg="sm")
#'## Spatial example using variogram, from Clement Dombry
#'#Variogram gamma(h) = scale*||h||^alpha
#'scale <- 0.5; alpha <- 1
#'vario <- function(x) scale*sqrt(sum(x^2))^alpha
#'#grid specification
#'grid.loc <- as.matrix(expand.grid(runif(4), runif(4)))
#'rmev(n=100, vario=vario,loc=grid.loc, model="br")
#'#Example with a grid (generating an array)
#'rmev(n=10, sigma=cbind(c(2,1),c(1,3)), loc=cbind(runif(4),runif(4)),model="smith", grid=TRUE)
#'## Example with Dirichlet mixture
#'alpha.mat <- cbind(c(2,1,1),c(1,2,1),c(1,1,2))
#'rmev(n=100, param=alpha.mat, weights=rep(1/3,3), model="dirmix")
rmev <-function(n, d, param, asy, sigma,
                model= c("log","alog","neglog","aneglog","bilog","negbilog","hr","br",
                         "xstud","smith","schlather","ct","dir","negdir","dirmix"),
                alg=c("ef","sm"), weights, vario, loc, grid=FALSE, ...){
  #Create gridded values if specification is for random field discretization
  if(!missing(loc)){
    if(ncol(loc)==1) grid=FALSE
    if(grid==TRUE){
      if(all(sapply(1:ncol(loc), function(i){length(unique(loc[,i]))==nrow(loc)}))){
      loc <- matrix(unlist(expand.grid(apply(loc, 2, as.list))),ncol=ncol(loc))
      } else{
       stop("Duplicate values in `loc' using `grid=TRUE' not allowed");
      }
    }
  }
  alg <- match.arg(alg)
  model <- match.arg(model)
  if(!model %in% c("alog","aneglog")){
    if(!missing(asy)){
      warning("Asymmetry parameter ignored")
    } else{
      asym <- matrix(TRUE, ncol=1,nrow=1)
    }
  }
	if(model == "schlather"){
		if(!missing(param)) warning("Parameter value (degrees of freedom) set to one for Schlather model")
		param <- 1;
		model <- "xstud"
	}
  #Define model families
  m1 <- c("log","neglog")
  m2 <- c("bilog","negbilog")
  m3 <- c("br","xstud","smith")
  m4 <- c("ct","dir","negdir")
  #Check that parameters are provided
  if(model %in% c(m1,m2, m4) && (!missing(param) && mode(param) != "numeric")){
    stop("Invalid parameter")
  }
  #Check spatial requirements
  if((!missing(loc) || !missing(vario)) && ! model %in% m3){
    if(model == "hr"){
      warning("Obsolete. Please use `model=br' instead of `hr' for spatial models.");
      model <- "br"
    } else{
    warning("Unused arguments `vario' or `loc'; only implemented for Extremal student or Brown-Resnick process");
    }
  }
  #One-parameter families
  if(model %in% m1){
    d <- as.integer(d)
    sigma = cbind(0)
    if(missing(param) || param < 0 || d < 1){
      stop("Invalid parameter value")
    }
    if(length(param)!=1){
      warning("Only first entry of param vector considered")
      param <- param[1]
    }
    if(model=="log"){
      if(param < 1.0){
        param <- 1.0/param
      }
    }
  } else if(model %in% m2){
    d <- as.integer(d)
    sigma = cbind(0)
    #if(model %in% c("bilog","negbilog")){
      if(missing(param) || length(param)!=d) stop("Invalid parameter value")
      #Check whether arguments are valid
      if(model=="bilog" && all(param>=1)){ param <- 1/param}
      if(model=="negbilog" && all(param>=0)){ param <- -param}
      if(any(param>1.0)) stop("Invalid param vector for bilogistic or negative bilogistic")
      if(any(param<0) && model=="bilog") warning("Negative parameter values in bilogistic");
      if(any(param>0) && model=="negbilog") warning("Positive parameter values in negative bilogistic");
    } else if(model %in% m4){
      sigma = cbind(0)
      if(missing(param)){
        stop("Invalid parameter value")
      }
      if(model=="ct"){
        if(length(param)!=d){
          if(length(param)==(d+1)){
            warning("Use `dir' model for the scaled extremal Dirichlet model.");
          } else{
           stop("Invalid arguments for the Coles and Tawn (extremal Dirichlet) model.")
          }
        }
      } else{
        if(length(param)!=(d+1)){
          warning("Use `dir' model for the scaled extremal Dirichlet");
        }
      }
    } else if(model %in% m3){
    if(missing(sigma) && !missing(vario) && !missing(loc)){
      if(!is.matrix(loc)) loc <- matrix(loc, ncol=1)
      stopifnot(is.function(vario))
      sigma <- sapply(1:nrow(loc), function(i) sapply(1:nrow(loc), function(j)
        vario(loc[i,], ...) + vario(loc[j,], ...) - vario(loc[i,]-loc[j,], ...)))
    }
		if(model=="xstud"){
			sigma <- cov2cor(sigma)
		}
    if(missing(sigma) || ncol(sigma)!=nrow(sigma)) stop("Invalid covariance matrix")
    if(any(diag(sigma)<=0)) stop("Degenerate covariance matrix; negative or zero entries found")
    if(model=="xstud" && any(diag(sigma)!=1)) warning("Extremal student requires correlation matrix")
    if(model=="xstud" && (missing(param) || length(param)!=1)) stop("Degrees of freedom argument missing or invalid")
    if(model=="smith" && missing(loc)) stop("Location should be provided for the Smith model")
    if(model=="smith" && ncol(as.matrix(loc))!=ncol(sigma)){
      stop("Covariance matrix of the Smith model should be of the same dimension as dimension of location vector")
    }
    d <- switch(model,
                xstud = ncol(sigma),
                br    = ncol(sigma),
                smith = nrow(loc)
    )
    if(model=="br"){
      param = 0
    } else if(model=="smith"){
      param <- 0
    }
  } else if(model %in% c("alog","aneglog")){
    #Sigma will be index of logistic sub-mixtures
    #param is vector of dependence parameters
    if(any(param<0)) stop("Parameter vector must be positive")
    param <- rep(param, length.out = 2^d - 1 - d) #if vector too short, recycle dep arguments
    if(model == "alog"){
      if(isTRUE(all.equal(param >= 1,rep(TRUE, length(param))))){
        param <- 1/param #sampler works with
      }
      sigma <- .mvasym.check(asy, param, d = d, model = "alog") #check parameters constraints
      param <- 1/param
    } else{
      sigma <- .mvasym.check(asy, param, d = d, model = "aneglog") #check parameters constraints
    }
    #Transform list to matrix, with correct correspondance
    asym <- sigma>0
    #Shed output to remove zero weight combinations
    if(d==2){
      param <- c(sigma[1,1]!=0, sigma[2,2]!=0, param) # not a matrix for d=2
    } else{
    zero_line <- which(rowSums(sigma[-(1:d),])==0)+d #Possibly empty,
    param <- c(!1:d %in% zero_line, param) #set dummies for point masses on edges
    if(length(zero_line)>0){
      sigma <- sigma[-zero_line,]
      asym <- asym[-zero_line,]
      param <- param[-zero_line]
    }
    }
    stopifnot(isTRUE(all.equal(nrow(sigma), nrow(asym),length(param))),
              isTRUE(all.equal(colSums(sigma), rep(1, ncol(sigma)))))
    ncompo <- rowSums(asym)
    } else if(model=="dirmix"){
    if(any(missing(param),
           length(weights)!=ncol(param) && ncol(param)!=1,
           any(param<0))){
      stop("Invalid arguments for the Dirichlet mixture")
    }
    if(!missing(weights)){
      if(any(weights<0))stop("Negative weights provided")
      if(sum(weights)!=1) warning("Weights do not sum to one")
      weights <- weights/sum(weights)
    }
    if(missing(d)){ d <- nrow(param)
    } else if(d != nrow(param)){
      stop("Dimension of `d' and provided `param' do not match")
    }
    #Checking for the mean constraints
    mar_mean <- colSums(t(param)/ colSums(param)*weights)
    if(! isTRUE(all.equal(mar_mean, rep(1/d,d),
          tolerance = .Machine$double.eps ^ 0.5))){
        	stop("Invalid mixture components")
    }
    #Switching parameters around to pass them to Rcpp function
    sigma <- param
    param <- weights
    } else if(model=="hr"){
      if(any(c(sigma<0,diag(sigma)!=rep(0,ncol(sigma)), ncol(sigma)!=nrow(sigma), ncol(sigma)!=d))){
        stop("Invalid parameter matrix for the Husler-Reiss model.
           Note: for Brown-Resnick model, please use model=`br' instead.");
      }
      if(!.is.CNSD(sigma)){
       stop("Parameter matrix for the Husler-Reiss model is not conditionally negative definite")
      }
     param <- 0
    }
  mod <- switch(model, log=1, neglog=2, dirmix=3, bilog=4, negbilog=4, xstud=5, br=6,
                  ct=7, dir=7, smith=8, hr=9, negdir=10)
  if(model %in% c("alog","aneglog")){
    .rmevasy(n=n, d=d, para=param, asym=asym, ncompo=ncompo, Sigma=sigma, model=mod)
  } else{
    if(model != "smith"){
      locat <- cbind(0)
    } else{
      locat <- loc
    }
    if(model %in% m3 && grid==TRUE){
      npdim <- d^(1/ncol(loc))
      if(!all.equal(npdim, as.integer(npdim))){
       stop("The dimension of the input grid does not match (square) covariance matrix")
      }
      ncompo <- c(1)
      array(t(switch(alg,
             ef=.rmevA2(n=n, d=d, para=param, model=mod,  Sigma=sigma, locat),
             sm=.rmevA1(n=n, d=d, para=param, model=mod,  Sigma=sigma, locat)
        )), dim=c(rep(npdim,ncol(loc)),n)
      )
    } else{
      ncompo <- c(1)
      switch(alg,
             ef=.rmevA2(n=n, d=d, para=param, model=mod,  Sigma=sigma, locat),
             sm=.rmevA1(n=n, d=d, para=param, model=mod,  Sigma=sigma, locat)
      )
    }
  }
}




#' Random samples from spectral distributions of multivariate extreme value models.
#'
#' Generate from \eqn{Q_i}{Qi}, the spectral measure of a given multivariate extreme value model based on the L1 norm.
#'
#' @section Note:
#'  This functionality can be useful to generate for example Pareto processes (by multiplying by a standard Pareto variable the output). Other functionals are not currently implemented.
#'
#' @inheritParams rmev
#'
#' @author Leo Belzile
#' @details The vector param differs depending on the model
#' \itemize{
#'  \item \code{log}: one dimensional parameter greater than 1
#'  \item \code{neglog}: one dimensional positive parameter
#'  \item \code{bilog}: \code{d}-dimensional vector of parameters in \eqn{[0,1]}
#'  \item \code{negbilog}: \code{d}-dimensional vector of negative parameters
#'  \item \code{ct}, \code{dir}, \code{negdir}: \code{d}-dimensional vector of positive (a)symmetry parameters. Alternatively, a \eqn{d+1}
#'  vector consisting of the \code{d} Dirichlet parameters and the last entry is an index of regular variation in \eqn{(0, 1]} treated as scale
#'  \item \code{xstud}: one dimensional parameter corresponding to degrees of freedom \code{alpha}
#'  \item \code{dirmix}: \code{d} by \code{m}-dimensional matrix of positive (a)symmetry parameters
#' }
#' @return an \code{n} by \code{d} exact sample from the corresponding multivariate extreme value model
#'
#' @references
#'Dombry, Engelke and Oesting (2016). Exact simulation of max-stable processes, \emph{Biometrika}, \bold{103}(2), 303--317.
#' @references Boldi (2009). A note on the representation of parametric models for multivariate extremes.
#' \emph{Extremes} \bold{12}, 211--218.
#'
#' @examples
#' set.seed(1)
#' rmevspec(n=100, d=3, param=2.5, model="log")
#' rmevspec(n=100, d=3, param=2.5, model="neglog")
#' rmevspec(n=100, d=4, param=c(0.2,0.1,0.9,0.5), model="bilog")
#' rmevspec(n=100, d=2, param=c(0.8,1.2), model="ct") #Dirichlet model
#' rmevspec(n=100, d=2, param=c(0.8,1.2,0.5), model="dir") #with additional scale parameter
#'#Variogram gamma(h) = scale*||h||^alpha
#' scale <- 0.5; alpha <- 1
#' vario <- function(x) scale*sqrt(sum(x^2))^alpha
#' #grid specification
#' grid.loc <- as.matrix(expand.grid(runif(4), runif(4)))
#' rmevspec(n=100, vario=vario,loc=grid.loc, model="br")
#' ## Example with Dirichlet mixture
#' alpha.mat <- cbind(c(2,1,1),c(1,2,1),c(1,1,2))
#' rmevspec(n=100, param=alpha.mat, weights=rep(1/3,3), model="dirmix")
#' @export
rmevspec <-function(n, d, param, sigma,
                    model=c("log","neglog","bilog","negbilog","hr","br","xstud","ct","dir","negdir","dirmix"),
                    weights, vario, loc,...){
  model <- match.arg(model)
  m1 <- c("log","neglog")
  m2 <- c("bilog","negbilog")
  m3 <- c("br","xstud")
  m4 <- c("ct","dir","negdir")
  if(model %in% c(m1,m2, m4) && (!missing(param) && mode(param) != "numeric")){
    stop("Invalid parameter")
  }

  if(model %in% m1){
    d <- as.integer(d)
    sigma = cbind(0)
    if(missing(param) || param < 0 || d < 1){
      stop("Invalid parameter value")
    }
    if(length(param)!=1){
      warning("Only first entry of param vector considered")
      param <- param[1]
    }
    if(model=="log"){
      if(param < 1.0){
        param <- 1.0/param
      }
    }
  } else if(model %in% m2){
    d <- as.integer(d)
    sigma = cbind(0)
    if(missing(param) || length(param)!=d) stop("Invalid parameter value")
      #Check whether arguments are valid
      if(model=="bilog" && all(param>=1)){ param <- 1/param}
      if(model=="negbilog" && all(param>=0)){ param <- -param}
      if(any(param>1.0)) stop("Invalid param vector for bilogistic or negative bilogistic")
      if(any(param<0) && model=="bilog") warning("Negative parameter values in bilogistic");
      if(any(param>0) && model=="negbilog") warning("Positive parameter values in negative bilogistic");
  } else if(model %in% m4){
    sigma = cbind(0)
      if(missing(param)){
        stop("Invalid parameter value")
      }
      if(model=="ct"){
        if(length(param)!=d){
          if(length(param)==(d+1)){
            warning("Use `dir' model for the scaled extremal Dirichlet");
          } else{
            stop("Invalid arguments for the Coles and Tawn extremal Dirichlet model")
          }
        }
      } else{ #model 'dir' or 'idir'
        if(length(param)!=(d+1)){
          warning("Use `dir' model for the scaled extremal Dirichlet");
        }
      }
    } else if(model %in% m3){
      if(missing(sigma) && !missing(vario) && !missing(loc)){
        if(!is.matrix(loc)) loc <- matrix(loc, ncol=1)
        stopifnot(is.function(vario))
        sigma <- sapply(1:nrow(loc), function(i) sapply(1:nrow(loc), function(j)
          vario(loc[i,],...) + vario(loc[j,],...) - vario(loc[i,]-loc[j,],...)))
      }
    	if(model=="xstud"){
  			sigma <- cov2cor(sigma)
  			if(any(diag(sigma)!=1)){warning("Extremal student requires correlation matrix") }
  			if(missing(param) || length(param)!=1){stop("Degrees of freedom argument missing or invalid")}
  		}
      d <- ncol(sigma)
      if(missing(sigma) || ncol(sigma)!=nrow(sigma)) stop("Invalid covariance matrix")
      if(any(diag(sigma)<=0)) stop("Degenerate covariance matrix; negative or zero entries found")
      if(model=="br"){
        param = 0
      }
    } else if(model=="dirmix"){
    if(any(missing(param),
           length(weights)!=ncol(param) && ncol(param)!=1,
           any(param<0))){
      stop("Invalid arguments for the Dirichlet mixture")
    }
    if(!missing(weights)){
      if(any(weights<0))stop("Negative weights provided")
      if(sum(weights)!=1) warning("weights do not sum to one")
      weights <- weights/sum(weights)
    }
    if(missing(d)){ d <- nrow(param)
    } else if(d != nrow(param)){
      stop("Dimension of d and provided param do not match")
    }
    #Checking for the mean constraints
    mar_mean <- colSums(t(param)/ colSums(param)*weights)
    if(! isTRUE(all.equal(mar_mean, rep(1/d,d),
          tolerance = .Machine$double.eps ^ 0.5))){
        	stop("Invalid mixture components")
    }
    #Switching parameters around to pass them to Rcpp function
    sigma <- param
    param <- weights
    } else if(model=="hr"){
     param = 0
    }
  if(!model=="smith"){
    loc <- cbind(0)
  }
  #Model
  mod <- switch(model, log=1, neglog=2, dirmix=3, bilog=4, negbilog=4, xstud=5, br=6,
                ct=7, dir=7, smith=8, hr=9, negdir=10)

  #Generate from spectral measure
  .rmevspec_cpp(n=n, d=d, para=param, model=mod, Sigma=sigma, loc=loc)
}


#' Internal function
#'
#' Takes a list of asymmetry parameters with an associated dependence vector and returns
#' the corresponding asymmetry matrix for the asymmetric logistic and asymmetric negative logistic models
#'
#' This function is extracted from the evd package and modified
#' (C) Alec Stephenson
#'
#' @param asy a list of \eqn{2^d-1} assymetry components, as in Stephenson bvevd functions
#' @param vector of \eqn{2^d-d-1} values for the dependence parameter
#' @param d dimension of the model
#' @param model, either \code{alog} for the asymmetric logistic or \code{aneglog}
#' for the asymmetric negative logistic
#'
#' @return a matrix of asymmetry components, enumerating all possible \eqn{2^d-1} subsets of
#' the power set
.mvasym.check <- function (asy, dep, d, model=c("alog","aneglog"))
{
  #Function subset is an internal function from the evd package
  subsets <- function(d){
      x <- 1:d
      k <- NULL
      for (m in x) k <- rbind(cbind(TRUE, k), cbind(FALSE, k))
      pset <- apply(k, 1, function(z) x[z])
      pset[sort.list(unlist(lapply(pset, length)))[-1]]
    }
    if(model=="alog"){
    if (mode(dep) != "numeric" || any(dep <= 0) || any(dep >
        1))
        stop("invalid argument for `dep'")
    } else{
      if (mode(dep) != "numeric" || any(dep <= 0))
        stop("invalid argument for `dep'")
    }
    nb <- 2^d - 1
    if (mode(asy) != "list" || length(asy) != nb)
        stop(paste("`asy' should be a list of length", nb))
    tasy <- function(theta, b) {
        trans <- matrix(0, nrow = nb, ncol = d)
        for (i in 1:nb) trans[i, (1:d %in% b[[i]])] <- theta[[i]]
        trans
    }
    b <- subsets(d)
    if (any(sapply(asy, length) != sapply(b, length)))
        stop("`asy' is not of the correct form")
    asy <- tasy(asy, b)
    if (!is.matrix(asy) || mode(asy) != "numeric")
        stop("`asy' is not of the correct form")
    if (min(asy) < 0 || max(asy) > 1)
        stop("`asy' must contain parameters in [0,1]")
    if (any(apply(asy, 2, sum) != 1) || any(asy[c(rep(FALSE,
        d), dep == 1), ] != 0) || any(apply(asy[-(1:d), , drop = FALSE],
        1, function(x) sum(x != 0)) == 1))
        stop("`asy' does not satisfy the appropriate constraints")
      asy
}

#' Is the matrix conditionally negative semi-definite?
#' Function adapted from "is.CNSD" in the CEGO package, v 2.1.0
#' @author Martin Zaefferer
#'
.is.CNSD <- function (X, tol = 1e-08){
    n <- nrow(X)
    P <- diag(n)
    if(n>2){
      diag(P[,-1]) <- -1
    } else if(n==2){#error with one dimensional case...
      P[1,-1]	<- -1
    }
    Xhat <- P %*% X %*% t(P)
    eigs <- eigen(Xhat[-n, -n], TRUE, TRUE)$values
    !eigs[1] > tol
}
