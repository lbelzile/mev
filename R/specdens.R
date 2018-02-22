.returnAng <- function(ang, R, Rnorm=c("l1","l2","linf"), wgt=c("Empirical","Euclidean"), region=c("sum","min","max")){
  #Other cases supported
  ang <- as.matrix(ang)
  if(region=="sum"){
    if(wgt=="Euclidean"){
      if(Rnorm=="l1"){
        return(list(ang=ang,rad=R,
                    wts=as.vector(.EuclideanWeights(ang, rep(1/(ncol(ang)+1), ncol(ang))))))
      } else{
        return(list(ang=ang,rad=R,
                    wts=as.vector(.EuclideanWeights(ang-cbind(ang[,-1],ang[,1]),rep(0, ncol(ang))))))
        #as.vector(.EuclideanWeights(ang,rep(1/(ncol(ang)+1), ncol(ang))))))
      }
    } else if(wgt=="Empirical"){
      if(Rnorm=="l1"){
        scel.fit <- .emplik_intern(z=ang, thresh=1e-10, mu=rep(1/(ncol(ang)+1), ncol(ang)), lam=rep(0,ncol(ang)), eps=1/nrow(ang), itermax = 10000L)
      } else{
        scel.fit <- .emplik_intern(z=ang-cbind(ang[,-1],ang[,1]), thresh=1e-10,
                            mu=rep(0, ncol(ang)), lam=rep(0,ncol(ang)), eps=1/nrow(ang), itermax = 5000L)
      }
      if(scel.fit$conv){
        return(list(ang=ang,rad=R, wts=as.vector(scel.fit$wts)))
      } else{
        warning("Self-concordant empirical likelihood for the mean did not converge.")
        return(list(ang=ang,rad=R, wts=rep(NA, nrow(ang))))
      }
    }
  } else{
    aw <- switch(region,
                 "min" = apply(cbind(ang, 1-rowSums(ang)),1,max),
                 "max" = apply(cbind(ang, 1-rowSums(ang)),1,min),
                 "sum" = rep(1, nrow(ang))
    )
    #b <- (ang-1/(ncol(ang)+1))/aw
    b <- (ang-cbind(ang[,-1],ang[,1]))/aw
    scel.fit <- .emplik_intern(z=b,mu=rep(0, ncol(ang)), lam=rep(0,ncol(ang)), eps=1/nrow(ang))
    if(scel.fit$conv){
      su <- 1/(sum((1/aw)*scel.fit$wts)) #mean of p_ia_i
      wts <-  su/aw*(scel.fit$wts)
      return(list(ang=ang,rad=R, wts=wts))
    } else{
      return(list(ang=ang, rad=R))
    }
  }
}


#' Rank-based transformation to angular measure
#'
#' The method uses the pseudo-polar transformation for suitable norms, transforming
#' the data to pseudo-observations, than marginally to unit Frechet or unit Pareto.
#' Empirical or Euclidean weights are computed and returned alongside with the angular and
#' radial sample for values above threshold(s) \code{th}, specified in terms of quantiles
#' of the radial component \code{R} or marginal quantiles. Only complete tuples are kept.
#'
#' The empirical likelihood weighted mean problem is implemented for all thresholds,
#'  while the Euclidean likelihood is only supported for diagonal thresholds specified
#'  via \code{region=sum}.
#'
#' @param x an \code{n} by \code{d} sample matrix
#' @param Rnorm  character string indicating the norm for the radial component.
#' @param Anorm character string indicating the norm for the angular component. \code{arctan} is only implemented for \eqn{d=2}
#' @param marg character string indicating choice of marginal transformation, either to Frechet or Pareto scale
#' @param wgt character string indicating weighting function for the equation. Can be based on Euclidean or empirical likelihood for the mean
#' @param th threshold of length 1 for \code{"sum"}, or \code{d} marginal thresholds otherwise.
#' @param region character string specifying which observations to consider (and weight). \code{"sum"} corresponds to a radial threshold
#' \eqn{\sum x_i > }\code{th}, \code{"min"} to \eqn{\min x_i >}\code{th} and \code{"max"} to \eqn{\max x_i >}\code{th}.
#' @param is.angle logical indicating whether observations are already angle with respect to \code{region}. Default to \code{FALSE}.
#' @return a list with arguments \code{ang} for the \eqn{d-1} pseudo-angular sample, \code{rad} with the radial component
#' and possibly \code{wts} if \code{Rnorm="l1"} and the empirical likelihood algorithm converged. The Euclidean algorithm always returns weights even if some of these are negative.
#' @author Leo Belzile
#' @references Einmahl, J.H.J. and J. Segers (2009). Maximum empirical likelihood estimation of the spectral measure of an extreme-value distribution, \emph{Annals of Statistics}, \bold{37}(5B), 2953--2989.
#' @references de Carvalho, M. and B. Oumow and J. Segers and M. Warchol (2013). A Euclidean likelihood estimator for bivariate tail dependence, \emph{Comm. Statist. Theory Methods}, \bold{42}(7), 1176--1192.
#' @references Owen, A.B. (2001). \emph{Empirical Likelihood}, CRC Press, 304p.
#' @export
#' @importFrom "graphics" "plot.new" "plot.window"
#' @return a list with components
#' \itemize{
#' \item \code{ang} matrix of pseudo-angular observations
#' \item \code{rad} vector of radial contributions
#' \item \code{wts} empirical or Euclidean likelihood weights for angular observations
#' }
#'
#' @examples
#' x <- rmev(n=25, d=3, param=0.5, model="log")
#' wts <- angmeas(x=x, th=0, Rnorm="l1", Anorm="l1", marg="Frechet", wgt="Empirical")
#' wts2 <- angmeas(x=x, Rnorm="l2", Anorm="l2", marg="Pareto", th=0)
angmeas <- function(x, th, Rnorm=c("l1","l2","linf"),
                    Anorm=c("l1","l2","linf","arctan"),
                    marg=c("Frechet","Pareto"),
                    wgt=c("Empirical","Euclidean"),
                    region=c("sum","min","max"),
                    is.angle = FALSE){
  if (!is.matrix(x)){ x <- rbind(x, deparse.level = 0L)}
  if(missing(th)){
    warning("Threshold set to zero. Using all the data")
    th <- 0
  } else if(any(th >= 1,th < 0)){
   stop("Threshold must be specified by a probability in [0,1)")
  }
    #Match arguments
    Rnorm <- match.arg(Rnorm[1],c("l1","l2","linf"))
    Anorm <- match.arg(Anorm[1],c("l1","l2","linf","arctan"))
    marg <- match.arg(marg[1],c("Frechet","Pareto"))
    wgt  <- match.arg(wgt[1],c("Euclidean","Empirical"))
    region <- match.arg(region[1],c("sum","min","max"))
    if(!is.angle){
    #Use only complete cases
    #x <- na.omit(x)
    #Margins are transformed to unit Frechet/Pareto (PIT)
    S <- switch(marg,
      Frechet= -1/log(na.omit(apply(x, 2, rank, na.last = "keep", ties.method = "random")/(nrow(x) +  1))),
      Pareto =  1/(1-na.omit(apply(x, 2, rank, na.last = "keep", ties.method = "random")/(nrow(x) +  1)))
    )

    #Obtain radius
    R <- switch(Rnorm,
      l1   = rowSums(S),
      l2   = apply(S,1,function(x){sqrt(sum(x^2))}),
      linf = apply(S,1,max)
    )
    #Ordering observations: this is not strictly necessary
    #but of some use for later inspection + for defn of `above` when region="sum"
    #Ordering is same regardless of above norm
    ordR <- order(R)
    S <- S[ordR,]
    #Order radius
    R <- R[ordR]

    #Check conditions for different regions
    if(region!="sum"){
     if(ncol(S)!=length(th)){
       if(length(th)==1){
        warning("The length of the threshold provided does not match the dimension of the problem; assuming identical thresholds")
        th <- rep(th[1],ncol(S))
       } else{
        stop("Invalid argument `th`: the length of the threshold provided does not match the dimension of the problem.")
       }
     }
     if(wgt=="Euclidean"){
       warning("Currently not implemented; switching to empirical likelihood")
       wgt <- "Empirical"
     }
      if(region=="min"){
        above <- Reduce(intersect,sapply(1:ncol(S), function(i){which(S[,i]>quantile(S[,i],th[i]))}, simplify=FALSE))
      } else if(region=="max"){
        above <- Reduce(union,sapply(1:ncol(S), function(i){which(S[,i]>quantile(S[,i],th[i]))}, simplify=FALSE))
      }
       #list if of unequal length, a matrix otherwise!
    } else{ #if region=="sum"
     if(length(th)!=1){
       warning("Only using the first entry of `th'")
       th <- th[1]
     }
      above <- floor(th*length(R)):length(R) #works because observations are ordered
    }



    if(length(above) < 5){
      warning("Less than five observations above the chosen threshold")
    }
    if(Rnorm==Anorm){
      ang <- S[,-ncol(S),drop=FALSE]/R
    } else if(Anorm=="arctan"){
      if(ncol(S)!=2){
        stop("Invalid norm for sample, arctan transformation only for bivariate samples")
      }
      ang <- atan(S[,2]/S[,1])
    } else{
      ang <- switch(Anorm,
                    l1   = S[,-ncol(S),drop=FALSE]/rowSums(S),
                    l2   = S[,-ncol(S),drop=FALSE]/apply(S,1,function(x){sqrt(sum(x^2))}),
                    linf = S[,-ncol(S),drop=FALSE]/apply(S,1,max)
      )
    }
    #Cast ang to a matrix for the bivariate case
    ang <- as.matrix(ang[above,])
    R <- as.vector(R[above])
    } else{
     R <- NULL
    }
    rownames(ang) <- NULL #remove names for time series
    .returnAng(ang=ang, R=R, Rnorm=Rnorm,  wgt=wgt, region=region)
}

#' Weighted empirical distribution function
#' @param x observations
#' @param w weights
#' @param ord logical indicating whether x values are already ordered. Default to \code{FALSE}
.wecdf <- function (x, w){
  #Adapted from spatstat (c) Adrian Baddeley and Rolf Turner
  stopifnot(length(x) == length(w)) #also returns error if x is multidimensional
  nbg <- is.na(x)
  x <- x[!nbg]
  w <- w[!nbg]
  n <- length(x)
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")
    ox <- order(x)
    x <- x[ox]
    w <- w[ox]
    vals <- sort(unique(x))
  xmatch <- factor(match(x, vals), levels = seq_along(vals))
  wmatch <- tapply(w, xmatch, sum)
  wmatch[is.na(wmatch)] <- 0
  rval <- approxfun(vals, cumsum(wmatch), method = "constant", yleft = 0,
                    yright = sum(wmatch), f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  invisible(rval)
}


.pickands.emp <- function(emp, plot=FALSE, add=FALSE, tikz=FALSE, ...){
  pos <- seq(0,1,length=1000)
  if(is.null(emp$wts)){
    warning("Invalid input; no `wts' arguments, likely due to convergence failure")
    pickands.pos <- rep(NA,length(pos))
  } else{
    pickands.pos <- .Pickands_emp(pos,emp$ang,emp$wts)
  }
  if(isTRUE(plot) || isTRUE(add)){
    if(plot){
      plot.new()
      plot.window(c(0,1), c(0.5,1))
      axis(side=2, at=seq(0.5,1,by=0.1), pos=0,las=2,tck=0.01)
      axis(side=1, at=seq(0,1,by=0.1), pos=0.5,las=0,tck=0.01)
      #title("Pickands dependence function")
      #Empirical bounds
      lines(c(0,0.5),c(1,0.5),lty=3,col="gray")
      lines(c(0.5,1),c(0.5,1),lty=3,col="gray")
      lines(c(0,1),c(1,1),lty=3,col="gray")
      if(tikz==TRUE){
        mtext("$t$", side=1, line=2)
        mtext("$\\mathrm{A}(t)$", side=2, line=2)
      } else{
        mtext(expression(t), side=1, line=2)
        mtext(expression(A(t)), side=2, line=2)
      }
    }
    lines(pos, pickands.pos,...)
  }
  return(invisible(cbind(pos,pickands.pos)))
  }


#' Self-concordant empirical likelihood for a vector mean
#'
#' @param dat \code{n} by \code{d} matrix of \code{d}-variate observations
#' @param mu  \code{d} vector of hypothesized mean of \code{dat}
#' @param lam  starting values for Lagrange multiplier vector, default to zero vector
#' @param eps  lower cutoff for \eqn{-\log}{-log}, with default \code{1/nrow(dat)}
#' @param M upper cutoff for \eqn{-\log}{-log}.
#' @param thresh  convergence threshold for log likelihood (default of \code{1e-30} is agressive)
#' @param itermax  upper bound on number of Newton steps.
#' @export
#' @author Art Owen, \code{C++} port by Leo Belzile
#' @references Owen, A.B. (2013). Self-concordance for empirical likelihood, \emph{Canadian Journal of Statistics}, \bold{41}(3), 387--397.
#' @return a list with components
#' \itemize{
#'  \item \code{logelr} log empirical likelihood ratio.
#'  \item \code{lam} Lagrange multiplier (vector of length \code{d}).
#'  \item \code{wts} \code{n} vector of observation weights (probabilities).
#'  \item \code{conv} boolean indicating convergence.
#'  \item \code{niter} number of iteration until convergence.
#'  \item \code{ndec} Newton decrement.
#'  \item \code{gradnorm} norm of gradient of log empirical likelihood.
#' }
emplik <- function(dat, mu=rep(0, ncol(dat)), lam = rep(0, ncol(dat)), eps = 1/nrow(dat), M=1e30, thresh=1e-30, itermax=100){
	if(is.infinite(M)){	M = 1e30}
	.emplik_intern(z=dat, mu=mu, lam, eps=eps, M=M, thresh=thresh, itermax=itermax)
}

#' Dirichlet mixture model for the spectral density
#'
#' This function computes the empirical or Euclidean likelihood
#' estimates of the spectral measure and uses the points returned from a call to \code{angmeas} to compute the Dirichlet
#' mixture smoothing of de Carvalho, Warchol and Segers (2012), placing a Dirichlet kernel at each observation.
#'
#' @details The cross-validation
#' bandwidth is the solution of
#'  \deqn{\max_{\nu} \sum_{i=1}^n \log \left\{ \sum_{k=1,k \neq i}^n p_{k, -i} f(\mathbf{w}_i; \nu \mathbf{w}_k)\right\},}
#' where \eqn{f} is the density of the Dirichlet distribution, \eqn{p_{k, -i}} is the Euclidean weight
#' obtained from estimating the Euclidean likelihood problem without observation \eqn{i}.
#' @return an invisible list with components
#' \itemize{
#'  \item \code{nu} bandwidth parameter obtained by cross-validation;
#'  \item \code{dirparmat} \code{n} by \code{d} matrix of Dirichlet parameters for the mixtures;
#'  \item \code{wts} mixture weights.
#' }
#' @inheritParams angmeas
#' @export
#' @examples
#' set.seed(123)
#' x <- rmev(n=250, d=2, param=0.5, model="log")
#' out <- angmeasdir(x=x, th=0, Rnorm="l1", Anorm="l1", marg="Frechet", wgt="Empirical")
angmeasdir <- function(x, th, Rnorm=c("l1","l2","linf"),
                       Anorm=c("l1","l2","linf","arctan"),
                       marg=c("Frechet","Pareto"),
                       wgt=c("Empirical","Euclidean"),
                       region=c("sum","min","max"),
                       is.angle = FALSE){
  Rnorm <- match.arg(Rnorm, c("l1","l2","linf"))
  Anorm <- match.arg(Anorm, c("l1","l2","linf","arctan"))
  marg <- match.arg(marg, c("Frechet","Pareto"))
  wgt <- match.arg(wgt, c("Empirical","Euclidean"))
  region <- match.arg(region, c("sum","min","max"))[1]
  #Obtain angles and weights
  angmeasure <- angmeas(x=x, th=th, Rnorm=Rnorm, Anorm=Anorm, marg=marg, wgt=wgt, region=region, is.angle = is.angle)
  n <- length(angmeasure$wts); #d <- ncol(angmeasure$ang) + 1
  loowts <- matrix(0, ncol=n, nrow=n)
  for(i in 1:n){
   loowts[i,-i] <- .returnAng(ang = angmeasure$ang[-i,, drop = FALSE], R = NULL, Rnorm=Rnorm, wgt=wgt, region=region)$wts
  }
  angles <- cbind(angmeasure$ang, 1-rowSums(angmeasure$ang))
  optnu <- optimize(f=.loocvdens, ang = angles, wts = angmeasure$wts,
                    loowts = loowts, lower = 1e-8, upper = 2000, maximum = FALSE, tol = 1e-5)
  return(invisible(list(nu = optnu$minimum, dirparmat = optnu$minimum*angles, wts = angmeasure$wts)))
}

