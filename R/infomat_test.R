#' Information matrix test statistic and MLE for the extremal index
#'
#' The Information Matrix Test (IMT), proposed by Suveges and Davison (2010) based on White (1982), is based
#' on the difference between the expected quadratic score and the second derivative of
#' the log-likelihood. The asymptotic distribution for each threshold \code{u} and gap \code{K}
#' is asymptotically \eqn{\chi^2}{chi-square} with one degree of freedom. The approximation is good for
#' \eqn{N >80} and conservative for smaller sample sizes. The test assumes independence between gaps.
#'
#' The procedure proposed in Suveges & Davison (2010) was corrected for erratas.
#' The maximum likelihood is based on the limiting mixture distribution of
#' the intervals between exceedances (an exponential with a point mass at zero).
#' The condition \eqn{D^{(K)}(u_n)}{D^K(u)} should be checked by the user.
#'
#' Fukutome et al. (2015) propose an ad hoc automated procedure
#' \enumerate{
#' \item Calculate the interexceedance times for each K-gap and each threshold, along with the number of clusters
#' \item  Select the (\code{u}, \code{K}) pairs for which IMT <0.05 (corresponding to a P-value of 0.82)
#' \item Among those, select the pair (\code{u}, \code{K}) for which the number of clusters is the largest
#' }
#' @author Leo Belzile
#'@references Fukutome, Liniger and Suveges (2015), Automatic threshold and run parameter selection: a climatology for extreme hourly precipitation in Switzerland. \emph{Theoretical and Applied Climatology}, \strong{120}(3), 403-416.
#' @references Suveges and Davison (2010), Model misspecification in peaks over threshold analysis. \emph{Annals of Applied Statistics}, \strong{4}(1), 203-221.
#' @references White (1982), Maximum Likelihood Estimation of Misspecified Models. \emph{Econometrica}, \strong{50}(1), 1-25.
#' @param x data vector
#' @param q vector of thresholds
#' @param K int specifying the largest K-gap
#' @param plot boolean: should the graphical diagnostic be plotted?
#'
#' @return an invisible list of matrices containing
#' \itemize{
#' \item \code{IMT} a matrix of test statistics
#' \item \code{pvals} a matrix of approximate p-values (corresponding to probabilities under a \eqn{\chi^2_1}{chis-square(1)} distribution)
#' \item \code{mle} a matrix of maximum likelihood estimates for each given pair (\code{u}, \code{K})
#' \item \code{loglik} a matrix of log-likelihood values at MLE for each given pair (\code{u}, \code{K})
#' \item \code{threshold} a vector of thresholds based on empirical quantiles at supplied levels.
#' \item \code{q} the vector \code{q} supplied by the user
#' }
#' @examples
#' library(evd)
#' infomat.test(x <- rgpd(n=1000), q = seq(0.9,0.995,length=10), K<- 3)
#' @export
infomat.test <- function(x, q, K, plot=TRUE){
	K <- as.integer(K)
	if(K<1){stop("Invalid K-gap specification")}
	if(isTRUE(all(q>0, q<1))){
		threshold <- quantile(x, q)
		} else{
		stop("Invalid vector of probabilities specified")
		}
	#Define functions

	ll.gap <- function(par, Ck, N0) {
		if(missing(N0)){ N0 <- sum(Ck==0) }
		ifelse(par>=1,0,N0*log(1-par))+2*(length(Ck)-N0)*log(par)-par*sum(Ck)
	}
	# ll.gap.score <- function(par, Ck, N0){
	# 	if(missing(N0)){N0 <- sum(Ck==0)}
	# 	ifelse(par>=1,0,-N0/(1-par))+2*(length(Ck)-N0)/par - sum(Ck)
	# }
	ll.gap.score.i <- function(par, Ck){
		ifelse(par>=1,0,-I(Ck==0)/(1-par))+2*I(Ck!=0)/par - Ck
	}
	Jn <- function(par, Ck, N0) {
			if(missing(N0)){N0 <- sum(Ck==0)}
		(ifelse(par>=1,0,N0/(1-par)^2)+4*(length(Ck)-N0)/par^2+sum(Ck^2-4*Ck/par))/length(Ck)
	}
	In <- function(par, Ck, N0) {
			if(missing(N0)){N0 <- sum(Ck==0)}
		(ifelse(par>=1,0,N0/(1-par)^2)+2*(length(Ck)-N0)/par^2)/length(Ck)
	}
	Dpn <- function(par, Ck, N0) {
		(-4*(length(Ck)-N0)/par^3+4*sum(Ck)/par^2)/length(Ck)
	}
	Vn <- function(par, Ck, N0) {
		mean((2*I(Ck!=0)/par^2+Ck^2-4*Ck/par-Dpn(par, Ck, N0)/In(par, Ck, N0)*ll.gap.score.i(par, Ck))^2)
		#Typo in Suveges and Davison (2010) and Suveges(2015) - should be squared
		#Otherwise the second part of the expression vanishes since it equals cst * score equation
	}
	Tn <-function(par, Ck, N0){
		(length(Ck)+1)*(Jn(par, Ck, N0)-In(par, Ck, N0))^2/Vn(par, Ck, N0)
	}

	#Interexceedance time
	#Define containers
	IMT <- llval <- mles <- matrix(0, ncol=K, nrow=length(q))
	#Exceedance values
	if(length(threshold)>1){
	  exceeds <- lapply(sapply(threshold, function(thresh){which(x>thresh)}), diff)
	} else if(length(threshold)==1){
	  exceeds <- list(diff(which(x>threshold)))
	} else{ stop("Invalid threshold")}
	for(k in 1:K){
		for(ind in 1:length(q)){
			#Interexceedance times, scaled by frequency
		Ck <- (length(exceeds[[ind]])/length(x))*pmax(exceeds[[ind]]-k,0)
		N0 <- sum(Ck==0) #Interclusters
		schi <- sum(Ck);
		N <- length(exceeds[[ind]])+1; Nc = N-1-N0
		mle <- (schi + N - 1 + Nc - sqrt((schi+N-1+Nc)^2-8*Nc*schi))/(2*schi)
		#   optimize(f=ll.gap,interval=c(0,1),lower = 0, upper=1, tol = .Machine$double.eps^0.5,
		# 	maximum=TRUE, Ck=Ck, N0=N0)
		mles[ind, k] <- mle
		llval[ind, k] <- ll.gap(mle,Ck,N0)
		#Information matrix test
		IMT[ind, k] <- Tn(par=mle, Ck=Ck, N0=N0)
		}
	pvals <- pchisq(IMT,df=1)
  }

	if(plot==TRUE){
		cols <- colorRampPalette(colors=c("red","white","blue"), bias=4.5, space="rgb")(100)
		image(x=q,y=seq(0.5,K+0.5,by=1),z=ceiling(100*(1-pvals)), breaks=1:101,col=cols,ylab="K",
		      main="Information Matrix Test",xlab="quantile")
		allpos <- expand.grid(1:length(q), 1:K)
		text(x=q[allpos[,1]], y=allpos[,2], labels=round(IMT,2),
		     col=c(rep("black",50),rep("white",50))[ceiling(100*(1-pvals))]) #
		text(x=q[allpos[,1]], y=allpos[,2], labels=paste0("(",round(pvals,2),")"), cex=0.75, pos=1,
		     col=c(rep("black",50),rep("white",50))[ceiling(100*(1-pvals))]) #
		mtext(side=3,line = 0,text="Test statistic (p-value)",adj=1)
	}
  invisible(list(IMT = IMT, pvals=pvals, loglik=llval, mle=mles, threshold=threshold, q=q))
}


#' Extremal index estimators based on interexceedance time and gap of exceedances
#'
#' The function implements the maximum likelihood estimator and iteratively reweigthed least
#' square estimators of Suveges (2007)  as well as the intervals estimator. The implementation
#' differs from the presentation of the paper in that an iteration limit is enforced to make sure
#' the iterative procedure terminates. sMultiple thresholds can be supplied
#'
#' The iteratively reweighted least square is a procedure based on the gaps of exceedances \eqn{S_n=T_n-1}{Tn-1}
#' The model is first fitted to non-zero gaps, which are rescaled to have unit exponential scale. The slope
#' between the theoretical quantiles and the normalized gap of exceedances is \eqn{b=1/\theta}{b=1/\theta},
#' with intercept \eqn{a=\log(\theta)/\theta}{a=log(\theta)/\theta}.
#' As such, the estimate of the extremal index is based on \eqn{\hat{\theta}=\exp(\hat{a}/\hat{b})}{\theta=exp(a/b)}.
#' The weights are chosen in such a way as to reduce the influence of the smallest values.
#' The estimator exploits the dual role of \eqn{\theta}{theta} as the parameter of the mean for
#' the interexceedance time as well as the mixture proportion for the non-zero component.
#'
#' The maximum likelihood is based on an independence likelihood for the rescaled gap of exceedances,
#' namely \eqn{\bar{F}(u_n)S(u_n)}{(1-F(u))*S(u)}. The score equation is equivalent to a quadratic equation in
#' \eqn{\theta}{theta} and the maximum likelihood estimate is available in closed form.
#' Its validity requires however condition \eqn{D^{(2)}(u_n)}{D2(u)} to apply;
#' this should be checked by the user beforehand.
#'
#'A warning is emitted if the effective sample size is less than 50 observations.
#'
#' @author Leo Belzile
#' @references Ferro and Segers (2003). Inference for clusters of extreme values,
#' JRSS: Series B, \strong{65}(2), 545-556.
#' @references Suveges (2007) Likelihood estimation of the extremal index. \emph{Extremes},
#'  \strong{10}(1), 41-55.
#' @references Suveges and Davison (2010), Model misspecification in peaks over threshold analysis. \emph{Annals of Applied Statistics}, \strong{4}(1), 203-221.
#'@references Fukutome, Liniger and Suveges (2015), Automatic threshold and run parameter selection: a climatology
#' for extreme hourly precipitation in Switzerland. \emph{Theoretical and Applied Climatology},
#' \strong{120}(3), 403-416.
#'
#'
#' @param x a vector containing of data points
#' @param q a vector of quantile levels in (0,1). Defaults to 0.95
#' @param method a string specifying the chosen method. Must be either \code{wls}
#' for weighted least squares, \code{mle} for maximum likelihood estimation or \code{intervals}
#' for the intervals estimator of Ferro and Segers (2003). Partial match is allowed.
#' @param plot a boolean specifying whether to plot the extremal index as a function of \code{q}
#' @return a vector or matrix of estimated extremal index of dimension \code{length(method)} by \code{length(q)}.
#' @examples
#' set.seed(234)
#' #Moving maxima model with theta=0.5
#' a <- 1; theta <-  1/(1+a)
#' sim <- evd::rgev(10001, loc=1/(1+a),scale=1/(1+a),shape=1)
#' x <- pmax(sim[-length(sim)]*a,sim[-1])
#' q <- seq(0.9,0.99,by=0.01)
#' ext.index(x=x,q=q,method=c("wls","mle"))
#' @export
ext.index <- function(x, q=0.95, method=c("wls","mle","intervals"), plot=FALSE){
  method <- match.arg(method, c("wls","mle","intervals"),several.ok = TRUE)
  stopifnot(all(q < 1), all(q>0))
  threshold <- quantile(x, q)
  q <- sort(q)
  #Main function, to be called recursively
  extmethods <- function(x=x, threshold=threshold[1], method=method){
    #The gap of exceedances
    exceeds <- c(diff(sapply(threshold, function(thresh){which(x>thresh)})))-1
    if(length(exceeds) < 50){
      warning("Small sample size - estimates may be unreliable")
    }
    N <- length(exceeds)
    chi <- sort(exceeds[which(exceeds>0)]); #Rescaled interexceedance time
    chi <- chi*length(exceeds)/length(x)
    Nc <- length(chi)
    if(method=="wls"){ #Suveges (2007); see Ferro (2003) for explanations
      xs <- -log((Nc:1)/(N+1)) #theoretical quantiles
      w <- cumsum(1/((N-1):(N-Nc))^2) #regression weights
      coefs <- lm(chi ~ xs, weights=w)$coef #slope is 1/th, intercept log(th)/th,
      #plot(y=chi, x=xs,xlim=c(0,7),ylim=c(0,11));abline(coefs)
      theta <- min(exp(coefs[1]/coefs[2]),1)
      Nc.new <- Nc
      trials <- 0
      while(floor(theta*(N-1)) != Nc.new && trials < 30){
      	trials <- trials +1
        Nc.new <- floor(theta*(N-1))
        if(Nc.new > Nc){ #pad estimates from original fit
          xs.new <- c(-log((Nc.new:(Nc+1))/(N+1)),xs)
          chi.new <- c(rep(0,Nc.new-Nc),chi)
          w.new <- c(cumsum( 1/((N-Nc.new):(N-Nc-1))^2),w)
        } else if(Nc.new < Nc){ #remove observations
          xs.new <- xs[-(1:(Nc-Nc.new))]
          chi.new <- chi[-(1:(Nc-Nc.new))]
          w.new <- w[-(1:(Nc-Nc.new))]
        }
        coefs <- lm(chi.new~xs.new, weights=w.new)$coef
        theta <- min(exp(coefs[1]/coefs[2]),1)
      }
    } else if(method=="mle"){ #Suveges (2007) and Suveges and Davison (2010)
      schi <- sum(chi)
      theta <- (schi+N-1+Nc-sqrt((schi+N-1+Nc)^2-8*Nc*schi))/(2*schi)
    } else if(method=="intervals"){ #Ferro and Segers (2003)
      if(max(exceeds-1) <= 2){
        theta <- min(1,2*(sum(exceeds+1)^2)/((N-1)*sum((exceeds+1)^2)))
      } else{
        theta <- min(1,2*((sum(exceeds))^2)/((N-1)*sum(exceeds*(exceeds-1))))
      }
    }
    return(theta)
  }
  #Provide a vector if multiple thresholds are supplied
  theta <- sapply(method, function(metho){
                  sapply(threshold, function(thresh){
                    extmethods(x=x, threshold=thresh, method=metho)})})
  if(plot){
    matplot(q, theta, type="l", xlab="q",ylim=c(0,1),ylab=expression(theta),
         main="Extremal index",lwd=2)
    legend(x="bottomright",legend=method,bty="n",lwd=2,lty=1:3,col=1:3)
  }
  return(t(theta))
}

