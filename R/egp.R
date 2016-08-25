### Code by Raphael Huser with contributions from P. Naveau for the methods described in
##@references Naveau, P., R. Huser, P. Ribereau, and A. Hannart (2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, \emph{Water Resour. Res.}, 52, 2753-2769, \code{doi:10.1002/2015WR018552}.


#' Fit an extended generalized Pareto distribution of Naveau et al.
#'
#' This is a wrapper function to obtain PWM or MLE estimates for
#' the extended GP models of Naveau et al. (2016) for rainfall intensities. The function calculates confidence intervals
#' by means of nonparametric percentile bootstrap and returns histograms and QQ plots of
#' the fitted distributions. The function handles both censoring and rounding.
#'
#' @details
#' The different models include the following transformations:
#' \itemize{
#' \item \code{model} 0 corresponds to uniform carrier, \eqn{G(u)=u}.
#' \item \code{model} 1 corresponds to a three parameters family, with carrier \eqn{G(u)=u^\kappa}.
#' \item \code{model} 2 corresponds to a three parameters family, with carrier \eqn{G(u)=1-V_\delta((1-u)^\delta)}.
#' \item \code{model} 3 corresponds to a four parameters family, with carrier \deqn{G(u)=1-V_\delta((1-u)^\delta))^{\kappa/2}}.
#' \item \code{model} 4 corresponds to a five parameter model (a mixture of \code{type} 2, with \eqn{G(u)=pu^\kappa + (1-p)*u^\delta}
#' }
#'
#' @author Raphael Huser and Philippe Naveau
#' @param data data vector.
#' @param model integer ranging from 0 to 4 indicating the model to select (see \code{\link{egp2}}).
#' @param method string; either \code{"mle"} for maximum likelihood, or \code{"pwm"} for probability weighted moments, or both.
#' @param init vector of initial values, comprising of \eqn{p}, \eqn{\kappa}, \eqn{\delta},\eqn{\sigma},\eqn{\xi} (in that order) for the optimization. All parameters may not appear depending on \code{model}.
#' @param censoring numeric vector of length 2 containing the lower and upper bound for censoring; \code{censoring=c(0,Inf)} is equivalent to no censoring.
#' @param rounded numeric giving the instrumental precision (and rounding of the data), with default of 0.
#' @param CI logical; should confidence interval be returned (percentile bootstrap).
#' @param R integer; number of bootstrap replications.
#' @param ncpus integer; number of CPUs for parallel calculations (default: 1).
#' @param plots logical; whether to produce histogram and density plots.
#' @export
#' @importFrom boot boot boot.ci
#' @importFrom grDevices rgb
#' @importFrom graphics hist layout
#' @importFrom evd dgpd pgpd qgpd
#' @importFrom gmm gmm
#' @seealso \code{\link{egp.fit}}, \code{\link{egp}}, \code{\link{egp2}}
#' @references Naveau, P., R. Huser, P. Ribereau, and A. Hannart (2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, \emph{Water Resour. Res.}, 52, 2753-2769, \code{doi:10.1002/2015WR018552}.
#' @examples
#' \donttest{
#' library(ismev)
#' data(rain)
#' egp2.fit(rain[rain>0], model=1, method="mle", init=c(0.9, gp.fit(rain,0, show=FALSE)$est),
#'  rounded=0.1,CI=TRUE, R=20)
#' }
egp2.fit <- function(data,model=1,method=c("mle","pwm"),init,censoring=c(0,Inf),
                        rounded=0,CI=FALSE,R=1000,ncpus=1,plots=TRUE){
  #Sanity checks
  initsize <- switch(model, 3, 3, 4, 5)
  if(length(init)!=initsize){stop("Invalid starting values in `init'; incorrect length.");}
  data = data[data>0]
  if(model==3 && "pwm" %in% method){stop("No explicit formula for PWMs for this model...")}

	x <- seq(0,max(data,na.rm=TRUE),by=0.05);
	degp2s <- c();
	qegp2s <- c();

	if(any(method=="pwm")){
		if(model==1){
			fit.pwm <- egp2.fit.pwm(x=data,type=1,kappa0=init[1],sigma0=init[2],xi0=init[3],censoring=censoring);
			if(CI){
				fit.pwm.boot <- boot::boot(data=data,statistic=egp2.fit.pwm.boot,R=R,
				                      type=1,kappa0=init[1],sigma0=init[2],ncpus=ncpus,
				                     xi0=init[3],censoring=censoring,parallel="multicore");
				CI.pwm.kappa <- boot::boot.ci(boot.out=fit.pwm.boot,index=1,type="perc")$perc[4:5];
				CI.pwm.sigma <- boot::boot.ci(boot.out=fit.pwm.boot,index=2,type="perc")$perc[4:5];
				CI.pwm.xi <- boot::boot.ci(boot.out=fit.pwm.boot,index=3,type="perc")$perc[4:5];
				CIs.pwm <- cbind(CI.pwm.kappa,CI.pwm.sigma,CI.pwm.xi);
			}
			degp2.pwm <- degp2(x=x,type=1,kappa=fit.pwm[1],sigma=fit.pwm[2],xi=fit.pwm[3]);
			degp2s <- c(degp2s,degp2.pwm);
			if(plots){
				qegp2.pwm <- qegp2(p=c(1:length(data))/(length(data)+1),type=1,kappa=fit.pwm[1],
				                   sigma=fit.pwm[2],xi=fit.pwm[3]);
				qegp2s <- c(qegp2s,qegp2.pwm);
				if(CI){
					q.pwm.boot <- mapply(FUN=qegp2,p=list(c(1:length(data))/(length(data)+1)),type=list(1),kappa=as.list(fit.pwm.boot$t[,1]),sigma=as.list(fit.pwm.boot$t[,2]),xi=as.list(fit.pwm.boot$t[,3]));
					q.pwm.L <- apply(q.pwm.boot,1,quantile,0.025,na.rm=TRUE);
					q.pwm.U <- apply(q.pwm.boot,1,quantile,0.975,na.rm=TRUE);
				}
			}
		} else if(model==2){
			fit.pwm <- egp2.fit.pwm(x=data,type=2,delta0=init[1],sigma0=init[2],xi0=init[3],censoring=censoring);
			if(CI){
				fit.pwm.boot <- boot::boot(data=data,statistic=egp2.fit.pwm.boot,R=R,type=2,delta0=init[1],sigma0=init[2],xi0=init[3],censoring=censoring,parallel="multicore",ncpus=ncpus);
				CI.pwm.delta <- boot::boot.ci(boot.out=fit.pwm.boot,index=1,type="perc")$perc[4:5];
				CI.pwm.sigma <- boot::boot.ci(boot.out=fit.pwm.boot,index=2,type="perc")$perc[4:5];
				CI.pwm.xi <- boot::boot.ci(boot.out=fit.pwm.boot,index=3,type="perc")$perc[4:5];
				CIs.pwm <- cbind(CI.pwm.delta,CI.pwm.sigma,CI.pwm.xi)
			}
			degp2.pwm <- degp2(x=x,type=2,delta=fit.pwm[1],sigma=fit.pwm[2],xi=fit.pwm[3]);
			degp2s <- c(degp2s,degp2.pwm);
			if(plots){
				qegp2.pwm <- qegp2(p=c(1:length(data))/(length(data)+1),type=2,delta=fit.pwm[1],sigma=fit.pwm[2],xi=fit.pwm[3]);
				qegp2s <- c(qegp2s,qegp2.pwm);
				if(CI){
					q.pwm.boot <- mapply(FUN=qegp2,p=list(c(1:length(data))/(length(data)+1)),type=list(2),delta=as.list(fit.pwm.boot$t[,1]),sigma=as.list(fit.pwm.boot$t[,2]),xi=as.list(fit.pwm.boot$t[,3]));
					q.pwm.L <- apply(q.pwm.boot,1,quantile,0.025,na.rm=TRUE);
					q.pwm.U <- apply(q.pwm.boot,1,quantile,0.975,na.rm=TRUE);
				}
			}
		} else if(model==3){
			fit.pwm <- egp2.fit.pwm(x=data,type=3,kappa0=init[1],delta0=init[2],sigma0=init[3],xi0=init[4],censoring=censoring);
			if(CI){
				fit.pwm.boot <- boot::boot(data=data,statistic=egp2.fit.pwm.boot,R=R,type=3,kappa0=init[1],delta0=init[2],sigma0=init[3],xi0=init[4],censoring=censoring,parallel="multicore",ncpus=ncpus);
				CI.pwm.kappa <- boot::boot.ci(boot.out=fit.pwm.boot,index=1,type="perc")$perc[4:5];
				CI.pwm.delta <- boot::boot.ci(boot.out=fit.pwm.boot,index=2,type="perc")$perc[4:5];
				CI.pwm.sigma <- boot::boot.ci(boot.out=fit.pwm.boot,index=3,type="perc")$perc[4:5];
				CI.pwm.xi <- boot::boot.ci(boot.out=fit.pwm.boot,index=4,type="perc")$perc[4:5];
				CIs.pwm <- cbind(CI.pwm.kappa,CI.pwm.delta,CI.pwm.sigma,CI.pwm.xi);
			}
			degp2.pwm <- degp2(x=x,type=3,kappa=fit.pwm[1],delta=fit.pwm[2],sigma=fit.pwm[3],xi=fit.pwm[4]);
			degp2s <- c(degp2s,degp2.pwm);
			if(plots){
				qegp2.pwm <- qegp2(p=c(1:length(data))/(length(data)+1),type=3,kappa=fit.pwm[1],delta=fit.pwm[2],sigma=fit.pwm[3],xi=fit.pwm[4]);
				qegp2s <- c(qegp2s,qegp2.pwm);
				if(CI){
					q.pwm.boot <- mapply(FUN=qegp2,p=list(c(1:length(data))/(length(data)+1)),type=list(3),kappa=as.list(fit.pwm.boot$t[,1]),delta=as.list(fit.pwm.boot$t[,2]),sigma=as.list(fit.pwm.boot$t[,3]),xi=as.list(fit.pwm.boot$t[,4]));
					q.pwm.L <- apply(q.pwm.boot,1,quantile,0.025,na.rm=TRUE);
					q.pwm.U <- apply(q.pwm.boot,1,quantile,0.975,na.rm=TRUE);
				}
			}
		} else if(model==4){
			fit.pwm <- egp2.fit.pwm(x=data,type=4,prob0=init[1],kappa0=init[2],delta0=init[3],sigma0=init[4],xi0=init[5],censoring=censoring);
			if(CI){
				fit.pwm.boot <- boot::boot(data=data,statistic=egp2.fit.pwm.boot,R=R,type=4,prob0=init[1],kappa0=init[2],delta0=init[3],sigma0=init[4],xi0=init[5],censoring=censoring,parallel="multicore",ncpus=ncpus);
				CI.pwm.prob <- boot::boot.ci(boot.out=fit.pwm.boot,index=1,type="perc")$perc[4:5];
				CI.pwm.kappa <- boot::boot.ci(boot.out=fit.pwm.boot,index=2,type="perc")$perc[4:5];
				CI.pwm.delta <- boot::boot.ci(boot.out=fit.pwm.boot,index=3,type="perc")$perc[4:5];
				CI.pwm.sigma <- boot::boot.ci(boot.out=fit.pwm.boot,index=4,type="perc")$perc[4:5];
				CI.pwm.xi <- boot::boot.ci(boot.out=fit.pwm.boot,index=5,type="perc")$perc[4:5];
				CIs.pwm <- cbind(CI.pwm.prob,CI.pwm.kappa,CI.pwm.delta,CI.pwm.sigma,CI.pwm.xi);
			}
			degp2.pwm <- degp2(x=x,type=4,prob=fit.pwm[1],kappa=fit.pwm[2],delta=fit.pwm[3],sigma=fit.pwm[4],xi=fit.pwm[5]);
			degp2s <- c(degp2s,degp2.pwm);
			if(plots){
				qegp2.pwm <- qegp2(p=c(1:length(data))/(length(data)+1),type=4,prob=fit.pwm[1],kappa=fit.pwm[2],delta=fit.pwm[3],sigma=fit.pwm[4],xi=fit.pwm[5]);
				qegp2s <- c(qegp2s,qegp2.pwm);
				if(CI){
					q.pwm.boot <- mapply(FUN=qegp2,p=list(c(1:length(data))/(length(data)+1)),type=list(4),prob=as.list(fit.pwm.boot$t[,1]),kappa=as.list(fit.pwm.boot$t[,2]),delta=as.list(fit.pwm.boot$t[,3]),sigma=as.list(fit.pwm.boot$t[,4]),xi=as.list(fit.pwm.boot$t[,5]));
					q.pwm.L <- apply(q.pwm.boot,1,quantile,0.025,na.rm=TRUE);
					q.pwm.U <- apply(q.pwm.boot,1,quantile,0.975,na.rm=TRUE);
				}
			}
		}
	}

	if(any(method=="mle")){
		if(model==1){
			fit.mle <- egp2.fit.ml(x=data,type=1,kappa0=init[1],sigma0=init[2],xi0=init[3],censoring=censoring,rounded=rounded);
			if(CI){
				fit.mle.boot <- boot::boot(data=data,statistic=egp2.fit.ml.boot,R=R,type=1,kappa0=init[1],sigma0=init[2],xi0=init[3],censoring=censoring,rounded=rounded,parallel="multicore",ncpus=ncpus);
				CI.mle.kappa <- boot::boot.ci(boot.out=fit.mle.boot,index=1,type="perc")$perc[4:5];
				CI.mle.sigma <- boot::boot.ci(boot.out=fit.mle.boot,index=2,type="perc")$perc[4:5];
				CI.mle.xi <- boot::boot.ci(boot.out=fit.mle.boot,index=3,type="perc")$perc[4:5];
				CIs.mle <- cbind(CI.mle.kappa,CI.mle.sigma,CI.mle.xi);
			}
			degp2.mle <- degp2(x=x,type=1,kappa=fit.mle[1],sigma=fit.mle[2],xi=fit.mle[3]);
			degp2s <- c(degp2s,degp2.mle);
			if(plots){
				qegp2.mle <- qegp2(p=c(1:length(data))/(length(data)+1),type=1,kappa=fit.mle[1],sigma=fit.mle[2],xi=fit.mle[3]);
				qegp2s <- c(qegp2s,qegp2.mle);
				if(CI){
					q.mle.boot <- mapply(FUN=qegp2,p=list(c(1:length(data))/(length(data)+1)),type=list(1),kappa=as.list(fit.mle.boot$t[,1]),sigma=as.list(fit.mle.boot$t[,2]),xi=as.list(fit.mle.boot$t[,3]));
					q.mle.L <- apply(q.mle.boot,1,quantile,0.025,na.rm=TRUE);
					q.mle.U <- apply(q.mle.boot,1,quantile,0.975,na.rm=TRUE);
				}
			}
		} else if(model==2){
			fit.mle <- egp2.fit.ml(x=data,type=2,delta0=init[1],sigma0=init[2],xi0=init[3],censoring=censoring,rounded=rounded);
			if(CI){
				fit.mle.boot <- boot::boot(data=data,statistic=egp2.fit.ml.boot,R=R,type=2,delta0=init[1],sigma0=init[2],xi0=init[3],censoring=censoring,rounded=rounded,parallel="multicore",ncpus=ncpus);
				CI.mle.delta <- boot::boot.ci(boot.out=fit.mle.boot,index=1,type="perc")$perc[4:5];
				CI.mle.sigma <- boot::boot.ci(boot.out=fit.mle.boot,index=2,type="perc")$perc[4:5];
				CI.mle.xi <- boot::boot.ci(boot.out=fit.mle.boot,index=3,type="perc")$perc[4:5];
				CIs.mle <- cbind(CI.mle.delta,CI.mle.sigma,CI.mle.xi);
			}
			degp2.mle <- degp2(x=x,type=2,delta=fit.mle[1],sigma=fit.mle[2],xi=fit.mle[3]);
			degp2s <- c(degp2s,degp2.mle);
			if(plots){
				qegp2.mle <- qegp2(p=c(1:length(data))/(length(data)+1),type=2,delta=fit.mle[1],sigma=fit.mle[2],xi=fit.mle[3]);
				qegp2s <- c(qegp2s,qegp2.mle);
				if(CI){
					q.mle.boot <- mapply(FUN=qegp2,p=list(c(1:length(data))/(length(data)+1)),type=list(2),delta=as.list(fit.mle.boot$t[,1]),sigma=as.list(fit.mle.boot$t[,2]),xi=as.list(fit.mle.boot$t[,3]));
					q.mle.L <- apply(q.mle.boot,1,quantile,0.025,na.rm=TRUE);
					q.mle.U <- apply(q.mle.boot,1,quantile,0.975,na.rm=TRUE);
				}
			}
		} else if(model==3){
			fit.mle <- egp2.fit.ml(x=data,type=3,kappa0=init[1],delta0=init[2],sigma0=init[3],xi0=init[4],censoring=censoring,rounded=rounded);
			if(CI){
				fit.mle.boot <- boot::boot(data=data,statistic=egp2.fit.ml.boot,R=R,type=3,kappa0=init[1],delta0=init[2],sigma0=init[3],xi0=init[4],censoring=censoring,rounded=rounded,parallel="multicore",ncpus=ncpus);
				CI.mle.kappa <- boot::boot.ci(boot.out=fit.mle.boot,index=1,type="perc")$perc[4:5];
				CI.mle.delta <- boot::boot.ci(boot.out=fit.mle.boot,index=2,type="perc")$perc[4:5];
				CI.mle.sigma <- boot::boot.ci(boot.out=fit.mle.boot,index=3,type="perc")$perc[4:5];
				CI.mle.xi <- boot::boot.ci(boot.out=fit.mle.boot,index=4,type="perc")$perc[4:5];
				CIs.mle <- cbind(CI.mle.kappa,CI.mle.delta,CI.mle.sigma,CI.mle.xi);
			}
			degp2.mle <- degp2(x=x,type=3,kappa=fit.mle[1],delta=fit.mle[2],sigma=fit.mle[3],xi=fit.mle[4]);
			degp2s <- c(degp2s,degp2.mle);
			if(plots){
				qegp2.mle <- qegp2(p=c(1:length(data))/(length(data)+1),type=3,kappa=fit.mle[1],delta=fit.mle[2],sigma=fit.mle[3],xi=fit.mle[4]);
				qegp2s <- c(qegp2s,qegp2.mle);
				if(CI){
					q.mle.boot <- mapply(FUN=qegp2,p=list(c(1:length(data))/(length(data)+1)),type=list(3),kappa=as.list(fit.mle.boot$t[,1]),delta=as.list(fit.mle.boot$t[,2]),sigma=as.list(fit.mle.boot$t[,3]),xi=as.list(fit.mle.boot$t[,4]));
					q.mle.L <- apply(q.mle.boot,1,quantile,0.025,na.rm=TRUE);
					q.mle.U <- apply(q.mle.boot,1,quantile,0.975,na.rm=TRUE);
				}
			}
		} else if(model==4){
			fit.mle <- egp2.fit.ml(x=data,type=4,prob0=init[1],kappa0=init[2],delta0=init[3],sigma0=init[4],xi0=init[5],censoring=censoring,rounded=rounded);
			if(CI){
				fit.mle.boot <- boot::boot(data=data,statistic=egp2.fit.ml.boot,R=R,type=4,prob0=init[1],kappa0=init[2],delta0=init[3],sigma0=init[4],xi0=init[5],censoring=censoring,rounded=rounded,parallel="multicore",ncpus=ncpus);
				CI.mle.prob <- boot::boot.ci(boot.out=fit.mle.boot,index=1,type="perc")$perc[4:5];
				CI.mle.kappa <- boot::boot.ci(boot.out=fit.mle.boot,index=2,type="perc")$perc[4:5];
				CI.mle.delta <- boot::boot.ci(boot.out=fit.mle.boot,index=3,type="perc")$perc[4:5];
				CI.mle.sigma <- boot::boot.ci(boot.out=fit.mle.boot,index=4,type="perc")$perc[4:5];
				CI.mle.xi <- boot::boot.ci(boot.out=fit.mle.boot,index=5,type="perc")$perc[4:5];
				CIs.mle <- cbind(CI.mle.prob,CI.mle.kappa,CI.mle.delta,CI.mle.sigma,CI.mle.xi);
			}
			degp2.mle <- degp2(x=x,type=4,prob=fit.mle[1],kappa=fit.mle[2],delta=fit.mle[3],sigma=fit.mle[4],xi=fit.mle[5]);
			degp2s <- c(degp2s,degp2.mle);
			if(plots){
				qegp2.mle <- qegp2(p=c(1:length(data))/(length(data)+1),type=4,prob=fit.mle[1],kappa=fit.mle[2],delta=fit.mle[3],sigma=fit.mle[4],xi=fit.mle[5]);
				qegp2s <- c(qegp2s,qegp2.mle);
				if(CI){
					q.mle.boot <- mapply(FUN=qegp2,p=list(c(1:length(data))/(length(data)+1)),type=list(4),prob=as.list(fit.mle.boot$t[,1]),kappa=as.list(fit.mle.boot$t[,2]),delta=as.list(fit.mle.boot$t[,3]),sigma=as.list(fit.mle.boot$t[,4]),xi=as.list(fit.mle.boot$t[,5]));
					q.mle.L <- apply(q.mle.boot,1,quantile,0.025,na.rm=TRUE);
					q.mle.U <- apply(q.mle.boot,1,quantile,0.975,na.rm=TRUE);
				}
			}
		}
	}

	if(any(method=="pwm")){
		if(any(method=="mle")){
			fits <- list(pwm=fit.pwm,mle=fit.mle)
			if(CI){
				CIs <- list(pwm=CIs.pwm,mle=CIs.mle)
			} else{
				CIs <- NULL
			}
		} else{
			fits <- list(pwm=fit.pwm)
			if(CI & plots){
				CIs <- list(pwm=CIs.pwm)
			} else{
				CIs <- NULL
			}
		}
	} else{
		if(any(method=="mle")){
			fits <- list(mle=fit.mle)
			if(CI){
				CIs <- list(mle=CIs.mle)
			} else{
				CIs <- NULL
			}
		} else{
			fits <- NULL
			CIs <- NULL
		}
	}

	if(plots){
		mat <- matrix(c(1:(1+plots)),nrow=1,ncol=1+plots,byrow=TRUE);
		layout(mat);
		par(mar=c(4,4,1,1))
		hist(data,breaks=c(0:30)*max(data+10^-1)/30,freq=FALSE,xlim=range(x),xlab="Rainfall amounts [mm]",ylab="Density",main="",col="lightgrey");
		dens <- density(data,from=0); lines(dens$x,dens$y,col="black");
		abline(v=censoring,col="lightgrey");

		if(any(method=="pwm")){
			lines(x,degp2.pwm,col="red",lty=2);
		}
		if(any(method=="mle")){
			lines(x,degp2.mle,col="blue");
		}

		if(any(method=="pwm")){
			plot(data,qegp2.pwm,asp=1,xlab="Empirical quantiles",ylab="Fitted quantiles",
			     ylim=range(qegp2s,na.rm=TRUE),type="n");
		} else{
			if(any(method=="mle")){
				plot(data,qegp2.mle,asp=1,xlab="Empirical quantiles",ylab="Fitted quantiles",
				     ylim=range(qegp2s,na.rm=TRUE),type="n");
			}
		}
		if(CI){
			if(any(method=="pwm")){
				polygon(x=c(sort(data),sort(data,decreasing=TRUE)),
				        y=c(q.pwm.L,q.pwm.U[length(q.pwm.U):1]),col=rgb(1,0,0,alpha=0.1),
				        lty=2,border=rgb(1,0,0,alpha=0.5));
			}
			if(any(method=="mle")){
				polygon(x=c(sort(data),sort(data,decreasing=TRUE)),
				        y=c(q.mle.L,q.mle.U[length(q.mle.U):1]),col=rgb(0,0,1,alpha=0.1),
				        lty=1,border=rgb(0,0,1,alpha=0.5));
			}
		}
		if(any(method=="pwm")){
			lines(sort(data),qegp2.pwm,lty=2,type="b",pch=20,col="red");
		}
		if(any(method=="mle")){
			lines(sort(data),qegp2.mle,lty=1,type="b",pch=20,col="blue");
		}
		abline(0,1,col="lightgrey");
	}

	return(list(fit=fits,CI=CIs))
}


################################################################
###                                                          ###
### Distribution/Density/Quantile/Random generator functions ###
###                                                          ###
################################################################

# ### GPD
# qgpd <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
# {
#   inds <- p <=0 | p>=1
#   res <- rep(NA,length(p))
#   if (min(scale) < 0)
#     stop("invalid scale")
#   if (length(shape) != 1)
#     stop("invalid shape")
#   if (lower.tail)
#     p <- 1 - p
#   if (shape == 0)
#     res[!inds] <- loc - scale * log(p[!inds])
#   else res[!inds] <- loc + scale * (p[!inds]^(-shape) - 1)/shape
#
#   return(res)
# }
#
# qgpd.fullrange <- function (p, loc = 0, scale = 1, shape = 0, prob.loc=0.95)
# {
#   res <- loc + qgpd((p-prob.loc)/(1-prob.loc),scale=scale,shape=shape)
#   return(res)
# }



### different types of distribution G(u):
### type=0: G(u)=u
### type=1: G(u)=u^kappa
### type=2: G(u)=1-V_delta((1-u)^delta)
### type=3: G(u)=(1-V_delta((1-u)^delta))^(kappa/2)
### type=4: G(u)=p*u^kappa + (1-p)*u^delta



#' @title Carrier distribution for the extended GP distributions of Naveau et al.
#' @description Density, distribution function, quantile function and random number
#' generation for the carrier distributions of the extended Generalized Pareto distributions.
#' @name egp2.G
#' @seealso \code{\link{egp2}}
#' @param u vector of observations (\code{degp2.G}), probabilities (\code{qegp2.G}) or quantiles (\code{pegp2.G}), in \eqn{[0,1]}
#' @param prob mixture probability for model \code{type} \code{4}
#' @param kappa shape parameter for \code{type} \code{1}, \code{3} and \code{4}
#' @param delta additional parameter for \code{type} \code{2}, \code{3} and \code{4}
#' @param type integer between 0 to 5 giving the model choice
#' @param log logical; should the log-density be returned (default to \code{FALSE})?
#' @param n sample size
#' @param unifsamp sample of uniform; if provided, the data will be used in place of new uniform random variates
#' @param censoring numeric vector of length 2 containing the lower and upper bound for censoring
#' @param direct logical; which method to use for sampling in model of \code{type} \code{4}?
#' @author  Raphael Huser and Philippe Naveau
#'
#' @section Usage: \code{pegp2.G(u, type=1,prob, kappa, delta)}
#' @section Usage: \code{degp2.G(u,type=1, prob=NA, kappa=NA, delta=NA, log=FALSE)}
#' @section Usage: \code{qegp2.G(u,type=1, prob=NA,kappa=NA,delta=NA)}
#' @section Usage: \code{regp2.G(n,prob=NA,kappa=NA,delta=NA,type=1,unifsamp=NULL,direct=FALSE,censoring=c(0,1))}
NULL

#' @rdname egp2-functions
#' @export
#' @keywords internal
pegp2.G <- function(u, type=1, prob, kappa, delta){
  if(! type %in% 1:5){stop("Invalid `type' argument")}
  type <- type[1]
  if(type %in% c(1,3,4) && missing(kappa)){stop("Argument `kappa' missing.")}
  if(type %in% c(2,3,4) && missing(delta)){stop("Argument `delta' missing.")}
  if(type == 4 && missing(prob)){stop("Argument `prob' missing.")}
  if(type==0){
    return(u)
  } else if(type==1){
    return(u^kappa)
  } else if(type==2){
    return( 1 - pbeta((1-u)^delta,1/delta,2) )
  } else if(type==3){
    return( (1 - pbeta((1-u)^delta,1/delta,2))^(kappa/2) )
  } else if(type==4){
    return( prob*u^kappa + (1-prob)*u^delta )
  }
}

#' @rdname egp2-functions
#' @export
#' @keywords internal
degp2.G <- function(u,type=1, prob=NA, kappa=NA, delta=NA, log=FALSE){
  if(! type %in% 1:5){stop("Invalid `type' argument")}
  type <- type[1]
  if(type %in% c(1,3,4) && missing(kappa)){stop("Argument `kappa' missing.")}
  if(type %in% c(2,3,4) && missing(delta)){stop("Argument `delta' missing.")}
  if(type == 4 && missing(prob)){stop("Argument `prob' missing.")}
  if(log==FALSE){
    if(type==0){
      return(1)
    } else if(type==1){
      return(kappa*u^(kappa-1))
    } else if(type==2){
      return( dbeta((1-u)^delta,1/delta,2)*delta*(1-u)^(delta-1) )
    } else if(type==3){
      return( (kappa/2)*(1 - pbeta((1-u)^delta,1/delta,2))^(kappa/2-1)*dbeta((1-u)^delta,1/delta,2)*
                delta*(1-u)^(delta-1) )
    } else if(type==4){
      return( prob*kappa*u^(kappa-1) + (1-prob)*delta*u^(delta-1) )
    }
  } else{
    if(type==0){
      return(0)
    } else if(type==1){
      return(log(kappa) + (kappa-1)*log(u))
    } else if(type==2){
      return( dbeta((1-u)^delta,1/delta,2,log=TRUE) + log(delta) + (delta-1)*log(1-u) )
    } else if(type==3){
      return( log(kappa/2) + (kappa/2-1)*log(1 - pbeta((1-u)^delta,1/delta,2)) +
                dbeta((1-u)^delta,1/delta,2,log=TRUE) + log(delta) + (delta-1)*log(1-u) )
    } else if(type==4){
      return( log(prob*kappa*u^(kappa-1) + (1-prob)*delta*u^(delta-1)) )
    }
  }
}

#' @rdname egp2-functions
#' @export
#' @keywords internal
qegp2.G <- function(u,type=1, prob=NA,kappa=NA,delta=NA){
  if(! type %in% 1:5){stop("Invalid `type' argument")}
  type <- type[1]
  if(type %in% c(1,3,4) && missing(kappa)){stop("Argument `kappa' missing.")}
  if(type %in% c(2,3,4) && missing(delta)){stop("Argument `delta' missing.")}
  if(type == 4 && missing(prob)){stop("Argument `prob' missing.")}
  if(type==0){
    return(u)
  } else if(type==1){
    return(u^(1/kappa))
  } else if(type==2){
    return( 1-qbeta(1-u,1/delta,2)^(1/delta) )
  } else if(type==3){
    return( 1-qbeta(1-u^(2/kappa),1/delta,2)^(1/delta) )
  } else if(type==4){
    dummy.func <- function(u,p,prob=NA,kappa=NA,delta=NA){
      return( pegp2.G(u=u,prob=prob,kappa=kappa,delta=delta,type=4)-p )
    }
    find.root <- function(u,prob=NA,kappa=NA,delta=NA){
      return( uniroot(dummy.func,interval=c(0,1),p=u,prob=prob,kappa=kappa,delta=delta)$root )
    }
    return( sapply(u,FUN=find.root,prob=prob,kappa=kappa,delta=delta) )
  }
}

#' @rdname egp2-functions
#' @export
#' @keywords internal
regp2.G <- function(n,prob=NA,kappa=NA,delta=NA,type=1,unifsamp=NULL,direct=FALSE,censoring=c(0,1)){
  if(! type %in% 1:5){stop("Invalid `type' argument")}
  type <- type[1]
  if(type %in% c(1,3,4) && missing(kappa)){stop("Argument `kappa' missing.")}
  if(type %in% c(2,3,4) && missing(delta)){stop("Argument `delta' missing.")}
  if(type == 4 && missing(prob)){stop("Argument `prob' missing.")}
  if(is.null(unifsamp)){
    unifsamp <- runif(n)
  }
  if(type!=4 | (type==4 & direct)){
    if(censoring[1]==0 & censoring[2]==1){
      return( qegp2.G(unifsamp,prob,kappa,delta,type) )
    } else{
      x.L <- pegp2.G(censoring[1],prob=prob,kappa=kappa,delta=delta,type=type)
      x.U <- pegp2.G(censoring[2],prob=prob,kappa=kappa,delta=delta,type=type)
      return( qegp2.G(x.L+(x.U-x.L)*unifsamp,prob=prob,kappa=kappa,delta=delta,type=type) )
    }
  } else if(type==4 & !direct){
    if(censoring[1]==0 & censoring[2]==1){
      components <- sample(x=c(1,2),size=n,replace=TRUE,prob=c(prob,1-prob))
      res <- c()
      res[components==1] <- qegp2.G(unifsamp[components==1],prob=NA,kappa=kappa,delta=NA,type=1)
      res[components==2] <- qegp2.G(unifsamp[components==2],prob=NA,kappa=delta,delta=NA,type=1)
      return(res)
    } else{
      x.L <- pegp2.G(censoring[1],prob=prob,kappa=kappa,delta=delta,type=type)
      x.U <- pegp2.G(censoring[2],prob=prob,kappa=kappa,delta=delta,type=type)
      prop1 <- prob*(pegp2.G(censoring[2],prob=NA,kappa=kappa,delta=NA,type=1)-
                       pegp2.G(censoring[1],prob=NA,kappa=kappa,delta=NA,type=1))
      prop2 <- (1-prob)*(pegp2.G(censoring[2],prob=NA,kappa=delta,delta=NA,type=1)-
                           pegp2.G(censoring[1],prob=NA,kappa=delta,delta=NA,type=1))
      new.prob <- prop1/(prop1+prop2)
      components <- sample(x=c(1,2),size=n,replace=TRUE,prob=c(new.prob,1-new.prob))
      res <- c()
      res[components==1] <- qegp2.G(x.L+(x.U-x.L)*unifsamp[components==1],prob=NA,kappa=kappa,delta=NA,type=1)
      res[components==2] <- qegp2.G(x.L+(x.U-x.L)*unifsamp[components==2],prob=NA,kappa=delta,delta=NA,type=1)
      return(res)
    }
  }
}


### different types of extended GP
### type=0: exact GP ---> 2 parameters (sigma,xi)
### type=1: egp2 with G(u)=u^kappa ---> 3 parameters (kappa,sigma,xi)
### type=2: egp2 with G(u)=1-V_delta((1-u)^delta) ---> 3 parameters (delta,sigma,xi)
### type=3: egp2 with G(u)=(1-V_delta((1-u)^delta))^(kappa/2) ---> 4 parameters (kappa,delta,sigma,xi)
### type=4: egp2 with mixture G(u)=p*u^kappa + (1-p)*u^delta ---> 5 parameters (p,kappa,delta,sigma,xi)

#' @title Extended generalised Pareto families of Naveau et al. (2016)
#'
#' @description Density function, distribution function, quantile function and random generation for the extended generalized Pareto distribution (GPD) with scale and shape parameters.
#'
#' @details The extended generalized Pareto families proposed in Naveau \emph{et al.} (2016)
#'retain the tail index of the distribution while being compliant with the theoretical behavior of extreme
#'low rainfall. There are five proposals, the first one being equivalent to the GP distribution.
#'\itemize{
#' \item \code{type} 0 corresponds to uniform carrier, \eqn{G(u)=u}.
#' \item \code{type} 1 corresponds to a three parameters family, with carrier \eqn{G(u)=u^\kappa}.
#' \item \code{type} 2 corresponds to a three parameters family, with carrier \eqn{G(u)=1-V_\delta((1-u)^\delta)}.
#' \item \code{type} 3 corresponds to a four parameters family, with carrier \deqn{G(u)=1-V_\delta((1-u)^\delta))^{\kappa/2}}.
#' \item \code{type} 4 corresponds to a five parameter model (a mixture of \code{type} 2, with \eqn{G(u)=pu^\kappa + (1-p)*u^\delta}
#' }
#' @name egp2
#' @param q vector of quantiles
#' @param x vector of observations
#' @param p vector of probabilities
#' @param n sample size
#' @param prob mixture probability for model \code{type} \code{4}
#' @param kappa shape parameter for \code{type} \code{1}, \code{3} and \code{4}
#' @param delta additional parameter for \code{type} \code{2}, \code{3} and \code{4}
#' @param sigma scale parameter
#' @param xi shape parameter
#' @param type integer between 0 to 5 giving the model choice
#' @param log logical; should the log-density be returned (default to \code{FALSE})?
#' @param unifsamp sample of uniform; if provided, the data will be used in place of new uniform random variates
#' @param censoring numeric vector of length 2 containing the lower and upper bound for censoring
#'
#' @section Usage: \code{pegp2(q, prob=NA, kappa=NA, delta=NA, sigma=NA, xi=NA, type=1)}
#' @section Usage: \code{degp2(x, prob=NA, kappa=NA, delta=NA, sigma=NA, xi=NA, type=1, log=FALSE)}
#' @section Usage: \code{qegp2(p, prob=NA, kappa=NA, delta=NA, sigma=NA, xi=NA, type=1)}
#' @section Usage: \code{regp2(n, prob=NA, kappa=NA, delta=NA, sigma=NA, xi=NA, type=1, unifsamp=NULL, censoring=c(0,Inf))}
#' @author  Raphael Huser and Philippe Naveau
#' @references Naveau, P., R. Huser, P. Ribereau, and A. Hannart (2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, \emph{Water Resour. Res.}, 52, 2753-2769, \code{doi:10.1002/2015WR018552}.
NULL


#' EGP functions
#'
#' These functions are documented in \code{\link{egp2}} and in \code{link{egp2.G}} for the carrier distributions supported in the unit interval.
#' @name egp2-functions
#' @export
#' @seealso \code{\link{egp2}}, \code{\link{egp2.G}}
#' @keywords internal
pegp2 <- function(q,prob=NA,kappa=NA,delta=NA,sigma=NA,xi=NA,type=1){
  return( pegp2.G(pgpd(q,scale=sigma,shape=xi),prob=prob,kappa=kappa,delta=delta,type=type) )
}

#' @rdname egp2-functions
#' @export
#' @keywords internal
degp2 <- function(x,prob=NA,kappa=NA,delta=NA,sigma=NA,xi=NA,type=1,log=FALSE){
  if(log==FALSE){
    return( degp2.G(pgpd(x,scale=sigma,shape=xi),prob=prob,kappa=kappa,delta=delta,type=type)*dgpd(x,scale=sigma,shape=xi) )
  } else{
    return( degp2.G(pgpd(x,scale=sigma,shape=xi),prob=prob,kappa=kappa,delta=delta,type=type,log=TRUE) + dgpd(x,scale=sigma,shape=xi,log=TRUE) )
  }
}

#' @rdname egp2-functions
#' @export
#' @keywords internal
qegp2 <- function(p,prob=NA,kappa=NA,delta=NA,sigma=NA,xi=NA,type=1){
  return( qgpd(qegp2.G(p,prob=prob,kappa=kappa,delta=delta,type=type),scale=sigma,shape=xi) )
}

#' @rdname egp2-functions
#' @export
#' @keywords internal
regp2 <- function(n,prob=NA,kappa=NA,delta=NA,sigma=NA,xi=NA,type=1,unifsamp=NULL,censoring=c(0,Inf)){
  if(censoring[1]==0 & censoring[2]==Inf){
    return( qgpd(regp2.G(n,prob=prob,kappa=kappa,delta=delta,type=type,unifsamp,censoring=c(0,1)),scale=sigma,shape=xi) )
  } else{
    H.L <- pgpd(censoring[1],loc=0,scale=sigma,shape=xi)
    H.U <- pgpd(censoring[2],loc=0,scale=sigma,shape=xi)
    return( qgpd(regp2.G(n,prob=prob,kappa=kappa,delta=delta,type=type,unifsamp,censoring=c(H.L,H.U)),scale=sigma,shape=xi) )
  }
}


#########################
###                   ###
### Inference via PWM ###
###                   ###
#########################

egp2.pwm <- function(orders,prob=NA,kappa=NA,delta=NA,sigma=NA,xi=NA,type=1,
                    censoring=c(0,Inf),empiric=FALSE,unifsamp=NULL,NbSamples=10^4,N=200){
  H.L <- ifelse(censoring[1]>0,pgpd(censoring[1],loc=0,scale=sigma,shape=xi),0)
  H.U <- ifelse(censoring[2]<Inf,pgpd(censoring[2],loc=0,scale=sigma,shape=xi),1)

  F.L <- ifelse(censoring[1]>0,pegp2(censoring[1],prob=prob,kappa=kappa,delta=delta,sigma=sigma,xi=xi,type=type),0)
  F.U <- ifelse(censoring[2]<Inf,pegp2(censoring[2],prob=prob,kappa=kappa,delta=delta,sigma=sigma,xi=xi,type=type),1)
  prob.LU <- F.U-F.L

  if(!empiric){
    E2 <- ((1-F.L)^(orders+1)-(1-F.U)^(orders+1))/(prob.LU*(orders+1))

    if(type==0){
      E1 <- (H.U-H.L)^(orders-xi)/(orders-xi+1)
      return( (sigma/xi)*(E1-E2) )
    } else if(type==1){
      beta.coefs <- beta(c(1:(max(orders)+1))*kappa,1-xi)
      probs.beta <- pbeta(H.U,shape1=c(1:(max(orders)+1))*kappa,shape2=1-xi)-pbeta(H.L,shape1=c(1:(max(orders)+1))*kappa,shape2=1-xi)
      E1 <- c()
      for(i in 1:length(orders)){
        E1[i] <- (kappa/prob.LU)*sum((-1)^c(0:orders[i])*choose(orders[i],c(0:orders[i]))*beta.coefs[1:(orders[i]+1)]*probs.beta[1:(orders[i]+1)])
      }
      return( (sigma/xi)*(E1-E2) )
    } else if(type==2){
      E1 <- c()
      for(i in 1:length(orders)){
        expos1 <- orders[i]-xi+delta*c(0:orders[i])+1
        expos2 <- orders[i]-xi+delta*(c(0:orders[i])+1)+1

        numer.j <- expos1*((1-H.U)^expos2-(1-H.L)^expos2) + expos2*((1-H.L)^expos1-(1-H.U)^expos1)
        denom.j <- expos1*expos2

        E1[i] <- (1/(delta^(orders[i]+1)*prob.LU))*sum( (-1)^c(0:orders[i])*choose(orders[i],c(0:orders[i]))*(1+delta)^(orders[i]-c(0:orders[i])+1)*numer.j/denom.j )
      }
      return( (sigma/xi)*(E1-E2) )
    } else if(type==3){
      print("No explicit formula for PWMs for this model...")
      return( NA )
    } else if(type==4){
      E1 <- c()
      for(i in 1:length(orders)){
        l.vals <- matrix(c(0:orders[i]),nrow=orders[i]+1,ncol=orders[i]+1,byrow=TRUE)
        m.vals <- matrix(c(0:orders[i]),nrow=orders[i]+1,ncol=orders[i]+1)
        inds <- l.vals<m.vals
        l.vals[inds] <- m.vals[inds] <- NA
        const <- choose(orders[i],l.vals)*choose(l.vals,m.vals)*(-1)^l.vals*prob^m.vals*(1-prob)^(l.vals-m.vals)
        beta.coefs1 <- beta(kappa*(m.vals+1)+delta*(l.vals-m.vals),1-xi)
        beta.coefs2 <- beta(kappa*m.vals+delta*(l.vals-m.vals+1),1-xi)
        probs.beta1 <- pbeta(H.U,shape1=kappa*(m.vals+1)+delta*(l.vals-m.vals),shape2=1-xi)-
          pbeta(H.L,shape1=kappa*(m.vals+1)+delta*(l.vals-m.vals),shape2=1-xi)
        probs.beta2 <- pbeta(H.U,shape1=kappa*m.vals+delta*(l.vals-m.vals+1),shape2=1-xi)-
          pbeta(H.L,shape1=kappa*m.vals+delta*(l.vals-m.vals+1),shape2=1-xi)
        E1[i] <- (1/prob.LU)*sum(const*(prob*kappa*beta.coefs1*probs.beta1 +
                                          (1-prob)*delta*beta.coefs2*probs.beta2),na.rm=TRUE)
      }
      return( (sigma/xi)*(E1-E2) )
    }
  } else{
    if(is.null(unifsamp)){
      unifsamp <- runif(NbSamples)
    }
    if(censoring[1]==0 & censoring[2]==Inf){
      V <- regp2.G(length(unifsamp),prob,kappa,delta,type,unifsamp,direct=TRUE,censoring=c(0,1))
    } else{
      V <- regp2.G(length(unifsamp),prob,kappa,delta,type,unifsamp,direct=TRUE,censoring=c(H.L,H.U))
    }
    X <- qgpd(V,scale=sigma,shape=xi)
    res <- c()
    for(i in 1:length(orders)){
      res[i] <- mean(X*(1-(F.L + (F.U-F.L)*unifsamp))^orders[i],na.rm=TRUE)
    }
    return( res )
  }
}

egp2.fit.pwm <- function(x,type=1,prob0=NA,kappa0=NA,delta0=NA,sigma0=NA,xi0=NA,censoring=c(0,Inf),empiric=FALSE,unifsamp=NULL,NbSamples=10^4){
  Fn = ecdf(x)
  inds <- x>censoring[1] & x<censoring[2]
  mu0hat = mean(x[inds])
  mu1hat = mean(x[inds]*(1-Fn(x[inds])))
  mu2hat = mean(x[inds]*(1-Fn(x[inds]))^2)
  mu3hat = mean(x[inds]*(1-Fn(x[inds]))^3)
  mu4hat = mean(x[inds]*(1-Fn(x[inds]))^4)

  if(type==0){
    if(censoring[1]==0 & censoring[2]==Inf){
      xihat <- (mu0hat-4*mu1hat)/(mu0hat-2*mu1hat)
      sigmahat <- mu0hat*(1-xihat)
      thetahat <- c(sigmahat,xihat)
    } else{
      fct <- function(theta,x){
        pwm.theor <- egp2.pwm(orders=c(0:1),sigma=theta[1],xi=theta[2],type=0,censoring=censoring,empiric=empiric,unifsamp=unifsamp,NbSamples=NbSamples)
        pwm.empir <- c(mu0hat,mu1hat)
        return(matrix(pwm.theor - pwm.empir,ncol=2))
      }
      theta0 <- c(sigma0,xi0)
      res <- gmm::gmm(fct, x, theta0, optfct = "nlminb", lower = c(0.0001, 10^(-6)),  upper = c(Inf, .99),vcov="iid")
      thetahat <- res$coefficients
    }
    names(thetahat) <- c("sigma","xi")
    return(thetahat)
  } else if(type==1){
    fct <- function(theta,x){
      pwm.theor <- egp2.pwm(orders=c(0:2),kappa=theta[1],sigma=theta[2],xi=theta[3],type=1,censoring=censoring,empiric=empiric,unifsamp=unifsamp,NbSamples=NbSamples)
      pwm.empir <- c(mu0hat,mu1hat,mu2hat)
      return(matrix(pwm.theor - pwm.empir,ncol=3))
    }
    theta0 <- c(kappa0,sigma0,xi0)
    res <- gmm::gmm(fct, x, theta0, optfct = "nlminb", lower = c(0.0001, 0.0001, 10^(-6)),  upper = c(Inf,Inf, .99),vcov="iid")
    thetahat <- res$coefficients
    names(thetahat) <- c("kappa","sigma","xi")
    return(thetahat)
  } else if(type==2){
    fct <- function(theta,x){
      pwm.theor <- egp2.pwm(orders=c(0:2),delta=theta[1],sigma=theta[2],xi=theta[3],
                           type=2,censoring=censoring,empiric=empiric,unifsamp=unifsamp,NbSamples=NbSamples)
      pwm.empir <- c(mu0hat,mu1hat,mu2hat)
      return(matrix(pwm.theor - pwm.empir,ncol=3))
    }
    theta0 <- c(delta0,sigma0,xi0)
    res <- gmm::gmm(fct, x, theta0, optfct = "nlminb", lower = c(0.0001, 0.0001, 10^(-6)),
               upper = c(100,Inf, .99),vcov="iid")
    thetahat <- res$coefficients
    names(thetahat) <- c("delta","sigma","xi")
    return(thetahat)
  } else if(type==3){
    fct <- function(theta,x){
      pwm.theor <- egp2.pwm(orders=c(0:3),kappa=theta[1],delta=theta[2],sigma=theta[3],
                           xi=theta[4],type=3,censoring=censoring,empiric=empiric,unifsamp=unifsamp,NbSamples=NbSamples)
      pwm.empir <- c(mu0hat,mu1hat,mu2hat,mu3hat)
      return(matrix(pwm.theor - pwm.empir,ncol=4))
    }
    theta0 <- c(kappa0,delta0,sigma0,xi0)
    res <- gmm::gmm(fct, x, theta0, optfct = "nlminb", lower = c(0.0001, 0.0001, 0.0001, 10^(-6)),
               upper = c(Inf,100,Inf, .99),vcov="iid")
    thetahat <- res$coefficients
    names(thetahat) <- c("kappa","delta","sigma","xi")
    return(thetahat)
  } else if(type==4){
    fct <- function(theta,x){
      pwm.theor <- egp2.pwm(orders=c(0:4),prob=theta[1],kappa=theta[2],delta=theta[2]+theta[3],sigma=theta[4],xi=theta[5],type=4,censoring=censoring,empiric=empiric,unifsamp=unifsamp,NbSamples=NbSamples)
      pwm.empir <- c(mu0hat,mu1hat,mu2hat,mu3hat,mu4hat)
      return(matrix(pwm.theor - pwm.empir,ncol=5))
    }
    res0 <- egp2.fit.pwm(x[x>quantile(x,0.9)]-quantile(x,0.9),type=0)
    if(is.na(prob0)){ prob0 <- 0.5 }
    if(is.na(kappa0)){ kappa0 <- mean(x[x<quantile(x,0.1)])/(quantile(x,0.1)-mean(x[x<quantile(x,0.1)])) }
    if(is.na(delta0)){ delta0 <- kappa0+0.01}
    Ddelta0 <- delta0-kappa0
    if(is.na(sigma0)){ sigma0 <- res0[1] }
    if(is.na(xi0)){ xi0 <- res0[2] }
    theta0 <- c(prob0,kappa0,Ddelta0,sigma0,xi0)
    names(theta0) <- c("prob","kappa","Ddelta","sigma","xi")
    #print(theta0)
    res <- gmm::gmm(fct, x, theta0, optfct = "nlminb", lower = c(0, 0.0001, 0, 0.0001, 10^(-6)),
               upper = c(1, Inf, Inf,Inf, .99),vcov="iid")
    thetahat <- res$coefficients; thetahat[3] <- thetahat[2] + thetahat[3]
    names(thetahat) <- c("prob","kappa","delta","sigma","xi")
    return(thetahat)
  }
}


egp2.fit.pwm.boot <- function(data,i,type=1,prob0=NA,kappa0=NA,delta0=NA,sigma0=NA,xi0=NA,
                             censoring=c(0,Inf),empiric=FALSE,unifsamp=NULL,NbSamples=10^4){
  return( egp2.fit.pwm(data[i],type=type,prob0=prob0,kappa0=kappa0,delta0=delta0,
                      sigma0=sigma0,xi0=xi0,censoring=censoring,empiric=empiric,
                      unifsamp=unifsamp,NbSamples=NbSamples) )
}


########################
###                  ###
### Inference via ML ###
###                  ###
########################

egp2.nll <- function(theta,x,censoring,rounded,type){
  x.cens1 <- x[x<censoring[1]]
  x.cens2 <- x[x>censoring[2]]
  x.not.cens <- x[x>=censoring[1] & x<=censoring[2]]
  if(type==0){
    if(theta[1]<=0 | theta[2]<=10^(-6) | theta[2]>0.99 ){
      return(Inf)
    } else{
      censor1 <- ifelse(censoring[1]>0,pegp2(censoring[1],sigma=theta[1],xi=theta[2],type=0),0);
      censor2 <- ifelse(censoring[2]<Inf,pegp2(censoring[2],sigma=theta[1],xi=theta[2],type=0),1);
      contrib.cens1 <- ifelse(length(x.cens1)>0,length(x.cens1)*log(censor1),0);
      contrib.cens2 <- ifelse(length(x.cens2)>0,length(x.cens2)*log(1-censor2),0);
      contrib.not.cens <- ifelse(rounded==0,sum(degp2(x.not.cens,sigma=theta[1],xi=theta[2],type=0,log=TRUE),na.rm=TRUE),sum(log(pegp2(x.not.cens+rounded,sigma=theta[1],xi=theta[2],type=0)-pegp2(x.not.cens,sigma=theta[1],xi=theta[2],type=0))));

      return( -(contrib.cens1+contrib.not.cens+contrib.cens2) )
    }
  } else if(type==1){
    if(theta[1]<=0 | theta[2]<=0  | theta[3]<=10^(-6) | theta[3]>0.99){
      return(Inf)
    } else{
      censor1 <- ifelse(censoring[1]>0,pegp2(censoring[1],kappa=theta[1],sigma=theta[2],xi=theta[3],type=1),0);
      censor2 <- ifelse(censoring[2]<Inf,pegp2(censoring[2],kappa=theta[1],sigma=theta[2],xi=theta[3],type=1),1);
      contrib.cens1 <- ifelse(length(x.cens1)>0,length(x.cens1)*log(censor1),0);
      contrib.cens2 <- ifelse(length(x.cens2)>0,length(x.cens2)*log(1-censor2),0);
      contrib.not.cens <- ifelse(rounded==0,sum(degp2(x.not.cens,kappa=theta[1],sigma=theta[2],xi=theta[3],type=1,log=TRUE),na.rm=TRUE),sum(log(pegp2(x.not.cens+rounded,kappa=theta[1],sigma=theta[2],xi=theta[3],type=1)-pegp2(x.not.cens,kappa=theta[1],sigma=theta[2],xi=theta[3],type=1))));

      return( -(contrib.cens1+contrib.not.cens+contrib.cens2) )
    }
  } else if(type==2){
    if(theta[1]<=0 | theta[1]>100 | theta[2]<=0 | theta[3]<=10^(-6) | theta[3]>0.99){
      return(Inf)
    } else{
      censor1 <- ifelse(censoring[1]>0,pegp2(censoring[1],delta=theta[1],sigma=theta[2],xi=theta[3],type=2),0);
      censor2 <- ifelse(censoring[2]<Inf,pegp2(censoring[2],delta=theta[1],sigma=theta[2],xi=theta[3],type=2),1);
      contrib.cens1 <- ifelse(length(x.cens1)>0,length(x.cens1)*log(censor1),0);
      contrib.cens2 <- ifelse(length(x.cens2)>0,length(x.cens2)*log(1-censor2),0);
      contrib.not.cens <- ifelse(rounded==0,
                                 sum(degp2(x.not.cens,delta=theta[1],sigma=theta[2],
                                          xi=theta[3],type=2,log=TRUE),na.rm=TRUE),
                                 sum(log(pegp2(x.not.cens+rounded,delta=theta[1],sigma=theta[2],
                                              xi=theta[3],type=2)-
                                           pegp2(x.not.cens,delta=theta[1],sigma=theta[2],
                                                xi=theta[3],type=2))));

      return( -(contrib.cens1+contrib.not.cens+contrib.cens2) )
    }
  } else if(type==3){
    if(theta[1]<=0 | theta[2]<=0 | theta[2]>100 | theta[3]<=0 | theta[4]<=10^(-6) | theta[4]>0.99){
      return(Inf)
    } else{
      censor1 <- ifelse(censoring[1]>0,pegp2(censoring[1],kappa=theta[1],delta=theta[2],sigma=theta[3],xi=theta[4],type=3),0);
      censor2 <- ifelse(censoring[2]<Inf,pegp2(censoring[2],kappa=theta[1],delta=theta[2],sigma=theta[3],xi=theta[4],type=3),1);
      contrib.cens1 <- ifelse(length(x.cens1)>0,length(x.cens1)*log(censor1),0);
      contrib.cens2 <- ifelse(length(x.cens2)>0,length(x.cens2)*log(1-censor2),0);
      contrib.not.cens <- ifelse(rounded==0,sum(degp2(x.not.cens,kappa=theta[1],delta=theta[2],
                                                     sigma=theta[3],xi=theta[4],type=3,log=TRUE),na.rm=TRUE),
                                 sum(log(pegp2(x.not.cens+rounded,kappa=theta[1],delta=theta[2],
                                              sigma=theta[3],xi=theta[4],type=3)-
                                           pegp2(x.not.cens,kappa=theta[1],delta=theta[2],
                                                sigma=theta[3],xi=theta[4],type=3))));

      return( -(contrib.cens1+contrib.not.cens+contrib.cens2) )
    }
  } else if(type==4){
    if(theta[1]<0 | theta[1]>1 | theta[2]<=0 | theta[3]<=0 | theta[4]<=0 | theta[5]<=10^(-6) | theta[5]>0.99 | theta[3]<theta[2]){
      return(Inf)
    } else{
      censor1 <- ifelse(censoring[1]>0,pegp2(censoring[1],prob=theta[1],kappa=theta[2],
                                            delta=theta[3],sigma=theta[4],xi=theta[5],type=4),0);
      censor2 <- ifelse(censoring[2]<Inf,pegp2(censoring[2],prob=theta[1],kappa=theta[2],
                                              delta=theta[3],sigma=theta[4],xi=theta[5],type=4),1);
      contrib.cens1 <- ifelse(length(x.cens1)>0,length(x.cens1)*log(censor1),0);
      contrib.cens2 <- ifelse(length(x.cens2)>0,length(x.cens2)*log(1-censor2),0);
      contrib.not.cens <- ifelse(rounded==0,
                                 sum(degp2(x.not.cens,prob=theta[1],kappa=theta[2],
                                          delta=theta[3],sigma=theta[4],xi=theta[5],type=4,log=TRUE),na.rm=TRUE),
                                 sum(log(pegp2(x.not.cens+rounded,prob=theta[1],kappa=theta[2],
                                              delta=theta[3],sigma=theta[4],xi=theta[5],type=4)-
                                           pegp2(x.not.cens,prob=theta[1],kappa=theta[2],delta=theta[3],
                                                sigma=theta[4],xi=theta[5],type=4))));

      return( -(contrib.cens1+contrib.not.cens+contrib.cens2) )
    }
  }
}

egp2.fit.ml <- function(x,type=1,prob0=NA,kappa0=NA,delta0=NA,sigma0=NA,xi0=NA,
                       censoring=c(0,Inf),rounded=0.1,print=FALSE){
  if(type==0){
    theta0 <- c(sigma0,xi0)
    opt <- optim(par=theta0,fn=egp2.nll,x=x,censoring=censoring,rounded=rounded,
                 type=type,method="Nelder-Mead",control=list(maxit=1000),hessian=FALSE)
    if(print){
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("sigma","xi")
    return(thetahat)
  } else if(type==1){
    theta0 <- c(kappa0,sigma0,xi0)
    opt <- optim(par=theta0,fn=egp2.nll,x=x,censoring=censoring,rounded=rounded,
                 type=type,method="Nelder-Mead",control=list(maxit=1000),hessian=FALSE)
    if(print){
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("kappa","sigma","xi")
    return(thetahat)
  } else if(type==2){
    theta0 <- c(delta0,sigma0,xi0)
    opt <- optim(par=theta0,fn=egp2.nll,x=x,censoring=censoring,rounded=rounded,
                 type=type,method="Nelder-Mead",control=list(maxit=1000),hessian=FALSE)
    if(print){
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("delta","sigma","xi")
    return(thetahat)
  } else if(type==3){
    theta0 <- c(kappa0,delta0,sigma0,xi0)
    opt <- optim(par=theta0,fn=egp2.nll,x=x,censoring=censoring,rounded=rounded,
                 type=type,method="Nelder-Mead",control=list(maxit=1000),hessian=FALSE)
    if(print){
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("kappa","delta","sigma","xi")
    return(thetahat)
  } else if(type==4){
    theta0 <- c(prob0,kappa0,delta0,sigma0,xi0)
    opt <- optim(par=theta0,fn=egp2.nll,x=x,censoring=censoring,rounded=rounded,
                 type=type,method="Nelder-Mead",control=list(maxit=1000),hessian=FALSE)
    if(print){
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("prob","kappa","delta","sigma","xi")
    return(thetahat)
  }
}

egp2.fit.ml.boot <- function(data,i,type=1,prob0=NA,kappa0=NA,delta0=NA,sigma0=NA,
                            xi0=NA,censoring=c(0,Inf),rounded=0.1,print=FALSE){
  return( egp2.fit.ml(data[i],type=type,prob0=prob0,kappa0=kappa0,delta0=delta0,
                     sigma0=sigma0,xi0=xi0,censoring=censoring,rounded=rounded,print=print) )
}

