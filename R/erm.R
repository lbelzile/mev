#' Exponential regression estimator
#'
#'
erm <- function(xdat, k, method = c("bdgm","fh")){
 method <- match.arg(method)
 k <- sort(as.integer(k))
 xdat <- sort(xdat[is.finite(xdat)], decreasing = TRUE)
 n <- length(xdat)
 kmax <- k[length(k)]
 xdat <- as.numeric(xdat[1:(kmax + 1)])
 stopifnot(xdat[(kmax+1)] > 0)
 logdata <- log(xdat)
 shape <- rho <- b <- numeric(length = length(k))
 if(method == "bdgm"){
   expreg <- function(par, z){
     k <- length(z)
     shape <- par[1]
     bn <- par[2]
     rho <- par[3]
     scale <- shape + bn*((1:k)/(k+1))^(-rho)
     if(isTRUE(any(scale < 0))){
       return(1e10)
     }
     -sum(dexp(z, rate = 1/scale, log = TRUE))
   }
   ineqfn <- function(par, z){
     k <- length(z)
     par[1] + par[2]*(c(1,k)/(k+1))^(-par[3])
   }
 for(i in seq_along(k)){
   ks <- k[i]
   Z <- 1:ks * (logdata[1:ks] - logdata[2:(ks+1)])
   opt <- Rsolnp::solnp(pars = c(0.25,1,-0.2),
                 fun = expreg,
                 LB = c(0, 0, -1.5),
                 UB = c(2, 10, 0),
                 ineqfun = ineqfn,
                 ineqLB = rep(0, 2),
                 ineqUB = rep(1e10, 2),
                 z = Z,
                 control = list(trace = 0))
   b[i] <- opt$par[2]
   shape[i] <- opt$par[1]
   rho[i] <- opt$par[3]
 }
 } else if(method == "fh"){
     expreg <- function(par, z){
       k <- length(z)
       shape <- par[1]
       bn <- par[2]
       rho <- par[3]
       scale <- shape*exp(bn*((1:k)/(k+1))^(-rho))
       -sum(dexp(z, rate = 1/scale, log = TRUE))
     }
     for(i in seq_along(k)){
       start <- c(0.25,1,-0.1)
       # if(i > 1){
       #   start <- opt$pars
       #   if(start[1] < 1e-4){
       #     start[1] <- 1e-1
       #   }
       #   if(start[3] < -0.98 | start[3] > -1e-2){
       #     start[3] <- -0.1
       #   }
       # }
       ks <- k[i]
       Z <- 1:ks * (logdata[1:ks] - logdata[2:(ks+1)])
       opt <- Rsolnp::solnp(pars = start,
                            fun = expreg,
                            LB = c(0, 0, -1),
                            UB = c(2, 10, 0),
                            z = Z,
                            control = list(trace = 0))
       b[i] <- opt$par[2]/opt$par[1]
       shape[i] <- opt$par[1]
       rho[i] <- opt$par[3]
     }
   }
   data.frame(k = k, shape = shape, rho = rho, b = b)
}
