#' Schlather's estimator of the extremal coefficient
#'
#' This functions returns a matrix with the Schlather and Tawn estimator of extremal coefficient
#' based on the extremal index of pairs of max-stable vectors. The likelihood of the naive estimator for a pair of two sites \eqn{A} is
#' \deqn{ \mathrm{card}\left\{ j: \max_{i \in A} X_i^{(j)}\bar{X}_i)>z \right\} \log(\theta_A) - \theta_A \sum_{j=1}^n \left[ \max \left\{z, \max_{i \in A} (X_i^{(j)}\bar{X}_i)\right\}\right]^{-1},}
#'where \eqn{\bar{X}_i = n^{-1} \sum_{j=1}^n 1/X_i^{(j)}} is the harmonic mean and \eqn{z} is a threshold on the unit Frechet scale.
#' The search for the maximum likelihood estimate for every pair \eqn{A} is restricted to the interval \eqn{[1,2]}.
#' The Schlather estimator is not self-consistent. A binned version of the extremal coefficient cloud is also reported and superimposed on the graph.
#'
#' The F-madogram estimator is
#' \deqn{\nu(h) = \frac{1}{2} \mathsf{E} |F(Z(x+h))-F(Z(x))|.}
#' We can write the extremal coefficient as
#' \deqn{\theta(h) = \frac{1+2\nu(h)}{1-2\nu(h)},}
#' leading to a non-parametric pairwise estimator. The empirical cloud of estimates is obtained by using a rank-based
#' approximation to the distribution function \eqn{F}.
#'
#' The implementation only uses complete pairs to calculate the relative ranks. Both estimators are coded in plain R and
#' the estimation time can be significant.
#'
#' @references Schlather, M. and J. Tawn (2003). A dependence measure for multivariate and spatial extremes, \emph{Biometrika}, \bold{90}(1), pp.139--156.
#' @references Cooley, D., P. Naveau and P. Poncet (2006). Variograms for spatial max-stable random fields,  In: Bertail P., Soulier P., Doukhan P. (eds) \emph{Dependence in Probability and Statistics}. Lecture Notes in Statistics, vol. 187. Springer, New York, NY
#' @param dat a matrix of unit Fréchet observations
#' @param thresh threshold parameter
#' @param prob probability of not exceeding threshold \code{thresh}
#' @param estimator string indicating which model estimates to compute, one of \code{schlather} or \code{fmado}.
#' @param which.plot string indicating which model estimates to plot, one of \code{schlather} or \code{fmado}. Can be abbreviated.
#' @param tikz logical indicating whether output be formatted to use with \code{tikzDevice} package? Default to \code{FALSE}.
#'
#' @examples
#' loc <- 10*cbind(runif(50), runif(50))
#' di <- as.matrix(dist(loc))
#' dat <- mev::rmev(n = 1000, d = 100, param = 3, sigma = exp(-di/2), model = "xstud")
#' extcoefschlather(dat = dat, loc = loc)
extcoefschlather <- function(dat, loc, thresh, estimator = c("schlather","fmado"),
                             prob = 0, which.plot = c("schlather","fmado"), tikz = FALSE){

  estimator <- match.arg(estimator, c("schlather", "fmado"), several.ok = TRUE)
  which.plot <- match.arg(which.plot, c("schlather", "fmado"), several.ok = TRUE)
  isSchlather <- "schlather" %in% estimator
  isFmado <- "fmado" %in% estimator
  plotind <- c()
  if("schlather" %in% which.plot && isSchlather){plotind <- c(plotind, 1)}
  if("fmado" %in% which.plot && isFmado){plotind <- c(plotind, 2)}
  if(missing(thresh)){
    thresh <- -1/log(prob) #threshold on unit Frechet scale
  }
   fr <- -1/log(apply(dat, 2, rank, ties.method = "random",
                          na.last = "keep")/(nrow(dat) + 1))
   #Negative log-likelihood function
   nll <- function(theta, X, thresh){
    - sum(X > thresh) * log(theta) + theta * sum(1/pmax(thresh, X))
  }

  harmo_mean <- apply(1/fr, 2, mean, na.rm = TRUE)
  transfo_fr <- t(t(fr) * harmo_mean)
  theta_mat <- matrix(NA, nrow = ncol(dat)*(ncol(dat)-1)/2, ncol = 4)
  k <- 0L
  for(i in 1:(ncol(dat)-1)){
    for(j in (i+1):ncol(dat)){
      X <- as.vector(apply((na.omit(transfo_fr[,c(i,j)])),1,max))
      if(length(X) > 1){
        k <- k + 1L
        theta_mat[k, 1] <- dist(loc[c(i,j),])
        if(isSchlather){
        theta_mat[k, 2] <- optim(par = 1.5, fn = nll, method = "Brent",
                                 X = X, lower = 1, upper = 3, thresh = thresh)$par
        theta_mat[k, 3] <- sqrt(length(X)) #sqrt of size of overlapping set
        }
        if(isFmado){
        #F-madogram estimator
        Y <- apply(na.omit(dat[,c(i,j)]), 2, rank, ties.method = "random")
        nu <- sum(abs(Y[,1]-Y[,2]))/(2*nrow(Y)^2)
        theta_mat[k, 4] <- (1+2*nu)/(1-2*nu)
        }
      }
    }
  }
  #Reorder theta observations by distance (sorted)
  theta_mat <- theta_mat[order(theta_mat[,1]),]
  #na.omit will strip columns with NA if either method not calculated
  d_max <- 0.75*theta_mat[nrow(theta_mat),1]
  #Keep fraction only of data and create bins
  del <- d_max/sqrt(nrow(theta_mat))
  h <- seq(del/2, d_max, by = del)
  if(isSchlather){
  #Bin thetas by interval, with halfwidth eps = del/2
  bin <- findInterval(theta_mat[,1], c(0, h+del/2))
  ubin <- unique(bin)
  binned_theta <- as.vector(by(data = theta_mat, INDICES = bin,
                     FUN = function(xmat){
                       weighted.mean(x = xmat[,2], w = xmat[,3])
                     }, simplify = TRUE))
  }
  #Plot results
  if(all(c(1,2) %in% plotind)){
    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(1,2))
  }
  if(1 %in% plotind){
    graphics::plot(theta_mat[,1], theta_mat[,2], xlab = ifelse(tikz, "$h$", "h"),
                   ylab = ifelse(tikz, "$\\theta$", expression(theta(h))),
                   bty = "l", cex = 0.4, col = scales::alpha("black", alpha = 0.5),
                   ylim=c(1,2.1), yaxs = "i", xaxs = "i", main = "Schlather")
    abline(h=2, col="grey")
    lines(h[ubin], binned_theta, col = 2, lwd = 2)
  }
  if(2 %in% plotind){
    graphics::plot(theta_mat[,1], theta_mat[,4], xlab = ifelse(tikz, "$h$", "h"),
                   ylab = ifelse(tikz, "$\\theta$", expression(theta(h))),
                   bty = "l", cex = 0.4, col = scales::alpha("black", alpha = 0.5),
                   ylim=c(1,2.1), yaxs = "i", xaxs = "i", main = "F-madogram")
    abline(h=2, col="grey")
  }
  if(all(c(1,2) %in% plotind)){
   par(old.par)
  }
  colnames(theta_mat[,1:2]) <- c("dist", "ext.coeff")
  #Return invisible list with coefs and binned
  reslist <- list(dist = theta_mat[,1])
  if(isSchlather){
    reslist$schlather <- theta_mat[,2]
    reslist$binned <- cbind(h = h[ubin], binned.ext.coef = binned_theta)
  }
  if(isFmado){
    reslist$fmado <- theta_mat[,4]
  }
  return(invisible(reslist))
}


#' Extremal coefficient function for the extremal Student process
#' @param h vector of pairwise distances
#' @param corrfun correlation function
#' @param ... additional parameters passed to \code{corrfun}
#' @examples
#' extcoef <- Xstud.extcoeffun(seq(0, 10, by = 0.1), nu = 3,
#'     corrfun = function(x){exp(-x/2)})
#' plot(extcoef, type = 'l', ylim = c(1,2))
Xstud.extcoeffun <- function(h, nu, corrfun, ...){
  if(!is.function(corrfun)){
    stop("Invalid function `corrfun`.")
  }
  h <- unique(as.vector(h))
  rhoh <- sapply(h, corrfun, ...)
  cbind(h = h, extcoef = 2*pt(sqrt((nu+1)*(1-rhoh)/(1+rhoh)), nu+1))
}

#' Extremal coefficient function for the Brown-Resnick process
#' @param h pairwise distance matrix
#' @param vario semivariogram function
#' @param ... additional parameters passed to \code{vario}
#' @examples
#' extcoef <- BR.extcoeffun(seq(0, 20, by = 0.25), vario = function(x){2*x^0.7})
#' plot(extcoef, type = 'l', ylim = c(1,2))
BR.extcoeffun <- function(h, vario, ...){
  if(!is.function(vario)){
    stop("Invalid function `vario`.")
  }
  h <- unique(as.vector(h))
  gammah <- sapply(h, vario, ...)
  cbind(h = h, extcoef = 2*pnorm(sqrt(gammah/4)))
}
