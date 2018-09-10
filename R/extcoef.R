#' Estimators of the extremal coefficient
#'
#' These functions estimate the extremal coefficient using an approxiate sample from the Frechet distribution.
#'
#' @details This function returns a matrix with the Schlather and Tawn estimator of extremal coefficient
#' based on the extremal index of pairs of max-stable vectors. The likelihood of the naive estimator for a pair of two sites \eqn{A} is
#' \deqn{ \mathrm{card}\left\{ j: \max_{i \in A} X_i^{(j)}\bar{X}_i)>z \right\} \log(\theta_A) - \theta_A \sum_{j=1}^n \left[ \max \left\{z, \max_{i \in A} (X_i^{(j)}\bar{X}_i)\right\}\right]^{-1},}
#'where \eqn{\bar{X}_i = n^{-1} \sum_{j=1}^n 1/X_i^{(j)}} is the harmonic mean and \eqn{z} is a threshold on the unit Frechet scale.
#' The search for the maximum likelihood estimate for every pair \eqn{A} is restricted to the interval \eqn{[1,2]}.
#' The Schlather estimator is not self-consistent. A binned version of the extremal coefficient cloud is also reported and superimposed on the graph.
#'
#' Additionally, a poor man's version of the F-madogram estimator is returned alongside with the Schlather and Tawn estimates.
#' If \code{estimator = c('schlather', 'fmado')} and \code{plot = TRUE}, both data clouds are plotted alongside one another.
#'
#' The F-madogram estimator is
#' \deqn{\nu(h) = \frac{1}{2} \mathsf{E} |F(Z(x+h))-F(Z(x))|.}
#' We can write the extremal coefficient as
#' \deqn{\theta(h) = \frac{1+2\nu(h)}{1-2\nu(h)},}
#' leading to a non-parametric pairwise estimator. The empirical cloud of estimates is obtained by using a rank-based
#' approximation to the distribution function \eqn{F}.
#'
#' The implementation only uses complete pairs to calculate the relative ranks. Both estimators are coded in plain R and computations are not optimized. The estimation time can therefore be significant for large data sets.
#' @importFrom grDevices rgb
#' @references Schlather, M. and J. Tawn (2003). A dependence measure for multivariate and spatial extremes, \emph{Biometrika}, \bold{90}(1), pp.139--156.
#' @references Cooley, D., P. Naveau and P. Poncet (2006). Variograms for spatial max-stable random fields,  In: Bertail P., Soulier P., Doukhan P. (eds) \emph{Dependence in Probability and Statistics}. Lecture Notes in Statistics, vol. 187. Springer, New York, NY
#' @param dat an \code{n} by \code{D} matrix of unit Frechet observations
#' @param thresh threshold parameter (dafault is to keep all data if \code{prob = 0}.
#' @param loc \code{d} by \code{D} matrix of location coordinates
#' @param prob probability of not exceeding threshold \code{thresh}
#' @param estimator string indicating which model estimates to compute, one of \code{schlather} or \code{fmado}.
#' @param which.plot string indicating which model estimates to plot, one of \code{schlather} or \code{fmado}. Can be abbreviated.
#' @param tikz logical indicating whether output be formatted to use with \code{tikzDevice} package? Default to \code{FALSE}.
#' @return an invisible list with components \code{dist}, \code{schlather}, \code{fmado} and \code{binned}. The first three are \code{n} vectors, while \code{binned} is a matrix with 2 columns containing the binned distance \code{h} and the binned extremal coefficient.
#' @examples
#' \dontrun{
#' loc <- 10*cbind(runif(50), runif(50))
#' di <- as.matrix(dist(loc))
#' dat <- rmev(n = 1000, d = 100, param = 3, sigma = exp(-di/2), model = 'xstud')
#' res <- extcoef(dat = dat, loc = loc)
#' Extremal Student extremal coefficient function
#'
#' XT.extcoeffun <- function(h, nu, corrfun, ...){
#'   if(!is.function(corrfun)){
#'     stop('Invalid function `corrfun`.')
#'   }
#'   h <- unique(as.vector(h))
#'   rhoh <- sapply(h, corrfun, ...)
#'   cbind(h = h, extcoef = 2*pt(sqrt((nu+1)*(1-rhoh)/(1+rhoh)), nu+1))
#' }
#' #This time, only one graph with theoretical extremal coef
#' plot(res$dist, res$schlather, ylim = c(1,2), pch = 20); abline(v = 2, col = 'gray')
#' extcoefxt <- XT.extcoeffun(seq(0, 10, by = 0.1), nu = 3,
#'                             corrfun = function(x){exp(-x/2)})
#' lines(extcoefxt[,'h'], extcoefxt[,'extcoef'], type = 'l', col = 'blue', lwd = 2)
#' # Brown--Resnick extremal coefficient function
#' BR.extcoeffun <- function(h, vario, ...){
#'   if(!is.function(vario)){
#'     stop('Invalid function `vario`.')
#'   }
#'   h <- unique(as.vector(h))
#'   gammah <- sapply(h, vario, ...)
#'   cbind(h = h, extcoef = 2*pnorm(sqrt(gammah/4)))
#' }
#' extcoefbr<- BR.extcoeffun(seq(0, 20, by = 0.25), vario = function(x){2*x^0.7})
#' lines(extcoefbr[,'h'], extcoefbr[,'extcoef'], type = 'l', col = 'orange', lwd = 2)
#' }
extcoef <- function(dat, loc, thresh, estimator = c("schlather", "fmado"), prob = 0,
                    which.plot = c("schlather", "fmado"), tikz = FALSE) {

    estimator <- match.arg(estimator, c("schlather", "fmado"), several.ok = TRUE)
    which.plot <- match.arg(which.plot, c("schlather", "fmado"), several.ok = TRUE)
    isSchlather <- "schlather" %in% estimator
    isFmado <- "fmado" %in% estimator
    plotind <- c()
    if ("schlather" %in% which.plot && isSchlather) {
        plotind <- c(plotind, 1)
    }
    if ("fmado" %in% which.plot && isFmado) {
        plotind <- c(plotind, 2)
    }
    if (missing(thresh)) {
        thresh <- -1/log(prob)  #threshold on unit Frechet scale
    }
    fr <- -1/log(apply(dat, 2, rank, ties.method = "random", na.last = "keep")/(nrow(dat) + 1))
    # Negative log-likelihood function
    nll <- function(theta, X, thresh) {
        -sum(X > thresh) * log(theta) + theta * sum(1/pmax(thresh, X))
    }

    harmo_mean <- apply(1/fr, 2, mean, na.rm = TRUE)
    transfo_fr <- t(t(fr) * harmo_mean)
    theta_mat <- matrix(NA, nrow = ncol(dat) * (ncol(dat) - 1)/2, ncol = 4)
    k <- 0L
    for (i in 1:(ncol(dat) - 1)) {
        for (j in (i + 1):ncol(dat)) {
            X <- as.vector(apply((na.omit(transfo_fr[, c(i, j)])), 1, max))
            if (length(X) > 1) {
                k <- k + 1L
                theta_mat[k, 1] <- dist(loc[c(i, j), ])
                if (isSchlather) {
                  theta_mat[k, 2] <- optim(par = 1.5, fn = nll, method = "Brent", X = X, lower = 1, upper = 3, thresh = thresh)$par
                  theta_mat[k, 3] <- sqrt(length(X))  #sqrt of size of overlapping set
                }
                if (isFmado) {
                  # F-madogram estimator
                  Y <- apply(na.omit(dat[, c(i, j)]), 2, rank, ties.method = "random")
                  nu <- sum(abs(Y[, 1] - Y[, 2]))/(2 * nrow(Y)^2)
                  theta_mat[k, 4] <- (1 + 2 * nu)/(1 - 2 * nu)
                }
            }
        }
    }
    # Reorder theta observations by distance (sorted)
    theta_mat <- theta_mat[order(theta_mat[, 1]), ]
    # na.omit will strip columns with NA if either method not calculated
    d_max <- 0.75 * theta_mat[nrow(theta_mat), 1]
    # Keep fraction only of data and create bins
    del <- d_max/sqrt(nrow(theta_mat))
    h <- seq(del/2, d_max, by = del)
    if (isSchlather) {
        # Bin thetas by interval, with halfwidth eps = del/2
        bin <- findInterval(theta_mat[, 1], c(0, h + del/2))
        ubin <- unique(bin)
        binned_theta <- as.vector(by(data = theta_mat, INDICES = bin, FUN = function(xmat) {
            weighted.mean(x = xmat[, 2], w = xmat[, 3])
        }, simplify = TRUE))
    }
    # Plot results
    if (all(c(1, 2) %in% plotind)) {
        old.par <- par(no.readonly = TRUE)
        par(mfrow = c(1, 2))
    }
    if (1 %in% plotind) {
        graphics::plot(theta_mat[, 1], theta_mat[, 2], xlab = ifelse(tikz, "$h$", "h"), ylab = ifelse(tikz, "$\\theta$", expression(theta(h))),
            bty = "l", cex = 0.4, col = grDevices::rgb(0, 0, 0, alpha = 0.5), ylim = c(1, 2.1), yaxs = "i", xaxs = "i", main = "Schlather")
        abline(h = 2, col = "grey")
        lines(h[ubin], binned_theta, col = 2, lwd = 2)
    }
    if (2 %in% plotind) {
        graphics::plot(theta_mat[, 1], theta_mat[, 4], xlab = ifelse(tikz, "$h$", "h"), ylab = ifelse(tikz, "$\\theta$", expression(theta(h))),
            bty = "l", cex = 0.4, col = grDevices::rgb(0, 0, 0, alpha = 0.5), ylim = c(1, 2.1), yaxs = "i", xaxs = "i", main = "F-madogram")
        abline(h = 2, col = "grey")
    }
    if (all(c(1, 2) %in% plotind)) {
        par(old.par)
    }
    colnames(theta_mat[, 1:2]) <- c("dist", "ext.coeff")
    # Return invisible list with coefs and binned
    reslist <- list(dist = theta_mat[, 1])
    if (isSchlather) {
        reslist$schlather <- theta_mat[, 2]
        reslist$binned <- cbind(h = h[ubin], binned.ext.coef = binned_theta)
    }
    if (isFmado) {
        reslist$fmado <- theta_mat[, 4]
    }
    return(invisible(reslist))
}


#' @details Consider the pairwise extremal coefficient estimate of Smith. Suppose \eqn{Z(x)} is  a max-stable vector with unit Frechet distribution. Then
#'  \eqn{1/Z} is unit exponential and \eqn{1/\max(Z(s_1), Z(s_2))} is exponential
#'  with rate \eqn{\theta = \max(Z(s_1), Z(s_2))}.
#'  The extremal index for the pair can therefore be calculated using the reciprocal mean.
#'  If \code{method = "parametric"}, a parametric GEV model is fitted to each column of \code{dat} using maximum likelihood
#'  estimation and transformed back using the probability integral transform. If \code{method = "nonparametric"},
#'  using the empirical distribution function. The latter is the default, as it is appreciatly faster.
#' @param dat an \code{n} by \code{d} matrix of maximas for \code{d} variables
#' @param loc \code{d} by \code{D} matrix of coordinates; default to \code{NULL}
#' @param method string indicating which method to use to transform the margins. See \bold{Details}
#' @param standardize logical; should observations be transformed to unit Frechet scale? Default is to transform
#' @references R. J. Erhardt, R. L. Smith (2012), Approximate Bayesian computing for spatial extremes, \emph{Computational Statistics and Data Analysis}, \bold{56}, pp.1468--1481.
#' @rdname extcoef
#' @export
#' @examples
#' loc <- 10*cbind(runif(20), runif(20))
#' di <- as.matrix(dist(loc))
#' dat <- rmev(n = 1000, d = 20, param = 3, sigma = exp(-di/2), model = 'xstud')
#' res <- extcoef.smith(dat = dat, loc = loc)
extcoef.smith <- function(dat, loc = NULL, standardize = TRUE, method = c("nonparametric", "parametric")){
  stopifnot(is.matrix(dat), ncol(dat) >= 2)
  #Preprocessing, transform observations to unit Frechet
  if(standardize){
    method <- match.arg(method[1], choices = c("parametric", "nonparametric"))
    fre_dat <- matrix(0, ncol = ncol(dat), nrow = nrow(dat))
    if(method == "nonparametric"){
      for(j in 1:ncol(dat)){
        nj <- sum(!is.na(dat[, 1])) + 1L
        fre_dat[,j] <- -1/log(rank(dat[,j], na.last = "keep", ties.method = "random")/nj)
      }
    } else{ #parametric
      for(j in 1:ncol(dat)){
        fre_dat[,j] <- -1/log(do.call(what = evd::pgev, c(list(q = dat[,j]), fit.gev(dat[,j])$estimate)))
      }
    }
  }
  #Two methods: multivariate with all pairs or else as function of Euclidean distance
    if(!is.null(loc)){
      extcoef <- matrix(0, nrow = ncol(dat)*(ncol(dat)-1)/2, ncol = 2)
      acc <- 0
      for(i in 1:(ncol(dat)-1)){
        for(j in (i+1):ncol(dat)){
          acc <- acc + 1L
          extcoef[acc,] <- c(dist(loc[c(i,j),]), 1/mean(1/apply(dat[,c(i, j)], 1, max)))
        }
      }
      colnames(extcoef) <- c("dist", "ext.coeff")
      return(extcoef)
    } else{
      extcoef <- matrix(0, nrow = ncol(dat)*(ncol(dat)-1)/2, ncol = 3)
      acc <- 0
    for(i in 1:(ncol(dat)-1)){
      for(j in (i+1):ncol(dat)){
        acc <- acc + 1L
        extcoef[acc,] <- c(1/mean(1/apply(dat[,c(i, j)], 1, max)), i, j)
      }
    }
      colnames(extcoef) <- c("dist", "index1", "index2")
      return(extcoef)
    }
}
