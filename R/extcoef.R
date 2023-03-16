#' Estimators of the extremal coefficient
#'
#' These functions estimate the extremal coefficient using an approximate sample
#' from the Frechet distribution.
#'
#' The \bold{Smith} estimator: suppose \eqn{Z(x)} is simple max-stable vector
#' (i.e., with unit Frechet marginals). Then
#'  \eqn{1/Z} is unit exponential and \eqn{1/\max(Z(s_1), Z(s_2))} is exponential
#'  with rate \eqn{\theta = \max\{Z(s_1), Z(s_2)\}}.
#'  The extremal index for the pair can therefore be calculated using the reciprocal mean.
#'
#' The \bold{Schlather and Tawn} estimator: the likelihood of the naive estimator for a pair
#' of two sites \eqn{A} is
#' \deqn{ \mathrm{card}\left\{ j: \max_{i \in A} X_i^{(j)}\bar{X}_i)>z \right\}
#' \log(\theta_A) - \theta_A \sum_{j=1}^n \left[ \max \left\{z, \max_{i \in A}
#' (X_i^{(j)}\bar{X}_i)\right\}\right]^{-1},}
#' where \eqn{\bar{X}_i = n^{-1} \sum_{j=1}^n 1/X_i^{(j)}} is the harmonic mean and \eqn{z}
#' is a threshold on the unit Frechet scale.
#' The search for the maximum likelihood estimate for every pair \eqn{A}
#'  is restricted to the interval \eqn{[1,3]}. A binned version of the extremal coefficient cloud is also returned.
#' The Schlather estimator is not self-consistent. The Schlather and Tawn estimator includes as special case
#' the Smith estimator if we do not censor the data (\code{p = 0}) and do not standardize observations by their harmonic mean.
#'
#'
#' The \bold{F-madogram} estimator is a non-parametric estimate based on a stationary process
#' \eqn{Z}; the extremal coefficient satisfies
#' \deqn{\theta(h)=\frac{1+2\nu(h)}{1-2\nu(h)},}
#' where
#' \deqn{\nu(h) = \frac{1}{2} \mathsf{E}[|F(Z(s+h)-F(Z(s))|]}
#' The implementation only uses complete pairs to calculate the relative ranks.
#'
#' All estimators are coded in plain R and computations are not optimized. The estimation
#' time can therefore be significant for large data sets. If there are no missing observations,
#' the routine \code{fmadogram} from the \code{SpatialExtremes} package should be prefered as it is
#' noticeably faster.
#'
#' The data will typically consist of max-stable vectors or block maxima.
#' Both of the Smith and the Schlather--Tawn estimators require unit Frechet margins; the margins will be standardized
#' to the unit Frechet scale, either parametrically or nonparametrically unless \code{standardize = FALSE}.
#' If \code{method = "parametric"}, a parametric GEV model is fitted to each column of \code{dat} using maximum likelihood
#'  estimation and transformed back using the probability integral transform. If \code{method = "nonparametric"},
#'  using the empirical distribution function. The latter is the default, as it is appreciably faster.
#' @export
#' @importFrom grDevices rgb
#' @references Schlather, M. and J. Tawn (2003). A dependence measure for multivariate and spatial extremes, \emph{Biometrika}, \bold{90}(1), pp.139--156.
#' @references Cooley, D., P. Naveau and P. Poncet (2006). Variograms for spatial max-stable random fields,  In: Bertail P., Soulier P., Doukhan P. (eds) \emph{Dependence in Probability and Statistics}. Lecture Notes in Statistics, vol. 187. Springer, New York, NY
#' @references R. J. Erhardt, R. L. Smith (2012), Approximate Bayesian computing for spatial extremes, \emph{Computational Statistics and Data Analysis}, \bold{56}, pp.1468--1481.
#' @param dat an \code{n} by \code{D} matrix of unit Frechet observations
#' @param coord an optional \code{d} by \code{D} matrix of location coordinates
#' @param estimator string indicating which model estimates to compute, one of \code{smith}, \code{schlather} or \code{fmado}.
#' @param method string indicating which method to use to transform the margins. See \bold{Details}
#' @param standardize logical; should observations be transformed to unit Frechet scale? Default is to transform
#' @param prob probability of not exceeding threshold \code{thresh}
#' @param thresh threshold parameter (default is to keep all data if \code{prob = 0}).
#' @param plot logical; should cloud of pairwise empirical estimates be plotted? Default to \code{TRUE}.
#' @param ... additional parameters passed to the function, currently ignored.
#' @return an invisible list with vectors \code{dist} if \code{coord} is non-null or else a matrix of pairwise indices \code{ind},
#'  \code{extcoef} and the supplied \code{estimator}, \code{fmado} and \code{binned}. If \code{estimator == "schlather"}, an additional matrix with 2 columns containing the binned distance \code{binned} with the \code{h} and the binned extremal coefficient.
#' @examples
#' \dontrun{
#' coord <- 10*cbind(runif(50), runif(50))
#' di <- as.matrix(dist(coord))
#' dat <- rmev(n = 1000, d = 100, param = 3, sigma = exp(-di/2), model = 'xstud')
#' res <- extcoef(dat = dat, coord = coord)
#' # Extremal Student extremal coefficient function
#'
#' XT.extcoeffun <- function(h, nu, corrfun, ...){
#'   if(!is.function(corrfun)){
#'     stop('Invalid function \"corrfun\".')
#'   }
#'   h <- unique(as.vector(h))
#'   rhoh <- sapply(h, corrfun, ...)
#'   cbind(h = h, extcoef = 2*pt(sqrt((nu+1)*(1-rhoh)/(1+rhoh)), nu+1))
#' }
#' #This time, only one graph with theoretical extremal coef
#' plot(res$dist, res$extcoef, ylim = c(1,2), pch = 20); abline(v = 2, col = 'gray')
#' extcoefxt <- XT.extcoeffun(seq(0, 10, by = 0.1), nu = 3,
#'                             corrfun = function(x){exp(-x/2)})
#' lines(extcoefxt[,'h'], extcoefxt[,'extcoef'], type = 'l', col = 'blue', lwd = 2)
#' # Brown--Resnick extremal coefficient function
#' BR.extcoeffun <- function(h, vario, ...){
#'   if(!is.function(vario)){
#'     stop('Invalid function \"vario\".')
#'   }
#'   h <- unique(as.vector(h))
#'   gammah <- sapply(h, vario, ...)
#'   cbind(h = h, extcoef = 2*pnorm(sqrt(gammah/4)))
#' }
#' extcoefbr<- BR.extcoeffun(seq(0, 20, by = 0.25), vario = function(x){2*x^0.7})
#' lines(extcoefbr[,'h'], extcoefbr[,'extcoef'], type = 'l', col = 'orange', lwd = 2)
#'
#' coord <- 10*cbind(runif(20), runif(20))
#' di <- as.matrix(dist(coord))
#' dat <- rmev(n = 1000, d = 20, param = 3, sigma = exp(-di/2), model = 'xstud')
#' res <- extcoef(dat = dat, coord = coord, estimator = "smith")
#' }
extcoef <- function(dat,
                    coord = NULL,
                    thresh = NULL,
                    estimator = c("schlather", "smith", "fmado"),
                    standardize = TRUE,
                    method = c("nonparametric", "parametric"),
                    prob = 0,
                    plot = TRUE,
                    ...) {
  if(is.null(coord) && !is.null(list(...)$loc)){
   coord <-  list(...)$loc
  }
    stopifnot(is.matrix(dat),
              ncol(dat) >= 2,
              nrow(coord) == ncol(dat))
    estimator <- match.arg(estimator)
    fr <- dat
    #Transform margins to unit Frechet
    if(standardize && estimator != "fmado"){
      method <- match.arg(method)
      if(method == "nonparametric"){
        for(j in 1:ncol(dat)){
          nj <- sum(!is.na(dat[, j])) + 1L
          fr[,j] <- -1/log(rank(dat[,j],
                                na.last = "keep",
                                ties.method = "random")/nj)
        }
      } else{ #parametric
        for(j in 1:ncol(dat)){
          fr[,j] <- -1/log(
            do.call(what = mev::pgev,
                    c(list(q = dat[,j]),
                      fit.gev(dat[,j])$estimate)))
        }
      }
    }
    if(estimator == "schlather"){
      if (is.null(thresh)) {
        thresh <- -1/log(prob)  #threshold on unit Frechet scale
      }
    # Negative log-likelihood function
      nll <- function(theta, X, thresh) {
        -sum(X > thresh) * log(theta) + theta * sum(1/pmax(thresh, X))
      }
    #Compute harmonic mean, reweight estimators
    harmo_mean <- apply(1/fr, 2, mean, na.rm = TRUE)
    transfo_fr <- sweep(x = fr,
                        MARGIN = 2,
                        STATS = harmo_mean,
                        FUN = "/")
    }
    N <- ncol(dat) * (ncol(dat) - 1)/2
    theta_vals <- rep(0, N)
    if(!is.null(coord)){
      dist_vals <- rep(0, N)
    } else{
     ind_vals <- matrix(0L, nrow = N, ncol = 2)
    }
    if(estimator == "schlather"){
     overlap_size <-  rep(0, N)
    }
    k <- 0L
    for (i in 1L:(ncol(dat) - 1L)) {
        for (j in (i + 1L):ncol(dat)) {
          k <- k + 1L
              if(!is.null(coord)){
                 dist_vals[k] <- dist(coord[c(i, j), ])
                } else{
                  ind_vals[k,] <- c(i, j)
                }
                if (estimator == "schlather") {
                  X <- as.vector(apply((na.omit(transfo_fr[, c(i, j)])), 1, max))
                  if (length(X) > 1) {
                    theta_vals[k] <- sum(X > thresh)/sum(1/pmax(thresh, X))
                  } else{
                  theta_vals[k] <- NA
                  }
                  overlap_size[k] <- sqrt(length(X))  #sqrt of size of overlapping set
                } else if (estimator == "fmado") {
                  # F-madogram estimator
                  Y <- apply(na.omit(dat[, c(i, j)]), 2,
                             rank, ties.method = "random")
                  if(is.null(dim(Y))){
                    next()
                  } else{
                  nu <- sum(abs(Y[, 1] - Y[, 2]))/(2 * nrow(Y)^2)
                  theta_vals[k] <- (1 + 2 * nu)/(1 - 2 * nu)
                  }
                } else if(estimator == "smith"){
                  theta_vals[k] <- 1/mean(1/apply(na.omit(fr[,c(i, j)]), 1, max))
                }
            }
       }
    reslist <- list()
    # Reorder theta observations by distance (sorted)
    if(!is.null(coord)){
      theta_vals <- theta_vals[order(dist_vals)]
      if(estimator == "schlather"){
        if(any(theta_vals > 3)){
         theta_vals[theta_vals > 3] <- 3
        } else if(any(theta_vals < 1)){
          theta_vals[theta_vals < 1] <- 1
        }
        overlap_size <- overlap_size[order(dist_vals)]
      }
      dist_vals <- sort(dist_vals)
    if (estimator == "schlather") {
        # na.omit will strip columns with NA if either method not calculated
        d_max <- 0.75 * dist_vals[length(dist_vals)]
        # Keep fraction only of data and create bins
        del <- d_max/sqrt(N)
        h <- seq(del/2, d_max, by = del)
        # Bin thetas by interval, with halfwidth eps = del/2
        bin <- findInterval(dist_vals, c(0, h + del/2))
        ubin <- unique(bin)
        binned_theta <- as.vector(by(data = cbind(theta_vals, overlap_size),
                                     INDICES = bin, FUN = function(xmat) {
            weighted.mean(x = xmat[, 1], w = xmat[, 2])
        }, simplify = TRUE))
        reslist$binned <- na.omit(cbind(h = h[ubin], binned.extcoef = binned_theta))

    }
    reslist$dist = dist_vals
    class(reslist) <- "mev_extcoef"
    } else{
      reslist$ind <- ind_vals
    }
    reslist$estimator <- estimator
    reslist$extcoef <- theta_vals
    if(plot && !is.null(coord)){
     plot(reslist)
    }
    return(invisible(reslist))
}

#' @export
plot.mev_extcoef <- function(x, ...){
    ellipsis <- list(...)
    nellips <- names(ellipsis)
    if("tikz" %in% nellips){
      tikz <- ellipsis$tikz
      stopifnot(is.logical(tikz))
    } else{
      tikz <- FALSE
    }
    ymax <- 1.1*min(c(max(c(2.1,x$extcoef),
                        na.rm = TRUE),
                    3), na.rm = TRUE)
    xmax <- 1.02*max(x$dist, na.rm = TRUE)
    if (x$estimator == "schlather"){
        graphics::plot(x = x$dist,
                       y = x$extcoef,
                       xlab = ifelse(tikz, "$h$", "h"),
                       ylab = ifelse(tikz, "$\\theta$", expression(theta(h))),
                       bty = "l",
                       cex = 0.8,
                       col = grDevices::rgb(0, 0, 0, alpha = 0.5),
                       ylim = c(1, ymax),
                       yaxs = "i",
                       xaxs = "i",
                       xlim = c(0, xmax))
          lines(x$binned[, 1], x$binned[, 2], col = 2, lwd = 2)
    } else{
      graphics::plot(x = x$dist,
                     y = x$extcoef,
                     xlab = ifelse(tikz, "$h$", "h"),
                     ylab = ifelse(tikz, "$\\theta$", expression(theta(h))),
                     bty = "l",
                     cex = 0.8,
                     col = grDevices::rgb(0, 0, 0, alpha = 0.5),
                     ylim = c(1, ymax),
                     yaxs = "i", xaxs = "i",
                     xlim = c(0, xmax))

    }
    abline(h = 2, col = "grey")
    title(main = "Extremal coefficient", sub = switch(x$estimator,
                       schlather = "Schlather and Tawn estimator",
                       fmado = "F-madogram estimator",
                       smith = "Smith estimator"))

}
