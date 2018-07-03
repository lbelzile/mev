#' Semi-parametric marginal transformation to uniform
#' 
#' The function \code{spunif} transforms a matrix or vector of data \code{x}
#' to the pseudo-uniform scale using a semiparametric transform. Data below the threshold
#' are transformed to pseudo-uniforms using a rank transform, while data above the threshold
#' are assumed to follow a generalized Pareto distribution. The parameters of the latter are
#' estimated using maximum likelihood if either \code{scale = NULL} or \code{shape = NULL}.
#' 
#' @param x matrix or vector of data
#' @param thresh vector of marginal thresholds
#' @param scale vector of marginal scale parameters for the generalized Pareto
#' @param shape vector of marginal shape parameters for the generalized Pareto
#' @return a matrix or vector of the same dimension as \code{x}, with pseudo-uniform observations
#' @export
#' @author Leo Belzile
#' @examples 
#' x <- rmev(1000, d = 3, param = 2, model = 'log')
#' thresh <- apply(x, 2, quantile, 0.95)
#' spunif(x, thresh)
spunif <- function(x, thresh, scale = NULL, shape = NULL) {
    # Routine for marginal transformation
    spunif_univ <- function(x, thresh, scale = NULL, shape = NULL) {
        if (any(is.null(scale), is.null(shape))) {
            # Routine to fit GP
            gp <- gp.fit(xdat = x, threshold = thresh, method = "Grimshaw")$est
            scale <- gp["scale"]
            shape <- gp["shape"]
        } else {
            if (scale < 0 || shape < -1) {
                stop("Invalid scale or shape parameter provided.")
            }
        }
        n <- length(x)
        # copy x
        un <- x
        below <- (x < thresh)
        zeta <- 1 - sum(below)/n  #proportion above
        if (zeta == 0) {
            stop("Threshold is too high, no data above")
        }
        un[below] <- rank(x, ties.method = "random", na.last = "keep")[below]/(n + 1)
        un[!below] <- 1 - zeta * (1 + shape * (x[!below] - thresh)/scale)^(-1/shape)
        return(un)
    }
    
    # Vector input
    if (is.vector(x)) {
        if (length(thresh) > 1) {
            stop("Invalid threshold in `spunif`")
        }
        spunif_univ(x = x, thresh = thresh, scale = scale, shape = shape)
    } else {
        # Matrix input
        if (length(thresh) != ncol(x)) {
            stop("Invalid input: `thresh` must be of the same length as `ncol(x)`")
        }
        if (length(scale) != length(shape)) {
            stop("Invalid `scale` or `shape` input in `spunif`")
        }
        sapply(1:ncol(x), function(j) {
            spunif_univ(x = x[, j], thresh = thresh[j], scale = scale[j], shape = shape[j])
        })
    }
}


