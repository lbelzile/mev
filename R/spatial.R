
#' Power variogram model
#'
#' The power variogram model is
#' \deqn{\gamma(h) = (\|h\|/\lambda)^\alpha, \quad \lambda>0, \alpha \in [0,2).}
#'
#' @param h vector or matrix of pairwise distances
#' @param scale scale parameter
#' @param alpha smoothness parameter
#' @return a vector or matrix of variogram values of the same length as \code{h}
#' @keywords internal
#' @export
power.vario <- function(h, alpha, scale = 1) {
  stopifnot(length(alpha) == 1, length(scale) == 1)
  if (scale < 0) {
    stop("Invalid scale parameter in \"power.vario\".")
  }
  if (alpha < 0 || alpha > 2 - 1e-10) {
    stop("Invalid shape parameter in \"power.vario\".")
  }
  (h / scale)^alpha
}


#' Power exponential correlation model
#'
#' The power correlation model is
#' \deqn{\rho(h) = \exp\{-(\|h\|/\lambda)^\alpha\}, \quad \lambda>0, \alpha \in [0,2).}
#'
#' @param h vector or matrix of pairwise distances
#' @param alpha smoothness parameter
#' @param scale scale parameter
#' @return a vector or matrix of correlations of the same dimensions as \code{h}
#' @export
#' @keywords internal
powerexp.cor <- function(h, alpha = 1, scale = 1) {
  stopifnot(length(alpha) == 1, length(scale) == 1)
  if (scale < 0) {
    stop("Invalid scale parameter in \"power.vario\"")
  }
  if (alpha < 0 || alpha > 2 - 1e-10) {
    stop("Invalid shape parameter in \"power.vario\"")
  }
  exp(-(h / scale)^alpha)
}

#' Variogram model of Schlather and Moreva
#'
#' The variogram model is
#' \deqn{\gamma(h) = \frac{[1+\{(\|h\|/\lambda\}^\alpha]^{\beta/\alpha}-1}{2^{\beta/\alpha}-1}, \quad 0 < \alpha \leq 2, \beta \leq 2.}
#' The model is defined at \eqn{\beta=0} by continuity.
#'
#' @param h vector or matrix of pairwise distances
#' @param alpha smoothness parameter
#' @param beta shape parameter, must be less than 2
#' @param scale scale parameter
#' @return a vector or matrix of variogram values of the same length as \code{h}
#' @export
#' @keywords internal
schlather.vario <- function(h, alpha, beta, scale = 1) {
  stopifnot(length(alpha) == 1, length(beta) == 1, length(scale) == 1)
  if (alpha <= 0 || alpha > 2 - 1e-10) {
    stop("Invalid shape parameter in \"schlather.vario\"")
  }
  if (beta > 2 - 1e-10) {
    stop("Invalid smoothness parameter in \"schlather.vario\"")
  }
  if (scale <= 0) {
    stop("Invalid scale parameter in \"schlather.vario\"")
  }
  if (!(abs(beta) < 1e-5)) {
    pow <- beta / alpha
    ((1 + (h / scale)^alpha)^pow - 1) / (exp(pow * log(2)) - 1)
  } else {
    log(1 + (h / scale)^alpha) / log(2)
  }
  # If alpha is too close to zero, may return matrix with NaN
}


#' Transform variogram matrix to covariance of conditional random field
#'
#' The matrix \code{Lambda} is half the semivariogram matrix. The function
#' returns the conditional covariance with respect to entries in \code{co},
#' restricted to the \code{subA} rows and the \code{subB} columns.
#' @param Lambda Negative definite matrix for the Huesler--Reiss model
#' @param co vector of integer with conditioning sites
#' @param subA vector of integers with sub-entries (not in \code{co}) for rows
#' @param subB vector of integers with sub-entries (not in \code{co}) for columns. If missing, default to \code{subA}.
#' @keywords internal
#' @export
Lambda2cov <- function(Lambda, co, subA, subB) {
  if (length(co) > 1) {
    stop("Conditioning site should have length 1")
  } # only if integer, could also use logical or negative subsetting
  if (is.logical(subA)) {
    if (length(subA) != ncol(Lambda)) {
      stop("Invalid length")
    }
    if (subA[co]) {
      stop("Covariance matrix should be shed with respect to conditioning site")
    }
  } else if (co %in% subA) { # is numeric
    stop("Invalid conditioning variable for covariance matrix")
  }
  if (missing(subB)) {
    subB <- subA
  } else {
    if (is.logical(subB)) {
      if (length(subB) != ncol(Lambda)) {
        stop("Invalid length")
      }
      if (subB[co]) {
        stop("Covariance matrix should be shed with respect to conditioning site")
      }
    } else if (co %in% subB) { # is numeric
      stop("Invalid conditioning variable for covariance matrix")
    }
  }
  2 * (outer(Lambda[co, subA], Lambda[co, subB], "+") - Lambda[subA, subB])
}
