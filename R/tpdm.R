#' Tail pairwise dependence matrix
#'
#' Given a multivariate sample of observations and a radial threshold quantile level,
#' estimate the tail pairwise dependence matrix (TPDM) empirically through marginal
#' transformation of the margins to unit Fréchet scale. This is entirely equivalent to
#' taking the definition of Cooley and Thibaud (2019) with GEV(1,0.5,0.5) margins.
#' @param xdat matrix of observations
#' @param qlev quantile level of threshold for the radial component
#' @param ties.method method for ties; see \link[base]{rank} for more details
#' @param margtrans string; if \code{"emp"} (default), apply a rank transformation to map observations to unit Frechet margins
#' @param standardize logical; if \code{TRUE} (default), matrix is standardized to correlation matrix
#' @return a positive definite matrix
#' @references Larsson, M. and S.I. Resnick(2012). \emph{Extremal dependence measure and extremogram: the regularly varying case}, Extremes 15, 231–256. <doi:10.1007/s10687-011-0135-9>
#' @references Cooley, D. and E. Thibaud (2019). \emph{Decompositions of dependence for high-dimensional extremes}, Biometrika, 106(\bold{3}), 587-604. <doi:10.1093/biomet/asz028>
#' @references Kiriliouk, A. and C. Zhou (2024+) Estimating probabilities of multivariate failure sets based on pairwise tail dependence coefficients, arXiv, <doi:10.48550/arXiv.2210.12618>
#' @examples
#' d <- 4L
#' xdat <- rmev(n = 1000, d = d, param = 0.9)
#' xdep.tpdm(xdat = xdat, qlev = 0.5, margtrans = "none")
#' # Equicorrelation matrix
#' Sigma <- 0.9 * diag(d) + matrix(0.1, d, d)
#' xdat <- rmnorm(n = 10000, mu = rep(0, d), Sigma = Sigma)
#' xdep.tpdm(xdat = xdat, qlev = 0.99)
#' @export
xdep.tpdm <- function(
  xdat,
  qlev,
  ties.method = "random",
  margtrans = c("emp", "none"),
  standardize = TRUE
) {
  stopifnot(is.matrix(xdat), ncol(xdat) >= 2L)
  d <- ncol(xdat)
  n <- nrow(xdat)
  ties.method <- match.arg(
    arg = ties.method,
    choices = eval(formals(rank)$ties.method)
  )
  margtrans <- match.arg(margtrans)
  standardize <- isTRUE(standardize[1])
  stopifnot(length(qlev) == 1L, isTRUE(qlev < 1 - 10 / n), isTRUE(qlev >= 0))
  if (margtrans == "emp") {
    fdat <- mev::qgev(
      p = apply(xdat, 2, rank, ties.method = ties.method) / (n + 1),
      loc = 1,
      scale = 1,
      shape = 1
    )
    # fdat2 <- mev::qgev(
    #   p = apply(xdat, 2, rank, ties.method = ties.method) / (n + 1),
    #   loc = 1,
    #   scale = 0.5,
    #   shape = 0.5
    # )
  } else if (margtrans == "none") {
    fdat <- xdat
  }
  rad <- rowSums(fdat)
  # rad2 <- sqrt(rowSums(fdat2^2))
  thresh <- quantile(rad, probs = qlev)
  # thresh2 <- quantile(rad2, probs = qlev)
  is_exc <- rad > thresh & is.finite(rad)
  wdat <- fdat[is_exc, ] / rad[is_exc]
  # wdat2 <- fdat2[is_exc, ] / rad2[is_exc]
  nexc <- nrow(wdat)
  tpdm <- (d / nexc) * crossprod(sqrt(wdat))
  # tpdm2 <- (d / nexc) * crossprod(wdat2)
  if (standardize) {
    tpdm <- cov2cor(tpdm)
  }
  return(tpdm)
}
