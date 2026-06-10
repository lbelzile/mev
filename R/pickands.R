#' Estimator of the Pickands dependence function
#'
#' This function computes the nonparametric angular measure
#' of the multivariate observations using a empirical rank
#' transformation to unit Frechet margins, and from there
#' calculates the weights for the self-concordant empirical
#' likelihood imposing the mean constraint on the angular measure
#'
#' \deqn{ \hat{A}(\boldsymbol{t}) = \sum_{i=1}^n p_i \max_{j=1}^d (w_{i,j}i t_j), \quad t \in \mathbb{S}_d}
#' @param xdat \code{n} by \code{d} matrix of observations
#' @param w \code{m} by \code{d} matrix of angles on the unit simplex
#' @param qlev quantile level for the threshold
#' @param region risk region for determining angles
#' @return a vector of length \code{m} with Pickands dependence function
#' @export
#' @examples
#' ang <- seq(0,1, by = 0.01)
#' xdat <- rmev(n = 1000, d = 2, param = 0.4)
#' pickands <- xdep.pickands(xdat, w = ang, qlev = 0)
#' plot(type = "n", x = 0.5,
#'      y = 1,
#'      xlim = c(0, 1),
#'      ylim = c(0.5, 1),
#'      xlab = "w",
#'      ylab = "Pickands dependence measure",
#'      bty = "l")
#' segments(x1 = 0.5, x0 = 0, y1 = 0.5, y0 = 1, col = "gray90")
#' segments(x1 = 0.5, x0 = 1, y1 = 0.5, y0 = 1, col = "gray90")
#' segments(x1 = 0, x0 = 1, y1 = 1, y0 = 1, col = "gray90")
#' lines(ang, pickands, lwd = 2)
xdep.pickands <- function(xdat, w, qlev, region = c('sum','max', 'min')){
  if(is.vector(xdat)){
    stop("Invalid input")
  }
  xdat <- as.matrix(xdat)
  d <- ncol(xdat)
  stopifnot(d >= 2)
  # Calculate angular measure
  if(region[1] %in% c("max","min")){
   qlev <- rep(qlev, d)
  }
  angm <- angmeas(xdat = xdat, thresh = qlev, region = region[1])
  angm$ang <- cbind(angm$ang, 1-rowSums(angm$ang))
  if(is.vector(w)){
    if(d == 2L){
     w <- cbind(w, 1-w)
    } else if(length(w) %in% c(d, d-1)){
      w <- matrix(w, nrow = 1)
    }
  }
  stopifnot(ncol(w) %in% c(d, d-1))
  if(ncol(w) == (d-1)){
    w <- cbind(w, 1-rowSums(w))
  } else if(ncol(w) == d){
    stopifnot(isTRUE(abs(sum(w[1,]) - 1) < 1e-12))
  }
  pickands.fn <- function(wpt){
    sum(angm$wts * d * apply(t(angm$ang) * wpt, 2, max))
  }
  pickands <- apply(w, 1, pickands.fn)
  return(pickands)
}


# Pickands estimator of Eastoe, Heffernan and Tawn
#
# @param w angles on the D-simplex at which to evaluate Pickands dependence function
# @param xdat data, transformed to unit Frechet margins through a rank transform if \code{standardize = TRUE}.
# @param ties.method method for ties; see \link[base]{rank} for more details
# @param margtrans string
# xdep.pickands <- function(
#     w,
#     xdat,
#     qlev,
#     ties.method = eval(formals(rank)$ties.method),
#     method = "eht",
#     standardize = c("emp","none"){
#   stopifnot(is.matrix(xdat), ncol(xdat) == 2L)
#   d <- ncol(xdat)
#   n <- nrow(xdat)
#   stopifnot(is.matrix(w), ncol(w) %in% c(d-1, d), min(w) < 0)
#   if(ncol(w) == (d-1)){
#    w <- cbind(w, 1-rowSums(w))
#   }
#   emp_ang_j <- function(fdat, wdat, t, j, u){
#    mean(apply(wdat[fdat[,j] > u,], 1, function(
#     w){ isTRUE(all(w < t))}))
#     sum(fdat[,j] >  u)
#   }
#   ties.method <- match.arg(ties.method)
#     stopifnot(length(qlev) %in% c(1L, d), isTRUE(all(qlev < 1, qlev > 0)))
#   qlev <- rep(qlev, d)
#   u <- mev::qgev(qlev, loc = 1, scale = 1, shape = 1)
#  fdat <- mev::qgev(apply(xdat, 2, rank, ties.method = ties.method) / (nrow(xdat) + 1), loc = 1, scale = 1, shape = 1)
#  ang <- fdat[,-d] / rowSums(fdat)
#  pick <- numeric(nrow(w))
#
#  for(pt in seq_along(pp)){
#  pick[pt] <- mean(((fdat[,1] > u[1]) / (1-qlev[1]) + (fdat[,2] > u[2]) /  (1-qlev[2]))* pmax(pp[pt] * ang, (1 - pp[pt]) * (1 - ang)))
#  }
# }
