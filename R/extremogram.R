#' Pairwise extremogram for max-risk functional
#'
#' The function computes the pairwise \eqn{chi} estimates and plots them as a function of the distance between sites.
#' @param dat data matrix
#' @param margp marginal probability above which to threshold observations
#' @param coord matrix of coordinates (one site per row)
#' @param scale geometric anisotropy scale parameter
#' @param rho geometric anisotropy angle parameter
#' @param plot logical; should a graph of the pairwise estimates against distance? Default to \code{FALSE}
#' @param ... additional arguments passed to plot
#' @return an invisible matrix with pairwise estimates of chi along with distance (unsorted)
#' @export
#' @examples
#' \dontrun{
#' lon <- seq(650, 720, length = 10)
#' lat <- seq(215, 290, length = 10)
#' # Create a grid
#' grid <- expand.grid(lon,lat)
#' coord <- as.matrix(grid)
#' dianiso <- distg(coord, 1.5, 0.5)
#' sgrid <- scale(grid, scale = FALSE)
#' # Specify marginal parameters `loc` and `scale` over grid
#' eta <- 26 + 0.05*sgrid[,1] - 0.16*sgrid[,2]
#' tau <- 9 + 0.05*sgrid[,1] - 0.04*sgrid[,2]
#' # Parameter matrix of Huesler--Reiss
#' # associated to power variogram
#' Lambda <- ((dianiso/30)^0.7)/4
#' # Regular Euclidean distance between sites
#' di <- distg(coord, 1, 0)
#' # Simulate generalized max-Pareto field
#' set.seed(345)
#' simu1 <- rgparp(n = 1000, thresh = 50, shape = 0.1, riskf = "max",
#'                 scale = tau, loc = eta, sigma = Lambda, model = "hr")
#' extdat <- extremo(dat = simu1, margp = 0.98, coord = coord,
#'                   scale = 1.5, rho = 0.5, plot = TRUE)
#'
#' # Constrained optimization
#' # Minimize distance between extremal coefficient from fitted variogram
#' mindistpvario <- function(par, emp, coord){
#' alpha <- par[1]; if(!isTRUE(all(alpha > 0, alpha < 2))){return(1e10)}
#' scale <- par[2]; if(scale <= 0){return(1e10)}
#' a <- par[3]; if(a<1){return(1e10)}
#' rho <- par[4]; if(abs(rho) >= pi/2){return(1e10)}
#' semivariomat <- power.vario(distg(coord, a, rho), alpha = alpha, scale = scale)
#'   sum((2*(1-pnorm(sqrt(semivariomat[lower.tri(semivariomat)]/2))) - emp)^2)
#' }
#'
#' hin <- function(par, ...){
#'   c(1.99-par[1], -1e-5 + par[1],
#'     -1e-5 + par[2],
#'     par[3]-1,
#'     pi/2 - par[4],
#'     par[4]+pi/2)
#'   }
#' opt <- alabama::auglag(par = c(0.7, 30, 1, 0),
#'                        hin = hin,
#'                         fn = function(par){
#'                           mindistpvario(par, emp = extdat[,'prob'], coord = coord)})
#' stopifnot(opt$kkt1, opt$kkt2)
#' # Plotting the extremogram in the deformed space
#' distfa <- distg(loc = coord, opt$par[3], opt$par[4])
#' plot(
#'  x = c(distfa[lower.tri(distfa)]), 
#'  y = extdat[,2], 
#'  pch = 20,
#'  yaxs = "i", 
#'  xaxs = "i", 
#'  bty = 'l',
#'  xlab = "distance", 
#'  ylab= "cond. prob. of exceedance", 
#'  ylim = c(0,1))
#' lines(
#'   x = (distvec <- seq(0,200, length = 1000)), 
#'   col = 2, lwd = 2,
#'   y = 2*(1-pnorm(sqrt(power.vario(distvec, alpha = opt$par[1],
#'                                scale = opt$par[2])/2))))
#' }
extremo <- function(dat, margp, coord, scale = 1, rho = 0, plot = FALSE, ...){
  stopifnot(isTRUE(all(margp >= 0, margp < 1, length(margp) == 1, nrow(coord) == ncol(dat))))
  # Local quantile - threshold data
  dat <- as.matrix(dat)
  margthresh <- apply(dat, 2, quantile, margp, na.rm = TRUE)
  dat <- t(t(dat)-margthresh)
  # Keep only instances where there is at least one exceedance
  dimat <- distg(coord, scale = scale, rho = rho)
  excind <- apply(dat, 1, function(x){isTRUE(max(x, na.rm = TRUE) > 0)}) #avoid NA
  dat <- dat[excind,]
  res <- matrix(0, ncol = 4, nrow = choose(ncol(dat),2))
  b <- 0L
  for(i in 1:(ncol(dat)-1)){
    for(j in (i+1):ncol(dat)){
      b <- b + 1L
      subdat <- na.omit(dat[,c(i,j)])
      if(length(subdat) > 0){
        res[b, ] <- c(i,j, nrow(subdat), mean(I(subdat[,2] > 0) * I(subdat[,1] > 0))/(0.5*mean(I(subdat[,1] > 0)) + 0.5*mean(I(subdat[,2] > 0))) )
      } else{
       res[b, ] <- NA
      }
    }
  }
  res <- na.omit(res)
  ellips <- list(...)
  ellips$y <- res[,4]
  ellips$x <- apply(res, 1, function(x){dimat[x[1],x[2]]})
  if(plot){
    if(is.null(ellips$xlab)){ellips$xlab <- "distance"}
    if(is.null(ellips$ylab)){ellips$ylab <- "conditional probability of exceedance"}
    if(is.null(ellips$pch)){ ellips$pch <- 20}
    if(is.null(ellips$yaxs)){ ellips$yaxs <- "i"}
    if(is.null(ellips$xlim)){ ellips$xlim <- c(0, max(ellips$x)*1.02)}
    if(is.null(ellips$xaxs)){ ellips$xaxs <- "i"}
    if(is.null(ellips$bty)){ ellips$bty <- "l"}
    if(is.null(ellips$ylim)){ ellips$ylim <- c(0,1)}
    if(is.null(ellips$pch)){ ellips$pch <- 20}
    if(is.null(ellips$col)){ ellips$col <- grDevices::rgb(0, 0, 0, alpha = 0.25)}
  do.call(what = graphics::plot, args = ellips)
  }
  return(invisible(cbind(dist = ellips$x, prob = ellips$y)))
}
