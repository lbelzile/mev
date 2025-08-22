#' Kernel-based threshold selection of Goegebeur, Beirlant and de Wet (2008)
#' @param xdat vector of observations
#' @param kmax maximum number of exceedances considered
#' @param kernel [string] kernel choice, one of \code{Jackson} or \code{Lewis}
#' @param rho string for the estimator of the second order regular variation. Can also be a negative scalar
thselect.gbw <- function(
  xdat,
  kmax,
  kernel = c("Jackson", "Lewis"),
  rho = c("gbw08", "ghp", "fagh", "cm", "gbw10"),
  ...
) {
  args <- list(...)
  kmax <- as.integer(kmax)
  args$k <- 10:kmax
  xdat <- sort(
    x = xdat[is.finite(xdat) & xdat > 0],
    decreasing = TRUE
  )
  n <- length(xdat)
  stopifnot(kmax <= n)
  k <- sort(as.integer(k))
  logdata <- log(xdat)
  Z <- 1:(kmax - 1) * (logdata[1:(kmax - 1)] - logdata[2:kmax])
  args$xdat <- xdat
  if (is.numeric(rho)) {
    stopifnot(length(rho) == 1L, rho < 0)
    erho <- rep(rho, length.out = length(k))
  } else {
    # String for estimator
    rho <- match.arg(rho)
    erho <- do.call(what = paste0("rho.", rho), args = args)$rho
  }
  kernel <- match.arg(kernel)
  if (isTRUE(any(erho == -1)) & kernel == "Lewis") {
    stop("Lewis kernel does not allow for rho=-1.")
  }
  kernel_jack <- function(u) {
    -1 - log(u)
  }

  kernfun <- switch(
    EXPR = kernel,
    Jackson = kernel_jack,
    Lewis = kernel_lewis
  )
  if(kernel == "Lewis"){
    kernfun <- function(u) {
      u - 0.5
    }
    Tstat <- sapply(k, function(ks){
      mean(kernfun((1:ks)/(ks+1))*Z[1:ks])/mean(Z[1:ks])
    })
  gof <- 1/k + (2*(2-erho)/(|abs(erho)|)*Tstat)^2
  } else if(kernel == "Jackson"){
    Kint <- erho/(1-erho)^2
    gof <- 1/k +
      sapply(seq_along(k), function(i){
        mean(kernfun((1:k[i])/(k[i]+1))*Z[1:k[i]])/(mean(Z[1:k[i]])*Kint[i]*(1-erho[i]))
      })^2
  }
  k0 <- k[which.min(gof)]
  return(
    list(k0,
         shape = mean(Z[1:k0]),
         rho = erho[which.min(gof)])
}
