#' Kernel-based threshold selection of Goegebeur, Beirlant and de Wet (2008)
#'
#' @param xdat [vector] sample exceedances
#' @param kmax [int] maximum number of exceedances considered
#' @param kernel [string] kernel choice, one of \code{Jackson} or \code{Lewis}
#' @param rho string for the estimator of the second order regular variation. Can also be a negative scalar
#' @param ... additional arguments, for backward compatibility purposes
#' @return a list with elements
#' \itemize{
#' \item{\code{k0}:} number of exceedances
#' \item{\code{shape}:} Hill's shape estimate
#' \item{\code{rho}:} second-order regular variation parameter estimate
#' \item{\code{gof}:} goodness-of-fit statistic for the chosen threshold.
#' }
#' @references Goegebeur , Y., Beirlant , J., and de Wet , T. (2008). Linking Pareto-Tail Kernel Goodness-of-fit Statistics with Tail Index at Optimal Threshold and Second Order Estimation. REVSTAT-Statistical Journal, 6(\bold{1}), 51–69. <doi:10.57805/revstat.v6i1.57>
#' @export
#' @examples
#' xdat <- rgp(n = 1000, scale = 2, shape = 0.5)
#' (thselect.gbw(xdat, kmax = 500))
thselect.gbw <- function(
  xdat,
  kmax,
  kernel = c("Jackson", "Lewis"),
  rho = c("gbw", "ghp", "fagh", "dk"),
  ...
) {
  args <- list(...)
  n <- length(xdat)
  if (missing(kmax)) {
    kmax <- n - 1L
  } else {
    kmax <- min(as.integer(kmax), n - 1L)
  }
  k <- args$k <- 10:kmax
  xdat <- sort(
    x = xdat[is.finite(xdat) & xdat > 0],
    decreasing = TRUE
  )
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
  if (kernel == "Lewis") {
    # print("Lewis")
    kernfun <- function(u) {
      u - 0.5
    }
    Tstat <- sapply(k, function(ks) {
      mean(kernfun((1:ks) / (ks + 1)) * Z[1:ks]) / mean(Z[1:ks])
    })
    gof <- 1 / k + (2 * (2 - erho) / (abs(erho)) * Tstat)^2
  } else if (kernel == "Jackson") {
    # print("Jackson")
    kernfun <- function(u) {
      -1 - log(u)
    }
    Kint <- erho / (1 - erho)^2
    gof <- 1 /
      k +
      sapply(seq_along(k), function(i) {
        mean(kernfun((1:k[i]) / (k[i] + 1)) * Z[1:k[i]]) /
          (mean(Z[1:k[i]]) * Kint[i] * (1 - erho[i]))
      })^2
  }
  k0 <- k[which.min(gof)]
  res <- list(
    k0 = k0,
    thresh0 = xdat[k0],
    shape = mean(Z[1:k0]),
    rho = erho[which.min(gof)],
    method = ifelse(is.numeric(rho), "user-supplied", rho),
    kernel = kernel,
    gof = gof[k0]
  )
  class(res) <- "mev_thselect_gbw"
  return(invisible(res))
}

#' @export
print.mev_thselect_gbw <- function(
  x,
  digits = min(3, getOption("digits") - 4),
  ...
) {
  cat("Threshold selection method:", x$kernel, "kernel", "\n")
  cat("Goegebeur, Beirlant and de Wet (2008)\n")
  cat(
    paste0("Second-order regular variation index (", x$method),
    "estimator): ",
    round(x$rho, digits),
    "\n"
  )
  cat("Number of exceedances:", x$k0, "\n")
  cat("Selected threshold:", round(x$thresh0, digits), "\n")
  cat("Shape estimate:", round(x$shape, digits), "\n")
  return(invisible(NULL))
}
