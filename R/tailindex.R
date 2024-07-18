#' Extreme U-statistic Pickands estimator
#'
#' Given a random sample of \code{n} exceedances, the estimator
#' returns an estimator of the shape parameter or extreme
#' value index using a kernel of order 3, based on
#' \code{k} largest exceedances of \code{xdat}. Oorschot et al. (2023) parametrize the model in terms of \eqn{m=n-k}, so that \eqn{m=n} corresponds to using only the three largest observations.
#'
#' The calculations are based on the recursions provided in Lemma 4.3 of Oorschot et al.
#' @param xdat vector of observations of length \eqn{n}
#' @param k number of largest order statistics \eqn{3 \leq k \leq n}.
#' @references Oorschot, J, J. Segers and C. Zhou (2023), Tail inference using extreme U-statistics,  Electron. J. Statist. 17(1): 1113-1159. \doi{10.1214/23-EJS2129}
#' @export
#' @examples
#' xdat <- rgp(n = 1000, shape = 0.2)
#' shape.pickandsxu(xdat, k = 10)
shape.pickandsxu <- function(xdat, k, ...){
  args <- list(...)
  xdat <- na.omit(as.numeric(xdat))
  xdat <- sort(xdat,
               decreasing = TRUE)
  n <- length(xdat)
  if(!is.null(args$m) & missing(k)){
   m <- as.integer(args$m)
  } else{
   m <- n - as.integer(k)
  }
  k <- n - m
  m <- sort(m, decreasing = TRUE)
  stopifnot(m[length(m)] >= 3,
            is.vector(xdat),
            length(xdat) >= m[1])

  shape <- numeric(length(m))
  for(i in seq_along(m)){
    ms <- m[i]
  # Initial recursion j=2
  shape0 <- (2*(n-1)/(ms-2)-2)*log(xdat[1] - xdat[2])
  mcst <- 1
  for(j in seq.int(
    from = 3L,
    to = n - ms + 3L,
    by = 1L)){
    mcst <- mcst * (n - j - ms + 4L) / (n - j + 1L)
     shape0 <- shape0 + mcst * (2 * (n - j + 1) / (ms - 2) - j) * sum(log(xdat[seq_len(j-1)] - xdat[j]))
  }
   shape[i] <- shape0 * ms * (ms - 1) * (ms - 2) / (n * (n - 1) * (n - ms + 1))
  }
  if(length(m) == 1L){
  return(as.numeric(shape))
  } else{
    return(data.frame(k = as.integer(k), shape = as.numeric(shape)))
  }
}

#' Extreme U-statistic Pickands estimator
#'
#' [Deprecated function]
#' Given a random sample of \code{n} exceedances, the estimator
#' returns an estimator of the shape parameter or extreme
#' value index using a kernel of order 3, based on
#' \code{k} largest exceedances of \code{xdat}. Oorschot et al. (2023) parametrize the model in terms of \eqn{m=n-k}, so that \eqn{m=n} corresponds to using only the three largest observations.
#'
#' The calculations are based on the recursions provided in Lemma 4.3 of Oorschot et al.
#' @param xdat vector of observations of length \eqn{n}
#' @param k number of largest order statistics \eqn{3 \leq k \leq n}.
#' @references Oorschot, J, J. Segers and C. Zhou (2023), Tail inference using extreme U-statistics,  Electron. J. Statist. 17(1): 1113-1159. \doi{10.1214/23-EJS2129}
#' @export
#' @keywords internal
PickandsXU <- function(xdat, m){
  .Deprecated(new = "fit.pickandsxu",
              package = "mev")
 shape.pickandsxu(xdat = xdat, m = m)
}

#' Hill's estimator of the shape parameter
#'
#' Given a sample of positive observations, calculate the tail index or
#' shape parameter. The shape estimate returned is positive.
#'
#' @param xdat vector of positive observations
#' @param kmin minimum number of upper order statistics (exceedances) for the estimator
#' @param kmax maximum number of upper order statistics (exceedances) for the estimator
#' @return a data frame with the number of order statistics \code{k} and the shape parameter estimate \code{shape}, or a single numeric value if \code{kmin} equals \code{kmax}.
#' @references Hill, B.M. (1975). \emph{A simple general approach to inference about the tail of a distribution.} Annals of Statistics, \bold{3}, 1163-1173.
#' @export
#' @examples
#' xdat <- mev::rgp(n = 200, loc = 1, scale = 0.5, shape = 0.5)
#' shape.hill(xdat)
shape.hill <- function(xdat, k) {
  logdata <- log(xdat[xdat>0 & is.finite(xdat)])
  logdata <- sort(logdata, decreasing = TRUE)
  if(missing(k)){
    k <- seq(10, length(logdata) - 1)
  }
  n <- length(logdata)
  k <- as.integer(sort(k))
  kmin <- k[1]
  kmax <- min(k[length(k)], n)
  n <- min(kmax + 1, n)
  ks <- kmin:(n - 1)
  cumlogdat <- cumsum(logdata[1:n])
  shape <- cumlogdat[ks] / ks - logdata[ks + 1]
  if(length(k) > 1){
    return(
      data.frame(
        k = k,
        shape = shape[k %in% ks]
    )
    )
  } else{
    return(as.numeric(cumsum(logdata[1:k])[k] / k - logdata[k + 1]))
  }
}

#' Beirlant et al. generalized quantile shape estimator
#'
#' This estimator estimates the real shape parameter based on generalized quantile plots based on mean excess functions, generalized median excesses or trimmed mean excesses.
#'
#' @references Beirlant, J., Vynckier P. and J.L. Teugels (1996). \emph{Excess functions and estimation of the extreme-value index.} Bernoulli, 2(\bold{4}), 293-318.
#' @export
#' @param xdat \code{n} vector of observations
#' @param k number of upper order statistics
#' @param type string indicating the estimator choice, one of \code{genmean}, \code{genmed} and \code{trimmean}.
#' @param weight weight a kernel function on \eqn{[0,1]}
shape.genquant <- function(xdat,
                       k,
                       type = c("genmean","genmed","trimmean"),
                       weight,
                       p = 0.5){
  type <- match.arg(type)
  if(type %in% c("genmed", "trimmean")){
    stopifnot(isTRUE(all(length(p) == 1L, p > 0, p < 1)))
  }
  xdat <- sort(xdat[is.finite(xdat) & xdat > 0],
               decreasing = TRUE)
  n <- length(xdat)
  k <- as.integer(sort(k))
  kmax <- k[length(k)]
  if(k[1] < 5){
    stop("Invalid argument: \"k\" must be larger than 5 to reliably estimate the shape parameter.")
  }
  stopifnot(kmax < n)
  weight_fun <- FALSE
  if(missing(weight)){
   we <- 1
  }  else{
    if(inherits(weight, "function")){
      weight_fun <- TRUE
    } else{
    if(length(k) == 1L & is.numeric(weight)){
      stopifnot(length(weight) %in% c(1, k))
      we <- weight
    } else{
      stop("Invalid \"weight\" function.")
    }
    }
  }
  if(type == "genmean"){
    hill <- shape.hill(k = seq_len(kmax+1), xdat = xdat)
    ES <- log(xdat[hill$k]*hill$shape)
  } else if(type == "genmed"){
    ES <- log(xdat[floor(p*seq_len(kmax+1))+1] - xdat[seq_len(kmax+1)+1])
  } else if(type == "trimmean"){
    ES <- log(sapply(seq_len(kmax+1), FUN = function(ks){
      mean(xdat[(floor(p*ks)+1):ks])
    }) - xdat[seq_len(kmax+1)+1])
  }
    shape <- length(k)
    for(j in seq_along(k)){
      ks <- k[j]
      if(weight_fun){
       we <- weight_fun((1:ks)/(ks+1))
      }
      shape[j] <- sum(we * (log(ks + 1) - log(1:ks))*(ES[1:ks] - ES[ks+1]))/ sum(we * (log(ks + 1) - log(1:ks))^2)
    }
    # Return estimator
  if(length(k) == 1L){
    return(as.numeric(shape))
  } else{
    return(data.frame(k = k, shape = shape))
  }


}



#' Pickand's shape estimator
#'
#' Given a sample of size \code{n} of positive exceedances, compute
#' the real shape parameter \eqn{\xi} based on the \code{k} largest order statistics.
#' @export
#' @references Pickands, III, J. (1975). \emph{Statistical inference using extreme order statistics.} Annals of Statistics, \bold{3}, 119-131.
#' @param xdat vector of positive observations of length \eqn{n}
#' @param k number of largest order statistics
#' @return a data frame with the number of order statistics \code{k} and the shape parameter estimate \code{shape}, or a single numeric value if \code{k} is a scalar.
shape.pickands <- function(xdat, k){
  k <- sort(as.integer(k))
  xdat <- sort(xdat[is.finite(xdat) & xdat > 0],
               decreasing = FALSE)
  n <- length(xdat)
  stopifnot(isTRUE(all(n >= k, k >= 5)))
  shape <- numeric(length = length(k))
  for(i in seq_along(k)){
    shape[i] <- as.numeric(
      (log(xdat[n - floor(k[i]/4)] - xdat[n - floor(k[i]/2)]) -
         log(xdat[n-floor(k[i]/2)] - xdat[n-k[i]]))
      )
  }
  shape <- shape / log(2)
  if(length(k) > 1){
   return(
     data.frame(k, shape)
     )
  } else{
    return(shape)
  }
}


#' Dekkers and de Haan moment estimator for the shape
#'
#' Given a sample of exceedances, compute the moment estimator of the real shape parameter.
#'
#' @references Dekkers, A.L.M. and de Haan, L. (1989). \emph{On the estimation of the extreme-value index and large quantile estimation.}, Annals of Statistics, \bold{17}, 1795-1833.
#' @export
#' @inheritParams shape.pickands
#' @return a data frame with the number of order statistics \code{k} and the shape parameter estimate \code{shape}, or a single numeric value if \code{k} is a scalar.
shape.moment <- function(xdat, k){
  k <- sort(as.integer(k))
  xdat <- sort(xdat[is.finite(xdat) & xdat > 0],
               decreasing = TRUE)
  logdata <- as.numeric(log(xdat))
  n <- min(k[length(k)], length(logdata))
  stopifnot(k[1] >= 5,
            k[length(k)] < length(logdata))
  cumlogdat <- cumsum(logdata[1:n])
  M1 <- cumlogdat[k] / k - logdata[k + 1]
  shape <- numeric(length = length(k))
  for(i in seq_along(k)){
   M2 <- mean((logdata[1:k[i]] - logdata[k[i]+1])^2)
   shape[i] <- M1[i]+1-0.5/(1-M1[i]^2/M2)
  }
  if(length(k) > 1){
    return(
      data.frame(k, shape)
    )
  } else{
    return(shape)
  }
}

#' Shape parameter estimates
#'
#' Wrapper to estimate the tail index or shape parameter of an extreme value distribution. Each function has similar sets of arguments, a vector or scalar number of order statistics \code{k} and
#' a vector of positive observations \code{xdat}. The \code{method} argument allows users to choose between different indicators, including the Hill estimator (\code{hill}, for positive observations and shape only), the moment estimator of Dekkers and de Haan (\code{mom} or \code{dekkers}), the Beirlant, Vynckier and Teugels generalized quantile estimator (\code{bvt} or \code{genquant}), the Pickands estimator (\code{pickands}), the extreme \eqn{U}-statistics estimator of Oorschot, Segers and Zhou (\code{osz}, or \code{pickandsxu}).
#' @export
#' @inheritParams shape.pickands
#' @param method estimation method.
#' @param ... additional parameters passed to functions
#' @return a data frame with the number of order statistics \code{k} and the shape parameter estimate \code{shape}, or a single numeric value if \code{k} is a scalar.
fit.shape <- function(xdat,
                      k,
                      method = c("hill","pickandsxu","osz","mom","dekkers","bvt","genquant","pickands"),
                      ...){
  method <- match.arg(method)
  if(method == "hill"){
    return(shape.hill(xdat = xdat, k = k))
  } else if(method %in% c("pickandsxu","osz")){
    return(shape.pickandsxu(xdat = xdat, k = k))
  } else if(method %in% c("mom","dekkers")){
   return(shape.moment(xdat = xdat, k = k))
  } else if(method %in% c("beirlant","genquant")){
   return(shape.genquant(xdat = xdat, k = k, ...))
  } else if(method == "pickands"){
    return(shape.pickands(xdat = xdat, k = k))
  }
}
