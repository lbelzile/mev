#' Generalized extreme value distribution
#'
#' Density function, distribution function, quantile function and
#' random number generation for the generalized extreme value
#' distribution.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n scalar number of observations
#' @param loc scalar or vector of location parameters whose length matches that of the input
#' @param scale scalar or vector of positive scale parameters whose length matches that of the input
#' @param shape scalar shape parameter
#' @param log,log.p logical; if \code{TRUE}, probabilities \eqn{p} are given as
#'   \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), returns the distribution function, otherwise the survival function
#' @details The distribution function of a GEV distribution with parameters
#'  \code{loc} = \eqn{\mu}, \code{scale} = \eqn{\sigma} and
#'  \code{shape} = \eqn{\xi} is
#'  \deqn{F(x) = \exp\{-[1 + \xi (x - \mu) / \sigma] ^ {-1/\xi} \}}{%
#'        F(x) = exp{ -[1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)} }
#'  for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  If \eqn{\xi = 0} the
#'  distribution function is defined as the limit as \eqn{\xi} tends to zero.
#'
#'  The quantile function, when evaluated at zero or one,
#'  returns the lower and upper endpoint, whether the latter is finite or not.
#'
#' @references Jenkinson, A. F. (1955) The frequency distribution of the
#'   annual maximum (or minimum) of meteorological elements.
#'   \emph{Quart. J. R. Met. Soc.}, \strong{81}, 158-171.
#'   Chapter 3: \doi{10.1002/qj.49708134804}
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   \doi{10.1007/978-1-4471-3675-0_3}
#' @name gevdist
#' @importFrom utils bibentry
NULL
## NULL

#' @rdname gevdist
#' @export
qgev <- function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE){
  n <- length(p)
  if (n == 0L) {
    return(numeric(0))
  }
  stopifnot(length(loc) %in% c(1L, n),
            length(scale) %in% c(1L, n),
            length(shape) == 1L,
            isTRUE(all(is.finite(scale))),
            is.finite(shape),
            isTRUE(all(is.finite(loc))),
            isTRUE(min(scale, na.rm = TRUE) > 0),
            is.logical(lower.tail),
            length(lower.tail) == 1L,
            is.logical(log.p),
            length(log.p) == 1L)
  log.p <- as.logical(log.p)
   if(isTRUE(log.p)){
     stopifnot(max(p, na.rm = TRUE) <= 0)
   } else{
    stopifnot(min(p, na.rm = TRUE) >= 0,
             max(p, na.rm = TRUE) <= 1)
   }
     if (!lower.tail){
      if(log.p){
       log.p <- FALSE
       p <- exp(p)
      }
        p <- 1 - p
     }
    if(!log.p){
     logp <- log(p)
    } else{
      logp <- p
    }
    if (isTRUE(all.equal(shape, 0, check.attributes = FALSE, tolerance = 1e-10))){
      return(loc - scale * log(-logp))
    } else if (abs(shape) < 1e-6){
      return(loc - scale * ifelse(
        is.infinite(logp),
        ifelse(shape > 0,
               1/shape,
               Inf),
        ifelse(logp == 0,
               ifelse(shape < 0,
                1/shape,
                -Inf),
               log(-logp) * (1 - log(-logp) * shape/2)))
      )
    } else{
    return(loc + scale * ((-logp)^(-shape) - 1)/shape)
    }
}

#' @rdname gevdist
#' @export
rgev <- function(n, loc = 0, scale = 1, shape = 0){
  n <- as.integer(n)
  stopifnot(is.finite(n),
	    n > 0)
  qgev(p = -rexp(n),
       loc = loc,
       scale = scale,
       shape = shape,
       log.p = TRUE)
}

#' @rdname gevdist
#' @export
dgev <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE) {
  n <- length(x)
    if (n == 0L) {
    return(numeric(0))
  }
  stopifnot(length(loc) %in% c(1L, n),
            length(scale) %in% c(1L, n),
            length(shape) == 1L,
            isTRUE(all(is.finite(loc))),
            isTRUE(all(is.finite(scale))),
            isTRUE(min(scale, na.rm = TRUE) > 0),
            is.finite(shape),
            is.logical(log),
            length(log) == 1L
            )
  x <- (x - loc) / scale
  if(isTRUE(all.equal(shape, 0, check.attributes = FALSE, tolerance = 1e-10))){
  # Gumbel subcase
  d <- -log(scale) - x - exp(-x)
  } else{
  xx <- pmax(-1, shape * x)
  # Negative numbers get mapped to zero, as well as points outside the support
   d <- -log(scale) - (1 + xx)^(-1/shape) - (1/shape + 1) * log1p(xx)
   d[(xx <= -1) | is.infinite(xx)] <- -Inf # density for NA/Inf is zero
   }
    if (!log){
     return(exp(d))
    } else{
     return(d)
    }
  }

# ----------------------------- pgev ---------------------------------

#' @rdname gevdist
#' @export
#' @author Leo Belzile, with code adapted from Paul Northrop
pgev <- function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                 log.p = FALSE) {
  n <- length(q)
  stopifnot(length(loc) %in% c(1L, n),
            length(scale) %in% c(1L, n),
            length(shape) == 1L,
            isTRUE(all(is.finite(loc))),
            isTRUE(all(is.finite(scale))),
            isTRUE(min(scale, na.rm = TRUE) > 0),
            is.finite(shape),
            is.logical(lower.tail),
            length(lower.tail) == 1L,
            is.logical(log.p),
            length(log.p) == 1L
            )
  if (n == 0L) {
    return(numeric(0))
  }
  q <- (q - loc) / scale

  if(abs(shape) > 1e-6){
    p <-  -pmax(1 + shape * q, 0)^(-1 / shape)
  } else {
     p <- ifelse(is.infinite(q),
              log((1 + sign(q)) / 2), # 0 or 1
             -exp(-q + shape * q^2 / 2))
   }
  if (lower.tail) {
    if (!log.p) {
      p <- exp(p)
    }
  } else {
    p <- -expm1(p)
    if (log.p) {
      p <- log(p)
    }
  }
  return(p)
}

#' Generalized Pareto distribution
#'
#' Density function, distribution function, quantile function and
#' random number generation for the generalized Pareto
#' distribution.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n scalar number of observations
#' @param loc location parameter.
#' @param scale scale parameter, strictly positive.
#' @param shape shape parameter.
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log.p,log logical; if \code{FALSE} (default), values are returned on the probability scale.
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   \doi{10.1007/978-1-4471-3675-0_3}
#' @name gpdist

#' @rdname gpdist
#' @export
pgp <- function(q,
                 loc = 0,
                 scale = 1,
                 shape = 0,
                 lower.tail = TRUE,
                 log.p = FALSE) {
  stopifnot(
    "\"loc\" must be a vector of length 1." = length(loc) == 1L,
    "\"loc\" must be finite" = isTRUE(is.finite(loc)),
    "\"scale\" must be a vector of the same length as  \"q\" or a scalar." = length(scale) %in% c(1L, length(q)),
    "\"scale\" must be finite." = isTRUE(all(is.finite(scale))),
    "\"scale\" must be positive." = isTRUE(min(scale, na.rm = TRUE) > 0),
    "\"shape\" must be a vector of length 1" = length(shape) == 1L,
    "\"shape\" must be finite." =   is.finite(shape),
    is.logical(lower.tail),
    length(lower.tail) == 1L,
    is.logical(log.p),
    length(log.p) == 1L
    )
  q <- pmax(q - loc, 0) / scale
  if (shape == 0) {
    p <- 1 - exp(-q)
  } else {
    p <- 1 - exp((-1 / shape) * log1p(pmax(-1, shape * q)))
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  return(p)
}

#' @rdname gpdist
#' @export
dgp <- function(x,
                  loc = 0,
                  scale = 1,
                  shape = 0,
                  log = FALSE) {
  stopifnot(
    "\"loc\" must be a vector of length 1." = length(loc) == 1L,
    "\"loc\" must be finite" = isTRUE(is.finite(loc)),
    "\"scale\" must be a vector of the same length as  \"x\" or a scalar." = length(scale) %in% c(1L, length(x)),
    "\"scale\" must be finite." = isTRUE(all(is.finite(scale))),
    "\"scale\" must be positive." = isTRUE(min(scale, na.rm = TRUE) > 0),
    "\"shape\" must be a vector of length 1" = length(shape) == 1L,
    "\"shape\" must be finite." =   is.finite(shape)
  )

  if (isTRUE(all.equal(shape, 0, check.attributes = FALSE, tolerance = 1e-10))) {
    return(dexp(
      x = x - loc,
      rate = 1 / scale,
      log = log
    ))
  }
  d <- (x - loc) / scale
  index <- (d >= 0 & ((1 + shape * d) >= 0)) | is.na(d)
  d[index] <- log(1 / scale) -
    (1 / shape + 1) * log1p(shape * d[index])
  d[!index & !is.na(d)] <- -Inf
  if (!log) {
    d <- exp(d)
  }
  return(d)
}

#' @rdname gpdist
#' @export
qgp <- function(p,
                 loc = 0,
                 scale = 1,
                 shape = 0,
                 lower.tail = TRUE) {
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1)
    stop("\"p\" must contain probabilities in (0,1)")
  stopifnot(
    "\"loc\" must be a vector of length 1." = length(loc) == 1L,
    "\"loc\" must be finite" = isTRUE(is.finite(loc)),
    "\"scale\" must be a vector of the same length as  \"p\" or a scalar." = length(scale) %in% c(1L, length(p)),
    "\"scale\" must be finite." = isTRUE(all(is.finite(scale))),
    "\"scale\" must be positive." = isTRUE(min(scale, na.rm = TRUE) > 0),
    "\"shape\" must be a vector of length 1" = length(shape) == 1L,
    "\"shape\" must be finite." =   is.finite(shape)
  )

  if (lower.tail) {
    p <- 1 - p
  }
  if (shape == 0) {
    return(loc - scale * log(p))
  } else {
    return(loc + scale * (p ^ (-shape) - 1) / shape)
  }
}

#' @rdname gpdist
#' @export
rgp <- function(n,
                 loc = 0,
                 scale = 1,
                 shape = 0) {
  n <- as.integer(n)
  stopifnot(n > 1, length(n) == 1L)
  stopifnot(
    "\"loc\" must be a vector of length 1." = length(loc) == 1L,
    "\"loc\" must be finite" = isTRUE(is.finite(loc)),
    "\"scale\" must be finite." = isTRUE(all(is.finite(scale))),
    "\"scale\" must be positive." = isTRUE(min(scale, na.rm = TRUE) > 0),
    "\"shape\" must be a vector of length 1" = length(shape) == 1L,
    "\"shape\" must be finite." =   is.finite(shape)
  )
  qgp(p = runif(n), loc = loc, scale = scale, shape = shape)
}



