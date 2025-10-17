#' Extreme U-statistic Pickands estimator
#'
#' Given a random sample of \code{n} exceedances, the estimator
#' returns an estimator of the shape parameter or extreme
#' value index using a kernel of order 3, based on
#' \code{k} largest exceedances of \code{xdat}. Note that the method does not allow for ties.
#'
#' The calculations are based on the recursions provided in Lemma 4.3 of Oorschot et al.
#' @param xdat vector of observations of length \eqn{n}
#' @param k number of largest order statistics \eqn{3 \leq k < n}.
#' @param ... additional arguments for backward compatibility
#' @references Oorschot, J, J. Segers and C. Zhou (2023), Tail inference using extreme U-statistics,  Electronic Journal of Statistics 17(1): 1113-1159. \doi{10.1214/23-EJS2129}
#' @export
#' @examples
#' xdat <- rgp(n = 1000, shape = 0.2)
#' shape.osz(xdat, k = 10)
shape.osz <- function(xdat, k, ...) {
  args <- list(...)
  xdat <- na.omit(as.numeric(xdat))
  xdat <- sort(xdat, decreasing = TRUE)
  k <- sort(as.integer(k))
  n <- length(xdat)
  if (!missing(k)) {
    k <- k[k >= 5 & k < n]
    if (length(k) == 0L) {
      stop("Invalid value for k")
    }
  }
  if (!is.null(args$m) & missing(k)) {
    m <- as.integer(args$m)
    k <- m
  }
  k <- sort(k)
  shape <- numeric(length(k))
  for (i in seq_along(k)) {
    ms <- k[i]
    # Initial recursion j=2
    shape0 <- (2 * (n - 1) / (ms - 2) - 2) * log(xdat[1] - xdat[2])
    mcst <- 1
    for (j in seq.int(
      from = 3L,
      to = n - ms + 3L,
      by = 1L
    )) {
      mcst <- mcst * (n - j - ms + 4L) / (n - j + 1L)
      shape0 <- shape0 +
        mcst *
          (2 * (n - j + 1) / (ms - 2) - j) *
          sum(log(xdat[seq_len(j - 1)] - xdat[j]))
    }
    shape[i] <- shape0 * ms * (ms - 1) * (ms - 2) / (n * (n - 1) * (n - ms + 1))
  }
  if (length(k) == 1L) {
    return(as.numeric(shape))
  } else {
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
#' @param m number of largest order statistics \eqn{3 \leq k \leq n}.
#' @references Oorschot, J, J. Segers and C. Zhou (2023), Tail inference using extreme U-statistics,  \emph{Electronic Journal of Statistics}, 17(1): 1113-1159. \doi{10.1214/23-EJS2129}
#' @export
#' @keywords internal
PickandsXU <- function(xdat, m) {
  .Deprecated(new = "fit.pickandsxu", package = "mev")
  m <- as.integer(m)
  stopifnot(m >= 3, is.vector(xdat), length(xdat) >= m)
  xdat <- na.omit(as.numeric(xdat))
  n <- length(xdat)
  xdat <- sort(xdat, decreasing = TRUE)
  shape <- (2 * (n - 1) / (m - 2) - 2) * log(xdat[1] - xdat[2])
  mcst <- 1
  for (j in seq.int(from = 3L, to = n - m + 3L, by = 1L)) {
    mcst <- mcst * (n - j - m + 4L) / (n - j + 1L)
    shape <- shape +
      mcst *
        (2 * (n - j + 1) / (m - 2) - j) *
        sum(log(xdat[seq_len(j - 1)] - xdat[j]))
  }
  shape <- shape *
    m *
    (m - 1) *
    (m - 2) /
    (n *
      (n - 1) *
      (n -
        m +
        1))
  return(shape)
}

#' Hill's estimator of the shape parameter
#'
#' Given a sample of positive observations, calculate the tail index or
#' shape parameter. The shape estimate returned is positive.
#'
#' @param xdat vector of positive observations
#' @param k vector of order statistics; if missing, a vector going from 10 to sample size minus one.
#' @return a data frame with the number of order statistics \code{k} and the shape parameter estimate \code{shape}, or a single numeric value if \code{k} is a scalar integer.
#' @references Hill, B.M. (1975). \emph{A simple general approach to inference about the tail of a distribution.} Annals of Statistics, \bold{3}, 1163-1173.
#' @export
#' @examples
#' xdat <- mev::rgp(n = 200, loc = 1, scale = 0.5, shape = 0.5)
#' shape.hill(xdat)
shape.hill <- function(xdat, k) {
  logdata <- log(xdat[xdat > 0 & is.finite(xdat)])
  logdata <- sort(logdata, decreasing = TRUE)
  if (missing(k)) {
    k <- seq(10, length(logdata) - 1)
  }
  n <- length(logdata)
  k <- as.integer(sort(k))
  k <- k[k < n]
  kmin <- k[1]
  kmax <- min(k[length(k)], n - 1)
  n <- min(kmax + 1, n)
  ks <- kmin:(n - 1)
  cumlogdat <- cumsum(logdata[1:n])
  shape <- cumlogdat[ks] / ks - logdata[ks + 1]
  if (length(k) > 1) {
    res <- data.frame(
      k = ks,
      shape = shape
    )[ks %in% k, ]
    attr(res, which = "row.names") <- 1:length(k)
    return(res)
  } else {
    return(as.numeric(cumsum(logdata[1:k])[k] / k - logdata[k + 1]))
  }
}

#' Threshold stability plot for Hill estimator
#'
#' @param xdat [vector] sample exceedances
#' @param kmax [int] maximum number of order statistics
#' @param method [string] name of estimator for shape parameter. Default to \code{hill}.
#' @param ... additional arguments passed to \code{fit.shape} for certain methods.
#' @param log [logical] should the x-axis for the number of order statistics used for estimation be displayed on the log scale? Default to \code{TRUE}
#' @return a plot of shape estimates as a function of the number of exceedances
#' @examples
#' xdat <- rgp(n = 250, loc = 1, scale = 2, shape = 0.5)
#' tstab.hill(xdat)
#' @export
tstab.hill <- function(
  xdat,
  kmax,
  method = "hill",
  ...,
  log = TRUE
) {
  xdat <- xdat[is.finite(xdat) & xdat > 0]
  if (missing(kmax)) {
    kmax <- length(xdat)
  } else {
    kmax <- as.integer(kmax[1])
  }
  method <- match.arg(
    method,
    choices = unlist(
      strsplit(
        as.character(
          formals(fit.shape)$method
        ),
        ","
      )
    )[-1],
    several.ok = FALSE
  )
  hill <- fit.shape(xdat = xdat, method = method, k = 10:kmax, ...)
  xp <- hill$k
  hill$lbound <- hill$shape + qnorm(0.025) * abs(hill$shape) / sqrt(hill$k)
  hill$ubound <- hill$shape + qnorm(0.975) * abs(hill$shape) / sqrt(hill$k)
  if (method %in% c("genquant", "osz", "mom", "dekkers", "pickands")) {
    lb <- min(hill$lbound)
  } else {
    lb <- 0
  }
  plot(
    x = hill$k,
    y = hill$shape,
    log = ifelse(isTRUE(log), "x", ""),
    ylim = c(lb, max(hill$ubound, na.rm = TRUE)),
    bty = "l",
    pch = 20,
    panel.first = {
      segments(
        x0 = xp,
        y0 = hill$lbound,
        y1 = hill$ubound,
        col = 'grey'
      )
    },
    xlab = "number of exceedances",
    ylab = "shape"
  )
  return(invisible(NULL))
}

#' Random block maxima shape estimator of Wager
#'
#' Computes the shape estimator for varying k up to sample size of maximum \code{kmax} largest observations
#' @param xdat [vector] sample exceedances
#' @param k [int] vector of integers for which to compute the estimator
#' @param ... additional parameters, currently ignored
#' @return a list with elements
#' \describe{
#'   \item{k}{number of exceedances}
#'   \item{shape}{tail index estimate, strictly positive}
#'   \item{risk}{empirical Bayes estimate of risk}
#'   \item{thresh}{threshold given by the smallest order statistic considered in the sample}
#' }
#' @references Wager, S. (2014). Subsampling extremes: From block maxima to smooth tail estimation, \emph{Journal of Multivariate Analysis}, 130, 335-353, \doi{10.1016/j.jmva.2014.06.010}
#' @export
shape.rbm = function(xdat, k = 10:floor(length(xdat) / 2), ...) {
  k <- as.integer(sort(k))
  xdat <- xdat[is.finite(xdat) & xdat > 0]
  n <- length(xdat)
  kmin <- k[1]
  kmax <- min(k[length(k)] + 1, n)
  ks <- kmin:kmax
  kw <- 2 * n / ks
  logxdat <- log(sort(xdat, decreasing = TRUE))
  M <- numeric(kmax - kmin + 2)
  j <- 0
  for (km in (kmin - 1):kmax) {
    j <- j + 1L
    M[j] <- sum(
      exp(lchoose(n - 1:(n - km + 1), km - 1) - lchoose(n, km)) *
        logxdat[1:(n - km + 1)]
    )
  }
  shape <- diff(M) * kmin:kmax
  empbayes <- c(
    (diff(shape) / diff(log(rev(kw))))^2 + 0.5 * shape[-1]^2 / kw[-1],
    NA
  )
  res <- data.frame(
    k = ks,
    shape = shape,
    risk = empbayes,
    thresh = xdat[ks]
  )[ks %in% k, ]
  class(res) <- c("mev_shape_rbm", "data.frame")
  attr(res, which = "row.names") <- 1:length(k)
  if (length(k) == 1L) {
    return(res$shape)
  } else {
    return(invisible(res))
  }
}

#' Plots for random block maximum estimator
#'
#' The function returns plot of the shape estimator along with the value (and 95\% Wald-based confidence interval) at the selected threshold, or a plot of the empirical Bayes risk.
#'
#' @param x object of class \code{mev_shape_rbm} returned by \code{shape.rbm}
#' @param type [string] type of plot, either \code{"shape"} for the tail index or \code{"risk"} for the empirical Bayes risk
#' @param log [logical] if \code{TRUE} (default), the x-axis for the number of exceedances is displayed on the log scale.
#' @param ... additional arguments, currently ignored
#' @return one or more plots
#' @export
plot.mev_shape_rbm <- function(x, type = c("shape", "risk"), log = TRUE, ...) {
  k0 <- which.min(x$risk)
  type <- match.arg(
    arg = type,
    choices = c("shape", "risk"),
    several.ok = TRUE
  )
  if ("shape" %in% type) {
    plot(
      x = x$k[x$k >= 10],
      y = x$shape[x$k >= 10],
      xlog = isTRUE(log),
      xlab = "log number of exceedances",
      ylab = "shape",
      ylim = c(0, max(x$shape[x$k >= 10])),
      bty = "l",
      type = "l"
    )
    segments(
      x0 = x$k[k0],
      y0 = x$shape[k0] + qnorm(0.025) * x$shape[k0] / sqrt(x$k[k0]),
      y1 = x$shape[k0] + qnorm(0.975) * x$shape[k0] / sqrt(x$k[k0]),
      col = "grey",
    )
    mtext(
      text = eval(bquote(expression(k[0] ~ "=" ~ .(x$k[k0])))),
      side = 3,
      adj = 1,
      line = 0
    )
    rug(x = x$k[k0], ticksize = -0.1, col = "grey50")
    rug(side = 2, x = x$shape[k0], ticksize = -0.1, col = "grey50")
  }
  if ("risk" %in% type) {
    plot(
      x = x$k[x$k >= 10],
      y = log(x$risk[x$k >= 10]),
      log = 'x',
      xlab = "log number of exceedances",
      ylab = "empirical Bayes log risk",
      bty = "l",
      type = "l"
    )
    mtext(
      text = eval(bquote(expression(k[0] ~ "=" ~ .(x$k[k0])))),
      side = 3,
      adj = 1,
      line = 0
    )
    rug(x = x$k[k0])
  }
}

#' Threshold selection for the random block maxima method
#' @inheritParams shape.rbm
#' @param kmax maximum number of exceedances to consider.
#' @return a list with elements
#' @export
thselect.rbm <- function(xdat, kmax = length(xdat)) {
  xdat <- sort(xdat, decreasing = TRUE)
  rbm <- shape.rbm(xdat = xdat, kmax = min(length(xdat) - 1L, kmax))
  ind <- which.min(rbm$risk)
  res <- list(
    k0 = rbm$k[ind],
    thresh0 = rbm$thresh[ind],
    shape = rbm$shape[ind]
  )
  class(res) <- "mev_thselect_rbm"
  return(res)
}

#' @export
print.mev_thselect_rbm <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  cat(
    "Threshold selection method: random block maxima of Wager (2014)\n"
  )
  cat("Selected threshold:", round(x$thresh0, digits), "\n")
  cat("Number of exceedances:", round(x$k0, digits), "\n")
  cat("Shape estimate:", round(x$shape, digits), "\n")
  return(invisible(NULL))
}


#' Beirlant et al. generalized quantile shape estimator
#'
#' This estimator estimates the real shape parameter based on generalized quantile plots based on mean excess functions, generalized median excesses or trimmed mean excesses.
#'
#' @references Beirlant, J., Vynckier P. and J.L. Teugels (1996). \emph{Excess functions and estimation of the extreme-value index.} Bernoulli, 2(\bold{4}), 293-318.
#' @export
#' @param xdat \code{n} vector of observations
#' @param k number of upper order statistics
#' @param p number between zero and one giving the proportion of order statistics for the second threshold
#' @param type string indicating the estimator choice, one of \code{genmean}, \code{genmed} and \code{trimmean}.
#' @param weight weight a kernel function on \eqn{[0,1]}
shape.genquant <- function(
  xdat,
  k,
  type = c("genmean", "genmed", "trimmean"),
  weight,
  p = 0.5
) {
  type <- match.arg(type)
  if (type %in% c("genmed", "trimmean")) {
    stopifnot(isTRUE(all(length(p) == 1L, p > 0, p < 1)))
  }
  xdat <- sort(xdat[is.finite(xdat) & xdat > 0], decreasing = TRUE)
  n <- length(xdat)
  k <- as.integer(sort(k))
  k <- k[k >= 5 & k < n]
  if (length(k) == 0) {
    stop(
      "Invalid number of exceedances: must be five or more, or less than the sample size."
    )
  }
  kmax <- k[length(k)] + 1L
  weight_fun <- FALSE
  if (missing(weight)) {
    we <- 1
  } else {
    if (inherits(weight, "function")) {
      weight_fun <- TRUE
    } else {
      if (length(k) == 1L & is.numeric(weight)) {
        stopifnot(length(weight) %in% c(1, k))
        we <- weight
      } else {
        stop("Invalid \"weight\" function.")
      }
    }
  }
  if (type == "genmean") {
    hill <- shape.hill(k = seq_len(kmax), xdat = xdat)
    ES <- log(xdat[hill$k] * hill$shape)
  } else if (type == "genmed") {
    ES <- log(xdat[floor(p * seq_len(kmax)) + 1] - xdat[seq_len(kmax) + 1])
  } else if (type == "trimmean") {
    ES <- log(
      sapply(seq_len(kmax), FUN = function(ks) {
        mean(xdat[(floor(p * ks) + 1):ks])
      }) -
        xdat[seq_len(kmax) + 1]
    )
  }
  shape <- length(k)
  for (j in seq_along(k)) {
    ks <- k[j]
    if (weight_fun) {
      we <- weight_fun((1:ks) / (ks + 1))
    }
    shape[j] <- sum(we * (log(ks + 1) - log(1:ks)) * (ES[1:ks] - ES[ks + 1])) /
      sum(we * (log(ks + 1) - log(1:ks))^2)
  }
  # Return estimator
  if (length(k) == 1L) {
    return(as.numeric(shape))
  } else {
    return(data.frame(k = k, shape = shape))
  }
}

#' Generalized jackknife shape estimator
#' @inheritParams shape.hill
#' @export
#' @return a data frame with the number of order statistics \code{k} and the shape parameter estimate \code{shape}, or a single numeric value if \code{k} is a scalar integer.
#' @references Gomes, I.M., João Martins, M. and Neves, M.  (2000) Alternatives to a Semi-Parametric Estimator of Parameters of Rare Events-The Jackknife Methodology. \emph{Extremes}, 3, 207–229. \doi{10.1023/A:1011470010228}
shape.genjack <- function(xdat, k) {
  xdat <- xdat[xdat > 0 & is.finite(xdat)]
  n <- length(xdat)
  if (missing(k)) {
    k <- seq(10, n - 1)
  }

  k <- as.integer(sort(k))
  k <- k[k < n]
  # shape <- 2 *
  #   shape.vries(xdat = xdat, k = k)$shape -
  #   shape.hill(xdat = xdat, k = k)$shape
  # if (length(k) > 1) {
  #   return(
  #     data.frame(
  #       k = k,
  #       shape = shape
  #     )
  #   )
  # } else {
  #   return(as.numeric(shape))
  # }

  logdat <- sort(log(xdat), decreasing = TRUE)

  kmin <- k[1]
  kmax <- min(k[length(k)], n - 1)
  cumlogdat <- cumsum(logdat[1:(kmax + 1L)])
  # Hill estimator for k=2, ..., n
  gamma_hill <- cumlogdat[k] / k - logdat[k + 1]
  # Compute square log spacing
  Mn <- vector(mode = "numeric", length = length(gamma_hill))
  for (i in seq_along(Mn)) {
    Mn[i] <- mean((logdat[1:k[i]] - logdat[k[i] + 1])^2)
  }
  # Compute de Vries and generalized jackknife estimators
  gamma_v <- 0.5 * Mn / gamma_hill
  gamma_gj <- 2 * gamma_v - gamma_hill
  if (length(k) > 1) {
    return(
      data.frame(
        k = k,
        shape = gamma_gj
      )
    )
  } else {
    return(as.numeric(gamma_gj))
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
shape.pickands <- function(xdat, k) {
  k <- sort(as.integer(k))
  xdat <- sort(xdat[is.finite(xdat) & xdat > 0], decreasing = FALSE)
  n <- min(k[length(k)] + 1, length(xdat))
  k <- k[k >= 5 & k < n]
  if (length(k) == 0) {
    stop(
      "Invalid number of exceedances: must be five or more, or less than the sample size."
    )
  }
  shape <- numeric(length = length(k))
  for (i in seq_along(k)) {
    shape[i] <- as.numeric(
      (log(xdat[n - floor(k[i] / 4)] - xdat[n - floor(k[i] / 2)]) -
        log(xdat[n - floor(k[i] / 2)] - xdat[n - k[i]]))
    )
  }
  shape <- shape / log(2)
  if (length(k) > 1) {
    return(
      data.frame(k, shape)
    )
  } else {
    return(as.numeric(shape))
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
shape.moment <- function(xdat, k) {
  k <- sort(as.integer(k))
  xdat <- sort(xdat[is.finite(xdat) & xdat > 0], decreasing = TRUE)
  logdata <- as.numeric(log(xdat))
  n <- min(k[length(k)] + 1, length(logdata))
  k <- k[k >= 5 & k < n]
  if (length(k) == 0) {
    stop(
      "Invalid number of exceedances: must be five or more, or less than the sample size."
    )
  }
  cumlogdat <- cumsum(logdata[1:n])
  M1 <- cumlogdat[k] / k - logdata[k + 1]
  shape <- numeric(length = length(k))
  for (i in seq_along(k)) {
    M2 <- mean((logdata[1:k[i]] - logdata[k[i] + 1])^2)
    shape[i] <- M1[i] + 1 - 0.5 / (1 - M1[i]^2 / M2)
  }
  if (length(k) > 1) {
    return(
      data.frame(k, shape)
    )
  } else {
    return(as.numeric(shape))
  }
}

#' de Vries shape estimator
#'
#' Given a sample of exceedances, compute the moment estimator of the positive shape parameter using the ratio of log ratio of exceedance and it's square.
#'
#' @references de Haan, L. and Peng, L. (1998). \emph{Comparison of tail index estimators}, Statistica Neerlandica 52, 60-70.
#' @export
#' @inheritParams shape.pickands
#' @return a data frame with the number of order statistics \code{k} and the shape parameter estimate \code{shape}, or a single numeric value if \code{k} is a scalar.
#' @references de Haan, L. and Peng, L. (1998). Comparison of tail index estimators. \emph{Statistica Neerlandica}, 52: 60-70. \doi{10.1111/1467-9574.00068}
shape.vries <- function(xdat, k) {
  xdat <- sort(
    xdat[is.finite(xdat) & xdat > 0],
    decreasing = TRUE
  )
  n <- length(xdat)
  if (missing(k)) {
    k <- seq(10, n - 1)
  } else {
    k <- sort(as.integer(k))
  }
  logdata <- as.numeric(log(xdat))
  k <- k[k >= 5 & k < n]
  nm <- min(k[length(k)] + 1, n)

  if (length(k) == 0) {
    stop(
      "Invalid number of exceedances: must be five or more, or less than the sample size."
    )
  }
  cumlogdat <- cumsum(logdata[1:nm])
  M1 <- cumlogdat[k] / k - logdata[k + 1]
  shape <- numeric(length = length(k))
  M2 <- sapply(
    seq_along(k),
    function(i) {
      mean((logdata[1:k[i]] - logdata[k[i] + 1])^2)
    }
  )
  shape <- 0.5 * M2 / M1
  if (length(k) > 1) {
    return(
      data.frame(k = k, shape = shape)
    )
  } else {
    return(as.numeric(shape))
  }
}


#' Shape parameter estimates
#'
#' Wrapper to estimate the tail index or shape parameter of an extreme value distribution. Each function has similar sets of arguments, a vector or scalar number of order statistics \code{k} and
#' a vector of positive observations \code{xdat}. The \code{method} argument allows users to choose between different indicators, including the Hill estimator (\code{hill}, for positive observations and shape only), the moment estimator of Dekkers and de Haan (\code{mom} or \code{dekkers}), the de Vries estimator of de Haan and Peng (\code{vries}), the generalized jackknife estimator of Gomes et al. (\code{genjack}), the Beirlant, Vynckier and Teugels generalized quantile estimator (\code{bvt} or \code{genquant}), the Pickands estimator (\code{pickands}), the extreme \eqn{U}-statistics estimator of Oorschot, Segers and Zhou (\code{osz}), or the exponential rgression model of Beirlant et al. (\code{erm}).
#' @export
#' @inheritParams shape.pickands
#' @param method estimation method.
#' @param ... additional parameters passed to functions
#' @return a data frame with the number of order statistics \code{k} and the shape parameter estimate \code{shape}, or a single numeric value if \code{k} is a scalar.
fit.shape <- function(
  xdat,
  k,
  method = c(
    "hill",
    "rbm",
    "osz",
    "vries",
    "genjack",
    "mom",
    "dekkers",
    "genquant",
    "pickands",
    "erm"
  ),
  ...
) {
  method <- match.arg(method)
  if (method == "hill") {
    return(shape.hill(xdat = xdat, k = k))
  } else if (method %in% c("pickandsxu", "osz")) {
    return(shape.osz(xdat = xdat, k = k))
  } else if (method %in% c("mom", "dekkers")) {
    return(shape.moment(xdat = xdat, k = k))
  } else if (method %in% c("beirlant", "genquant")) {
    return(shape.genquant(xdat = xdat, k = k, ...))
  } else if (method == "pickands") {
    return(shape.pickands(xdat = xdat, k = k))
  } else if (method == "vries") {
    return(shape.vries(xdat = xdat, k = k))
  } else if (method == "genjack") {
    return(shape.genjack(xdat = xdat, k = k))
  } else if (method == "rbm") {
    return(shape.rbm(xdat = xdat, kmax = max(k)))
  } else if (method == "erm") {
    res <- shape.erm(xdat = xdat, k = k, ...)
    if (length(k) == 1L) {
      return(res$shape)
    } else {
      return(data.frame(k = res$k, shape = res$shape))
    }
  }
}


#' Estimator of the second order tail index parameter
#'
#'
#' @param xdat vector of positive observations
#' @param k number of highest order statistics to use for estimation
#' @param method string for the estimator
#' @param ... additional arguments passed to individual routinescurrently ignored.
#' @export
#' @examples
#' # Example with rho = -0.2
#' n <- 1000
#' xdat <- mev::rgp(n = n, shape = 0.2)
#' kmin <- floor(n^0.995)
#' kmax <- ceiling(n^0.999)
#' rho_est <- fit.rho(
#'    xdat = xdat,
#'    k = n - kmin:kmax)
#' rho_med <- mean(rho_est$rho)
fit.rho <- function(
  xdat,
  k,
  method = c("fagh", "dk", "ghp", "gbw"),
  ...
) {
  method <- match.arg(method)
  if (method == "fagh") {
    rho.fagh(xdat = xdat, k = k, ...)
  } else if (method == "dk") {
    rho.dk(xdat = xdat, k = k, ...)
  } else if (method == "ghp") {
    rho.ghp(xdat = xdat, k = k, ...)
  } else if (method == "gbw") {
    rho.gbw(xdat = xdat, k = k, ...)
  }
}
#' Second order tail index estimator of Drees and Kaufmann
#'
#' Estimator of the second order regular variation parameter \eqn{\rho \leq 0} parameter for heavy-tailed data proposed by Drees and Kaufmann (1998)
#'
#' @references Drees, H. and E. Kaufmann (1998). Selecting the optimal sample fraction in univariate extreme value estimation, \emph{Stochastic Processes and their Applications}, 75(\bold{2}), 149-172, <doi:10.1016/S0304-4149(98)00017-9>.
#' @param xdat vector of positive observations
#' @param k number of highest order statistics to use for estimation
#' @param tau tuning parameter \eqn{\tau \in (0,1)}
#' @export
rho.dk <- function(xdat, k, tau = 0.5) {
  stopifnot(length(tau) == 1L, tau > 0, tau < 1)
  xdat <- as.numeric(xdat[is.finite(xdat) & xdat > 0])
  logdata <- log(sort(xdat, decreasing = TRUE))
  n <- length(logdata)
  k <- as.integer(sort(k))
  stopifnot(k[1] >= 5, k[length(k)] < n)
  tau <- as.numeric(tau)[1]
  rho <- numeric(length(k))
  for (j in seq_along(k)) {
    ks <- k[j]
    Hl2k <- mean(logdata[1:floor(tau^2 * ks)]) - logdata[floor(tau^2 * ks) + 1]
    Hlk <- mean(logdata[1:floor(tau * ks)]) - logdata[floor(tau * ks) + 1]
    Hk <- mean(logdata[1:ks]) - logdata[ks + 1]
    rho[j] <- min(0, -log(abs((Hl2k - Hlk) / (Hlk - Hk))) / log(tau))
  }
  if (length(k) == 1L) {
    return(as.numeric(rho))
  } else {
    data.frame(k = k, rho = as.numeric(rho))
  }
}

#' Second order tail index estimator of Fraga Alves et al.
#
#' Estimator of the second order regular variation parameter \eqn{\rho \leq 0} parameter for heavy-tailed data proposed by Fraga Alves et al. (2003)
#'
#' @references Fraga Alves, M.I., Gomes, M. Ivette, and de Haan, Laurens (2003). A new class of semi-parametric estimators of the second order parameter. \emph{Portugaliae Mathematica}. Nova Serie 60(\bold{2}), 193-213. <http://eudml.org/doc/50867>.
#' @param xdat vector of positive observations
#' @param k number of highest order statistics to use for estimation
#' @param tau scalar real tuning parameter. Default values is 0, which is typically chosen whenever \eqn{\rho \ge -1}. The choice \eqn{\tau=1} otherwise.
#' @export
#' @examples
#' # Example with rho = -0.2
#' n <- 1000
#' xdat <- mev::rgp(n = n, shape = 0.2)
#' kmin <- floor(n^0.995)
#' kmax <- ceiling(n^0.999)
#' rho_est <- rho.fagh(
#'    xdat = xdat,
#'    k = n - kmin:kmax)
#' rho_med <- mean(rho_est$rho)
rho.fagh <- function(xdat, k, tau = 0) {
  xdat <- as.numeric(xdat[is.finite(xdat) & xdat > 0])
  logdata <- log(sort(xdat, decreasing = TRUE))
  n <- length(logdata)
  k <- as.integer(sort(k))
  stopifnot(k[1] >= 5, k[length(k)] < n)
  tau <- as.numeric(tau)[1]
  rho <- numeric(length(k))
  for (j in seq_along(k)) {
    ks <- k[j]
    mk1 <- mean(logdata[1:ks] - logdata[ks + 1])
    mk2 <- mean((logdata[1:ks] - logdata[ks + 1])^2)
    mk3 <- mean((logdata[1:ks] - logdata[ks + 1])^3)
    if (!isTRUE(all.equal(tau, 0))) {
      w <- (mk1^tau - (0.5 * mk2)^(0.5 * tau)) /
        ((0.5 * mk2)^(0.5 * tau) - (mk3 / 6)^(tau / 3))
    } else {
      w <- (log(mk1) - log(0.5 * mk2) / 2) /
        (log(0.5 * mk2) / 2 - log(mk3 / 6) / 3)
    }
    rho[j] <- min(0, (3 * (w - 1) / (w - 3)))
  }
  if (length(k) == 1L) {
    return(as.numeric(rho))
  } else {
    data.frame(k = k, rho = as.numeric(rho))
  }
}

#' Second order tail index estimator of Gomes et al.
#'
#' Estimator of the second order regular variation parameter \eqn{\rho \leq 0} parameter for heavy-tailed data proposed by Gomes et al. (2003)
#'
#' @references Gomes, M.I., de Haan, L. & Peng, L. (2002). Semi-parametric Estimation of the Second Order Parameter in Statistics of Extremes. \emph{Extremes} 5, 387–414. <doi:10.1023/A:1025128326588>
#'@param xdat vector of positive observations
#'@param k number of highest order statistics to use for estimation
#'@param alpha positive scalar tuning parameter
#'@export
rho.ghp <- function(xdat, k, alpha = 2) {
  xdat <- as.numeric(xdat[is.finite(xdat) & xdat > 0])
  logdata <- log(sort(xdat, decreasing = TRUE))
  n <- length(logdata)
  k <- as.integer(sort(k))
  # alpha must be positive and different from 0.5 and 1
  stopifnot(
    alpha > 0,
    length(alpha) == 1L,
    !isTRUE(all.equal(alpha, 0.5)),
    !isTRUE(all.equal(alpha, 1))
  )
  stopifnot(k[1] >= 5, k[length(k)] < n)
  rho <- numeric(length(k))
  bounds <- sort(c(
    (2 * alpha - 1) / alpha^2,
    4 * (2 * alpha - 1) / (alpha * (alpha + 1)^2)
  ))
  sa_rho <- function(rho, sa, alpha) {
    sa -
      rho^2 *
        (1 -
          (1 - rho)^(2 * alpha) -
          2 * alpha * rho * (1 - rho)^(2 * alpha - 1)) /
        ((1 - (1 - rho)^(alpha + 1) - (alpha + 1) * rho * (1 - rho)^alpha)^2)
  }
  for (j in seq_along(k)) {
    ks <- k[j]
    mk1 <- mean(logdata[1:ks] - logdata[ks + 1])
    mk2 <- mean((logdata[1:ks] - logdata[ks + 1])^2)
    mkap1 <- mean((logdata[1:ks] - logdata[ks + 1])^(alpha + 1))
    mk2a <- mean((logdata[1:ks] - logdata[ks + 1])^(2 * alpha))
    Qa <- function(alpha) {
      mean(
        (logdata[1:ks] - logdata[ks + 1])^alpha - gamma(alpha + 1) * mk1^alpha
      ) /
        (mk2 - 2 * (mk1^2))
    }
    Sa <- exp(
      log(alpha) +
        2 * log(alpha + 1) +
        2 * lgamma(alpha) -
        log(4) -
        lgamma(2 * alpha)
    ) *
      Qa(2 * alpha) /
      (Qa(alpha + 1))^2

    if (Sa > bounds[1] & Sa < bounds[2]) {
      if (!isTRUE(all.equal(alpha, 2))) {
        rootsolve <- try(
          uniroot(
            f = sa_rho,
            lower = -25,
            upper = -1e-8,
            sa = Sa,
            alpha = alpha
          ),
          silent = TRUE
        )
        if (!inherits(rootsolve, "try-error")) {
          rho[j] <- rootsolve$root
        }
      } else {
        rho[j] <- -(2 * (3 * Sa - 2) + sqrt(3 * Sa - 2)) / (3 - 4 * Sa)
      }
    }
  }
  if (length(k) == 1L) {
    return(as.numeric(rho))
  } else {
    data.frame(k = k, rho = as.numeric(rho))
  }
}


#' Second order tail index estimator of Goegebeur et al. (2008)
#'
#' Estimator of the second order regular variation parameter \eqn{\rho \leq 0} parameter for heavy-tailed data based on ratio of kernel goodness-of-fit statistics.
#'
#' @references Goegebeur, Y., J. Beirlant and T. de Wet (2008). Linking Pareto-tail kernel goodness-of-fit statistics with tail index at optimal threshold and second order estimation.  \emph{REVSTAT-Statistical Journal}, 6(\bold{1}), 51-69. <doi:10.57805/revstat.v6i1.57>
#'@param xdat vector of positive observations
#'@param k number of highest order statistics to use for estimation
#'@export
rho.gbw <- function(
  xdat,
  k
) {
  k <- sort(as.integer(k))
  kernel_jack <- function(u) {
    -1 - log(u)
  }
  kernel_lewis <- function(u) {
    u - 0.5
  }
  xdat <- sort(xdat, decreasing = TRUE)
  kmax <- max(k)
  stopifnot(xdat[kmax + 1] > 0)
  Z <- (1:kmax) * diff(log(xdat[1:(kmax + 1)]))
  rho <- numeric(length(k))
  for (j in seq_along(k)) {
    T1 <- mean(kernel_jack((1:k[j]) / (k[j] + 1)) * Z[1:k[j]])
    T2 <- mean(kernel_lewis((1:k[j]) / (k[j] + 1)) * Z[1:k[j]])
    rho[j] <- (4 * T2 + T1) / (2 * T2 + T1)
  }
  rho <- pmin(0, rho)
  if (length(k) == 1L) {
    return(as.numeric(rho))
  } else {
    data.frame(k = k, rho = as.numeric(rho))
  }
}

## Ciuperca and Mercadier estimator of second-order regular variation
# rho.cm <- function(xdat, k) {
#   xdat <- as.numeric(xdat[is.finite(xdat)])
#   k <- sort(k)
#   kmax <- k[length(k)]
#   stopifnot(xdat[kmax] > 0)
#   logdata <- log(sort(xdat, decreasing = TRUE))
#   n <- length(logdata)
#   k <- as.integer(sort(k))
# }

#' Weissman's quantile estimator
#'
#' Given a small probability of exceedance \code{p},
#' the number of exceedances \code{k} out of \code{n} observation
#' above the threshold \eqn{u} (\code{thresh}) (corresponding typically to the (\eqn{k+1})th order statistic, compute the tail quantile at level \eqn{Q(1-p)} using the estimator of Weissman (1978) under the assumption of Pareto tail (positive shape \eqn{\xi}), viz.
#' \deqn{ Q(1-p) = u \left(\frac{k}{pn}\right)^{\xi}.}
#'
#' @param p tail probability, must be larger than the proportion of exceedances \code{k/n}.
#' @param k vector of the number of exceedances above \code{thresh}
#' @param n integer, total sample size
#' @param thresh vector of thresholds
#' @param shape vector of positive shape parameters
#' @references Weissman, I. (1978). Estimation of Parameters and Larger Quantiles Based on the \emph{k} Largest Observations. \emph{Journal of the American Statistical Association}, 73(\bold{364}), 812–815. <doi:10.2307/2286285>.
#' @return a vector of tail quantiles
#' @export
#' @examples
#' set.seed(2025)
#' p <- 1/100
#' xdat <- rgp(n = 1000, loc = 2, scale = 2, shape = 0.4)
#' hill <- shape.hill(xdat, k = seq(20L, 100L, by = 10L))
#' thresh <- sort(xdat, decreasing = TRUE)[hill$k+1]
#' qweissman(
#'    p = 1/100,
#'    k = hill$k,
#'    n = length(xdat),
#'    thresh = thresh,
#'    shape = hill$shape)
#' # Compare with true quantile
#' qgp(1/100, loc = 2, scale = 2, shape = 0.4, lower.tail = FALSE)
qweissman <- function(p, k, n, thresh, shape) {
  m <- length(k)
  stopifnot(
    "Vectors \"k\" and \"shape\" must be of the same length." = (m ==
      length(shape)),
    "Vectors \"k\" and \"thresh\" must be of the same length." = (m ==
      length(thresh)),
    "Sample size must be a scalar." = (length(n) == 1L),
    "Number of exceedances \"k\" must be less than sample size \"n\"." = isTRUE(all(
      k < n
    )),
    "Tail probability must be a scalar." = (length(p) == 1L) | (m == 1L),
    "Tail probability must be smaller than fraction of exceedances." = isTRUE(all(
      p < k / n
    )),
    "Shape parameter must be strictly positive." = isTRUE(all(shape > 0))
  )
  thresh * (k / n / p)^shape
}
