#' Coefficient of extremal asymmetry
#'
#' This function implements estimators of the bivariate
#' coefficient of extremal asymmetry proposed in Semadeni's (2021) PhD thesis.
#'  Two estimators are implemented: one based on empirical distributions, the second using empirical likelihood.
#' @details Let \code{U}, \code{V} be uniform random variables and define the partial extremal dependence coefficients
#' \deqn{\varphi_{+}(u) = \Pr(V > U \mid U > u, V > u),},
#' \deqn{\varphi_{-}(u) = \Pr(V < U \mid U > u, V > u),}
#' \deqn{\varphi_0(u) = \Pr(V = U \mid U > u, V > u).}
#' Define
#' \deqn{ \varphi(u) = \frac{\varphi_{+} - \varphi_{-}}{\varphi_{+} + \varphi_{-}}}
#' and the coefficient of extremal asymmetry as \eqn{\varphi = \lim_{u \to 1} \varphi(u)}.
#'
#' The empirical likelihood estimator, derived for max-stable vectors with unit Frechet margins, is
#' \deqn{\widehat{\varphi}_{\mathrm{el}} = \frac{\sum_i p_i \mathrm{I}(w_i \leq 0.5) - 0.5}{0.5 - 2\sum_i p_i(0.5-w_i) \mathrm{I}(w_i \leq 0.5)}}
#' where \eqn{p_i} is the empirical likelihood weight for observation \eqn{i}, \eqn{\mathrm{I}} is an indicator function and \eqn{w_i} is the pseudo-angle associated to the first coordinate, derived based on exceedances above \eqn{u}.
#' @param xdat an \code{n} by 2 matrix of observations
#' @param qlev vector of quantile levels at which to evaluate extremal asymmetry
#' @param nq integer; number of quantiles at which to evaluate the coefficient if \code{u} is \code{NULL}
#' @param qlim a vector of length 2 with the probability limits for the quantiles
#' @param estimator string indicating the estimation method, one of \code{emp} or empirical likelihood (\code{elik})
#' @param confint string for the method used to derive confidence intervals, either \code{none} (default) or a nonparametric \code{bootstrap}
#' @param level probability level for confidence intervals, default to 0.95 or bounds for the interval
#' @param B integer; number of bootstrap replicates (if applicable)
#' @param ties.method string; method for handling ties. See the documentation of \link[base]{rank} for available options.
#' @param plot logical; if \code{TRUE}, return a plot.
#' @param ... additional arguments for backward compatibility
#' @return an invisible data frame with columns
#' \describe{
#' \item{\code{qlev}}{quantile level of thresholds}
#' \item{\code{coef}}{extremal asymmetry coefficient estimates}
#' \item{\code{lower}}{either \code{NULL} or a vector containing the lower bound of the confidence interval}
#' \item{\code{upper}}{either \code{NULL} or a vector containing the lower bound of the confidence interval}
#' }
#' @references Semadeni, C. (2020). Inference on the Angular Distribution of Extremes, PhD thesis, EPFL, no. 8168.
#' @examples
#' \dontrun{
#' samp <- rmev(n = 1000,
#'              d = 2,
#'              param = 0.2,
#'              model = "log")
#' xdep.asym(samp, confint = "wald")
#' xdep.asym(samp, method = "emplik", confint = "none")
#' }
#' @export
xdep.asym <- function(
  xdat,
  qlev = NULL,
  nq = 40,
  qlim = c(0.8, 0.99),
  estimator = c("emp", "elik"),
  confint = c("none", "wald", "bootstrap"),
  level = 0.95,
  B = 999L,
  ties.method = "random",
  plot = TRUE,
  ...
) {
  # TODO support the multivariate case?
  stopifnot(
    "Coefficient of extremal asymmetry only defined for bivariate inputs." = ncol(
      xdat
    ) ==
      2L
  )
  estimator <- match.arg(estimator)
  xasym(
    xdat = xdat,
    qlev = qlev,
    nq = nq,
    qlim = qlim,
    estimator = switch(estimator, emp = "empirical", elik = "emplik"),
    confint = match.arg(confint),
    level = level,
    B = B,
    ties.method = ties.method,
    plot = plot,
    ...
  )
}

#' Coefficient of tail correlation
#'
#' The coefficient of tail correlation \eqn{\chi} is
#' \deqn{\chi = \lim_{u \to 1} \frac{\Pr(F_1(X_1)>u, \ldots, F_D(X_D)>u)}{1-u}.}
#' Asymptotically independent vectors have \eqn{\chi = 0}. The estimator uses an estimator of the survival copula
#' @export
#' @param xdat an \eqn{n} by \eqn{d} matrix of multivariate observations
#' @param qlev vector of percentiles between 0 and 1
#' @param nq number of quantiles of the structural variable at which to form a grid; only used if \code{u = NULL}.
#' @param qlim limits for the sequence \code{u} of the structural variable
#' @param confint string indicating the type of confidence interval, one of \code{"wald"} or \code{"lrt"}
#' @param level the confidence level required (default to 0.95).
#' @param ties.method string indicating the type of method for \code{rank}; see \code{\link[base]{rank}} for a list of options. Default to \code{"random"}
#' @param margtrans string giving the marginal transformation, one of \code{emp} for rank-based transformation or \code{none} if data are already on the uniform scale
#' @param estimator string giving estimator to employ
#' @param plot logical; if \code{TRUE}, return a plot
#' @param ... additional arguments to \code{taildep}, currently ignored
#' @return a data frame
#' \itemize{
#' \item \code{qlev}: quantile level of estimates
#' \item \code{coef}: point estimates
#' \item \code{lower}: lower bound of confidence interval
#' \item \code{upper}: lower bound of confidence interval
#' }
#' @examples
#' \dontrun{
#' set.seed(765)
#' # Max-stable model
#' dat <- rmev(n = 1000, d = 2, param = 0.7, model = "log")
#' xdep.chi(dat, confint = 'wald')
#' }
xdep.chi <- function(
  xdat,
  qlev = NULL,
  nq = 40,
  qlim = c(0.8, 0.99),
  estimator = c("emp", "betacop", "gpd", "hill"),
  confint = c("wald", "lrt"),
  level = 0.95,
  margtrans = c("emp", "none"),
  ties.method = "random",
  plot = TRUE,
  ...
) {
  out <- taildep(
    xdat = xdat,
    qlev = qlev,
    nq = nq,
    qlim = qlim,
    depmeas = "chi",
    estimator = list("chi" = match.arg(estimator)),
    confint = match.arg(confint),
    level = level,
    trunc = TRUE,
    margtrans = match.arg(margtrans),
    ties.method = match.arg(ties.method),
    plot = FALSE,
    ...
  )
  res <- data.frame(
    qlev = out$qlev,
    coef = out$chi[, "coef"],
    lower = out$chi[, "lowerci"],
    upper = out$chi[, "upperci"]
  )
  attr(res, "estimator") <- out$chi_method
  attr(res, "confint") <- out$chi_confint_method
  attr(res, "measure") <- "chi"
  class(res) <- c("mev_xdep", "data.frame")
  if (isTRUE(plot)) {
    plot(res)
  }
  return(invisible(res))
}


#' Coefficient of tail dependence
#'
#' For data with unit Pareto margins, the coefficient of tail dependence \eqn{\eta} is defined  via \deqn{\Pr(\min(X) > x) = L(x)x^{-1/\eta},}
#' where \eqn{L(x)} is a slowly varying function. Ignoring the latter, several estimators of \eqn{\eta} can be defined. In unit Pareto margins, \eqn{\eta} is a nonnegative shape parameter that can be estimated by fitting a generalized Pareto distribution above a high threshold. In exponential margins, \eqn{\eta} is a scale parameter and the maximum likelihood estimator of the latter is the Hill estimator. Both methods are based on peaks-over-threshold and the user can choose between pointwise confidence obtained through a likelihood ratio test statistic (\code{"lrt"}) or the Wald statistic (\code{"wald"}).
#'
#' The most common approach for estimation is the empirical survival copula, by evaluating the proportion of sample minima with uniform margins that exceed a given \eqn{x}. An alternative estimator uses a smoothed estimator of the survival copula using Bernstein polynomial, resulting in the so-called \code{betacop} estimator. Approximate pointwise confidence intervals for the latter are obtained by assuming the proportion of points is binomial.
#' @export
#' @param xdat an \eqn{n} by \eqn{d} matrix of multivariate observations
#' @param qlev vector of percentiles between 0 and 1
#' @param nq number of quantiles of the structural variable at which to form a grid; only used if \code{u = NULL}.
#' @param qlim limits for the sequence \code{u} of the structural variable
#' @param confint string indicating the type of confidence interval, one of \code{"wald"} or \code{"lrt"}
#' @param level the confidence level required (default to 0.95).
#' @param ties.method string indicating the type of method for \code{rank}; see \code{\link[base]{rank}} for a list of options. Default to \code{"random"}
#' @param margtrans string giving the marginal transformation, one of \code{emp} for rank-based transformation or \code{none} if data are already on the uniform scale
#' @param estimator string giving estimator to employ
#' @param plot logical; if \code{TRUE}, return a plot
#' @param mqlev marginal quantile levels for semiparametric estimation for estimator \code{kj}; data above this are modelled using a generalized Pareto distribution. If missing, empirical estimation is used throughout
#' @param ... additional arguments to \code{taildep}, currently ignored
#' @return a data frame
#' \itemize{
#' \item \code{qlev}: quantile level of estimates
#' \item \code{coef}: point estimates
#' \item \code{lower}: lower bound of confidence interval
#' \item \code{upper}: lower bound of confidence interval
#' }
#' @references Ledford, A.W. and J. A. Tawn (1996), Statistics for near independence in multivariate extreme values. \emph{Biometrika}, \bold{83}(1), 169--187.
#' @examples
#' \dontrun{
#' set.seed(765)
#' # Max-stable model
#' dat <- rmev(n = 1000, d = 2, param = 0.7, model = "log")
#' xdep.eta(dat, confint = 'wald')
#' }
xdep.eta <- function(
  xdat,
  qlev = NULL,
  nq = 40,
  qlim = c(0.8, 0.99),
  estimator = c("emp", "betacop", "gpd", "hill", "kj"),
  confint = c("wald", "lrt"),
  level = 0.95,
  margtrans = c("emp", "sp", "none"),
  ties.method = "random",
  plot = TRUE,
  mqlev = NULL,
  ...
) {
  if (estimator != "kj") {
    out <- taildep(
      xdat = xdat,
      qlev = qlev,
      nq = nq,
      qlim = qlim,
      depmeas = "eta",
      estimator = list(eta = match.arg(estimator)),
      confint = match.arg(confint),
      level = level,
      trunc = TRUE,
      margtrans = match.arg(margtrans),
      ties.method = ties.method,
      plot = FALSE,
      ...
    )
    res <- data.frame(
      qlev = out$qlev,
      coef = out$eta[, "coef"],
      lower = out$eta[, "lowerci"],
      upper = out$eta[, "upperci"]
    )
    attr(res, "estimator") <- out$eta_method
    attr(res, "confint") <- out$eta_confint_method
    attr(res, "measure") <- "eta"
  } else {
    if (is.null(qlev)) {
      if (length(qlim) != 2L) {
        stop("\"qlim\" must be a bivariate numeric vector.")
      }
      if (qlim[1] < 0 || qlim[2] >= 1) {
        stop("\"qlim\" must contain probabilities.")
      }
      qlev <- seq(qlim[1], qlim[2], length = nq)
    }
    out <- kjtail(
      xdat = xdat,
      qlev = qlev,
      ptail = NULL,
      mqu = mqlev,
      ties.method = ties.method
    )
    res <- data.frame(
      qlev = out$p,
      coef = out$eta[, "coef"],
      lower = pmax(
        0,
        out$eta[, "coef"] + qnorm(0.5 - level / 2) * out$eta[, "sd"]
      ),
      upper = pmax(
        0,
        out$eta[, "coef"] + qnorm(0.5 + level / 2) * out$eta[, "sd"]
      )
    )
    attr(res, "estimator") <- "kj"
    attr(res, "confint") <- "wald"
    attr(res, "measure") <- "eta"
  }
  class(res) <- c("mev_xdep", "data.frame")
  if (isTRUE(plot)) {
    plot(res)
  }
  return(invisible(res))
}

#' Coefficient chi-bar
#'
#' For data with unit Pareto margins, the coefficient \eqn{\bar{\chi} = 2\eta-1} is defined  via \deqn{\Pr(\min(X) > x) = L(x)x^{-1/\eta},}
#' where \eqn{L(x)} is a slowly varying function. Ignoring the latter, several estimators of \eqn{\eta} can be defined. In unit Pareto margins, \eqn{\eta} is a nonnegative shape parameter that can be estimated by fitting a generalized Pareto distribution above a high threshold. In exponential margins, \eqn{\eta} is a scale parameter and the maximum likelihood estimator of the latter is the Hill estimator. Both methods are based on peaks-over-threshold and the user can choose between pointwise confidence obtained through a likelihood ratio test statistic (\code{"lrt"}) or the Wald statistic (\code{"wald"}).
#'
#' The most common approach for estimation is the empirical survival copula, by evaluating the proportion of sample minima with uniform margins that exceed a given \eqn{x}. An alternative estimator uses a smoothed estimator of the survival copula using Bernstein polynomial, resulting in the so-called \code{betacop} estimator. Approximate pointwise confidence intervals for the latter are obtained by assuming the proportion of points is binomial.
#' @export
#' @param xdat an \eqn{n} by \eqn{d} matrix of multivariate observations
#' @param qlev vector of percentiles between 0 and 1
#' @param nq number of quantiles of the structural variable at which to form a grid; only used if \code{u = NULL}.
#' @param qlim limits for the sequence \code{u} of the structural variable
#' @param confint string indicating the type of confidence interval, one of \code{"wald"} or \code{"lrt"}
#' @param level the confidence level required (default to 0.95).
#' @param ties.method string indicating the type of method for \code{rank}; see \code{\link[base]{rank}} for a list of options. Default to \code{"random"}
#' @param margtrans string giving the marginal transformation, one of \code{emp} for rank-based transformation or \code{none} if data are already on the uniform scale
#' @param estimator string giving estimator to employ
#' @param plot logical; if \code{TRUE}, return a plot
#' @param ... additional arguments to \code{taildep}, currently ignored
#' @return a data frame
#' \itemize{
#' \item \code{qlev}: quantile level of estimates
#' \item \code{coef}: point estimates
#' \item \code{lower}: lower bound of confidence interval
#' \item \code{upper}: lower bound of confidence interval
#' }
#' @references Ledford, A.W. and J. A. Tawn (1996), Statistics for near independence in multivariate extreme values. \emph{Biometrika}, \bold{83}(1), 169--187.
#' @examples
#' \dontrun{
#' set.seed(765)
#' # Max-stable model
#' dat <- rmev(n = 1000, d = 2, param = 0.7, model = "log")
#' xdep.chibar(dat, confint = 'wald')
#' }
#' @references Ledford, A.W. and J. A. Tawn (1996), Statistics for near independence in multivariate extreme values. \emph{Biometrika}, \bold{83}(1), 169--187.
xdep.chibar <- function(
  xdat,
  qlev = NULL,
  nq = 40,
  qlim = c(0.8, 0.99),
  estimator = c("emp", "betacop"),
  confint = c("wald", "lrt"),
  level = 0.95,
  margtrans = c("emp", "none"),
  ties.method = "random",
  plot = TRUE,
  ...
) {
  stopifnot(
    "Measure chi-bar only defined for bivariate inputs." = ncol(xdat) == 2L
  )
  out <- taildep(
    xdat = xdat,
    qlev = qlev,
    nq = nq,
    qlim = qlim,
    depmeas = "eta",
    estimator = list(eta = match.arg(estimator)),
    confint = match.arg(confint),
    level = level,
    trunc = TRUE,
    margtrans = match.arg(margtrans),
    ties.method = ties.method,
    plot = FALSE,
    ...
  )
  res <- data.frame(
    qlev = out$qlev,
    coef = 2 * out$eta[, "coef"] - 1,
    lower = 2 * out$eta[, "lowerci"] - 1,
    upper = 2 * out$eta[, "upperci"] - 1
  )
  attr(res, "estimator") <- out$eta_method
  attr(res, "confint") <- out$eta_confint_method
  attr(res, "measure") <- "chibar"
  class(res) <- c("mev_xdep", "data.frame")
  if (isTRUE(plot)) {
    plot(res)
  }
  return(invisible(res))
}


#' @export
plot.mev_xdep <- function(x, ...) {
  measure <- attributes(x)$measure
  ellips <- list(...)
  if (is.null(ellips$bty)) {
    ellips$bty <- 'l'
  }
  if (is.null(ellips$xlab)) {
    ellips$xlab <- "quantile level"
  }
  if (is.null(ellips$ylab)) {
    ellips$ylab <- switch(
      measure,
      xasym = expression("extremal asymmetry" ~ phi),
      chi = expression("tail correlation" ~ chi),
      eta = expression("tail dependence" ~ eta),
      chibar = expression(bar(chi)),
      xindex = "extremal index"
    )
  }
  if (is.null(ellips$pch)) {
    ellips$pch <- 20
  }
  if (is.null(ellips$ylim)) {
    ellips$ylim <- switch(
      measure,
      xasym = c(-1, 1),
      chi = c(0, 1),
      eta = c(0, 1),
      chibar = c(-1, 1),
      xindex = c(0, 1)
    )
  }
  if (is.null(ellips$xlim)) {
    ellips$xlim <- range(x$qlev)
    ellips$xlim <- pmax(0, pmin(ellips$xlim, 1))
  }
  if (is.null(ellips$yaxs)) {
    ellips$yaxs <- "i"
  }
  ellips$x <- x$qlev
  ellips$y <- x$coef
  do.call("plot", ellips)
  if (attr(x, "confint") != "none") {
    points(x$qlev, x$lower, pch = 95, col = "grey")
    points(x$qlev, x$upper, pch = 95, col = "grey")
  }
  return(invisible(NULL))
}
