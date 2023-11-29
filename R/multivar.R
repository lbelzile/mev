#' Parametric estimates of \eqn{\bar{\chi}}{chi bar}
#'
#' The function fits a generalized Pareto distribution to minima of Pareto variates,
#' using the representation \deqn{\Pr(\min(X) > x) = \frac{L(x)}{x^{1/\eta}},}
#' where \eqn{\bar{\chi}=2\eta-1}. The data are transformed to the unit Pareto scale and
#' a generalized Pareto variable is fitted to the minimum. The parameter \eqn{\eta} corresponds to the shape of the latter.
#' The confidence intervals can be based either on the delta-method, a profile likelihood or a tangent exponential model approximation.
#' @export
#' @param dat an \eqn{n} by \eqn{d} matrix of multivariate observations
#' @param qu percentile level at which to threshold. Default to all observations.
#' @param confint string indicating the type of confidence interval.
#' @param level the confidence level required
#' @return a named vector of length 3 containing the point estimate, the lower and the upper confidence intervals
#' @seealso \code{\link[evd]{chiplot}} for empirical estimates of \eqn{\chi}{chi} and \eqn{\bar{\chi}}{chibar}.
#' @keywords internal
#' @examples
#' \dontrun{
#' set.seed(765)
#' # Max-stable model, chibar = 1
#' dat <- rmev(n = 1000, model = "log", d = 2, param = 0.5)
#' chibar(dat, 'profile', qu = 0.5)
#' s <- seq(0.05,1, length = 30)
#' chibar_est <- t(sapply(s, function(keep){chibar(dat, 'delta', qu = keep)}))
#' matplot(s, chibar_est, type = 'l', col = c(1, 2, 2),  lty = c(1, 2, 2),
#'  ylab = expression(bar(chi)), xlab = 'p')
#' abline(h = 1, lty = 3, col = 'grey')
#' # Multivariate normal sample, chibar = 0 - strong asymptotic independence at penultimate level
#' dat <- mvrnorm(n = 1000, mu = c(0, 0), Sigma = cbind(c(1, 0.75), c(0.75, 1)))
#' chibar(dat, 'tem', q = 0.1)
#' chibar_est <- t(sapply(s, function(keep){chibar(dat, 'profile', qu = keep)}))
#' matplot(s, chibar_est, type = 'l', col = c(1, 2, 2),  lty = c(1, 2, 2),
#'  ylab = expression(bar(chi)), xlab = 'p')
#' abline(h = 1, lty = 3, col = 'grey')
#' }
chibar <- function(dat, confint = c("delta", "profile", "tem"), qu = 0, level = 0.95) {
  base::.Deprecated(new = ".chibar", package = "mev", msg = "The 'chibar' function is being depreceated. Use 'taildep' with 'depmeas = eta' to obtain estimates of 'chibar = 2*eta-1'")
    if (ncol(as.matrix(dat)) < 2) {
        stop("The method is valid for multivariate data only.")
    }
    confint <- match.arg(arg = confint)
    # Transform variables to standard Pareto margin
    sp <- apply(dat, 2, function(x) {
        1/(1 - rank(x, na.last = "keep", ties.method = "average")/(length(x) + 1))
    })
    sp <- apply(sp, 1, min) - 1
    qu <- quantile(sp, 1 - min(max(0, qu[1]), 1))
    sp <- sp[sp > qu] - qu
    if ("delta" == confint) {
        gpfit_min_par <- gp.fit(sp, threshold = 0)
        chibar_est <- as.vector((2 * gpfit_min_par$est[2] - 1))
        return(c(Estimate = chibar_est, `Lower CI` = chibar_est - qnorm(1 - (1 - level)/2) * 2 * as.vector(gpfit_min_par$std.err[2]),
            `Upper CI` = chibar_est + qnorm(1 - (1 - level)/2) * 2 * as.vector(gpfit_min_par$std.err[2])))
    } else {
        confint_profile <- 2 * confint(gpd.pll(param = "shape", psi = NA, dat = sp, mod = confint, plot = FALSE), level = level, print = FALSE) - 1
        return(switch(confint, profile = confint_profile, tem = confint_profile[,2]))
    }
}


#' Bivariate angular function for extrapolation based on rays
#'
#' The scale parameter \eqn{g(w)} in the Ledford and Tawn approach is estimated empirically for
#' \eqn{x} large as \deqn{\frac{\Pr(X_P>xw, Y_P>x(1-w))}{\Pr(X_P>x, Y_P>x)}}
#' where the sample (\eqn{X_P, Y_P}) are observations on a common unit Pareto scale.
#'  The coefficient \eqn{\eta} is estimated using maximum likelihood as the
#'  shape parameter of a generalized Pareto distribution on \eqn{\min(X_P, Y_P)}.
#'
#' @param dat an \eqn{n} by \eqn{2} matrix of multivariate observations
#' @param qu quantile level on uniform scale at which to threshold data. Default to 0.95
#' @param w vector of unique angles between 0 and 1 at which to evaluate scale empirically.
#' @return a list with elements
#' \itemize{
#' \item \code{w}: angles between zero and one
#' \item \code{g}: scale function at a given value of \code{w}
#' \item \code{eta}: Ledford and Tawn tail dependence coefficient
#' }
#' @export
#' @references Ledford, A.W. and J. A. Tawn (1996), Statistics for near independence in multivariate extreme values. \emph{Biometrika}, \bold{83}(1), 169--187.
#' @examples
#' angextrapo(rmev(n = 1000, model = 'log', d = 2, param = 0.5))
angextrapo <- function(dat, qu = 0.95, w = seq(0.05, 0.95, length = 20)) {
    if (ncol(dat) != 2) {
        stop("Only implemented in the bivariate case")
    }
    sp <- apply(dat, 2, function(x) {
        1/(1 - rank(x, na.last = "keep", ties.method = "average")/(length(x) + 1))
    })
    # Estimate of eta at w = 1/2
    x <- 1/(1 - qu)
    eta <- fit.gpd(apply(sp, 1, min), threshold = x)$estimate["shape"]
    # Angles
    if (any(c(w < 0, w > 1, length(unique(w)) != length(w)))) {
        stop("Invalid argument \"w\" to angextrapo")
    }
    g <- sapply(w, function(wi) {
        sum((sp[, 1] > wi * x) + (sp[, 2] > (1 - wi) * x))
    })/sum(rowSums(sp > x) == 2)
    # Return angles, empirical estimates of g(w) and eta coefficient
    return(list(w = w, g = g, eta = eta))
}

##################
#' Estimation of the bivariate lambda function of Wadsworth and Tawn (2013)
#'
#' @param dat an \eqn{n} by \eqn{2} matrix of multivariate observations
#' @param qu quantile level on uniform scale at which to threshold data. Default to 0.95
#' @param method string indicating the estimation method
#' @param plot logical indicating whether to return the graph of \code{lambda}
#'
#'
#' The confidence intervals are based on normal quantiles. The standard errors for the \code{hill}
#' are based on the asymptotic covariance and that of the \code{mle} derived using the delta-method.
#' Bayesian posterior predictive interval estimates are obtained using ratio-of-uniform sampling with flat priors:
#' the shape parameters are constrained to lie within the triangle, as are frequentist point estimates
#' which are adjusted post-inference.
#'
#' @importFrom utils capture.output
#' @importFrom graphics segments
#' @return a plot of the lambda function if \code{plot=TRUE}, plus an invisible list with components
#' \itemize{
#' \item \code{w} the sequence of angles in (0,1) at which the \code{lambda} values are evaluated
#' \item \code{lambda} point estimates of lambda
#' \item \code{lower.confint} 95% confidence interval for lambda (lower bound)
#' \item \code{upper.confint} 95% confidence interval for lambda (upper bound)
#' }
#' @examples
#' set.seed(12)
#' dat <- mev::rmev(n = 1000, d = 2, model = "log", param = 0.1)
#' lambdadep(dat, method = 'hill')
#' \dontrun{
#' lambdadep(dat, method = 'bayes')
#' lambdadep(dat, method = 'mle')
#' # With independent observations
#' dat <- matrix(runif(n = 2000), ncol = 2)
#' lambdadep(dat, method = 'hill')
#' }
#' @export
lambdadep <- function(dat, qu = 0.95, method = c("hill", "mle", "bayes"), plot = TRUE) {
    ## Hill estimator for fixed kth order statistic
    hill_thresh <- function(dat, qu = 0.95, thresh = quantile(dat, qu)) {
        dat <- as.numeric(dat)
        excess <- dat[dat > thresh]
        1/(mean(log(dat[dat > thresh])) - log(thresh))
    }
    if (method == "bayes") {
      if (!requireNamespace("revdbayes", quietly = TRUE)) {
        stop("Package \"revdbayes\" needed for this function to work. Please install it.",
             call. = FALSE)
      }
    }
    # Transform variables to the exponential scale
    Xexp <- t(apply(dat, 2, function(x) {
        -log(1 - rank(x, na.last = "keep", ties.method = "average")/(length(na.omit(x)) + 1))
    }))
    # Form a bivariate minima pair for a grid of values of w in Sd, the unit simplex
    v <- rep(1, 2)
    w_seq <- seq(0, 1, by = 0.02)
    lambda_seq <- sapply(w_seq, function(w) {
        ang_weighted_dat <- exp(apply(Xexp/(c(w, 1 - w)), 2, min))
        if (method == "mle") {
            fit <- mev::gp.fit(ang_weighted_dat, thresh = quantile(ang_weighted_dat, qu))
            return(c(1/fit$estimate["shape"], fit$std.err["shape"]/(fit$estimate["shape"]^2)))
        } else if (method == "hill") {
            hillest <- hill_thresh(dat = ang_weighted_dat, qu = qu)
            return(c(hillest, hillest/sqrt((1 - qu) * nrow(dat))))
        } else if (method == "bayes") {
            # If at endpoint, posterior shape is degenerate
            if (w == 0 || w == 1) {
                return(rep(1, 3))
            } else {
                # Try fitting a GP model, mode of posterior should not be too far
                pot_stval0 <- mev::gp.fit(ang_weighted_dat, thresh = quantile(ang_weighted_dat, qu))$estimate
                if (pot_stval0[2] < 1 || pot_stval0[2] > 1/max(c(1 - w, w))) {
                  # If values are not within the allowed interval, fit GP fixing the shape to a legit value
                  xi <- max(
                    1 + 0.001,
                    min(v[2], 1 / max(c(1 - w, w)) - 0.001,
                        na.rm = TRUE)
                  )
                  pot_stval <- try(
                    mev::fit.gpd(xdat = ang_weighted_dat,
                                 threshold = quantile(ang_weighted_dat,
                                                      qu),
                                 fpar = list(shape = xi)
                                 ))
                  # Make sure that the result is valid and optim converged
                  if (!inherits(pot_stval, what = "try-error")) {
                    start <- c(pot_stval$estimate[[1]], xi)
                  } else {
                    start <- NA
                  }
                } else {
                  start <- as.vector(pot_stval0)
                }
                # Generate independent samples from the posterior Catch and sink error messages, print statements and warnings - invalid input is
                # removed anyway and cast to NA if needs be
                invisible(utils::capture.output(postsamp <- try(suppressWarnings(revdbayes::rpost(n = 300, model = "gp", data = ang_weighted_dat,
                  thresh = quantile(ang_weighted_dat, qu), prior = revdbayes::set_prior(prior = "flat", model = "gp", min_xi = 1,
                    max_xi = 1/max(c(1 - w, w))), init_ests = start, trans = "BC")), silent = TRUE)))
                if (inherits(postsamp, what = "try-error")) {
                  # Try again if it failed, with different values
                  invisible(capture.output(postsamp <- try(suppressWarnings(revdbayes::rpost(n = 300, model = "gp", data = ang_weighted_dat,
                    thresh = quantile(ang_weighted_dat, qu), prior = revdbayes::set_prior(prior = "flat", model = "gp", min_xi = 1,
                      max_xi = 1/max(c(1 - w, w))), init_ests = start)), silent = TRUE)))
                }
                if (inherits(postsamp, what = "try-error")) {
                  return(rep(NA, 3))
                } else {
                  v <- as.vector(apply(postsamp$sim_vals, 2, median))
                  return(quantile(1/postsamp$sim_vals[, 2], c(0.025, 0.5, 0.975)))
                }
            }
        }
    })
    if (method %in% c("hill", "mle")) {
        lower <- pmax(c(1 - w_seq[w_seq < 0.5], w_seq[w_seq <= 0.5] + 0.5), pmin(1, lambda_seq[1, ] - qnorm(0.975) * lambda_seq[2,
            ]))
        upper <- pmax(c(1 - w_seq[w_seq < 0.5], w_seq[w_seq <= 0.5] + 0.5), pmin(1, lambda_seq[1, ] + qnorm(0.975) * lambda_seq[2,
            ]))
        pe <- pmax(c(1 - w_seq[w_seq < 0.5], w_seq[w_seq <= 0.5] + 0.5), pmin(1, lambda_seq[1, ]))
    } else if (method == "bayes") {
        lower <- lambda_seq[1, ]
        upper <- lambda_seq[3, ]
        pe <- lambda_seq[2, ]
    }
    if (plot) {
        plot(type = "n", x = 0.5, y = 1, xlim = c(0, 1), ylim = c(0.5, 1), xlab = expression(omega), ylab = expression(lambda(omega)),
            bty = "l")
        segments(x1 = 0.5, x0 = 0, y1 = 0.5, y0 = 1, col = "gray")
        segments(x1 = 0.5, x0 = 1, y1 = 0.5, y0 = 1, col = "gray")
        segments(x1 = 0, x0 = 1, y1 = 1, y0 = 1, col = "gray")
        lines(w_seq, lower, lty = 1, col = "red")
        lines(w_seq, upper, lty = 1, col = "red")
        lines(w_seq, pe, lwd = 2)
    }
    invisible(list(w = w_seq, lambda = pe, lower.confint = lower, upper.confint = upper))
}
