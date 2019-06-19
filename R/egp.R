### Code by Raphael Huser with contributions from P. Naveau for the methods described in @references Naveau, P., R. Huser, P.
### Ribereau, and A. Hannart (2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection,
### \emph{Water Resour. Res.}, 52, 2753-2769, \code{doi:10.1002/2015WR018552}.


#' Fit an extended generalized Pareto distribution of Naveau et al.
#'
#' This is a wrapper function to obtain PWM or MLE estimates for
#' the extended GP models of Naveau et al. (2016) for rainfall intensities. The function calculates confidence intervals
#' by means of nonparametric percentile bootstrap and returns histograms and QQ plots of
#' the fitted distributions. The function handles both censoring and rounding.
#'
#' @details
#' The different models include the following transformations:
#' \itemize{
#' \item \code{model} 0 corresponds to uniform carrier, \eqn{G(u)=u}.
#' \item \code{model} 1 corresponds to a three parameters family, with carrier \eqn{G(u)=u^\kappa}.
#' \item \code{model} 2 corresponds to a three parameters family, with carrier \eqn{G(u)=1-V_\delta((1-u)^\delta)}.
#' \item \code{model} 3 corresponds to a four parameters family, with carrier \deqn{G(u)=1-V_\delta((1-u)^\delta))^{\kappa/2}}.
#' \item \code{model} 4 corresponds to a five parameter model (a mixture of \code{type} 2, with \eqn{G(u)=pu^\kappa + (1-p)*u^\delta}
#' }
#'
#' @author Raphael Huser and Philippe Naveau
#' @param data data vector.
#' @param model integer ranging from 0 to 4 indicating the model to select (see \code{\link{extgp}}).
#' @param method string; either \code{'mle'} for maximum likelihood, or \code{'pwm'} for probability weighted moments, or both.
#' @param init vector of initial values, comprising of \eqn{p}, \eqn{\kappa}, \eqn{\delta},\eqn{\sigma},\eqn{\xi} (in that order) for the optimization. All parameters may not appear depending on \code{model}.
#' @param censoring numeric vector of length 2 containing the lower and upper bound for censoring; \code{censoring=c(0,Inf)} is equivalent to no censoring.
#' @param rounded numeric giving the instrumental precision (and rounding of the data), with default of 0.
#' @param confint logical; should confidence interval be returned (percentile bootstrap).
#' @param R integer; number of bootstrap replications.
#' @param ncpus integer; number of CPUs for parallel calculations (default: 1).
#' @param plots logical; whether to produce histogram and density plots.
#' @export
#' @importFrom boot boot boot.ci
#' @importFrom grDevices rgb
#' @importFrom graphics hist layout
#' @importFrom evd dgpd pgpd qgpd
#' @seealso \code{\link{egp.fit}}, \code{\link{egp}}, \code{\link{extgp}}
#' @references Naveau, P., R. Huser, P. Ribereau, and A. Hannart (2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, \emph{Water Resour. Res.}, 52, 2753-2769, \code{doi:10.1002/2015WR018552}.
#' @examples
#' \dontrun{
#' data(rain, package = "ismev")
#' fit.extgp(rain[rain>0], model=1, method = 'mle', init = c(0.9, gp.fit(rain, 0)$est),
#'  rounded = 0.1, confint = TRUE, R = 20)
#' }
fit.extgp <- function(data, model = 1, method = c("mle", "pwm"), init, censoring = c(0, Inf),
                      rounded = 0, confint = FALSE, R = 1000, ncpus = 1, plots = TRUE) {
    method <- match.arg(method, c("mle", "pwm"), several.ok = TRUE)
    if("pwm" %in% method){
      if (!requireNamespace("gmm", quietly = TRUE)) {
        stop("Package \"gmm\" needed for this function to work. Please install it.",
             call. = FALSE)
      }
    }
    # Sanity checks
    initsize <- switch(model, 3, 3, 4, 5)
    if (length(init) != initsize) {
        stop("Invalid starting values in `init'; incorrect length.")
    }
    data = data[data > 0]
    if (model == 3 && "pwm" %in% method) {
        stop("No explicit formula for PWMs for this model...")
    }

    x <- seq(0, max(data, na.rm = TRUE), by = 0.05)
    dextgps <- c()
    qextgps <- c()

    if (any(method == "pwm")) {
        if (model == 1) {
            fit.pwm <- extgp.fit.pwm(x = data, type = 1, kappa0 = init[1], sigma0 = init[2], xi0 = init[3], censoring = censoring)
            if (confint) {
                fit.pwm.boot <- boot::boot(data = data, statistic = extgp.fit.pwm.boot, R = R, type = 1, kappa0 = init[1], sigma0 = init[2],
                  ncpus = ncpus, xi0 = init[3], censoring = censoring, parallel = "multicore")
                CI.pwm.kappa <- boot::boot.ci(boot.out = fit.pwm.boot, index = 1, type = "perc")$perc[4:5]
                CI.pwm.sigma <- boot::boot.ci(boot.out = fit.pwm.boot, index = 2, type = "perc")$perc[4:5]
                CI.pwm.xi <- boot::boot.ci(boot.out = fit.pwm.boot, index = 3, type = "perc")$perc[4:5]
                CIs.pwm <- cbind(CI.pwm.kappa, CI.pwm.sigma, CI.pwm.xi)
            }
            dextgp.pwm <- dextgp(x = x, type = 1, kappa = fit.pwm[1], sigma = fit.pwm[2], xi = fit.pwm[3])
            dextgps <- c(dextgps, dextgp.pwm)
            if (plots) {
                qextgp.pwm <- qextgp(p = c(1:length(data))/(length(data) + 1), type = 1, kappa = fit.pwm[1], sigma = fit.pwm[2], xi = fit.pwm[3])
                qextgps <- c(qextgps, qextgp.pwm)
                if (confint) {
                  q.pwm.boot <- mapply(FUN = qextgp, p = list(c(1:length(data))/(length(data) + 1)), type = list(1),
                                       kappa = as.list(fit.pwm.boot$t[,
                    1]), sigma = as.list(fit.pwm.boot$t[, 2]), xi = as.list(fit.pwm.boot$t[, 3]))
                  q.pwm.L <- apply(q.pwm.boot, 1, quantile, 0.025, na.rm = TRUE)
                  q.pwm.U <- apply(q.pwm.boot, 1, quantile, 0.975, na.rm = TRUE)
                }
            }
        } else if (model == 2) {
            fit.pwm <- extgp.fit.pwm(x = data, type = 2, delta0 = init[1], sigma0 = init[2], xi0 = init[3], censoring = censoring)
            if (confint) {
                fit.pwm.boot <- boot::boot(data = data, statistic = extgp.fit.pwm.boot, R = R, type = 2, delta0 = init[1], sigma0 = init[2],
                  xi0 = init[3], censoring = censoring, parallel = "multicore", ncpus = ncpus)
                CI.pwm.delta <- boot::boot.ci(boot.out = fit.pwm.boot, index = 1, type = "perc")$perc[4:5]
                CI.pwm.sigma <- boot::boot.ci(boot.out = fit.pwm.boot, index = 2, type = "perc")$perc[4:5]
                CI.pwm.xi <- boot::boot.ci(boot.out = fit.pwm.boot, index = 3, type = "perc")$perc[4:5]
                CIs.pwm <- cbind(CI.pwm.delta, CI.pwm.sigma, CI.pwm.xi)
            }
            dextgp.pwm <- dextgp(x = x, type = 2, delta = fit.pwm[1], sigma = fit.pwm[2], xi = fit.pwm[3])
            dextgps <- c(dextgps, dextgp.pwm)
            if (plots) {
                qextgp.pwm <- qextgp(p = c(1:length(data))/(length(data) + 1), type = 2, delta = fit.pwm[1], sigma = fit.pwm[2], xi = fit.pwm[3])
                qextgps <- c(qextgps, qextgp.pwm)
                if (confint) {
                  q.pwm.boot <- mapply(FUN = qextgp, p = list(c(1:length(data))/(length(data) + 1)),
                                       type = list(2), delta = as.list(fit.pwm.boot$t[,
                    1]), sigma = as.list(fit.pwm.boot$t[, 2]), xi = as.list(fit.pwm.boot$t[, 3]))
                  q.pwm.L <- apply(q.pwm.boot, 1, quantile, 0.025, na.rm = TRUE)
                  q.pwm.U <- apply(q.pwm.boot, 1, quantile, 0.975, na.rm = TRUE)
                }
            }
        } else if (model == 3) {
            fit.pwm <- extgp.fit.pwm(x = data, type = 3, kappa0 = init[1], delta0 = init[2], sigma0 = init[3], xi0 = init[4], censoring = censoring)
            if (confint) {
                fit.pwm.boot <- boot::boot(data = data, statistic = extgp.fit.pwm.boot, R = R, type = 3, kappa0 = init[1], delta0 = init[2],
                  sigma0 = init[3], xi0 = init[4], censoring = censoring, parallel = "multicore", ncpus = ncpus)
                CI.pwm.kappa <- boot::boot.ci(boot.out = fit.pwm.boot, index = 1, type = "perc")$perc[4:5]
                CI.pwm.delta <- boot::boot.ci(boot.out = fit.pwm.boot, index = 2, type = "perc")$perc[4:5]
                CI.pwm.sigma <- boot::boot.ci(boot.out = fit.pwm.boot, index = 3, type = "perc")$perc[4:5]
                CI.pwm.xi <- boot::boot.ci(boot.out = fit.pwm.boot, index = 4, type = "perc")$perc[4:5]
                CIs.pwm <- cbind(CI.pwm.kappa, CI.pwm.delta, CI.pwm.sigma, CI.pwm.xi)
            }
            dextgp.pwm <- dextgp(x = x, type = 3, kappa = fit.pwm[1], delta = fit.pwm[2], sigma = fit.pwm[3], xi = fit.pwm[4])
            dextgps <- c(dextgps, dextgp.pwm)
            if (plots) {
                qextgp.pwm <- qextgp(p = c(1:length(data))/(length(data) + 1), type = 3, kappa = fit.pwm[1], delta = fit.pwm[2], sigma = fit.pwm[3],
                  xi = fit.pwm[4])
                qextgps <- c(qextgps, qextgp.pwm)
                if (confint) {
                  q.pwm.boot <- mapply(FUN = qextgp, p = list(c(1:length(data))/(length(data) + 1)), type = list(3), kappa = as.list(fit.pwm.boot$t[,
                    1]), delta = as.list(fit.pwm.boot$t[, 2]), sigma = as.list(fit.pwm.boot$t[, 3]), xi = as.list(fit.pwm.boot$t[,
                    4]))
                  q.pwm.L <- apply(q.pwm.boot, 1, quantile, 0.025, na.rm = TRUE)
                  q.pwm.U <- apply(q.pwm.boot, 1, quantile, 0.975, na.rm = TRUE)
                }
            }
        } else if (model == 4) {
            fit.pwm <- extgp.fit.pwm(x = data, type = 4, prob0 = init[1], kappa0 = init[2], delta0 = init[3], sigma0 = init[4], xi0 = init[5],
                censoring = censoring)
            if (confint) {
                fit.pwm.boot <- boot::boot(data = data, statistic = extgp.fit.pwm.boot, R = R, type = 4, prob0 = init[1], kappa0 = init[2],
                  delta0 = init[3], sigma0 = init[4], xi0 = init[5], censoring = censoring, parallel = "multicore", ncpus = ncpus)
                CI.pwm.prob <- boot::boot.ci(boot.out = fit.pwm.boot, index = 1, type = "perc")$perc[4:5]
                CI.pwm.kappa <- boot::boot.ci(boot.out = fit.pwm.boot, index = 2, type = "perc")$perc[4:5]
                CI.pwm.delta <- boot::boot.ci(boot.out = fit.pwm.boot, index = 3, type = "perc")$perc[4:5]
                CI.pwm.sigma <- boot::boot.ci(boot.out = fit.pwm.boot, index = 4, type = "perc")$perc[4:5]
                CI.pwm.xi <- boot::boot.ci(boot.out = fit.pwm.boot, index = 5, type = "perc")$perc[4:5]
                CIs.pwm <- cbind(CI.pwm.prob, CI.pwm.kappa, CI.pwm.delta, CI.pwm.sigma, CI.pwm.xi)
            }
            dextgp.pwm <- dextgp(x = x, type = 4, prob = fit.pwm[1], kappa = fit.pwm[2], delta = fit.pwm[3], sigma = fit.pwm[4], xi = fit.pwm[5])
            dextgps <- c(dextgps, dextgp.pwm)
            if (plots) {
                qextgp.pwm <- qextgp(p = c(1:length(data))/(length(data) + 1), type = 4, prob = fit.pwm[1], kappa = fit.pwm[2], delta = fit.pwm[3],
                  sigma = fit.pwm[4], xi = fit.pwm[5])
                qextgps <- c(qextgps, qextgp.pwm)
                if (confint) {
                  q.pwm.boot <- mapply(FUN = qextgp, p = list(c(1:length(data))/(length(data) + 1)), type = list(4), prob = as.list(fit.pwm.boot$t[,
                    1]), kappa = as.list(fit.pwm.boot$t[, 2]), delta = as.list(fit.pwm.boot$t[, 3]), sigma = as.list(fit.pwm.boot$t[,
                    4]), xi = as.list(fit.pwm.boot$t[, 5]))
                  q.pwm.L <- apply(q.pwm.boot, 1, quantile, 0.025, na.rm = TRUE)
                  q.pwm.U <- apply(q.pwm.boot, 1, quantile, 0.975, na.rm = TRUE)
                }
            }
        }
    }

    if (any(method == "mle")) {
        if (model == 1) {
            fit.mle <- extgp.fit.ml(x = data, type = 1, kappa0 = init[1], sigma0 = init[2], xi0 = init[3], censoring = censoring, rounded = rounded)
            if (confint) {
                fit.mle.boot <- boot::boot(data = data, statistic = extgp.fit.ml.boot, R = R, type = 1, kappa0 = init[1], sigma0 = init[2],
                  xi0 = init[3], censoring = censoring, rounded = rounded, parallel = "multicore", ncpus = ncpus)
                CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
                CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
                CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
                CIs.mle <- cbind(CI.mle.kappa, CI.mle.sigma, CI.mle.xi)
            }
            dextgp.mle <- dextgp(x = x, type = 1, kappa = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
            dextgps <- c(dextgps, dextgp.mle)
            if (plots) {
                qextgp.mle <- qextgp(p = c(1:length(data))/(length(data) + 1), type = 1, kappa = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
                qextgps <- c(qextgps, qextgp.mle)
                if (confint) {
                  q.mle.boot <- mapply(FUN = qextgp, p = list(c(1:length(data))/(length(data) + 1)), type = list(1), kappa = as.list(fit.mle.boot$t[,
                    1]), sigma = as.list(fit.mle.boot$t[, 2]), xi = as.list(fit.mle.boot$t[, 3]))
                  q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
                  q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
                }
            }
        } else if (model == 2) {
            fit.mle <- extgp.fit.ml(x = data, type = 2, delta0 = init[1], sigma0 = init[2], xi0 = init[3], censoring = censoring, rounded = rounded)
            if (confint) {
                fit.mle.boot <- boot::boot(data = data, statistic = extgp.fit.ml.boot, R = R, type = 2, delta0 = init[1], sigma0 = init[2],
                  xi0 = init[3], censoring = censoring, rounded = rounded, parallel = "multicore", ncpus = ncpus)
                CI.mle.delta <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
                CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
                CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
                CIs.mle <- cbind(CI.mle.delta, CI.mle.sigma, CI.mle.xi)
            }
            dextgp.mle <- dextgp(x = x, type = 2, delta = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
            dextgps <- c(dextgps, dextgp.mle)
            if (plots) {
                qextgp.mle <- qextgp(p = c(1:length(data))/(length(data) + 1), type = 2, delta = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
                qextgps <- c(qextgps, qextgp.mle)
                if (confint) {
                  q.mle.boot <- mapply(FUN = qextgp, p = list(c(1:length(data))/(length(data) + 1)), type = list(2), delta = as.list(fit.mle.boot$t[,
                    1]), sigma = as.list(fit.mle.boot$t[, 2]), xi = as.list(fit.mle.boot$t[, 3]))
                  q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
                  q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
                }
            }
        } else if (model == 3) {
            fit.mle <- extgp.fit.ml(x = data, type = 3, kappa0 = init[1], delta0 = init[2], sigma0 = init[3], xi0 = init[4], censoring = censoring,
                rounded = rounded)
            if (confint) {
                fit.mle.boot <- boot::boot(data = data, statistic = extgp.fit.ml.boot, R = R, type = 3, kappa0 = init[1], delta0 = init[2],
                  sigma0 = init[3], xi0 = init[4], censoring = censoring, rounded = rounded, parallel = "multicore", ncpus = ncpus)
                CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
                CI.mle.delta <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
                CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
                CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 4, type = "perc")$perc[4:5]
                CIs.mle <- cbind(CI.mle.kappa, CI.mle.delta, CI.mle.sigma, CI.mle.xi)
            }
            dextgp.mle <- dextgp(x = x, type = 3, kappa = fit.mle[1], delta = fit.mle[2], sigma = fit.mle[3], xi = fit.mle[4])
            dextgps <- c(dextgps, dextgp.mle)
            if (plots) {
                qextgp.mle <- qextgp(p = c(1:length(data))/(length(data) + 1), type = 3, kappa = fit.mle[1], delta = fit.mle[2], sigma = fit.mle[3],
                  xi = fit.mle[4])
                qextgps <- c(qextgps, qextgp.mle)
                if (confint) {
                  q.mle.boot <- mapply(FUN = qextgp, p = list(c(1:length(data))/(length(data) + 1)), type = list(3), kappa = as.list(fit.mle.boot$t[,
                    1]), delta = as.list(fit.mle.boot$t[, 2]), sigma = as.list(fit.mle.boot$t[, 3]), xi = as.list(fit.mle.boot$t[,
                    4]))
                  q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
                  q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
                }
            }
        } else if (model == 4) {
            fit.mle <- extgp.fit.ml(x = data, type = 4, prob0 = init[1], kappa0 = init[2], delta0 = init[3], sigma0 = init[4], xi0 = init[5],
                censoring = censoring, rounded = rounded)
            if (confint) {
                fit.mle.boot <- boot::boot(data = data, statistic = extgp.fit.ml.boot, R = R, type = 4, prob0 = init[1], kappa0 = init[2],
                  delta0 = init[3], sigma0 = init[4], xi0 = init[5], censoring = censoring, rounded = rounded, parallel = "multicore",
                  ncpus = ncpus)
                CI.mle.prob <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
                CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
                CI.mle.delta <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
                CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 4, type = "perc")$perc[4:5]
                CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 5, type = "perc")$perc[4:5]
                CIs.mle <- cbind(CI.mle.prob, CI.mle.kappa, CI.mle.delta, CI.mle.sigma, CI.mle.xi)
            }
            dextgp.mle <- dextgp(x = x, type = 4, prob = fit.mle[1], kappa = fit.mle[2], delta = fit.mle[3], sigma = fit.mle[4], xi = fit.mle[5])
            dextgps <- c(dextgps, dextgp.mle)
            if (plots) {
                qextgp.mle <- qextgp(p = c(1:length(data))/(length(data) + 1), type = 4, prob = fit.mle[1], kappa = fit.mle[2], delta = fit.mle[3],
                  sigma = fit.mle[4], xi = fit.mle[5])
                qextgps <- c(qextgps, qextgp.mle)
                if (confint) {
                  q.mle.boot <- mapply(FUN = qextgp, p = list(c(1:length(data))/(length(data) + 1)), type = list(4), prob = as.list(fit.mle.boot$t[,
                    1]), kappa = as.list(fit.mle.boot$t[, 2]), delta = as.list(fit.mle.boot$t[, 3]), sigma = as.list(fit.mle.boot$t[,
                    4]), xi = as.list(fit.mle.boot$t[, 5]))
                  q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
                  q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
                }
            }
        }
    }

    if (any(method == "pwm")) {
        if (any(method == "mle")) {
            fits <- list(pwm = fit.pwm, mle = fit.mle)
            if (confint) {
                CIs <- list(pwm = CIs.pwm, mle = CIs.mle)
            } else {
                CIs <- NULL
            }
        } else {
            fits <- list(pwm = fit.pwm)
            if (confint & plots) {
                CIs <- list(pwm = CIs.pwm)
            } else {
                CIs <- NULL
            }
        }
    } else {
        if (any(method == "mle")) {
            fits <- list(mle = fit.mle)
            if (confint) {
                CIs <- list(mle = CIs.mle)
            } else {
                CIs <- NULL
            }
        } else {
            fits <- NULL
            CIs <- NULL
        }
    }

    if (plots) {
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
        mat <- matrix(c(1:(1 + plots)), nrow = 1, ncol = 1 + plots, byrow = TRUE)
        layout(mat)
        par(mar = c(4, 4, 1, 1))
        hist(data, breaks = c(0:30) * max(data + 10^-1)/30, freq = FALSE, xlim = range(x),
             xlab = "Rainfall amounts [mm]", ylab = "Density",
            main = "", col = "lightgrey")
        dens <- density(data, from = 0)
        lines(dens$x, dens$y, col = "black")
        abline(v = censoring, col = "lightgrey")

        if (any(method == "pwm")) {
            lines(x, dextgp.pwm, col = "red", lty = 2)
        }
        if (any(method == "mle")) {
            lines(x, dextgp.mle, col = "blue")
        }

        if (any(method == "pwm")) {
            plot(data, qextgp.pwm, asp = 1, xlab = "Empirical quantiles", ylab = "Fitted quantiles", ylim = range(qextgps, na.rm = TRUE),
                type = "n")
        } else {
            if (any(method == "mle")) {
                plot(data, qextgp.mle, asp = 1, xlab = "Empirical quantiles", ylab = "Fitted quantiles", ylim = range(qextgps, na.rm = TRUE),
                  type = "n")
            }
        }
        if (confint) {
            if (any(method == "pwm")) {
                polygon(x = c(sort(data), sort(data, decreasing = TRUE)), y = c(q.pwm.L, q.pwm.U[length(q.pwm.U):1]), col = rgb(1,
                  0, 0, alpha = 0.1), lty = 2, border = rgb(1, 0, 0, alpha = 0.5))
            }
            if (any(method == "mle")) {
                polygon(x = c(sort(data), sort(data, decreasing = TRUE)), y = c(q.mle.L, q.mle.U[length(q.mle.U):1]), col = rgb(0,
                  0, 1, alpha = 0.1), lty = 1, border = rgb(0, 0, 1, alpha = 0.5))
            }
        }
        if (any(method == "pwm")) {
            lines(sort(data), qextgp.pwm, lty = 2, type = "b", pch = 20, col = "red")
        }
        if (any(method == "mle")) {
            lines(sort(data), qextgp.mle, lty = 1, type = "b", pch = 20, col = "blue")
        }
        abline(0, 1, col = "lightgrey")
    }

    return(list(fit = fits, confint = CIs))
}

#' Fit an extended generalized Pareto distribution of Naveau et al.
#'
#' Deprecated function name to fit an extended generalized Pareto family. The user should call \code{\link[mev]{fit.extgp}} instead.
#'
#' @seealso \code{\link[mev]{fit.extgp}}
#' @export
#' @keywords internal
#' @inheritParams fit.extgp
egp2.fit <- function(data, model = 1, method = c("mle", "pwm"), init, censoring = c(0, Inf), rounded = 0, CI = FALSE, R = 1000, ncpus = 1, plots = TRUE) {
  fit.extgp(data = data, model = model, method = method, init = init, censoring = censoring, rounded = rounded,
            confint = CI, R = R, ncpus = ncpus, plots = plots)
}
################################################################ ### Distribution/Density/Quantile/Random generator functions ### ###
### different types of distribution G(u):
### type=0: G(u)=u type=1: G(u)=u^kappa
### type=2: G(u)=1-V_delta((1-u)^delta)
### type=3: G(u)=(1-V_delta((1-u)^delta))^(kappa/2)
### type=4: G(u)=p*u^kappa + (1-p)*u^delta



#' @title Carrier distribution for the extended GP distributions of Naveau et al.
#' @description Density, distribution function, quantile function and random number
#' generation for the carrier distributions of the extended Generalized Pareto distributions.
#' @name extgp.G
#' @seealso \code{\link{extgp}}
#' @param u vector of observations (\code{dextgp.G}), probabilities (\code{qextgp.G}) or quantiles (\code{pextgp.G}), in \eqn{[0,1]}
#' @param prob mixture probability for model \code{type} \code{4}
#' @param kappa shape parameter for \code{type} \code{1}, \code{3} and \code{4}
#' @param delta additional parameter for \code{type} \code{2}, \code{3} and \code{4}
#' @param type integer between 0 to 5 giving the model choice
#' @param log logical; should the log-density be returned (default to \code{FALSE})?
#' @param n sample size
#' @param unifsamp sample of uniform; if provided, the data will be used in place of new uniform random variates
#' @param censoring numeric vector of length 2 containing the lower and upper bound for censoring
#' @param direct logical; which method to use for sampling in model of \code{type} \code{4}?
#' @author  Raphael Huser and Philippe Naveau
#'
#' @section Usage: \code{pextgp.G(u, type=1,prob, kappa, delta)}
#' @section Usage: \code{dextgp.G(u,type=1, prob=NA, kappa=NA, delta=NA, log=FALSE)}
#' @section Usage: \code{qextgp.G(u,type=1, prob=NA,kappa=NA,delta=NA)}
#' @section Usage: \code{rextgp.G(n,prob=NA,kappa=NA,delta=NA,type=1,unifsamp=NULL,direct=FALSE,censoring=c(0,1))}
NULL

#' @rdname extgp-functions
#' @export
#' @keywords internal
pextgp.G <- function(u, type = 1, prob, kappa, delta) {
    if (!type %in% 1:5) {
        stop("Invalid `type' argument")
    }
    type <- type[1]
    if (type %in% c(1, 3, 4) && missing(kappa)) {
        stop("Argument `kappa' missing.")
    }
    if (type %in% c(2, 3, 4) && missing(delta)) {
        stop("Argument `delta' missing.")
    }
    if (type == 4 && missing(prob)) {
        stop("Argument `prob' missing.")
    }
    if (type == 0) {
        return(u)
    } else if (type == 1) {
        return(u^kappa)
    } else if (type == 2) {
        return(1 - pbeta((1 - u)^delta, 1/delta, 2))
    } else if (type == 3) {
        return((1 - pbeta((1 - u)^delta, 1/delta, 2))^(kappa/2))
    } else if (type == 4) {
        return(prob * u^kappa + (1 - prob) * u^delta)
    }
}

#' @rdname extgp-functions
#' @export
#' @keywords internal
dextgp.G <- function(u, type = 1, prob = NA, kappa = NA, delta = NA, log = FALSE) {
    if (!type %in% 1:5) {
        stop("Invalid `type' argument")
    }
    type <- type[1]
    if (type %in% c(1, 3, 4) && missing(kappa)) {
        stop("Argument `kappa' missing.")
    }
    if (type %in% c(2, 3, 4) && missing(delta)) {
        stop("Argument `delta' missing.")
    }
    if (type == 4 && missing(prob)) {
        stop("Argument `prob' missing.")
    }
    if (log == FALSE) {
        if (type == 0) {
            return(1)
        } else if (type == 1) {
            return(kappa * u^(kappa - 1))
        } else if (type == 2) {
            return(dbeta((1 - u)^delta, 1/delta, 2) * delta * (1 - u)^(delta - 1))
        } else if (type == 3) {
            return((kappa/2) * (1 - pbeta((1 - u)^delta, 1/delta, 2))^(kappa/2 - 1) * dbeta((1 - u)^delta, 1/delta, 2) * delta * (1 -
                u)^(delta - 1))
        } else if (type == 4) {
            return(prob * kappa * u^(kappa - 1) + (1 - prob) * delta * u^(delta - 1))
        }
    } else {
        if (type == 0) {
            return(0)
        } else if (type == 1) {
            return(log(kappa) + (kappa - 1) * log(u))
        } else if (type == 2) {
            return(dbeta((1 - u)^delta, 1/delta, 2, log = TRUE) + log(delta) + (delta - 1) * log(1 - u))
        } else if (type == 3) {
            return(log(kappa/2) + (kappa/2 - 1) * log(1 - pbeta((1 - u)^delta, 1/delta, 2)) + dbeta((1 - u)^delta, 1/delta, 2, log = TRUE) +
                log(delta) + (delta - 1) * log(1 - u))
        } else if (type == 4) {
            return(log(prob * kappa * u^(kappa - 1) + (1 - prob) * delta * u^(delta - 1)))
        }
    }
}

#' @rdname extgp-functions
#' @export
#' @keywords internal
qextgp.G <- function(u, type = 1, prob = NA, kappa = NA, delta = NA) {
    if (!type %in% 1:5) {
        stop("Invalid `type' argument")
    }
    type <- type[1]
    if (type %in% c(1, 3, 4) && missing(kappa)) {
        stop("Argument `kappa' missing.")
    }
    if (type %in% c(2, 3, 4) && missing(delta)) {
        stop("Argument `delta' missing.")
    }
    if (type == 4 && missing(prob)) {
        stop("Argument `prob' missing.")
    }
    if (type == 0) {
        return(u)
    } else if (type == 1) {
        return(u^(1/kappa))
    } else if (type == 2) {
        return(1 - qbeta(1 - u, 1/delta, 2)^(1/delta))
    } else if (type == 3) {
        return(1 - qbeta(1 - u^(2/kappa), 1/delta, 2)^(1/delta))
    } else if (type == 4) {
        dummy.func <- function(u, p, prob = NA, kappa = NA, delta = NA) {
            return(pextgp.G(u = u, prob = prob, kappa = kappa, delta = delta, type = 4) - p)
        }
        find.root <- function(u, prob = NA, kappa = NA, delta = NA) {
            return(uniroot(dummy.func, interval = c(0, 1), p = u, prob = prob, kappa = kappa, delta = delta)$root)
        }
        return(sapply(u, FUN = find.root, prob = prob, kappa = kappa, delta = delta))
    }
}

#' @rdname extgp-functions
#' @export
#' @keywords internal
rextgp.G <- function(n, prob = NA, kappa = NA, delta = NA, type = 1, unifsamp = NULL, direct = FALSE, censoring = c(0, 1)) {
    if (!type %in% 1:5) {
        stop("Invalid `type' argument")
    }
    type <- type[1]
    if (type %in% c(1, 3, 4) && missing(kappa)) {
        stop("Argument `kappa' missing.")
    }
    if (type %in% c(2, 3, 4) && missing(delta)) {
        stop("Argument `delta' missing.")
    }
    if (type == 4 && missing(prob)) {
        stop("Argument `prob' missing.")
    }
    if (is.null(unifsamp)) {
        unifsamp <- runif(n)
    }
    if (type != 4 | (type == 4 & direct)) {
        if (censoring[1] == 0 & censoring[2] == 1) {
            return(qextgp.G(unifsamp, prob = prob, kappa = kappa, delta = delta, type = type))
        } else {
            x.L <- pextgp.G(censoring[1], prob = prob, kappa = kappa, delta = delta, type = type)
            x.U <- pextgp.G(censoring[2], prob = prob, kappa = kappa, delta = delta, type = type)
            return(qextgp.G(x.L + (x.U - x.L) * unifsamp, prob = prob, kappa = kappa, delta = delta, type = type))
        }
    } else if (type == 4 & !direct) {
        if (censoring[1] == 0 & censoring[2] == 1) {
            components <- sample(x = c(1, 2), size = n, replace = TRUE, prob = c(prob, 1 - prob))
            res <- c()
            res[components == 1] <- qextgp.G(unifsamp[components == 1], prob = NA, kappa = kappa, delta = NA, type = 1)
            res[components == 2] <- qextgp.G(unifsamp[components == 2], prob = NA, kappa = delta, delta = NA, type = 1)
            return(res)
        } else {
            x.L <- pextgp.G(censoring[1], prob = prob, kappa = kappa, delta = delta, type = type)
            x.U <- pextgp.G(censoring[2], prob = prob, kappa = kappa, delta = delta, type = type)
            prop1 <- prob * (pextgp.G(censoring[2], prob = NA, kappa = kappa, delta = NA, type = 1) - pextgp.G(censoring[1], prob = NA,
                kappa = kappa, delta = NA, type = 1))
            prop2 <- (1 - prob) * (pextgp.G(censoring[2], prob = NA, kappa = delta, delta = NA, type = 1) - pextgp.G(censoring[1], prob = NA,
                kappa = delta, delta = NA, type = 1))
            new.prob <- prop1/(prop1 + prop2)
            components <- sample(x = c(1, 2), size = n, replace = TRUE, prob = c(new.prob, 1 - new.prob))
            res <- c()
            res[components == 1] <- qextgp.G(x.L + (x.U - x.L) * unifsamp[components == 1], prob = NA, kappa = kappa, delta = NA, type = 1)
            res[components == 2] <- qextgp.G(x.L + (x.U - x.L) * unifsamp[components == 2], prob = NA, kappa = delta, delta = NA, type = 1)
            return(res)
        }
    }
}


### different types of extended GP type=0:
### exact GP ---> 2 parameters (sigma,xi)
### type=1: extgp with G(u)=u^kappa ---> 3 parameters(kappa,sigma,xi)
### type=2: extgp with G(u)=1-V_delta((1-u)^delta) ---> 3 parameters (delta,sigma,xi)
### type=3: extgp with G(u)=(1-V_delta((1-u)^delta))^(kappa/2) ---> 4 parameters (kappa,delta,sigma,xi)
### type=4: extgp with mixture G(u)=p*u^kappa + (1-p)*u^delta ---> 5 parameters (p,kappa,delta,sigma,xi)

#' @title Extended generalised Pareto families of Naveau et al. (2016)
#'
#' @description Density function, distribution function, quantile function and random generation for the extended generalized Pareto distribution (GPD) with scale and shape parameters.
#'
#' @details The extended generalized Pareto families proposed in Naveau \emph{et al.} (2016)
#'retain the tail index of the distribution while being compliant with the theoretical behavior of extreme
#'low rainfall. There are five proposals, the first one being equivalent to the GP distribution.
#'\itemize{
#' \item \code{type} 0 corresponds to uniform carrier, \eqn{G(u)=u}.
#' \item \code{type} 1 corresponds to a three parameters family, with carrier \eqn{G(u)=u^\kappa}.
#' \item \code{type} 2 corresponds to a three parameters family, with carrier \eqn{G(u)=1-V_\delta((1-u)^\delta)}.
#' \item \code{type} 3 corresponds to a four parameters family, with carrier \deqn{G(u)=1-V_\delta((1-u)^\delta))^{\kappa/2}}.
#' \item \code{type} 4 corresponds to a five parameter model (a mixture of \code{type} 2, with \eqn{G(u)=pu^\kappa + (1-p)*u^\delta}
#' }
#' @name extgp
#' @param q vector of quantiles
#' @param x vector of observations
#' @param p vector of probabilities
#' @param n sample size
#' @param prob mixture probability for model \code{type} \code{4}
#' @param kappa shape parameter for \code{type} \code{1}, \code{3} and \code{4}
#' @param delta additional parameter for \code{type} \code{2}, \code{3} and \code{4}
#' @param sigma scale parameter
#' @param xi shape parameter
#' @param type integer between 0 to 5 giving the model choice
#' @param log logical; should the log-density be returned (default to \code{FALSE})?
#' @param unifsamp sample of uniform; if provided, the data will be used in place of new uniform random variates
#' @param censoring numeric vector of length 2 containing the lower and upper bound for censoring
#'
#' @section Usage: \code{pextgp(q, prob=NA, kappa=NA, delta=NA, sigma=NA, xi=NA, type=1)}
#' @section Usage: \code{dextgp(x, prob=NA, kappa=NA, delta=NA, sigma=NA, xi=NA, type=1, log=FALSE)}
#' @section Usage: \code{qextgp(p, prob=NA, kappa=NA, delta=NA, sigma=NA, xi=NA, type=1)}
#' @section Usage: \code{rextgp(n, prob=NA, kappa=NA, delta=NA, sigma=NA, xi=NA, type=1, unifsamp=NULL, censoring=c(0,Inf))}
#' @author  Raphael Huser and Philippe Naveau
#' @references Naveau, P., R. Huser, P. Ribereau, and A. Hannart (2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection, \emph{Water Resour. Res.}, 52, 2753-2769, \code{doi:10.1002/2015WR018552}.
NULL


#' Extended GP functions
#'
#' These functions are documented in \code{\link{extgp}} and in \code{\link{extgp.G}} for the carrier distributions supported in the unit interval.
#' @name extgp-functions
#' @export
#' @seealso \code{\link{extgp}}, \code{\link{extgp.G}}
#' @keywords internal
pextgp <- function(q, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
    return(pextgp.G(pgpd(q, scale = sigma, shape = xi), prob = prob, kappa = kappa, delta = delta, type = type))
}

#' @rdname extgp-functions
#' @export
#' @keywords internal
dextgp <- function(x, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, log = FALSE) {
    if (log == FALSE) {
        return(dextgp.G(pgpd(x, scale = sigma, shape = xi), prob = prob, kappa = kappa, delta = delta, type = type) * dgpd(x, scale = sigma,
            shape = xi))
    } else {
        return(dextgp.G(pgpd(x, scale = sigma, shape = xi), prob = prob, kappa = kappa, delta = delta, type = type, log = TRUE) + dgpd(x,
            scale = sigma, shape = xi, log = TRUE))
    }
}

#' @rdname extgp-functions
#' @export
#' @keywords internal
qextgp <- function(p, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
    return(qgpd(qextgp.G(p, prob = prob, kappa = kappa, delta = delta, type = type), scale = sigma, shape = xi))
}

#' @rdname extgp-functions
#' @export
#' @keywords internal
rextgp <- function(n, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL, censoring = c(0, Inf)) {
    if (censoring[1] == 0 & censoring[2] == Inf) {
        return(qgpd(rextgp.G(n, prob = prob, kappa = kappa, delta = delta, type = type, unifsamp, censoring = c(0, 1)), scale = sigma,
            shape = xi))
    } else {
        H.L <- pgpd(censoring[1], loc = 0, scale = sigma, shape = xi)
        H.U <- pgpd(censoring[2], loc = 0, scale = sigma, shape = xi)
        return(qgpd(rextgp.G(n, prob = prob, kappa = kappa, delta = delta, type = type, unifsamp, censoring = c(H.L, H.U)), scale = sigma,
            shape = xi))
    }
}


######################### ### Inference via PWM ### ###

extgp.pwm <- function(orders, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, censoring = c(0, Inf), empiric = FALSE,
    unifsamp = NULL, NbSamples = 10^4, N = 200) {
    H.L <- ifelse(censoring[1] > 0, pgpd(censoring[1], loc = 0, scale = sigma, shape = xi), 0)
    H.U <- ifelse(censoring[2] < Inf, pgpd(censoring[2], loc = 0, scale = sigma, shape = xi), 1)

    F.L <- ifelse(censoring[1] > 0, pextgp(censoring[1], prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type),
        0)
    F.U <- ifelse(censoring[2] < Inf, pextgp(censoring[2], prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type),
        1)
    prob.LU <- F.U - F.L

    if (!empiric) {
        E2 <- ((1 - F.L)^(orders + 1) - (1 - F.U)^(orders + 1))/(prob.LU * (orders + 1))

        if (type == 0) {
            E1 <- (H.U - H.L)^(orders - xi)/(orders - xi + 1)
            return((sigma/xi) * (E1 - E2))
        } else if (type == 1) {
            beta.coefs <- beta(c(1:(max(orders) + 1)) * kappa, 1 - xi)
            probs.beta <- pbeta(H.U, shape1 = c(1:(max(orders) + 1)) * kappa, shape2 = 1 - xi) - pbeta(H.L, shape1 = c(1:(max(orders) +
                1)) * kappa, shape2 = 1 - xi)
            E1 <- c()
            for (i in 1:length(orders)) {
                E1[i] <- (kappa/prob.LU) * sum((-1)^c(0:orders[i]) * choose(orders[i], c(0:orders[i])) * beta.coefs[1:(orders[i] +
                  1)] * probs.beta[1:(orders[i] + 1)])
            }
            return((sigma/xi) * (E1 - E2))
        } else if (type == 2) {
            E1 <- c()
            for (i in 1:length(orders)) {
                expos1 <- orders[i] - xi + delta * c(0:orders[i]) + 1
                expos2 <- orders[i] - xi + delta * (c(0:orders[i]) + 1) + 1

                numer.j <- expos1 * ((1 - H.U)^expos2 - (1 - H.L)^expos2) + expos2 * ((1 - H.L)^expos1 - (1 - H.U)^expos1)
                denom.j <- expos1 * expos2

                E1[i] <- (1/(delta^(orders[i] + 1) * prob.LU)) * sum((-1)^c(0:orders[i]) * choose(orders[i], c(0:orders[i])) * (1 +
                  delta)^(orders[i] - c(0:orders[i]) + 1) * numer.j/denom.j)
            }
            return((sigma/xi) * (E1 - E2))
        } else if (type == 3) {
            print("No explicit formula for PWMs for this model...")
            return(NA)
        } else if (type == 4) {
            E1 <- c()
            for (i in 1:length(orders)) {
                l.vals <- matrix(c(0:orders[i]), nrow = orders[i] + 1, ncol = orders[i] + 1, byrow = TRUE)
                m.vals <- matrix(c(0:orders[i]), nrow = orders[i] + 1, ncol = orders[i] + 1)
                inds <- l.vals < m.vals
                l.vals[inds] <- m.vals[inds] <- NA
                const <- choose(orders[i], l.vals) * choose(l.vals, m.vals) * (-1)^l.vals * prob^m.vals * (1 - prob)^(l.vals - m.vals)
                beta.coefs1 <- beta(kappa * (m.vals + 1) + delta * (l.vals - m.vals), 1 - xi)
                beta.coefs2 <- beta(kappa * m.vals + delta * (l.vals - m.vals + 1), 1 - xi)
                probs.beta1 <- pbeta(H.U, shape1 = kappa * (m.vals + 1) + delta * (l.vals - m.vals), shape2 = 1 - xi) - pbeta(H.L,
                  shape1 = kappa * (m.vals + 1) + delta * (l.vals - m.vals), shape2 = 1 - xi)
                probs.beta2 <- pbeta(H.U, shape1 = kappa * m.vals + delta * (l.vals - m.vals + 1), shape2 = 1 - xi) - pbeta(H.L, shape1 = kappa *
                  m.vals + delta * (l.vals - m.vals + 1), shape2 = 1 - xi)
                E1[i] <- (1/prob.LU) * sum(const * (prob * kappa * beta.coefs1 * probs.beta1 + (1 - prob) * delta * beta.coefs2 *
                  probs.beta2), na.rm = TRUE)
            }
            return((sigma/xi) * (E1 - E2))
        }
    } else {
        if (is.null(unifsamp)) {
            unifsamp <- runif(NbSamples)
        }
        if (censoring[1] == 0 & censoring[2] == Inf) {
            V <- rextgp.G(length(unifsamp), prob, kappa, delta, type, unifsamp, direct = TRUE, censoring = c(0, 1))
        } else {
            V <- rextgp.G(length(unifsamp), prob, kappa, delta, type, unifsamp, direct = TRUE, censoring = c(H.L, H.U))
        }
        X <- qgpd(V, scale = sigma, shape = xi)
        res <- c()
        for (i in 1:length(orders)) {
            res[i] <- mean(X * (1 - (F.L + (F.U - F.L) * unifsamp))^orders[i], na.rm = TRUE)
        }
        return(res)
    }
}

extgp.fit.pwm <- function(x, type = 1, prob0 = NA, kappa0 = NA, delta0 = NA, sigma0 = NA, xi0 = NA, censoring = c(0, Inf), empiric = FALSE,
    unifsamp = NULL, NbSamples = 10^4) {
    Fn = ecdf(x)
    inds <- x > censoring[1] & x < censoring[2]
    mu0hat = mean(x[inds])
    mu1hat = mean(x[inds] * (1 - Fn(x[inds])))
    mu2hat = mean(x[inds] * (1 - Fn(x[inds]))^2)
    mu3hat = mean(x[inds] * (1 - Fn(x[inds]))^3)
    mu4hat = mean(x[inds] * (1 - Fn(x[inds]))^4)

    if (type == 0) {
        if (censoring[1] == 0 & censoring[2] == Inf) {
            xihat <- (mu0hat - 4 * mu1hat)/(mu0hat - 2 * mu1hat)
            sigmahat <- mu0hat * (1 - xihat)
            thetahat <- c(sigmahat, xihat)
        } else {
            fct <- function(theta, x) {
                pwm.theor <- extgp.pwm(orders = c(0:1), sigma = theta[1], xi = theta[2], type = 0, censoring = censoring, empiric = empiric,
                  unifsamp = unifsamp, NbSamples = NbSamples)
                pwm.empir <- c(mu0hat, mu1hat)
                return(matrix(pwm.theor - pwm.empir, ncol = 2))
            }
            theta0 <- c(sigma0, xi0)
            res <- gmm::gmm(fct, x, theta0, optfct = "nlminb", lower = c(1e-04, 10^(-6)), upper = c(Inf, 0.99), vcov = "iid")
            thetahat <- res$coefficients
        }
        names(thetahat) <- c("sigma", "xi")
        return(thetahat)
    } else if (type == 1) {
        fct <- function(theta, x) {
            pwm.theor <- extgp.pwm(orders = c(0:2), kappa = theta[1], sigma = theta[2], xi = theta[3], type = 1, censoring = censoring,
                empiric = empiric, unifsamp = unifsamp, NbSamples = NbSamples)
            pwm.empir <- c(mu0hat, mu1hat, mu2hat)
            return(matrix(pwm.theor - pwm.empir, ncol = 3))
        }
        theta0 <- c(kappa0, sigma0, xi0)
        res <- gmm::gmm(fct, x, theta0, optfct = "nlminb", lower = c(1e-04, 1e-04, 10^(-6)), upper = c(Inf, Inf, 0.99), vcov = "iid")
        thetahat <- res$coefficients
        names(thetahat) <- c("kappa", "sigma", "xi")
        return(thetahat)
    } else if (type == 2) {
        fct <- function(theta, x) {
            pwm.theor <- extgp.pwm(orders = c(0:2), delta = theta[1], sigma = theta[2], xi = theta[3], type = 2, censoring = censoring,
                empiric = empiric, unifsamp = unifsamp, NbSamples = NbSamples)
            pwm.empir <- c(mu0hat, mu1hat, mu2hat)
            return(matrix(pwm.theor - pwm.empir, ncol = 3))
        }
        theta0 <- c(delta0, sigma0, xi0)
        res <- gmm::gmm(fct, x, theta0, optfct = "nlminb", lower = c(1e-04, 1e-04, 10^(-6)), upper = c(100, Inf, 0.99), vcov = "iid")
        thetahat <- res$coefficients
        names(thetahat) <- c("delta", "sigma", "xi")
        return(thetahat)
    } else if (type == 3) {
        fct <- function(theta, x) {
            pwm.theor <- extgp.pwm(orders = c(0:3), kappa = theta[1], delta = theta[2], sigma = theta[3], xi = theta[4], type = 3,
                censoring = censoring, empiric = empiric, unifsamp = unifsamp, NbSamples = NbSamples)
            pwm.empir <- c(mu0hat, mu1hat, mu2hat, mu3hat)
            return(matrix(pwm.theor - pwm.empir, ncol = 4))
        }
        theta0 <- c(kappa0, delta0, sigma0, xi0)
        res <- gmm::gmm(fct, x, theta0, optfct = "nlminb", lower = c(1e-04, 1e-04, 1e-04, 10^(-6)), upper = c(Inf, 100, Inf, 0.99),
            vcov = "iid")
        thetahat <- res$coefficients
        names(thetahat) <- c("kappa", "delta", "sigma", "xi")
        return(thetahat)
    } else if (type == 4) {
        fct <- function(theta, x) {
            pwm.theor <- extgp.pwm(orders = c(0:4), prob = theta[1], kappa = theta[2], delta = theta[2] + theta[3], sigma = theta[4],
                xi = theta[5], type = 4, censoring = censoring, empiric = empiric, unifsamp = unifsamp, NbSamples = NbSamples)
            pwm.empir <- c(mu0hat, mu1hat, mu2hat, mu3hat, mu4hat)
            return(matrix(pwm.theor - pwm.empir, ncol = 5))
        }
        res0 <- extgp.fit.pwm(x[x > quantile(x, 0.9)] - quantile(x, 0.9), type = 0)
        if (is.na(prob0)) {
            prob0 <- 0.5
        }
        if (is.na(kappa0)) {
            kappa0 <- mean(x[x < quantile(x, 0.1)])/(quantile(x, 0.1) - mean(x[x < quantile(x, 0.1)]))
        }
        if (is.na(delta0)) {
            delta0 <- kappa0 + 0.01
        }
        Ddelta0 <- delta0 - kappa0
        if (is.na(sigma0)) {
            sigma0 <- res0[1]
        }
        if (is.na(xi0)) {
            xi0 <- res0[2]
        }
        theta0 <- c(prob0, kappa0, Ddelta0, sigma0, xi0)
        names(theta0) <- c("prob", "kappa", "Ddelta", "sigma", "xi")
        # print(theta0)
        res <- gmm::gmm(fct, x, theta0, optfct = "nlminb", lower = c(0, 1e-04, 0, 1e-04, 10^(-6)), upper = c(1, Inf, Inf, Inf, 0.99),
            vcov = "iid")
        thetahat <- res$coefficients
        thetahat[3] <- thetahat[2] + thetahat[3]
        names(thetahat) <- c("prob", "kappa", "delta", "sigma", "xi")
        return(thetahat)
    }
}


extgp.fit.pwm.boot <- function(data, i, type = 1, prob0 = NA, kappa0 = NA, delta0 = NA, sigma0 = NA, xi0 = NA, censoring = c(0, Inf),
    empiric = FALSE, unifsamp = NULL, NbSamples = 10^4) {
    return(extgp.fit.pwm(data[i], type = type, prob0 = prob0, kappa0 = kappa0, delta0 = delta0, sigma0 = sigma0, xi0 = xi0, censoring = censoring,
        empiric = empiric, unifsamp = unifsamp, NbSamples = NbSamples))
}


######################## ### Inference via ML ### ###

extgp.nll <- function(theta, x, censoring, rounded, type) {
    x.cens1 <- x[x < censoring[1]]
    x.cens2 <- x[x > censoring[2]]
    x.not.cens <- x[x >= censoring[1] & x <= censoring[2]]
    if (type == 0) {
        if (theta[1] <= 0 | theta[2] <= 10^(-6) | theta[2] > 0.99) {
            return(Inf)
        } else {
            censor1 <- ifelse(censoring[1] > 0, pextgp(censoring[1], sigma = theta[1], xi = theta[2], type = 0), 0)
            censor2 <- ifelse(censoring[2] < Inf, pextgp(censoring[2], sigma = theta[1], xi = theta[2], type = 0), 1)
            contrib.cens1 <- ifelse(length(x.cens1) > 0, length(x.cens1) * log(censor1), 0)
            contrib.cens2 <- ifelse(length(x.cens2) > 0, length(x.cens2) * log(1 - censor2), 0)
            contrib.not.cens <- ifelse(rounded == 0, sum(dextgp(x.not.cens, sigma = theta[1], xi = theta[2], type = 0, log = TRUE),
                na.rm = TRUE), sum(log(pextgp(x.not.cens + rounded, sigma = theta[1], xi = theta[2], type = 0) - pextgp(x.not.cens,
                sigma = theta[1], xi = theta[2], type = 0))))

            return(-(contrib.cens1 + contrib.not.cens + contrib.cens2))
        }
    } else if (type == 1) {
        if (theta[1] <= 0 | theta[2] <= 0 | theta[3] <= 10^(-6) | theta[3] > 0.99) {
            return(Inf)
        } else {
            censor1 <- ifelse(censoring[1] > 0, pextgp(censoring[1], kappa = theta[1], sigma = theta[2], xi = theta[3], type = 1),
                0)
            censor2 <- ifelse(censoring[2] < Inf, pextgp(censoring[2], kappa = theta[1], sigma = theta[2], xi = theta[3], type = 1),
                1)
            contrib.cens1 <- ifelse(length(x.cens1) > 0, length(x.cens1) * log(censor1), 0)
            contrib.cens2 <- ifelse(length(x.cens2) > 0, length(x.cens2) * log(1 - censor2), 0)
            contrib.not.cens <- ifelse(rounded == 0, sum(dextgp(x.not.cens, kappa = theta[1], sigma = theta[2], xi = theta[3], type = 1,
                log = TRUE), na.rm = TRUE), sum(log(pextgp(x.not.cens + rounded, kappa = theta[1], sigma = theta[2], xi = theta[3],
                type = 1) - pextgp(x.not.cens, kappa = theta[1], sigma = theta[2], xi = theta[3], type = 1))))

            return(-(contrib.cens1 + contrib.not.cens + contrib.cens2))
        }
    } else if (type == 2) {
        if (theta[1] <= 0 | theta[1] > 100 | theta[2] <= 0 | theta[3] <= 10^(-6) | theta[3] > 0.99) {
            return(Inf)
        } else {
            censor1 <- ifelse(censoring[1] > 0, pextgp(censoring[1], delta = theta[1], sigma = theta[2], xi = theta[3], type = 2),
                0)
            censor2 <- ifelse(censoring[2] < Inf, pextgp(censoring[2], delta = theta[1], sigma = theta[2], xi = theta[3], type = 2),
                1)
            contrib.cens1 <- ifelse(length(x.cens1) > 0, length(x.cens1) * log(censor1), 0)
            contrib.cens2 <- ifelse(length(x.cens2) > 0, length(x.cens2) * log(1 - censor2), 0)
            contrib.not.cens <- ifelse(rounded == 0, sum(dextgp(x.not.cens, delta = theta[1], sigma = theta[2], xi = theta[3], type = 2,
                log = TRUE), na.rm = TRUE), sum(log(pextgp(x.not.cens + rounded, delta = theta[1], sigma = theta[2], xi = theta[3],
                type = 2) - pextgp(x.not.cens, delta = theta[1], sigma = theta[2], xi = theta[3], type = 2))))

            return(-(contrib.cens1 + contrib.not.cens + contrib.cens2))
        }
    } else if (type == 3) {
        if (theta[1] <= 0 | theta[2] <= 0 | theta[2] > 100 | theta[3] <= 0 | theta[4] <= 10^(-6) | theta[4] > 0.99) {
            return(Inf)
        } else {
            censor1 <- ifelse(censoring[1] > 0, pextgp(censoring[1], kappa = theta[1], delta = theta[2], sigma = theta[3], xi = theta[4],
                type = 3), 0)
            censor2 <- ifelse(censoring[2] < Inf, pextgp(censoring[2], kappa = theta[1], delta = theta[2], sigma = theta[3], xi = theta[4],
                type = 3), 1)
            contrib.cens1 <- ifelse(length(x.cens1) > 0, length(x.cens1) * log(censor1), 0)
            contrib.cens2 <- ifelse(length(x.cens2) > 0, length(x.cens2) * log(1 - censor2), 0)
            contrib.not.cens <- ifelse(rounded == 0, sum(dextgp(x.not.cens, kappa = theta[1], delta = theta[2], sigma = theta[3], xi = theta[4],
                type = 3, log = TRUE), na.rm = TRUE), sum(log(pextgp(x.not.cens + rounded, kappa = theta[1], delta = theta[2], sigma = theta[3],
                xi = theta[4], type = 3) - pextgp(x.not.cens, kappa = theta[1], delta = theta[2], sigma = theta[3], xi = theta[4],
                type = 3))))

            return(-(contrib.cens1 + contrib.not.cens + contrib.cens2))
        }
    } else if (type == 4) {
        if (theta[1] < 0 | theta[1] > 1 | theta[2] <= 0 | theta[3] <= 0 | theta[4] <= 0 | theta[5] <= 10^(-6) | theta[5] > 0.99 |
            theta[3] < theta[2]) {
            return(Inf)
        } else {
            censor1 <- ifelse(censoring[1] > 0, pextgp(censoring[1], prob = theta[1], kappa = theta[2], delta = theta[3], sigma = theta[4],
                xi = theta[5], type = 4), 0)
            censor2 <- ifelse(censoring[2] < Inf, pextgp(censoring[2], prob = theta[1], kappa = theta[2], delta = theta[3], sigma = theta[4],
                xi = theta[5], type = 4), 1)
            contrib.cens1 <- ifelse(length(x.cens1) > 0, length(x.cens1) * log(censor1), 0)
            contrib.cens2 <- ifelse(length(x.cens2) > 0, length(x.cens2) * log(1 - censor2), 0)
            contrib.not.cens <- ifelse(rounded == 0, sum(dextgp(x.not.cens, prob = theta[1], kappa = theta[2], delta = theta[3], sigma = theta[4],
                xi = theta[5], type = 4, log = TRUE), na.rm = TRUE), sum(log(pextgp(x.not.cens + rounded, prob = theta[1], kappa = theta[2],
                delta = theta[3], sigma = theta[4], xi = theta[5], type = 4) - pextgp(x.not.cens, prob = theta[1], kappa = theta[2],
                delta = theta[3], sigma = theta[4], xi = theta[5], type = 4))))

            return(-(contrib.cens1 + contrib.not.cens + contrib.cens2))
        }
    }
}

extgp.fit.ml <- function(x, type = 1, prob0 = NA, kappa0 = NA, delta0 = NA, sigma0 = NA, xi0 = NA, censoring = c(0, Inf), rounded = 0.1,
    print = FALSE) {
    if (type == 0) {
        theta0 <- c(sigma0, xi0)
        opt <- optim(par = theta0, fn = extgp.nll, x = x, censoring = censoring, rounded = rounded, type = type, method = "Nelder-Mead",
            control = list(maxit = 1000), hessian = FALSE)
        if (print) {
            print(opt)
        }
        thetahat <- opt$par
        names(thetahat) <- c("sigma", "xi")
        return(thetahat)
    } else if (type == 1) {
        theta0 <- c(kappa0, sigma0, xi0)
        opt <- optim(par = theta0, fn = extgp.nll, x = x, censoring = censoring, rounded = rounded, type = type, method = "Nelder-Mead",
            control = list(maxit = 1000), hessian = FALSE)
        if (print) {
            print(opt)
        }
        thetahat <- opt$par
        names(thetahat) <- c("kappa", "sigma", "xi")
        return(thetahat)
    } else if (type == 2) {
        theta0 <- c(delta0, sigma0, xi0)
        opt <- optim(par = theta0, fn = extgp.nll, x = x, censoring = censoring, rounded = rounded, type = type, method = "Nelder-Mead",
            control = list(maxit = 1000), hessian = FALSE)
        if (print) {
            print(opt)
        }
        thetahat <- opt$par
        names(thetahat) <- c("delta", "sigma", "xi")
        return(thetahat)
    } else if (type == 3) {
        theta0 <- c(kappa0, delta0, sigma0, xi0)
        opt <- optim(par = theta0, fn = extgp.nll, x = x, censoring = censoring, rounded = rounded, type = type, method = "Nelder-Mead",
            control = list(maxit = 1000), hessian = FALSE)
        if (print) {
            print(opt)
        }
        thetahat <- opt$par
        names(thetahat) <- c("kappa", "delta", "sigma", "xi")
        return(thetahat)
    } else if (type == 4) {
        theta0 <- c(prob0, kappa0, delta0, sigma0, xi0)
        opt <- optim(par = theta0, fn = extgp.nll, x = x, censoring = censoring, rounded = rounded, type = type, method = "Nelder-Mead",
            control = list(maxit = 1000), hessian = FALSE)
        if (print) {
            print(opt)
        }
        thetahat <- opt$par
        names(thetahat) <- c("prob", "kappa", "delta", "sigma", "xi")
        return(thetahat)
    }
}

extgp.fit.ml.boot <- function(data, i, type = 1, prob0 = NA, kappa0 = NA, delta0 = NA, sigma0 = NA, xi0 = NA, censoring = c(0, Inf),
    rounded = 0.1, print = FALSE) {
    return(extgp.fit.ml(data[i], type = type, prob0 = prob0, kappa0 = kappa0, delta0 = delta0, sigma0 = sigma0, xi0 = xi0, censoring = censoring,
        rounded = rounded, print = print))
}

