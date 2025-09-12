# This function is added to ensure that the local solver can be fixed to something different
# than COBYLA - see https://github.com/jyypma/nloptr/pull/38 and
# https://github.com/stevengj/nlopt/issues/118
# .auglag <- function(x0, fn, gr = NULL, lower = NULL, upper = NULL, hin = NULL, hinjac = NULL,
#     heq = NULL, heqjac = NULL, localsolver = c("COBYLA"), localtol = 1e-06, ineq2local = FALSE,
#     nl.info = FALSE, control = list(), ...) {
#     if (ineq2local) {
#         stop("Inequalities to local solver: feature not yet implemented.")
#     }
#     localsolver <- toupper(localsolver)
#     if (localsolver %in% c("COBYLA", "BOBYQA")) {
#         # changed this line to add BOBYQA local solver
#         dfree <- TRUE
#         gsolver <- "NLOPT_LN_AUGLAG"
#         lsolver <- paste("NLOPT_LN_", localsolver, sep = "")
#     } else if (localsolver %in% c("LBFGS", "MMA", "SLSQP")) {
#         dfree <- FALSE
#         gsolver <- "NLOPT_LD_AUGLAG"
#         lsolver <- paste("NLOPT_LD_", localsolver, sep = "")
#     } else {
#         stop("Only local solvers allowed: BOBYQA, COBYLA, LBFGS, MMA, SLSQP.")
#     }
#     .fn <- match.fun(fn)
#     fn <- function(x) .fn(x, ...)
#     if (!dfree && is.null(gr)) {
#         gr <- function(x) nloptr::nl.grad(x, fn)
#     }
#     opts <- nloptr::nl.opts(control)
#     opts$algorithm <- gsolver
#     local_opts <- list(algorithm = lsolver, xtol_rel = localtol, eval_grad_f = if (!dfree) gr else NULL)
#     opts$local_opts <- local_opts
#     if (!is.null(hin)) {
#         .hin <- match.fun(hin)
#         hin <- function(x) (-1) * .hin(x)
#     }
#     if (!dfree) {
#         if (is.null(hinjac)) {
#             hinjac <- function(x) nloptr::nl.jacobian(x, hin)
#         } else {
#             .hinjac <- match.fun(hinjac)
#             hinjac <- function(x) (-1) * .hinjac(x)
#         }
#     }
#     if (!is.null(heq)) {
#         .heq <- match.fun(heq)
#         heq <- function(x) .heq(x)
#     }
#     if (!dfree) {
#         if (is.null(heqjac)) {
#             heqjac <- function(x) nloptr::nl.jacobian(x, heq)
#         } else {
#             .heqjac <- match.fun(heqjac)
#             heqjac <- function(x) .heqjac(x)
#         }
#     }
#     S0 <- nloptr::nloptr(x0, eval_f = fn, eval_grad_f = gr, lb = lower, ub = upper, eval_g_ineq = hin,
#         eval_jac_g_ineq = hinjac, eval_g_eq = heq, eval_jac_g_eq = heqjac, opts = opts)
#     if (nl.info)
#         print(S0)
#     S1 <- list(par = S0$solution, value = S0$objective, iter = S0$iterations, global_solver = gsolver,
#         local_solver = lsolver, convergence = S0$status, message = S0$message)
#     return(S1)
# }

#' @export
print.eprof <- function(x, ...) {
  confint.eprof(x, print = TRUE)
}

#' @export
summary.eprof <- function(object, ...) {
  confint.eprof(object, print = TRUE)
}

#' Confidence intervals for profile likelihood objects
#'
#' Computes confidence intervals for the parameter psi for profile likelihood objects.
#' This function uses spline interpolation to derive \code{level} confidence intervals
#'
#' @param object an object of class \code{eprof}, normally the output of \link{gpd.pll} or \link{gev.pll}.
#' @param parm a specification of which parameters are to be given confidence intervals,
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level	confidence level, with default value of 0.95
#' @param prob percentiles, with default giving symmetric 95\% confidence intervals
#' @param method string for the method, either \code{cobs} (constrained robust B-spline from eponym package) or \code{smooth.spline}
#' @param ... additional arguments passed to functions. Providing a logical \code{warn=FALSE} turns off warning messages when the lower or upper confidence interval for \code{psi} are extrapolated beyond the provided calculations.
#' @param print should a summary be printed. Default to \code{FALSE}.
#' @param boundary logical; if \code{TRUE}, the null distribution is assumed to be a mixture of a point mass and half a chi-square with one degree of freedom.
#' @return returns a 2 by 3 matrix containing point estimates, lower and upper confidence intervals based on the likelihood root and modified version thereof
#' @export
#' @importFrom Rsolnp solnp
#' @importFrom alabama auglag
#' @importFrom utils tail
confint.eprof <-
  function(
    object,
    parm,
    level = 0.95,
    prob = c((1 - level) / 2, 1 - (1 - level) / 2),
    print = FALSE,
    method = c("cobs", "smooth.spline"),
    boundary = FALSE,
    ...
  ) {
    if (!isTRUE(all.equal(diff(prob), level, check.attributes = FALSE))) {
      warning("Incompatible arguments: \"level\" does not match \"prob\".")
    }
    method <- match.arg(method[1], c("cobs", "smooth.spline"))
    args <- list(...)
    if ("warn" %in% names(args) && is.logical(args$warn)) {
      warn <- args$warn
    } else {
      warn <- TRUE
    }
    if (length(prob) != 2) {
      stop("\"prob\" must be a vector of size 2")
      prob <- sort(prob)
    }
    if (missing(parm)) {
      parm <- NULL
      ind <- args$ind
      if (!is.null(object$pll) || !is.null(object$r)) {
        parm <- c(parm, "profile")
        ind <- c(ind, 1)
      }
      if (!is.null(object$rstar)) {
        parm <- c(parm, "tem")
        ind <- c(ind, 2)
      }
      if (!is.null(object$tem.pll)) {
        parm <- c(parm, "modif.tem")
        ind <- c(ind, 3)
      }
      if (!is.null(object$empcov.pll)) {
        parm <- c(parm, "modif.empcov")
        ind <- c(ind, 4)
      }
    } else {
      if (is.numeric(parm)) {
        ind <- parm
        parm <-
          c("profile", "tem", "modif.tem", "modif.empcov")[ind]
      } else {
        parm <-
          match.arg(
            arg = parm,
            choices = c(
              "profile",
              "tem",
              "modif.tem",
              "modif.empcov",
              "r",
              "rstar"
            ),
            several.ok = TRUE
          )
        parm[parm %in% "r"] <- "profile"
        parm[parm %in% "rstar"] <- "tem"
        ind <-
          which(c("profile", "tem", "modif.tem", "modif.empcov") %in% parm)
      }
      parm <- unique(parm)
      ind <- unique(ind[ind %in% 1:4])
    }
    if (length(ind) == 0) {
      stop("Invalid \"parm\" argument.")
    }
    qulev <- qnorm(1 - prob)
    if (isTRUE(boundary[1])) {
      qulev <- sqrt(0.5) * qulev
    }
    conf <- matrix(ncol = 4, nrow = 3)
    for (i in ind) {
      if (i == 1) {
        if (is.null(object$pll) && is.null(object$r)) {
          break
        }
        if (is.null(object$r)) {
          # no r object, but must have pll + maxpll
          object$r <-
            suppressWarnings(
              sign(object$psi.max - object$psi) *
                sqrt(2 * (object$maxpll - object$pll))
            )
        } else {
          object$r[is.infinite(object$r)] <- NA
        }
        if (is.null(object$normal)) {
          object$normal <- c(object$psi.max, object$std.error)
        }
        if (
          method == "cobs" &&
            requireNamespace("cobs", quietly = TRUE)
        ) {
          fit.r <-
            cobs::cobs(
              x = object$r,
              y = object$psi,
              constraint = "decrease",
              lambda = 0,
              ic = "SIC",
              pointwise = cbind(0, 0, object$normal[1]),
              knots.add = TRUE,
              repeat.delete.add = TRUE,
              print.mesg = FALSE,
              print.warn = FALSE
            )
          pr <- predict(fit.r, c(0, qulev))[, 2]
        } else {
          fit.r <-
            stats::smooth.spline(
              x = na.omit(cbind(object$r, object$psi)),
              cv = FALSE
            )
          pr <- predict(fit.r, c(0, qulev))$y
          pr[1] <- object$normal[1]
        }
        conf[, i] <- pr
        if (warn) {
          if (!any(object$r > qnorm(prob[1]))) {
            warning(
              "Extrapolating the lower confidence interval for the profile likelihood ratio test"
            )
          }
          if (!any(object$r < qnorm(prob[2]))) {
            warning(
              "Extrapolating the upper confidence interval for the profile likelihood ratio test"
            )
          }
        }
      } else if (i == 2) {
        if (is.null(object$rstar)) {
          break
        }
        if (
          method == "cobs" &&
            requireNamespace("cobs", quietly = TRUE)
        ) {
          fit.rst <-
            cobs::cobs(
              x = object$rstar,
              y = object$psi,
              constraint = "decrease",
              lambda = 0,
              ic = "SIC",
              knots.add = TRUE,
              repeat.delete.add = TRUE,
              print.mesg = FALSE,
              print.warn = FALSE
            )
          prst <- predict(fit.rst, c(0, qulev))[, 2]
        } else {
          fit.rst <-
            stats::smooth.spline(
              x = na.omit(cbind(object$rstar, object$psi)),
              cv = FALSE
            )
          prst <- predict(fit.rst, c(0, qulev))$y
        }
        if (!is.null(object$tem.psimax)) {
          prst[1] <- object$tem.psimax
        } else {
          object$tem.psimax <- prst[1]
        }
        conf[, i] <- prst
        # lines(x=object$rstar,fit.rst$fitted,col=2,pch=19)
        if (warn) {
          if (!any(object$rstar > qulev[2])) {
            warning(
              "Extrapolating the adjusted lower confidence interval for rstar."
            )
          }
          if (!any(object$rstar < qulev[1])) {
            warning(
              "Extrapolating the adjusted upper confidence interval for rstar"
            )
          }
        }
      } else if (i == 3) {
        if (is.null(object$tem.pll)) {
          break
        }
        if (
          method == "cobs" &&
            requireNamespace("cobs", quietly = TRUE)
        ) {
          fit.mtem <-
            cobs::cobs(
              x = sign(object$tem.mle - object$psi) *
                suppressWarnings(sqrt(
                  -2 *
                    (object$tem.pll -
                      object$tem.maxpll)
                )),
              y = object$psi,
              constraint = "decrease",
              lambda = 0,
              ic = "SIC",
              knots.add = TRUE,
              repeat.delete.add = TRUE,
              print.mesg = FALSE,
              print.warn = FALSE
            )
          ptem <- predict(fit.mtem, c(0, qulev))[, 2]
        } else {
          fit.mtem <-
            stats::smooth.spline(
              x = na.omit(cbind(
                sign(object$tem.mle - object$psi) *
                  suppressWarnings(sqrt(
                    -2 * (object$tem.pll - object$tem.maxpll)
                  )),
                object$psi
              )),
              cv = FALSE
            )
          ptem <- predict(fit.mtem, c(0, qulev))$y
        }
        ptem[1] <- object$tem.mle
        conf[, i] <- ptem
        #TODO add warnings about extrapolation
      } else if (i == 4) {
        if (is.null(object$empcov.pll)) {
          break
        }
        if (
          method == "cobs" &&
            requireNamespace("cobs", quietly = TRUE)
        ) {
          fit.mempcov <-
            cobs::cobs(
              x = sign(object$empcov.mle - object$psi) *
                suppressWarnings(sqrt(
                  -2 *
                    (object$empcov.pll - object$empcov.maxpll)
                )),
              y = object$psi,
              constraint = "decrease",
              lambda = 0,
              ic = "SIC",
              knots.add = TRUE,
              repeat.delete.add = TRUE,
              print.mesg = FALSE,
              print.warn = FALSE
            )
          pempcov <- predict(fit.mempcov, c(0.5, qulev))[, 2]
        } else {
          fit.mempcov <-
            stats::smooth.spline(
              x = na.omit(cbind(
                sign(
                  object$empcov.mle -
                    object$psi
                ) *
                  suppressWarnings(sqrt(
                    -2 * (object$empcov.pll - object$empcov.maxpll)
                  )),
                object$psi
              )),
              cv = FALSE
            )
          pempcov <- predict(fit.mempcov, c(0.5, qulev))$y
        }
        pempcov[1] <- object$empcov.mle
        conf[, i] <- pempcov
      }
    }
    if (!is.null(conf)) {
      colnames(conf) <-
        c("Profile", "TEM", "Severini (TEM)", "Severini (emp. cov.)")
      rownames(conf) <- c("Estimate", "Lower CI", "Upper CI")
      #Check output makes sense - lower CI is lower than estimate
      #and similarly upper CI is above estimate
      wrong_below <- which(conf[2, ] > conf[1, ])
      if (length(wrong_below) > 0) {
        conf[2, ][wrong_below] <- NA
      }
      wrong_above <- which(conf[3, ] < conf[1, ])
      if (length(wrong_above) > 0) {
        conf[3, ][wrong_above] <- NA
      }
      if (print) {
        if (is.null(object$param)) {
          object$param <- "parameter psi"
        }
        cat(paste0("Point estimate for ", object$param, ":\n"))
        cat("Maximum likelihood          :", round(object$psi.max, 3), "\n")
        if (2 %in% ind) {
          cat(
            "Tangent exponential model   :",
            round(object$tem.psimax, 3),
            "\n"
          )
        }
        if (3 %in% ind) {
          cat("Severini's profile (TEM)    :", round(object$tem.mle, 3), "\n")
        }
        if (4 %in% ind) {
          cat(
            "Severini's profile (empcov) :",
            round(object$empcov.mle, 3),
            "\n"
          )
        }
        cat("\n")
        cat("Confidence intervals, levels :", prob, "\n")
        cat(
          "Wald intervals               :",
          round(object$psi.max + sort(qulev) * object$std.error, 3),
          "\n"
        )
        cat("Profile likelihood           :", round(conf[2:3, 1], 3), "\n")
        if (2 %in% ind) {
          cat("Tangent exponential model    :", round(conf[2:3, 2], 3), "\n")
        }
        if (3 %in% ind) {
          cat("Severini's profile (TEM)     :", round(conf[2:3, 3], 3), "\n")
        }
        if (4 %in% ind) {
          cat("Severini's profile (empcov)  :", round(conf[2:3, 4], 3), "\n")
        }
      }
      return(invisible(conf[, ind]))
    }
  }


#' Plot of (modified) profile likelihood
#'
#' The function plots the (modified) profile likelihood and the tangent exponential profile likelihood
#'
#' @param x an object of class \code{eprof} returned by \code{\link{gpd.pll}} or \code{\link{gev.pll}}.
#' @param ... further arguments to \code{plot}.
#' @return a graph of the (modified) profile likelihoods
#' @references Brazzale, A. R., Davison, A. C. and Reid, N. (2007). \emph{Applied Asymptotics: Case Studies in Small-Sample Statistics}. Cambridge University Press, Cambridge.
#' @references Severini, T. A. (2000). \emph{Likelihood Methods in Statistics}. Oxford University Press, Oxford.
#' @export
plot.eprof <- function(x, ...) {
  # plot the profile log-likelihoods
  #old.pars <- par(no.readonly = TRUE)
  args <- list(...)
  lik <- list()
  if (is.null(x$pll) && !is.null(x$r)) {
    lik$npll <- -x$r^2 / 2
  } else if (!is.null(x$pll)) {
    lik$npll <- x$pll - x$maxpll
  } else {
    stop("Invalid object provided")
  }
  if (!is.null(x$tem.pll)) {
    lik$tem.npll <- x$tem.pll - x$tem.maxpll
  }
  if (!is.null(x$empcov.pll)) {
    lik$empcov.npll <- x$empcov.pll - x$empcov.maxpll
  }

  if (is.null(args$ylim)) {
    ylim <- c(
      max(
        -8,
        min(unlist(
          lapply(lik, min, na.rm = TRUE)
        ))
      ),
      0
    )
  } else {
    ylim <- args$ylim
  }

  tikz <- FALSE
  level <- c(0.95, 0.99)
  if (!is.null(args$level)) {
    level <- args$level[1]
  }
  if (!is.null(args$tikz)) {
    if (isTRUE(args$tikz)) {
      tikz <- TRUE
    }
  }
  if (any(!is.null(args$which), !is.null(args$ind), !is.null(args$parm))) {
    if (!is.null(args$parm)) {
      parm <-
        match.arg(
          arg = args$parm,
          choices = c(
            "profile",
            "tem",
            "modif.tem",
            "modif.empcov",
            "r",
            "rstar"
          ),
          several.ok = TRUE
        )
      parm[parm %in% "r"] <- "profile"
      parm[parm %in% "rstar"] <- "tem"
      ind <-
        which(c("profile", "tem", "modif.tem", "modif.empcov") %in% parm)
    } else if (!is.null(args$ind)) {
      ind <- args$ind[args$ind %in% 1:4]
      parm <-
        c("profile", "tem", "modif.tem", "modif.empcov")[ind]
    } else if (!is.null(args$which)) {
      ind <- args$which[args$which %in% 1:4]
      parm <-
        c("profile", "tem", "modif.tem", "modif.empcov")[ind]
    }
    parm <- unique(parm)
    ind <- unique(ind[ind %in% 1:4])
  } else {
    ind <- 1:4
    parm <- c("profile", "tem", "modif.tem", "modif.empcov")
  }

  if (is.null(args$xlim)) {
    xlim <- c(min(x$psi), max(x$psi))
  } else {
    xlim <- args$xlim
  }
  if (is.null(args$xlab)) {
    xlab <- ifelse(!tikz, expression(psi), "$\\psi$")
  } else {
    xlab <- args$xlab
  }
  if (is.null(args$ylab)) {
    ylab <- "profile log likelihood"
  } else {
    ylab <- args$ylab
  }
  plot(
    NULL,
    type = "n",
    bty = "l",
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab
  )
  abline(h = -qchisq(level, 1) / 2, col = "gray")
  # Legend
  lcols <- NULL
  llty <- NULL
  llwd <- NULL
  llegend <- NULL

  if (4 %in% ind && !is.null(x$empcov.mle)) {
    abline(
      v = x$empcov.mle,
      lwd = 0.5,
      col = 4,
      lty = 4
    )
  }
  if (3 %in% ind && !is.null(x$tem.mle)) {
    abline(
      v = x$tem.mle,
      lwd = 0.5,
      col = 2,
      lty = 2
    )
  }
  if (2 %in% ind && !is.null(x$tem.psimax)) {
    abline(v = x$tem.psimax, lwd = 0.5, col = 3)
  }
  abline(v = x$psi.max, lwd = 0.5)
  if (4 %in% ind && !is.null(lik$empcov.npll)) {
    lines(
      x$psi,
      lik$empcov.npll,
      lty = 4,
      col = 4,
      lwd = 2
    )
    lcols <- c(lcols, 4)
    llty <- c(llty, 4)
    llwd <- c(llwd, 2)
    llegend <- c(llegend, "modif. emp. cov.")
  }
  if (3 %in% ind && !is.null(lik$tem.npll)) {
    lines(x$psi, lik$tem.npll, lty = 2, col = 2, lwd = 2)
    lcols <- c(lcols, 2)
    llty <- c(llty, 2)
    llwd <- c(llwd, 2)
    llegend <- c(llegend, "modif. tem.")
  }
  if (2 %in% ind && !is.null(x$rstar)) {
    lines(x$psi, -x$rstar^2 / 2, lwd = 2, col = 3, lty = 5)
    lcols <- c(lcols, 3)
    llty <- c(llty, 5)
    llwd <- c(llwd, 2)
    llegend <- c(llegend, "tem")
  }
  lcols <- c(lcols, 1)
  llty <- c(llty, 1)
  llwd <- c(llwd, 2)
  llegend <- c(llegend, "profile")
  lines(x$psi, lik$npll, lwd = 2)
  # add the legend in the top right corner
  if (!isTRUE(all.equal(llegend, "profile"))) {
    legend(
      x = "topright",
      legend = rev(llegend),
      lty = rev(llty),
      lwd = rev(llwd),
      col = rev(lcols),
      bty = "n",
      x.intersp = 0.2,
      seg.len = 0.5,
      cex = 0.9
    )
  }
  #par(old.pars)
}


#' Profile log-likelihood for the generalized extreme value distribution
#'
#' This function calculates the profile likelihood along with two small-sample corrections
#' based on Severini's (1999) empirical covariance and the Fraser and Reid tangent exponential
#' model approximation.
#'
#' @details The two additional \code{mod} available are \code{tem}, the tangent exponential model (TEM) approximation and
#' \code{modif} for the penalized profile likelihood based on \eqn{p^*} approximation proposed by Severini.
#' For the latter, the penalization is based on the TEM or an empirical covariance adjustment term.
#'
#' @param psi parameter vector over which to profile (unidimensional)
#' @param param string indicating the parameter to profile over
#' @param mod string indicating the model, one of \code{profile}, \code{tem} or \code{modif}.See \bold{Details}.
#' @param dat sample vector
#' @param N size of block over which to take maxima. Required only for \code{param} \code{Nmean} and \code{Nquant}.
#' @param p tail probability. Required only for \code{param} \code{quant}.
#' @param q probability level of quantile. Required only for \code{param} \code{Nquant}.
#' @param correction logical indicating whether to use \code{spline.corr} to smooth the tem approximation.
#' @param plot logical; should the profile likelihood be displayed? Default to \code{TRUE}
#' @param ... additional arguments such as output from call to \code{Vfun} if \code{mode='tem'}.
#'
#' @return a list with components
#' \itemize{
#' \item \code{mle}: maximum likelihood estimate
#' \item \code{psi.max}: maximum profile likelihood estimate
#' \item \code{param}: string indicating the parameter to profile over
#' \item \code{std.error}: standard error of \code{psi.max}
#' \item \code{psi}: vector of parameter \eqn{\psi} given in \code{psi}
#' \item \code{pll}: values of the profile log likelihood at \code{psi}
#' \item \code{maxpll}: value of maximum profile log likelihood
#' }
#'
#'
#' In addition, if \code{mod} includes \code{tem}
#' \itemize{
#' \item \code{normal}: maximum likelihood estimate and standard error of the interest parameter \eqn{\psi}
#' \item \code{r}: values of likelihood root corresponding to \eqn{\psi}
#' \item \code{q}: vector of likelihood modifications
#' \item \code{rstar}: modified likelihood root vector
#' \item \code{rstar.old}: uncorrected modified likelihood root vector
#' \item \code{tem.psimax}: maximum of the tangent exponential model likelihood
#' }
#' In addition, if \code{mod} includes \code{modif}
#' \itemize{
#' \item \code{tem.mle}: maximum of tangent exponential modified profile log likelihood
#' \item \code{tem.profll}: values of the modified profile log likelihood at \code{psi}
#' \item \code{tem.maxpll}: value of maximum modified profile log likelihood
#' \item \code{empcov.mle}: maximum of Severini's empirical covariance modified profile log likelihood
#' \item \code{empcov.profll}: values of the modified profile log likelihood at \code{psi}
#' \item \code{empcov.maxpll}: value of maximum modified profile log likelihood
#' }
#'
#' @references Fraser, D. A. S., Reid, N. and Wu, J. (1999), A simple general formula for tail probabilities for frequentist and Bayesian inference. \emph{Biometrika}, \bold{86}(2), 249--264.
#' @references Severini, T. (2000) Likelihood Methods in Statistics. Oxford University Press. ISBN 9780198506508.
#' @references Brazzale, A. R., Davison, A. C. and Reid, N. (2007) Applied asymptotics: case studies in small-sample statistics. Cambridge University Press, Cambridge. ISBN 978-0-521-84703-2
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' dat <- rgev(n = 100, loc = 0, scale = 2, shape = 0.3)
#' gev.pll(psi = seq(0,0.5, length = 50), param = 'shape', dat = dat)
#' gev.pll(psi = seq(-1.5, 1.5, length = 50), param = 'loc', dat = dat)
#' gev.pll(psi = seq(10, 40, length = 50), param = 'quant', dat = dat, p = 0.01)
#' gev.pll(psi = seq(12, 100, length = 50), param = 'Nmean', N = 100, dat = dat)
#' gev.pll(psi = seq(12, 90, length = 50), param = 'Nquant', N = 100, dat = dat, q = 0.5)
#' }
gev.pll <-
  function(
    psi,
    param = c("loc", "scale", "shape", "quant", "Nmean", "Nquant"),
    mod = "profile",
    dat,
    N = NULL,
    p = NULL,
    q = NULL,
    correction = TRUE,
    plot = TRUE,
    ...
  ) {
    param <- match.arg(param)
    mod <-
      match.arg(mod, c("profile", "tem", "modif"), several.ok = TRUE)
    # Parametrization profiling for quant over scale is more numerically stable
    if (param == "quant") {
      stopifnot(!is.null(p))
      q <- 1 - p
      N = 1
      param <- "Nquant"
    }
    oldpar <-
      match.arg(
        param,
        choices = c(
          "loc",
          "scale",
          "shape",
          "quant",
          "Nmean",
          "Nquant"
        ),
        several.ok = FALSE
      )
    # Arguments for parametrization of the log likelihood
    if (param %in% c("loc", "scale", "shape")) {
      args <- c("loc", "scale", "shape")
    } else if (param == "quant") {
      args <- c(param, "scale", "shape")
    } else {
      args <- c("loc", param, "shape")
    }
    # Sanity checks to ensure all arguments are provided
    if (is.null(N)) {
      if (param %in% c("Nmean", "Nquant")) {
        stop("Argument \"N\" missing. Procedure aborted")
      } else {
        N <- NA
      }
    }
    if (is.null(q)) {
      if (param == "Nquant") {
        stop("Argument \"q\" missing. Procedure aborted")
      } else {
        q <- NA
      }
    }
    if (is.null(p)) {
      if (param == "quant") {
        stop("Argument \"p\" missing. Procedure aborted")
      } else {
        p <- NA
      }
    }
    xmin <- min(dat)
    xmax <- max(dat)
    # Find maximum likelihood estimates
    mle <- gev.mle(
      xdat = dat,
      args = args,
      q = q,
      N = N,
      p = p
    )

    # if(missing(psi) || any(is.null(psi)) || any(is.na(psi))){ psi <- mle[param] } psi <- as.vector(psi)

    # Extract the components, notably V for model `tem`. Keep other components for optimization
    Vprovided <- FALSE
    extra.args <- list(...)
    if ("V" %in% names(extra.args)) {
      V <- extra.args$V
      extra.args$V <- NULL
      if (isTRUE(all.equal(dim(V), c(length(dat), 2)))) {
        Vprovided <- TRUE
      }
    }
    if (!Vprovided) {
      V <- switch(
        param,
        loc = gev.Vfun(par = mle, dat = dat),
        scale = gev.Vfun(par = mle, dat = dat),
        shape = gev.Vfun(par = mle, dat = dat),
        quant = gevr.Vfun(par = mle, dat = dat, p = p),
        Nmean = gevN.Vfun(
          par = mle,
          dat = dat,
          N = N,
          qty = "mean"
        ),
        Nquant = gevN.Vfun(
          par = mle,
          dat = dat,
          q = q,
          N = N,
          qty = "quantile"
        )
      )
    }

    # Obtained constrained maximum likelihood estimates for given value of psi
    if (param %in% c("loc", "scale", "shape")) {
      # Define observation-wise gradient
      gev.score.f <- function(par, dat) {
        dat <- as.vector(dat)
        mu = par[1]
        sigma = par[2]
        xi = as.vector(par[3])
        if (!isTRUE(all.equal(xi, 0, tolerance = 1e-8))) {
          cbind(
            -(-(mu - dat) * xi / sigma + 1)^(-1 / xi - 1) /
              sigma -
              xi *
                (1 / xi + 1) /
                (sigma *
                  ((mu - dat) * xi / sigma - 1)),
            -(dat - mu) *
              ((dat - mu) * xi / sigma + 1)^(-1 / xi - 1) /
              sigma^2 +
              (dat - mu) *
                xi *
                (1 / xi + 1) /
                (sigma^2 * ((dat - mu) * xi / sigma + 1)) -
              1 / sigma,
            -(mu - dat) *
              (1 / xi + 1) /
              (sigma * ((mu - dat) * xi / sigma - 1)) -
              (log1p(-(mu - dat) * xi / sigma) /
                xi^2 -
                (mu - dat) /
                  (sigma *
                    ((mu -
                      dat) *
                      xi /
                      sigma -
                      1) *
                    xi)) /
                (-(mu - dat) * xi / sigma + 1)^(1 / xi) +
              log1p(
                -(mu -
                  dat) *
                  xi /
                  sigma
              ) /
                xi^2
          )
        } else {
          cbind(
            -exp(mu / sigma - dat / sigma) / sigma + 1 / sigma,
            mu *
              exp(mu / sigma - dat / sigma) /
              sigma^2 -
              dat * exp(mu / sigma - dat / sigma) / sigma^2 -
              mu /
                sigma^2 -
              1 / sigma +
              dat / sigma^2,
            rep(0, length(dat))
          )
        }
      }
      ind <- switch(param, loc = 1, scale = 2, shape = 3)
      maxll <- gev.ll(mle, dat = dat)
      std.error <-
        sqrt(solve(gev.infomat(
          par = mle,
          dat = dat,
          method = "exp"
        ))[ind, ind])
      constr.mle.scale <- function(sigmat, dat = dat) {
        x0 = c(median(dat / sigmat), 0.05)
        if (is.nan(gev.ll(c(x0[1], sigmat, x0[2]), dat = dat))) {
          constr_fit <-
            try(
              mev::fit.gev(
                xdat = dat,
                fpar = list(scale = sigmat, shape = x0[2])
              ),
              silent = TRUE
            )
          if (!inherits(constr_fit, what = "try-error")) {
            if (constr_fit$convergence == "successful") {
              x0 <- as.vector(c(constr_fit$estimate["loc"], 0.05))
            } else {
              stop("Could not find starting values for optimization routine")
            }
          } else {
            stop("Could not find starting values for optimization routine")
          }
        }
        # opt <- suppressMessages(nloptr::sbplx(x0 = x0, fn = function(par) {
        #     -gev.ll(c(par[1], sigmat, par[2]), dat = dat)
        # }))
        # opt2 <- suppressMessages(nloptr::slsqp(x0 = opt$par, fn = function(par) {
        #     -gev.ll(c(par[1], sigmat, par[2]), dat = dat)
        # }, gr = function(par) {
        #     -gev.score(c(par[1], sigmat, par[2]), dat = dat)[-2]
        # }))
        opt <-
          try(
            suppressWarnings(suppressMessages(
              alabama::auglag(
                par = x0,
                fn = function(par) {
                  -gev.ll(c(par[1], sigmat, par[2]), dat = dat)
                },
                gr = function(par) {
                  -gev.score(c(par[1], sigmat, par[2]), dat = dat)[-2]
                },
                hin = function(par) {
                  ifelse(
                    par[2] <= 0,
                    sigmat + par[2] * (xmax - par[1]),
                    sigmat +
                      par[2] *
                        (xmin -
                          par[1])
                  )
                },
                control.outer = list(trace = FALSE, method = "nlminb")
              )
            )),
            silent = TRUE
          )
        if (!inherits(opt, what = "try-error")) {
          if (opt$convergence == 0 && !isTRUE(all.equal(opt$value, 1e+10))) {
            return(c(opt$par, opt$value))
          }
        }
        opt2 <- try(
          suppressWarnings(suppressMessages(Rsolnp::solnp(
            pars = x0,
            fun = function(par) {
              -gev.ll(c(par[1], sigmat, par[2]), dat = dat)
            },
            control = list(trace = 0)
          ))),
          silent = TRUE
        )
        if (!inherits(opt2, what = "try-error")) {
          if (
            opt2$convergence == 0 &&
              !isTRUE(all.equal(tail(opt2$values, 1), 1e+10))
          ) {
            return(c(opt2$par, tail(opt2$values, 1)))
          }
        }
        return(rep(NA, 3))
      }
      constr.mle.loc <- function(mut, dat = dat) {
        opt <- try(
          suppressMessages(suppressWarnings(
            alabama::auglag(
              par = c(mad(dat, constant = 1), 0.1),
              fn = function(par) {
                val <- -gev.ll(c(mut, par[1:2]), dat = dat)
                ifelse(is.infinite(val), 1e+10, val)
              },
              gr = function(par) {
                -gev.score(c(mut, par[1:2]), dat = dat)[-ind]
              },
              hin = function(par) {
                ifelse(
                  par[2] <= 0,
                  par[1] + par[2] * (xmax - mut),
                  par[1] +
                    par[2] *
                      (xmin -
                        mut)
                )
              },
              control.outer = list(method = "nlminb", trace = FALSE)
            )
          )),
          silent = TRUE
        )
        if (!inherits(opt, what = "try-error")) {
          if (opt$convergence == 0 && !isTRUE(all.equal(opt$value, 1e+10))) {
            return(c(opt$par, opt$value))
          }
        }
        opt2 <- try(
          suppressWarnings(suppressMessages(Rsolnp::solnp(
            pars = c(mad(dat, constant = 1), 0.1),
            fun = function(par) {
              -gev.ll(c(mut, par), dat = dat)
            },
            control = list(trace = 0)
          ))),
          silent = TRUE
        )
        if (!inherits(opt2, what = "try-error")) {
          if (
            opt2$convergence == 0 &&
              !isTRUE(all.equal(tail(opt2$values, 1), 1e+10))
          ) {
            return(c(opt2$par, tail(opt2$values, 1)))
          }
        }
        return(rep(NA, 3))
      }
      # constr.mle.shape <- function(xit, dat = dat){
      # start.scale <- max(1e-2,mad(dat, constant = 1)/median(abs(mev::rgev(10000, shape=xit)-mev::qgev(0.5, shape=xit))))
      # start.loc <- median(dat) - start.scale*ifelse(xit == 0, log(log(2)), (log(2)^(-xit)-1)/xit)
      # if(any(c(start.scale + xit*(xmax-start.loc) <= 0, start.scale + xit*(xmin-start.loc) <= 0)))
      # {
      #   if(xit < 0){ start.loc <- start.scale/xit + xmax + 1e-3
      # } else { start.loc <- start.scale/xit + xmin + 1e-3 } }
      # #Subplex - simplex type algorithm, more robust than Nelder-Mead
      # opt <- nloptr::sbplx(x0 = c(start.loc, start.scale),
      # fn = function(par){-gev.ll(c(par,xit), dat=dat)}, lower= c(-Inf, 0),
      # control = list(xtol_rel = 1e-8, maxeval = 2000, ftol_rel = 1e-10))
      # #If solution is not within the region
      # if(ifelse(xit < 0, opt$par[2] + xit*(xmax-opt$par[1]) <= 0, opt$par[2] + xit*(xmin-opt$par[1]) <= 0))
      # {
      #   opt <- nloptr::auglag(x0 = c(start.loc, start.scale), fn = function(par){
      #   val <- gev.ll.optim(c(par[1:2], xit), dat = dat);
      #   ifelse(is.infinite(val) || is.na(val), 1e10, val)},
      #   gr = function(par){ val <- -gev.score(c(par[1], exp(par[2]), xit), dat = dat)[-ind]},
      #   hin = function(par){c(ifelse(xit <= 0, exp(par[2]) + xit*(xmax-par[1]), exp(par[2]) + xit*(xmin-par[1])))} ) }
      #    if(!inherits(opt, what = "try-error")){
      # if(all(c(opt$convergence > 0, abs(gev.score(c(opt$par, xit), dat = dat)[1:2]) < 5e-4))){
      # return(c(opt$par, opt$value)) } else { #mev::fgev(start = list(loc = opt$par[1], scale =
      # opt$par[2]), shape = xit, # x = dat, method = 'BFGS', control=list(reltol=1e-10, abstol =
      # 1e-9)) opt2 <- suppressWarnings(nloptr::slsqp(x0 = opt$par, fn = function(par){ val <-
      # gev.ll.optim(c(par[1:2], xit), dat = dat); ifelse(is.infinite(val) || is.na(val), 1e10,
      # val)}, gr = function(par){ val <- -gev.score(c(par[1], exp(par[2]), xit), dat =
      # dat)[-ind]}, hin = function(par){c(ifelse(xit <= 0, exp(par[2]) + xit*(xmax-par[1]),
      # exp(par[2]) + xit*(xmin-par[1])))} )) opt2$par[2] <- exp(opt2$par[2])
      # if(all(c(opt2$convergence > 0, !isTRUE(all.equal(opt2$value, 1e10)),
      # abs(gev.score(c(opt2$par, xit), dat = dat)[1:2]) < 1e-4))){ return(c(opt2$par,
      # opt2$value)) } } } return(rep(NA, 3)) }
      constr.mle.shape <- function(xit, dat = dat) {
        if (abs(xit) < 1e-08) {
          xit <- 0
        } #because rgev does not handle this case!
        start.scale <-
          mad(dat, constant = 1) /
          median(abs(
            mev::rgev(2000, shape = xit) -
              mev::qgev(0.5, shape = xit)
          ))
        start.loc <-
          median(dat) -
          start.scale *
            ifelse(
              xit == 0,
              log(log(2)),
              (log(2)^(-xit) -
                1) /
                xit
            )
        if (start.scale + xit * (xmax - start.loc) <= 0) {
          if (xit < 0) {
            start.loc <- start.scale / xit + xmax + 0.001
          } else {
            start.loc <- start.scale / xit + xmin + 0.001
          }
        }
        opt <- try(
          suppressMessages(suppressWarnings(
            alabama::auglag(
              par = c(start.loc, start.scale),
              fn = function(par) {
                val <- -gev.ll(c(par[1:2], xit), dat = dat)
                ifelse(is.infinite(val) || is.na(val), 1e+10, val)
              },
              gr = function(par) {
                -gev.score(c(par[1:2], xit), dat = dat)[-ind]
              },
              hin = function(par) {
                ifelse(
                  xit <= 0,
                  par[2] + xit * (xmax - par[1]),
                  par[2] + xit * (xmin - par[1])
                )
              },
              control.outer = list(method = "nlminb", trace = FALSE)
            )
          )),
          silent = TRUE
        )
        if (!inherits(opt, what = "try-error")) {
          if (opt$convergence == 0 && !isTRUE(all.equal(opt$value, 1e+10))) {
            return(c(opt$par, opt$value))
          }
        }
        opt2 <- try(
          suppressMessages(suppressWarnings(
            Rsolnp::solnp(
              pars = c(start.loc, start.scale),
              fun = function(par) {
                val <- -gev.ll(c(par[1:2], xit), dat = dat)
                ifelse(
                  is.infinite(val) ||
                    is.na(val),
                  1e+10,
                  val
                )
              },
              ineqfun = function(par) {
                ifelse(
                  xit <= 0,
                  par[2] + xit * (xmax - par[1]),
                  par[2] +
                    xit *
                      (xmin -
                        par[1])
                )
              },
              ineqLB = 0,
              ineqUB = Inf,
              control = list(trace = 0)
            )
          )),
          silent = TRUE
        )
        if (!inherits(opt2, what = "try-error")) {
          if (
            isTRUE(opt2$convergence == 0) &&
              !isTRUE(all.equal(tail(opt2$values, 1), 1e+10))
          ) {
            return(c(opt2$par, tail(opt2$values, 1)))
          }
        }
        return(rep(NA, 3))
      }
      # Missing psi vector
      if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
        if (ind == 2) {
          psirangelow <-
            unique(min(
              1e-05,
              seq(-4, -1.5, length = 10) * std.error + mle[param]
            ))
        } else {
          psirangelow <- seq(-4, -1.5, length = 10) * std.error + mle[param]
        }
        lowvals <- -t(sapply(psirangelow, function(par) {
          switch(
            param,
            loc = constr.mle.loc(par, dat = dat),
            scale = constr.mle.scale(par, dat = dat),
            shape = constr.mle.shape(par, dat = dat)
          )
        }))[, 3] -
          maxll
        psirangehigh <-
          seq(1.5, 4, length = 10) * std.error + mle[param]
        highvals <- -t(sapply(psirangehigh, function(par) {
          switch(
            param,
            loc = constr.mle.loc(par, dat = dat),
            scale = constr.mle.scale(par, dat = dat),
            shape = constr.mle.shape(par, dat = dat)
          )
        }))[, 3] -
          maxll
        lo <- approx(x = lowvals, y = psirangelow, xout = -4)$y
        hi <-
          approx(x = highvals, y = psirangehigh, xout = -4)$y
        lo <-
          ifelse(
            is.na(lo),
            predict(
              lm(psirangelow ~ lowvals),
              newdata = data.frame(lowvals = -4)
            )[[1]],
            lo
          )
        hi <-
          ifelse(
            is.na(hi),
            predict(
              lm(psirangehigh ~ highvals),
              newdata = data.frame(highvals = -4)
            )[[1]],
            hi
          )
        psi <- seq(lo, hi, length = 55)
      }
      if (param == "shape" && min(psi) < -1) {
        warning(
          "psi includes shape parameters below -1. These will be automatically removed."
        )
        psi <- psi[psi > -1]
      }
      # Calculate profile likelihood at psi
      lambda <- t(sapply(psi, function(par) {
        switch(
          param,
          loc = constr.mle.loc(par, dat = dat),
          scale = constr.mle.scale(par, dat = dat),
          shape = constr.mle.shape(par, dat = dat)
        )
      }))
      pars <-
        switch(
          param,
          loc = cbind(psi, lambda[, 1:2, drop = FALSE]),
          scale = cbind(
            lambda[,
              1,
              drop = FALSE
            ],
            psi,
            lambda[, 2, drop = FALSE]
          ),
          shape = cbind(lambda[, 1:2, drop = FALSE], psi)
        )
      # Profile log likelihood values for psi
      profll <- -lambda[, 3]
      r <- sign(mle[param] - psi) * sqrt(2 * (maxll - profll))
      if ("tem" %in% mod) {
        # Tangent exponential model approximation of Fraser and Reid to the profile likelihood
        phi.mle <- gev.phi(par = mle, dat = dat, V = V)
        q2num <-
          ifelse(ind %% 2 == 0, -1, 1) *
          apply(pars, 1, function(par) {
            det(rbind(
              c(
                c(phi.mle) -
                  gev.phi(
                    par = par,
                    dat = dat,
                    V = V
                  )
              ),
              gev.dphi(
                par = par,
                dat = dat,
                V = V
              )[-ind, ]
            ))
          })
        if (isTRUE(any(sign(q2num) * sign(r) < 0, na.rm = TRUE))) {
          warning(
            "Correction factor and likelihood root are of opposite sign - check output"
          )
        }
        logq <- apply(pars, 1, function(par) {
          -0.5 *
            log(det(gev.infomat(
              par = par,
              dat = dat,
              method = "obs"
            )[-ind, -ind]))
        }) +
          log(abs(q2num))
        qmlecontrib <-
          -log(det(gev.dphi(
            par = mle,
            dat = dat,
            V = V
          ))) +
          0.5 *
            log(det(gev.infomat(
              par = mle,
              dat = dat,
              method = "obs"
            )))
        logq <- logq + qmlecontrib
        qcor <- sign(q2num) * exp(logq)
        rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r))) / r)
        ###

        tem.max.opt <- function(psi, dat = dat) {
          lambda <-
            switch(
              param,
              loc = constr.mle.loc(psi, dat = dat),
              scale = constr.mle.scale(psi, dat = dat),
              shape = constr.mle.shape(psi, dat = dat)
            )[1:2]
          para <-
            switch(
              ind,
              c(psi, lambda),
              c(lambda[1], psi, lambda[2]),
              c(lambda, psi)
            )
          pll <- gev.ll(par = para, dat = dat)
          rs <- 2 * (maxll - pll)
          logq <-
            -0.5 *
            log(det(gev.infomat(
              par = para,
              dat = dat,
              method = "obs"
            )[-ind, -ind])) +
            qmlecontrib +
            log(abs(det(rbind(
              c(
                c(phi.mle) -
                  gev.phi(
                    par = para,
                    dat = dat,
                    V = V
                  )
              ),
              gev.dphi(
                par = para,
                dat = dat,
                V = V
              )[-ind, ]
            ))))
          abs(rs + logq - log(sqrt(abs(rs))))
          # rs is r squared - finding the maximum by multiplying rstar by r
          # when rstar vanishes at zero, so does rstar*r
          # this ensures that we do not divide by r = zero
        }
        tem.max <-
          optim(
            par = mle[ind],
            fn = tem.max.opt,
            method = "Brent",
            dat = dat,
            lower = ifelse(
              ind == 2,
              max(1e-05, mle[ind] - std.error),
              mle[ind] - std.error
            ),
            upper = mle[ind] + std.error,
            control = list(abstol = 1e-10)
          )$par
        ###
      }
      # Modified profile likelihood based on p* approximation, two modifications due to Severini
      if ("modif" %in% mod) {
        tem.objfunc <- function(par, dat = dat) {
          0.5 *
            log(det(gev.infomat(
              par,
              dat = dat,
              method = "obs"
            )[-ind, -ind])) -
            log(abs(det(
              gev.dphi(
                par = par,
                dat = dat,
                V = V[, -ind]
              )[-ind, ]
            )))
        }
        optim.tem.fn <-
          function(psi, dat = dat, param = param) {
            theta.psi.opt <-
              switch(
                param,
                loc = constr.mle.loc(psi, dat = dat),
                scale = constr.mle.scale(psi, dat = dat),
                shape = constr.mle.shape(psi, dat = dat)
              )
            para <-
              switch(
                param,
                loc = c(psi, theta.psi.opt[1:2]),
                scale = c(theta.psi.opt[1], psi, theta.psi.opt[2]),
                shape = c(theta.psi.opt[1:2], psi)
              )
            ll <- -theta.psi.opt[3]
            ll + tem.objfunc(para, dat = dat)
          }
        # TEM profile log likelihood values for psi
        proflltem <-
          profll + apply(pars, 1, tem.objfunc, dat = dat)
        # Maximum objective function for TEM (line search in neighborhood of the MLE)
        tem.mle.opt <-
          optim(
            par = mle[ind],
            fn = optim.tem.fn,
            method = "Brent",
            dat = dat,
            param = param,
            lower = ifelse(
              param == "scale",
              max(1e-05, mle[ind] - std.error),
              mle[ind] - std.error
            ),
            upper = mle[ind] + std.error,
            control = list(fnscale = -1)
          )

        # Severini empirical covariance function adjustment to profile likelihood Score function -
        # observation-wise

        # Score at MLE (sums to zero)
        score.scale.mle <-
          gev.score.f(mle, dat)[, -ind] #keep s_lambda
        empcov.objfunc <- function(par, dat) {
          0.5 *
            log(det(gev.infomat(
              par = par,
              dat = dat,
              method = "obs"
            )[-ind, -ind])) -
            log(abs(sum(
              score.scale.mle * gev.score.f(par, dat)[, -ind]
            )))
        }
        profllempcov <-
          profll + apply(pars, 1, empcov.objfunc, dat = dat)
        optim.empcov.fn <-
          function(psi, param = param, dat = dat) {
            theta.psi.opt <-
              switch(
                param,
                loc = constr.mle.loc(psi, dat = dat),
                scale = constr.mle.scale(psi, dat = dat),
                shape = constr.mle.shape(psi, dat = dat)
              )
            para <-
              switch(
                param,
                loc = c(psi, theta.psi.opt[1:2]),
                scale = c(theta.psi.opt[1], psi, theta.psi.opt[2]),
                shape = c(theta.psi.opt[1:2], psi)
              )
            ll <- -theta.psi.opt[3]
            ll + empcov.objfunc(para, dat = dat)
          }

        empcov.mle.opt <-
          optim(
            par = mle[ind],
            fn = optim.empcov.fn,
            method = "Brent",
            dat = dat,
            param = param,
            lower = ifelse(
              param == "scale",
              max(
                1e-05,
                mle[ind] -
                  std.error
              ),
              mle[ind] - std.error
            ),
            upper = mle[ind] + std.error,
            control = list(fnscale = -1)
          )
      }

      # Return levels, quantiles or value-at-risk
    } else if (param %in% c("Nquant", "Nmean")) {
      qty <- switch(param, Nquant = "quantile", Nmean = "mean")
      maxll <- gevN.ll(
        mle,
        dat = dat,
        N = N,
        q = q,
        qty = qty
      )
      std.error <-
        sqrt(solve(
          gevN.infomat(
            par = mle,
            dat = dat,
            method = "exp",
            N = N,
            q = q,
            qty = qty
          )
        )[2, 2])
      constr.mle.N <- function(zt, dat = dat) {
        st_vals <- c(median(dat), 0.1)
        if (isTRUE(as.vector(mle["shape"] > 0))) {
          st_vals <- mle[c("loc", "shape")]
        }
        # COBYLA unfortunately hangs from time to time.  so it cannot be used in optimization
        # routine without risking hanging despite the algorithm begin more robust and faster for
        # this problem...
        opt <-
          try(
            suppressMessages(suppressWarnings(
              alabama::auglag(
                par = st_vals,
                fn = function(par) {
                  val <-
                    -gevN.ll(
                      par = c(par[1], zt, par[2]),
                      dat = dat,
                      q = q,
                      qty = qty,
                      N = N
                    )
                  ifelse(
                    isTRUE(any(
                      is.infinite(val),
                      is.na(val),
                      val <= -maxll
                    )),
                    1e+10,
                    val
                  )
                },
                hin = function(par) {
                  sigma <-
                    switch(
                      qty,
                      quantile = (zt - par[1]) *
                        par[2] /
                        (N^par[2] * (log(1 / q))^(-par[2]) - 1),
                      mean = (zt - par[1]) *
                        par[2] /
                        (N^par[2] * gamma(1 - par[2]) - 1)
                    )
                  c(
                    sigma,
                    sigma + par[2] * (xmax - par[1]),
                    sigma + par[2] * (xmin - par[1])
                  )
                },
                control.outer = list(method = "nlminb", trace = FALSE)
              )
            )),
            silent = TRUE
          )
        #With `alabama` package opt <- try(suppressWarnings( alabama::auglag(par = c(median(dat),
        # 0.1), fn = function(par){ val <- -gevN.ll(par = c(par[1], zt, par[2]), dat = dat, q = q,
        # qty = qty, N = N); ifelse(is.infinite(val) || is.na(val), 1e10, val)}, gr =
        # function(par){-gevN.score(par = c(par[1], zt, par[2]), dat = dat, q = q, qty = qty, N =
        # N)[-2]}, hin = function(par){ sigma <- switch(qty, quantile =
        # (zt-par[1])*par[2]/(N^par[2]*(log(1/q))^(-par[2])-1), mean =
        # (zt-par[1])*par[2]/(N^par[2]*gamma(1-par[2])-1)) c(sigma, sigma + par[2]*(xmax-par[1]),
        # sigma + par[2]*(xmin-par[1]))}, control.outer = list (trace = FALSE) )))
        if (!inherits(opt, what = "try-error")) {
          # if(opt$convergence == 0 && !isTRUE(all.equal(opt$value, 1e10))){
          if (
            isTRUE(opt$convergence == 0) &&
              !isTRUE(all.equal(opt$value, 1e+10))
          ) {
            return(c(opt$par, opt$value))
          } else {
            st_vals <- opt$par
          }
        }
        opt2 <- try(
          suppressMessages(suppressWarnings(
            Rsolnp::solnp(
              par = st_vals,
              fun = function(par) {
                val <-
                  -gevN.ll(
                    par = c(par[1], zt, par[2]),
                    dat = dat,
                    q = q,
                    qty = qty,
                    N = N
                  )
                ifelse(
                  is.infinite(val) ||
                    is.na(val),
                  1e+10,
                  val
                )
              },
              ineqfun = function(par) {
                sigma <-
                  switch(
                    qty,
                    quantile = (zt - par[1]) *
                      par[2] /
                      (N^par[2] * (log(1 / q))^(-par[2]) - 1),
                    mean = (zt - par[1]) *
                      par[2] /
                      (N^par[2] * gamma(1 - par[2]) - 1)
                  )
                c(
                  sigma,
                  sigma + par[2] * (xmax - par[1]),
                  sigma + par[2] * (xmin - par[1])
                )
              },
              ineqLB = rep(0, 3),
              ineqUB = rep(Inf, 3),
              control = list(trace = 0)
            )
          )),
          silent = TRUE
        )
        if (!inherits(opt2, what = "try-error")) {
          if (
            isTRUE(opt2$convergence == 0) &&
              !isTRUE(all.equal(tail(opt2$values, 1), 1e+10))
          ) {
            return(c(opt2$par, tail(opt2$values, 1)))
          }
        }
        return(rep(NA, 3))
      }
      # Missing psi vector
      if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
        psirangelow <-
          pmax(
            1e-05,
            seq(ifelse(mle[3] < 0, -3.5, -2.5), -0.5, length = 6) *
              std.error +
              mle[2]
          )
        lowvals <-
          -sapply(psirangelow, constr.mle.N, dat = dat)[3, ] - maxll
        psirangehigh <-
          seq(2.5, (1 + mle[3]) * 8, length = 10) * std.error + mle[2]
        highvals <-
          -sapply(psirangehigh, constr.mle.N, dat = dat)[3, ] - maxll

        lo <- approx(x = lowvals, y = psirangelow, xout = -4)$y
        hi <-
          approx(x = highvals, y = psirangehigh, xout = -4)$y
        lo <-
          ifelse(is.na(lo), lm(psirangelow ~ lowvals)$coef[2] * -4 + mle[2], lo)
        hi <-
          ifelse(
            is.na(hi),
            lm(psirangehigh ~ highvals)$coef[2] * -4 + mle[2],
            hi
          )
        psi <-
          c(seq(lo, mle[2], length = 20)[-20], seq(mle[2], hi, length = 55)[-1])
      }

      lambda <- t(sapply(psi, constr.mle.N, dat = dat))
      pars <-
        cbind(lambda[, 1, drop = FALSE], psi, lambda[, 2, drop = FALSE])
      # Profile log likelihood values for psi
      profll <- -lambda[, 3]
      r <- sign(mle[param] - psi) * sqrt(2 * (maxll - profll))
      if ("tem" %in% mod) {
        phi.mle <-
          gevN.phi(
            par = mle,
            dat = dat,
            V = V,
            N = N,
            q = q,
            qty = qty
          )
        q2num <- apply(pars, 1, function(par) {
          -det(rbind(
            c(
              c(phi.mle) -
                gevN.phi(
                  par = par,
                  dat = dat,
                  V = V,
                  N = N,
                  q = q,
                  qty = qty
                )
            ),
            gevN.dphi(
              par = par,
              dat = dat,
              V = V,
              N = N,
              q = q,
              qty = qty
            )[-2, ]
          ))
        })
        if (isTRUE(any(sign(q2num) * sign(r) < 0, na.rm = TRUE))) {
          warning(
            "Correction factor and likelihood root are of opposite sign - check output"
          )
        }

        logq <- apply(pars, 1, function(par) {
          -0.5 *
            log(det(
              gevN.infomat(
                par = par,
                dat = dat,
                method = "obs",
                N = N,
                q = q,
                qty = qty
              )[-2, -2]
            ))
        }) +
          log(abs(q2num))
        qmlecontrib <-
          -log(det(gevN.dphi(
            par = mle,
            dat = dat,
            V = V,
            N = N,
            q = q,
            qty = qty
          ))) +
          0.5 *
            log(det(
              gevN.infomat(
                par = mle,
                dat = dat,
                method = "obs",
                N = N,
                q = q,
                qty = qty
              )
            ))
        logq <- logq + qmlecontrib
        qcor <- sign(q2num) * exp(logq)
        rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r))) / r)

        tem.max.opt <- function(psi, dat = dat) {
          lam <- constr.mle.N(psi, dat = dat)
          para <- c(lam[1], psi, lam[2])
          pll <- -lam[3]
          rs <- 2 * (maxll - pll)
          logq <-
            -0.5 *
            log(det(
              gevN.infomat(
                par = para,
                dat = dat,
                method = "obs",
                N = N,
                q = q,
                qty = qty
              )[-2, -2]
            )) +
            qmlecontrib +
            log(abs(det(rbind(
              c(
                c(phi.mle) -
                  gevN.phi(
                    par = para,
                    dat = dat,
                    V = V,
                    q = q,
                    N = N,
                    qty = qty
                  )
              ),
              gevN.dphi(
                par = para,
                dat = dat,
                V = V,
                q = q,
                N = N,
                qty = qty
              )[-2, ]
            ))))
          abs(rs + logq - log(sqrt(abs(rs))))
        }
        tem.max <-
          optim(
            par = mle[2],
            fn = tem.max.opt,
            method = "Brent",
            dat = dat,
            lower = max(1e-05, mle[2] - std.error),
            upper = mle[2] + std.error,
            control = list(abstol = 1e-10)
          )$par
        names(mle)[2] <- oldpar
      }
      if ("modif" %in% mod) {
        # Tangent exponential model approximation of Fraser and Reid to the profile likelihood
        tem.objfunc.N <-
          function(par, dat = dat, N = N, qty = qty, q = q) {
            0.5 *
              log(det(
                gevN.infomat(
                  par = par,
                  dat = dat,
                  method = "obs",
                  N = N,
                  qty = qty,
                  q = q
                )[-2, -2]
              )) -
              log(abs(det(
                gevN.dphi(
                  par = par,
                  dat = dat,
                  N = N,
                  V = V[,
                    -2
                  ],
                  qty = qty,
                  q = q
                )[-2, ]
              )))
          }
        optim.tem.fn.N <-
          function(psi, dat = dat, q = q, qty = qty, N = N) {
            theta.psi.opt <- constr.mle.N(psi, dat = dat)
            param <- c(theta.psi.opt[1], psi, theta.psi.opt[2])
            ll <-
              gevN.ll(
                param,
                dat = dat,
                N = N,
                q = q,
                qty = qty
              )
            ll +
              tem.objfunc.N(
                param,
                dat = dat,
                N = N,
                qty = qty,
                q = q
              )
          }
        # TEM profile log likelihood values for psi
        proflltem <-
          profll +
          suppressWarnings(apply(
            pars,
            1,
            tem.objfunc.N,
            dat = dat,
            N = N,
            qty = qty,
            q = q
          ))
        # Maximum objective function for TEM
        tem.mle.opt <-
          optim(
            par = mle[2],
            fn = optim.tem.fn.N,
            method = "Brent",
            dat = dat,
            q = q,
            qty = qty,
            N = N,
            lower = max(min(dat), mle[2] - std.error),
            upper = mle[2] +
              std.error,
            control = list(fnscale = -1)
          )
        # Severini empirical covariance function adjustment to profile likelihood Score function -
        # observation-wise for xi
        gevN.score.f <-
          function(par, dat, N, q = q, qty = qty) {
            qty <- match.arg(qty, c("mean", "quantile"))
            mu = par[1]
            z = par[2]
            xi = par[3]
            Npxi <- N^xi
            log1q <- log(1 / q)
            logN <- log(N)
            if (qty == "quantile") {
              # quantiles at prob. q
              cbind(
                ((-(Npxi * log1q^(-xi) - 1) * (dat - mu) / (mu - z) + 1)^(-1 /
                  xi -
                  1) *
                  ((Npxi * log1q^(-xi) - 1) *
                    (dat - mu) /
                    (mu - z)^2 +
                    (Npxi * log1q^(-xi) - 1) / (mu - z)) /
                  xi +
                  ((Npxi * log1q^(-xi) - 1) *
                    (dat - mu) /
                    (mu - z)^2 +
                    (Npxi * log1q^(-xi) - 1) / (mu - z)) *
                    (1 / xi + 1) /
                    ((Npxi * log1q^(-xi) - 1) * (dat - mu) / (mu - z) - 1) -
                  1 / (mu - z)),
                (-(Npxi * log1q^(-xi) - 1) *
                  (dat - mu) *
                  (-(Npxi * log1q^(-xi) - 1) * (dat - mu) / (mu - z) + 1)^(-1 /
                    xi -
                    1) /
                  ((mu - z)^2 * xi) -
                  (Npxi * log1q^(-xi) - 1) *
                    (dat -
                      mu) *
                    (1 / xi + 1) /
                    (((Npxi * log1q^(-xi) - 1) * (dat - mu) / (mu - z) - 1) *
                      (mu - z)^2) +
                  1 / (mu - z)),
                ((-(Npxi * log1q^(-xi) - 1) *
                  (dat -
                    mu) /
                  (mu - z) +
                  1)^(-1 / xi) *
                  ((Npxi *
                    log1q^(-xi) *
                    logN -
                    Npxi * log1q^(-xi) * log(log1q)) *
                    (dat - mu) /
                    (((Npxi * log1q^(-xi) - 1) *
                      (dat - mu) /
                      (mu -
                        z) -
                      1) *
                      (mu - z) *
                      xi) -
                    log1p(
                      -(Npxi * log1q^(-xi) - 1) *
                        (dat - mu) /
                        (mu -
                          z)
                    ) /
                      xi^2) -
                  (Npxi *
                    log1q^(-xi) *
                    logN -
                    Npxi * log1q^(-xi) * log(log1q)) *
                    (dat - mu) *
                    (1 / xi + 1) /
                    (((Npxi * log1q^(-xi) - 1) *
                      (dat - mu) /
                      (mu - z) -
                      1) *
                      (mu - z)) +
                  (Npxi * log1q^(-xi) - 1) *
                    ((Npxi *
                      log1q^(-xi) *
                      logN -
                      Npxi * log1q^(-xi) * log(log1q)) *
                      (mu -
                        z) *
                      xi /
                      (Npxi * log1q^(-xi) - 1)^2 -
                      (mu - z) /
                        (Npxi * log1q^(-xi) - 1)) /
                    ((mu - z) * xi) +
                  log1p(
                    -(Npxi * log1q^(-xi) - 1) *
                      (dat - mu) /
                      (mu -
                        z)
                  ) /
                    xi^2)
              )
            } else {
              # Mean
              cbind(
                ((-(Npxi * gamma(-xi + 1) - 1) *
                  (dat - mu) /
                  (mu - z) +
                  1)^(-1 / xi - 1) *
                  ((Npxi * gamma(-xi + 1) - 1) *
                    (dat - mu) /
                    (mu - z)^2 +
                    (Npxi *
                      gamma(
                        -xi +
                          1
                      ) -
                      1) /
                      (mu - z)) /
                  xi +
                  ((Npxi * gamma(-xi + 1) - 1) *
                    (dat - mu) /
                    (mu -
                      z)^2 +
                    (Npxi * gamma(-xi + 1) - 1) / (mu - z)) *
                    (1 / xi + 1) /
                    ((Npxi *
                      gamma(
                        -xi +
                          1
                      ) -
                      1) *
                      (dat - mu) /
                      (mu - z) -
                      1) -
                  1 / (mu - z)),
                (-(Npxi *
                  gamma(
                    -xi +
                      1
                  ) -
                  1) *
                  (dat - mu) *
                  (-(Npxi * gamma(-xi + 1) - 1) *
                    (dat - mu) /
                    (mu -
                      z) +
                    1)^(-1 / xi - 1) /
                  ((mu - z)^2 * xi) -
                  (Npxi * gamma(-xi + 1) - 1) *
                    (dat -
                      mu) *
                    (1 / xi + 1) /
                    (((Npxi * gamma(-xi + 1) - 1) * (dat - mu) / (mu - z) - 1) *
                      (mu - z)^2) +
                  1 / (mu - z)),
                ((-(Npxi * gamma(-xi + 1) - 1) *
                  (dat -
                    mu) /
                  (mu - z) +
                  1)^(-1 / xi) *
                  ((Npxi *
                    logN *
                    gamma(-xi + 1) -
                    Npxi *
                      psigamma(
                        -xi +
                          1
                      ) *
                      gamma(-xi + 1)) *
                    (dat - mu) /
                    (((Npxi * gamma(-xi + 1) - 1) *
                      (dat -
                        mu) /
                      (mu - z) -
                      1) *
                      (mu - z) *
                      xi) -
                    log1p(
                      -(Npxi * gamma(-xi + 1) - 1) *
                        (dat - mu) /
                        (mu - z)
                    ) /
                      xi^2) -
                  (Npxi *
                    logN *
                    gamma(-xi + 1) -
                    Npxi *
                      psigamma(-xi + 1) *
                      gamma(-xi + 1)) *
                    (dat - mu) *
                    (1 / xi + 1) /
                    (((Npxi *
                      gamma(-xi + 1) -
                      1) *
                      (dat - mu) /
                      (mu - z) -
                      1) *
                      (mu - z)) +
                  (Npxi *
                    gamma(
                      -xi +
                        1
                    ) -
                    1) *
                    ((Npxi *
                      logN *
                      gamma(-xi + 1) -
                      Npxi * psigamma(-xi + 1) * gamma(-xi + 1)) *
                      (mu - z) *
                      xi /
                      (Npxi * gamma(-xi + 1) - 1)^2 -
                      (mu - z) /
                        (Npxi *
                          gamma(-xi + 1) -
                          1)) /
                    ((mu - z) * xi) +
                  log1p(
                    -(Npxi * gamma(-xi + 1) - 1) *
                      (dat - mu) /
                      (mu - z)
                  ) /
                    xi^2)
              )
            }
          }
        # Score at MLE (sums to zero)
        score.N.mle <-
          gevN.score.f(mle, dat, N, qty = qty, q = q)
        empcov.objfunc.N <-
          function(par, dat = dat, q = q, qty = qty, N = N) {
            0.5 *
              log(det(
                gevN.infomat(
                  par = par,
                  dat = dat,
                  method = "obs",
                  N = N,
                  q = q,
                  qty = qty
                )[-2, -2]
              )) -
              log(abs(sum(
                score.N.mle *
                  gevN.score.f(
                    par,
                    dat,
                    N = N,
                    qty = qty,
                    q = q
                  )
              )))
          }
        profllempcov <-
          profll +
          suppressWarnings(apply(
            pars,
            1,
            empcov.objfunc.N,
            N = N,
            q = q,
            dat = dat,
            qty = qty
          ))
        optim.empcov.fn.N <-
          function(psi, dat = dat, q = q, qty = qty, N = N) {
            theta.psi.opt <- constr.mle.N(psi, dat = dat)
            param <- c(theta.psi.opt[1], psi, theta.psi.opt[2])
            ll <-
              gevN.ll(
                param,
                dat = dat,
                N = N,
                q = q,
                qty = qty
              )
            ll +
              empcov.objfunc.N(
                param,
                dat = dat,
                q = q,
                qty = qty,
                N = N
              )
          }
        empcov.mle.opt <-
          optim(
            par = mle[2],
            fn = optim.empcov.fn.N,
            method = "Brent",
            dat = dat,
            qty = qty,
            q = q,
            N = N,
            lower = max(min(dat), mle[2] - std.error),
            upper = mle[2] + std.error,
            control = list(fnscale = -1)
          )
      }
    }
    # Return profile likelihood and quantities of interest (modified likelihoods)
    colnames(pars) <- names(mle)
    ans <-
      list(
        mle = mle,
        pars = pars,
        psi.max = as.vector(mle[oldpar]),
        param = oldpar,
        std.error = std.error,
        psi = psi,
        pll = profll,
        maxpll = maxll,
        r = r
      )
    if ("tem" %in% mod) {
      ans$q <- qcor
      ans$rstar <- rstar
      ans$normal <- c(ans$psi.max, ans$std.error)
      if (correction && length(psi) > 10) {
        ans <- spline.corr(ans)
      }
      ans$tem.psimax <- tem.max
    }
    if ("modif" %in% mod) {
      ans$tem.mle <- tem.mle.opt$par
      ans$tem.pll <- proflltem
      ans$tem.maxpll <- tem.mle.opt$value
      ans$empcov.mle <- empcov.mle.opt$par
      ans$empcov.pll <- profllempcov
      ans$empcov.maxpll <- empcov.mle.opt$value
    }

    if ("tem" %in% mod) {
      class(ans) <- c("eprof", "fr")
    } else {
      class(ans) <- "eprof"
    }
    ans$family <- "gev"
    if (plot) {
      plot(ans)
    }
    return(invisible(ans))
  }

#' Profile log-likelihood for the generalized Pareto distribution
#'
#' This function calculates the (modified) profile likelihood based on the \eqn{p^*} formula.
#' There are two small-sample corrections that use a proxy for
#' \eqn{\ell_{\lambda; \hat{\lambda}}}{the sample space derivative of the nuisance},
#' which are based on Severini's (1999) empirical covariance
#' and the Fraser and Reid tangent exponential model approximation.
#' @details The three \code{mod} available are \code{profile} (the default), \code{tem}, the tangent exponential model (TEM) approximation and
#' \code{modif} for the penalized profile likelihood based on \eqn{p^*} approximation proposed by Severini.
#' For the latter, the penalization is based on the TEM or an empirical covariance adjustment term.
#'
#' @param psi parameter vector over which to profile (unidimensional)
#' @param param string indicating the parameter to profile over
#' @param mod string indicating the model. See \bold{Details}.
#' @param mle maximum likelihood estimate in \eqn{(\psi, \xi)} parametrization if \eqn{\psi \neq \xi} and \eqn{(\sigma, \xi)} otherwise (optional).
#' @param dat sample vector of excesses, unless \code{thresh} is provided (in which case user provides original data)
#' @param m number of observations of interest for return levels. Required only for \code{args} values \code{'VaR'} or \code{'ES'}
#' @param N size of block over which to take maxima. Required only for \code{args} \code{Nmean} and \code{Nquant}.
#' @param p tail probability, equivalent to \eqn{1/m}. Required only for \code{args} \code{quant}.
#' @param q level of quantile for N-block maxima. Required only for \code{args} \code{Nquant}.
#' @param correction logical indicating whether to use \code{spline.corr} to smooth the tem approximation.
#' @param thresh numerical threshold above which to fit the generalized Pareto distribution
#' @param plot logical; should the profile likelihood be displayed? Default to \code{TRUE}
#' @param ... additional arguments such as output from call to \code{Vfun} if \code{mode='tem'}.
#'
#' @return a list with components
#' \itemize{
#' \item \code{mle}: maximum likelihood estimate
#' \item \code{psi.max}: maximum profile likelihood estimate
#' \item \code{param}: string indicating the parameter to profile over
#' \item \code{std.error}: standard error of \code{psi.max}
#' \item \code{psi}: vector of parameter \eqn{\psi} given in \code{psi}
#' \item \code{pll}: values of the profile log likelihood at \code{psi}
#' \item \code{maxpll}: value of maximum profile log likelihood
#' \item \code{family}: a string indicating "gpd"
#' \item \code{thresh}: value of the threshold, by default zero
#' }
#' In addition, if \code{mod} includes \code{tem}
#' \itemize{
#' \item \code{normal}: maximum likelihood estimate and standard error of the interest parameter \eqn{\psi}
#' \item \code{r}: values of likelihood root corresponding to \eqn{\psi}
#' \item \code{q}: vector of likelihood modifications
#' \item \code{rstar}: modified likelihood root vector
#' \item \code{rstar.old}: uncorrected modified likelihood root vector
#' \item \code{tem.psimax}: maximum of the tangent exponential model likelihood
#' }
#' In addition, if \code{mod} includes \code{modif}
#' \itemize{
#' \item \code{tem.mle}: maximum of tangent exponential modified profile log likelihood
#' \item \code{tem.profll}: values of the modified profile log likelihood at \code{psi}
#' \item \code{tem.maxpll}: value of maximum modified profile log likelihood
#' \item \code{empcov.mle}: maximum of Severini's empirical covariance modified profile log likelihood
#' \item \code{empcov.profll}: values of the modified profile log likelihood at \code{psi}
#' \item \code{empcov.maxpll}: value of maximum modified profile log likelihood
#' }
#' @export
#' @examples
#' \dontrun{
#' dat <- rgp(n = 100, scale = 2, shape = 0.3)
#' gpd.pll(psi = seq(-0.5, 1, by=0.01), param = 'shape', dat = dat)
#' gpd.pll(psi = seq(0.1, 5, by=0.1), param = 'scale', dat = dat)
#' gpd.pll(psi = seq(20, 35, by=0.1), param = 'quant', dat = dat, p = 0.01)
#' gpd.pll(psi = seq(20, 80, by=0.1), param = 'ES', dat = dat, m = 100)
#' gpd.pll(psi = seq(15, 100, by=1), param = 'Nmean', N = 100, dat = dat)
#' gpd.pll(psi = seq(15, 90, by=1), param = 'Nquant', N = 100, dat = dat, q = 0.5)
#' }
gpd.pll <-
  function(
    psi,
    param = c(
      "scale",
      "shape",
      "quant",
      "retlev",
      "VaR",
      "ES",
      "Nmean",
      "Nquant"
    ),
    mod = "profile",
    mle = NULL,
    dat,
    m = NULL,
    N = NULL,
    p = NULL,
    q = NULL,
    correction = TRUE,
    thresh = NULL,
    plot = TRUE,
    ...
  ) {
    param <- match.arg(param)
    mod <- match.arg(
      arg = mod,
      choices = c("profile", "tem", "modif"),
      several.ok = TRUE
    )
    #If there is a threshold
    rate <- 1
    if (!is.null(thresh)) {
      stopifnot(is.numeric(thresh), length(thresh) == 1)
      if (min(dat) < thresh) {
        ntot <- length(dat)
        dat <- dat[dat > thresh] - thresh
        nabove <- length(dat)
        rate <- nabove / ntot
      } else {
        dat <- dat - thresh
      }
    } else {
      thresh <- 0
    }

    # Sanity checks to ensure all arguments are provided
    if (is.null(N)) {
      if (param %in% c("Nmean", "Nquant")) {
        stop("Argument \"N\" missing. Procedure aborted")
      } else {
        N <- NA
      }
    }
    if (is.null(m)) {
      if (param %in% c("VaR", "ES")) {
        stop("Argument \"m\" missing. Procedure aborted")
      } else {
        m <- NA
      }
    }

    if (is.null(q)) {
      if (param == "Nquant") {
        stop("Argument \"q\" missing. Procedure aborted")
      } else {
        q <- NA
      }
    }
    if (is.null(p)) {
      if (param %in% c("quant", "retlev")) {
        stop("Argument \"p\" missing. Procedure aborted")
      } else {
        p <- NA
      }
    }
    if (param == "retlev") {
      p <- p / rate
      param <- "quant"
    }
    xmin <- min(dat)
    xmax <- max(dat)
    # Arguments for parametrization of the log likelihood
    if (param == "shape") {
      args <- c("scale", "shape")
    } else {
      args <- c(param, "shape")
    }
    shiftres <- param %in% c("Nmean", "Nquant", "VaR", "quant")
    # If maximum likelihood estimates are not provided, find them
    if (is.null(mle) || length(mle) != 2) {
      mle <- gpd.mle(
        xdat = dat,
        args = args,
        m = m,
        N = N,
        p = p,
        q = q
      )
    }
    # Extract the components, notably V for model `tem`. Keep other components for optimization
    Vprovided <- FALSE
    extra.args <- list(...)
    if ("V" %in% names(extra.args)) {
      V <- extra.args$V
      extra.args$V <- NULL
      if (isTRUE(all.equal(dim(V), c(length(dat), 1)))) {
        Vprovided <- TRUE
      }
    }
    if (!Vprovided) {
      V <-
        switch(
          param,
          scale = gpd.Vfun(par = mle, dat = dat),
          shape = gpd.Vfun(par = mle, dat = dat),
          quant = gpdr.Vfun(par = mle, dat = dat, m = 1 / p),
          VaR = gpdr.Vfun(par = mle, dat = dat, m = m),
          ES = gpde.Vfun(par = mle, dat = dat, m = m),
          Nmean = gpdN.Vfun(par = mle, dat = dat, N = N),
          Nquant = gpdr.Vfun(
            par = mle,
            dat = dat,
            m = 1 / (1 - q^(1 / N))
          )
        )
    }
    # Obtained constrained maximum likelihood estimates for given value of psi
    if (param == "scale") {
      maxll <- gpd.ll(mle, dat = dat)
      std.error <-
        sqrt(solve(gpd.infomat(
          par = mle,
          dat = dat,
          method = "exp"
        ))[1, 1])
      constr.mle.scale <- function(sigmat) {
        as.vector(suppressWarnings(
          optim(
            par = 0.01,
            fn = function(par, scale) {
              -gpd.ll(par = c(scale, par), dat = dat)
            },
            method = "Brent",
            lower = max(-1, -sigmat / xmax),
            upper = min(10, sigmat / xmin),
            scale = sigmat
          )$par
        ))
      }

      # Missing psi vector
      if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
        psirangelow <-
          pmax(1e-05, seq(-3, -1.5, length = 6) * std.error + mle[1])
        lowvals <- sapply(psirangelow, function(par) {
          gpd.ll(c(par, constr.mle.scale(par)), dat = dat)
        }) -
          maxll
        psirangehigh <-
          seq(2.5, 4, length = 6) * std.error + mle[1]
        highvals <- sapply(psirangehigh, function(par) {
          gpd.ll(c(par, constr.mle.scale(par)), dat = dat)
        }) -
          maxll

        lo <- approx(x = lowvals, y = psirangelow, xout = -4)$y
        hi <-
          approx(x = highvals, y = psirangehigh, xout = -4)$y
        lo <-
          ifelse(is.na(lo), lm(psirangelow ~ lowvals)$coef[2] * -4 + mle[1], lo)
        hi <-
          ifelse(
            is.na(hi),
            lm(psirangehigh ~ highvals)$coef[2] * -4 + mle[1],
            hi
          )
        psi <- seq(lo, hi, length = 55)
      }
      if (any(as.vector(psi) < 0)) {
        warning("Negative scale values provided.")
        psi <- psi[psi > 0]
        if (length(psi) == 0) {
          psi <- mle[1]
        }
      }

      pars <- cbind(psi, sapply(psi, constr.mle.scale))
      # Profile log likelihood values for psi
      profll <- apply(pars, 1, function(par) {
        gpd.ll(par = par, dat = dat)
      })
      r <- sign(mle[param] - psi) * sqrt(2 * (maxll - profll))
      if ("tem" %in% mod) {
        phi.mle <- gpd.phi(par = mle, dat = dat, V = V)
        q2num <- apply(pars, 1, function(par) {
          det(rbind(
            c(
              c(phi.mle) -
                gpd.phi(
                  par = par,
                  dat = dat,
                  V = V
                )
            ),
            gpd.dphi(
              par = par,
              dat = dat,
              V = V
            )[-1, ]
          ))
        })
        if (isTRUE(any(sign(q2num) * sign(r) < 0, na.rm = TRUE))) {
          warning(
            "Correction factor and likelihood root are of opposite sign - check output"
          )
        }
        logq <- apply(pars, 1, function(par) {
          -0.5 *
            log(gpd.infomat(
              par = par,
              dat = dat,
              method = "obs"
            )[-1, -1])
        }) +
          log(abs(q2num))
        qmlecontrib <-
          -log(det(gpd.dphi(
            par = mle,
            dat = dat,
            V = V
          ))) +
          0.5 *
            log(det(gpd.infomat(
              par = mle,
              dat = dat,
              method = "obs"
            )))
        logq <- logq + qmlecontrib
        qcor <- sign(q2num) * exp(logq)
        rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r))) / r)
        tem.max.opt <- function(psi, dat = dat) {
          para <- c(psi, constr.mle.scale(psi))
          pll <- gpd.ll(par = para, dat = dat)
          rs <- 2 * (maxll - pll)
          logq <-
            -0.5 *
            log(gpd.infomat(
              par = para,
              dat = dat,
              method = "obs"
            )[-1, -1]) +
            qmlecontrib +
            log(abs(det(rbind(
              c(
                c(phi.mle) -
                  gpd.phi(
                    par = para,
                    dat = dat,
                    V = V
                  )
              ),
              gpd.dphi(
                par = para,
                dat = dat,
                V = V
              )[-1, ]
            ))))
          rs + logq - log(sqrt(abs(rs)))
        }
        tem.max <-
          optim(
            par = mle[1],
            fn = tem.max.opt,
            method = "Brent",
            dat = dat,
            lower = max(1e-05, mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(abstol = 1e-10)
          )$par
      }

      if ("modif" %in% mod) {
        # Tangent exponential model approximation of Fraser and Reid to the profile likelihood
        tem.objfunc.scale <- function(par) {
          0.5 *
            log(gpd.infomat(
              par = par,
              dat = dat,
              method = "obs"
            )[2, 2]) -
            log(abs(gpd.dphi(
              par = par,
              dat = dat,
              V = V[, 2, drop = FALSE]
            )[2, 1]))
        }
        optim.tem.fn.scale <- function(psi) {
          theta.psi.opt <- constr.mle.scale(psi)
          param <-
            c(psi, theta.psi.opt) #ll <- -theta.psi.opt$nllh
          ll <- gpd.ll(param, dat = dat)
          ll + tem.objfunc.scale(param)
        }
        # TEM profile log likelihood values for psi
        proflltem <- profll + apply(pars, 1, tem.objfunc.scale)
        # Maximum objective function for TEM (line search in neighborhood of the MLE)

        tem.mle.opt <-
          optim(
            par = mle[1],
            fn = optim.tem.fn.scale,
            method = "Brent",
            lower = max(1e-05, mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(fnscale = -1)
          )
        tem.mle <-
          c(tem.mle.opt$par, constr.mle.scale(tem.mle.opt$par))

        # Severini empirical covariance function adjustment to profile likelihood Score function -
        # observation-wise
        gpd.score.f <- function(par, dat) {
          sigma = par[1]
          xi = par[2]
          if (!isTRUE(all.equal(0, xi))) {
            cbind(
              dat * (xi + 1) / (sigma^2 * (dat * xi / sigma + 1)) - 1 / sigma,
              -dat *
                (1 / xi + 1) /
                (sigma * (dat * xi / sigma + 1)) +
                log(pmax(dat * xi / sigma + 1, 0)) / xi^2
            )
          } else {
            cbind(
              (dat - sigma) / sigma^2,
              1 / 2 * (dat - 2 * sigma) * dat / sigma^2
            )
          }
        }
        # Score at MLE (sums to zero)
        score.scale.mle <-
          gpd.score.f(mle, dat)[, 2] #keep s_lambda
        empcov.objfunc.scale <- function(par) {
          0.5 *
            log(gpd.infomat(
              par = par,
              dat = dat,
              method = "obs"
            )[2, 2]) -
            log(abs(sum(
              score.scale.mle * gpd.score.f(par, dat)[, 2]
            )))
        }
        profllempcov <-
          profll + apply(pars, 1, empcov.objfunc.scale)
        optim.empcov.fn.scale <- function(psi) {
          theta.psi.opt <- constr.mle.scale(psi)
          param <- c(psi, theta.psi.opt)
          ll <- gpd.ll(param, dat = dat)
          ll + empcov.objfunc.scale(param)
        }
        empcov.mle.opt <-
          optim(
            par = mle[1],
            fn = optim.empcov.fn.scale,
            method = "Brent",
            lower = max(1e-05, mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(fnscale = -1)
          )
        empcov.mle <-
          c(empcov.mle.opt$par, constr.mle.scale(empcov.mle.opt$par))
      }
      # Shape parameter
    } else if (param == "shape") {
      maxll <- try(gpd.ll(mle, dat = dat))
      if (inherits(maxll, "try-error")) {
        stop(paste(
          "Could not find maximum likelihood estimation for the sample of size",
          length(dat)
        ))
      }
      if (isTRUE(mle['shape'] > -0.5)) {
        std.error <-
          sqrt(solve(gpd.infomat(
            par = mle,
            dat = dat,
            method = "exp"
          ))[2, 2])
      } else {
        std.error <- NA
      }
      constr.mle.shape <- function(xit) {
        as.vector(suppressWarnings(
          optim(
            par = 2 * abs(xit) * xmax,
            fn = function(par, shape) {
              -gpd.ll(par = c(par, shape), dat = dat)
            },
            method = "Brent",
            shape = xit,
            lower = ifelse(xit < 0, abs(xit) * xmax + 1e-05, 1e-05),
            upper = 1e+10
          )$par
        ))
      }
      # Missing psi vector
      if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
        if (isTRUE(is.finite(std.error))) {
          psirangelow <-
            seq(ifelse(mle[2] < 0, -7, -5), -1.5, length = 10) *
            std.error +
            mle[2]
          psirangelow <- psirangelow[psirangelow > -1]
          if (length(psirangelow) > 0L) {
            lowvals <- sapply(psirangelow, function(par) {
              gpd.ll(c(constr.mle.shape(par), par), dat = dat)
            }) -
              maxll
          }
          psirangehigh <-
            seq(1.5, 10, length = 10) * std.error + mle[2]
          highvals <- sapply(psirangehigh, function(par) {
            gpd.ll(c(constr.mle.shape(par), par), dat = dat)
          }) -
            maxll

          lo <- approx(x = lowvals, y = psirangelow, xout = -4)$y
          hi <-
            approx(x = highvals, y = psirangehigh, xout = -4)$y
          lo <-
            ifelse(
              is.na(lo),
              lm(psirangelow ~ lowvals)$coef[2] * -4 + mle[2],
              lo
            )
          hi <-
            ifelse(
              is.na(hi),
              lm(psirangehigh ~ highvals)$coef[2] * -4 + mle[2],
              hi
            )
          psi <- seq(lo, hi, length = 55)
          psi <- psi[psi > -1]
        } else {
          psi <- seq(-1, 0, length.out = 21)
        }
      }

      pars <- cbind(sapply(psi, constr.mle.shape), psi)
      # Profile log likelihood values for psi
      profll <- apply(pars, 1, function(par) {
        gpd.ll(par = par, dat = dat)
      })
      r <- sign(mle[param] - psi) * sqrt(2 * (maxll - profll))
      if ("tem" %in% mod) {
        phi.mle <- gpd.phi(par = mle, dat = dat, V = V)
        q2num <- apply(pars, 1, function(par) {
          det(rbind(
            gpd.dphi(
              par = par,
              dat = dat,
              V = V
            )[-2, ],
            c(
              c(phi.mle) -
                gpd.phi(
                  par = par,
                  dat = dat,
                  V = V
                )
            )
          ))
        })
        if (isTRUE(any(sign(q2num) * sign(r) < 0, na.rm = TRUE))) {
          warning(
            "Correction factor and likelihood root are of opposite sign - check output"
          )
        }
        logq <- apply(pars, 1, function(par) {
          -0.5 *
            log(gpd.infomat(
              par = par,
              dat = dat,
              method = "obs"
            )[-2, -2])
        }) +
          log(abs(q2num))
        qmlecontrib <-
          -log(det(gpd.dphi(
            par = mle,
            dat = dat,
            V = V
          ))) +
          0.5 *
            log(det(gpd.infomat(
              par = mle,
              dat = dat,
              method = "obs"
            )))
        logq <- logq + qmlecontrib
        qcor <- sign(q2num) * exp(logq)
        rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r))) / r)
        # Maximum of TEM likelihood - indirect estimation via rstar
        tem.max.opt <- function(psi, dat = dat) {
          para <- c(constr.mle.shape(psi), psi)
          pll <- gpd.ll(par = para, dat = dat)
          rs <- 2 * (maxll - pll)
          logq <-
            -0.5 *
            log(gpd.infomat(
              par = para,
              dat = dat,
              method = "obs"
            )[-2, -2]) +
            qmlecontrib +
            log(abs(det(rbind(
              c(
                c(phi.mle) -
                  gpd.phi(
                    par = para,
                    dat = dat,
                    V = V
                  )
              ),
              gpd.dphi(
                par = para,
                dat = dat,
                V = V
              )[-2, ]
            ))))
          rs + logq - log(sqrt(abs(rs)))
        }
        tem.max <-
          optim(
            par = mle[2],
            fn = tem.max.opt,
            method = "Brent",
            dat = dat,
            lower = mle[2] -
              std.error,
            upper = mle[2] + std.error,
            control = list(abstol = 1e-10)
          )$par
      }
      if ("modif" %in% mod) {
        # Tangent exponential model approximation of Fraser and Reid to the profile likelihood
        tem.objfunc.shape <- function(par) {
          0.5 *
            log(gpd.infomat(
              par = par,
              dat = dat,
              method = "obs"
            )[1, 1]) -
            log(abs(gpd.dphi(
              par = par,
              dat = dat,
              V = V[, 1, drop = FALSE]
            )[1, 1]))
        }
        optim.tem.fn.shape <- function(psi) {
          theta.psi.opt <- constr.mle.shape(psi)
          param <- c(theta.psi.opt, psi)
          ll <- gpd.ll(param, dat = dat)
          ll + tem.objfunc.shape(param)
        }
        # TEM profile log likelihood values for psi
        proflltem <- profll + apply(pars, 1, tem.objfunc.shape)
        # Maximum objective function for TEM (line search in neighborhood of the MLE)
        tem.mle.opt <-
          optim(
            par = mle[2],
            fn = optim.tem.fn.shape,
            method = "Brent",
            lower = mle[2] -
              std.error,
            upper = mle[2] + std.error,
            control = list(fnscale = -1)
          )
        tem.mle <-
          c(constr.mle.shape(tem.mle.opt$par), tem.mle.opt$par)

        # Severini empirical covariance function adjustment to profile likelihood Score function -
        # observation-wise
        gpd.score.f <- function(par, dat) {
          sigma = par[1]
          xi = par[2]
          if (!isTRUE(all.equal(0, xi))) {
            cbind(
              dat * (xi + 1) / (sigma^2 * (dat * xi / sigma + 1)) - 1 / sigma,
              -dat *
                (1 / xi + 1) /
                (sigma * (dat * xi / sigma + 1)) +
                log(pmax(dat * xi / sigma + 1, 0)) / xi^2
            )
          } else {
            cbind(
              (dat - sigma) / sigma^2,
              1 / 2 * (dat - 2 * sigma) * dat / sigma^2
            )
          }
        }
        # Score at MLE (sums to zero)
        score.shape.mle <-
          gpd.score.f(mle, dat)[, 1] #keep s_lambda
        empcov.objfunc.shape <- function(par) {
          0.5 *
            log(gpd.infomat(
              par = par,
              dat = dat,
              method = "obs"
            )[1, 1]) -
            log(abs(sum(
              score.shape.mle * gpd.score.f(par, dat)[, 1]
            )))
        }
        profllempcov <-
          profll + apply(pars, 1, empcov.objfunc.shape)
        optim.empcov.fn.shape <- function(psi) {
          theta.psi.opt <- constr.mle.shape(psi)
          param <- c(theta.psi.opt, psi)
          ll <- gpd.ll(param, dat = dat)
          ll + empcov.objfunc.shape(param)
        }
        empcov.mle.opt <-
          optim(
            par = mle[2],
            fn = optim.empcov.fn.shape,
            method = "Brent",
            lower = mle[2] - std.error,
            upper = mle[2] + std.error,
            control = list(fnscale = -1)
          )
        empcov.mle <-
          c(constr.mle.shape(empcov.mle.opt$par), empcov.mle.opt$par)
      }

      # Return levels, quantiles or value-at-risk
    } else if (param %in% c("quant", "VaR", "Nquant")) {
      if (param == "quant") {
        m <- 1 / p
      }
      if (param == "Nquant") {
        m <- 1 / (1 - q^(1 / N))
      }
      maxll <- gpdr.ll(mle, dat = dat, m = m)
      std.error <-
        sqrt(solve(gpdr.infomat(
          par = mle,
          dat = dat,
          method = "exp",
          m = m
        ))[1, 1])
      constr.mle.quant <- function(quant) {
        suppressWarnings(suppressMessages(
          Rsolnp::solnp(
            par = 0.01,
            function(lambda, psi, m) {
              -gpdr.ll(par = c(psi, lambda), dat = dat, m = m)
            },
            psi = quant,
            m = m,
            control = list(trace = 0)
          )$par
        ))
      }

      # Missing psi vector
      if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
        psirangelow <-
          unique(pmax(
            mean(dat),
            seq(-3, -0.5, length = 12) * std.error + mle[1]
          ))
        lowvals <- sapply(psirangelow, function(par) {
          gpdr.ll(c(par, constr.mle.quant(par)), m = m, dat = dat)
        }) -
          maxll
        psirangehigh <-
          seq(2, 10, length = 12) * std.error + mle[1]
        highvals <- sapply(psirangehigh, function(par) {
          gpdr.ll(c(par, constr.mle.quant(par)), m = m, dat = dat)
        }) -
          maxll

        lo <- approx(x = lowvals, y = psirangelow, xout = -4)$y
        hi <-
          approx(x = highvals, y = psirangehigh, xout = -4)$y
        lo <-
          ifelse(is.na(lo), lm(psirangelow ~ lowvals)$coef[2] * -4 + mle[1], lo)
        hi <-
          ifelse(
            is.na(hi),
            lm(psirangehigh ~ highvals)$coef[2] * -4 + mle[1],
            hi
          )
        psi <- seq(lo, hi, length = 55)
      } else {
        psi <- psi - thresh
        if (any(psi < 0)) {
          stop(
            "Invalid psi sequence: rescaled psi for quantiles must be positive"
          )
        }
      }

      pars <- cbind(psi, sapply(psi, constr.mle.quant))
      # Profile log likelihood values for psi
      profll <- apply(pars, 1, function(par) {
        gpdr.ll(par = par, dat = dat, m = m)
      })
      r <- sign(mle[param] - psi) * sqrt(2 * (maxll - profll))
      if ("tem" %in% mod) {
        phi.mle <- gpdr.phi(
          par = mle,
          dat = dat,
          m = m,
          V = V
        )
        q2num <- apply(pars, 1, function(par) {
          det(rbind(
            c(
              c(phi.mle) -
                gpdr.phi(
                  par = par,
                  dat = dat,
                  V = V,
                  m = m
                )
            ),
            gpdr.dphi(
              par = par,
              dat = dat,
              V = V,
              m = m
            )[-1, ]
          ))
        })
        if (isTRUE(any(sign(q2num) * sign(r) < 0, na.rm = TRUE))) {
          warning(
            "Correction factor and likelihood root are of opposite sign - check output"
          )
        }
        logq <- apply(pars, 1, function(par) {
          -0.5 *
            log(gpdr.infomat(
              par = par,
              dat = dat,
              method = "obs",
              m = m
            )[-1, -1])
        }) +
          log(abs(q2num))
        qmlecontrib <-
          -log(det(gpdr.dphi(
            par = mle,
            dat = dat,
            V = V,
            m = m
          ))) +
          0.5 *
            log(det(gpdr.infomat(
              par = mle,
              dat = dat,
              method = "obs",
              m = m
            )))
        logq <- logq + qmlecontrib
        qcor <- sign(q2num) * exp(logq)
        rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r))) / r)

        tem.max.opt <- function(psi, dat = dat) {
          para <- c(psi, constr.mle.quant(psi))
          pll <- gpdr.ll(par = para, dat = dat, m = m)
          rs <- 2 * (maxll - pll)
          logq <-
            -0.5 *
            log(gpdr.infomat(
              par = para,
              dat = dat,
              method = "obs",
              m = m
            )[-1, -1]) +
            qmlecontrib +
            log(abs(det(rbind(
              c(
                c(phi.mle) -
                  gpdr.phi(
                    par = para,
                    dat = dat,
                    V = V,
                    m = m
                  )
              ),
              gpdr.dphi(
                par = para,
                dat = dat,
                V = V,
                m = m
              )[-1, ]
            ))))
          rs + logq - log(sqrt(abs(rs)))
        }
        tem.max <-
          optim(
            par = mle[1],
            fn = tem.max.opt,
            method = "Brent",
            dat = dat,
            lower = max(1e-05, mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(abstol = 1e-10)
          )$par
      }
      if ("modif" %in% mod) {
        # Tangent exponential model approximation of Fraser and Reid to the profile likelihood
        tem.objfunc.quant <- function(par) {
          0.5 *
            log(gpdr.infomat(
              par = par,
              dat = dat,
              method = "obs",
              m = m
            )[2, 2]) -
            log(abs(gpdr.dphi(
              par = par,
              dat = dat,
              m = m,
              V = V[, 2, drop = FALSE]
            )[2, 1]))
        }
        optim.tem.fn.quant <- function(psi) {
          theta.psi.opt <- constr.mle.quant(psi)
          param <- c(psi, theta.psi.opt)
          ll <- gpdr.ll(param, dat = dat, m = m)
          ll + tem.objfunc.quant(param)
        }
        # TEM profile log likelihood values for psi
        proflltem <-
          profll + suppressWarnings(apply(pars, 1, tem.objfunc.quant))
        # Maximum objective function for TEM
        tem.mle.opt <-
          optim(
            par = mle[1],
            fn = optim.tem.fn.quant,
            method = "Brent",
            lower = max(1e-05, mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(fnscale = -1)
          )
        tem.mle <-
          c(tem.mle.opt$par, constr.mle.quant(tem.mle.opt$par))

        # Severini empirical covariance function adjustment to profile likelihood Score function -
        # observation-wise for xi
        gpdr.score.f <- function(par, dat, m) {
          xi = par[2]
          r = par[1]
          -dat *
            m^xi *
            (1 / xi + 1) *
            log(m) /
            ((dat * (m^xi - 1) / r + 1) * r) +
            (m^xi *
              r *
              xi *
              log(m) /
              (m^xi - 1)^2 -
              r / (m^xi - 1)) *
              (m^xi - 1) /
              (r * xi) +
            log(
              dat *
                (m^xi - 1) /
                r +
                1
            ) /
              xi^2
        }
        # Score at MLE (sums to zero)
        score.quant.mle <- gpdr.score.f(mle, dat, m)
        empcov.objfunc.quant <- function(par) {
          0.5 *
            log(gpdr.infomat(
              par = par,
              dat = dat,
              method = "obs",
              m = m
            )[2, 2]) -
            log(abs(sum(
              score.quant.mle * gpdr.score.f(par, dat, m = m)
            )))
        }
        profllempcov <-
          profll + suppressWarnings(apply(pars, 1, empcov.objfunc.quant))
        optim.empcov.fn.quant <- function(psi) {
          theta.psi.opt <- constr.mle.quant(psi)
          param <- c(psi, theta.psi.opt)
          ll <- gpdr.ll(param, dat = dat, m = m)
          ll + empcov.objfunc.quant(param)
        }
        empcov.mle.opt <-
          optim(
            par = mle[1],
            fn = optim.empcov.fn.quant,
            method = "Brent",
            lower = max(1e-05, mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(fnscale = -1)
          )
        empcov.mle <-
          c(empcov.mle.opt$par, constr.mle.quant(empcov.mle.opt$par))
      }
    } else if (param == "ES") {
      maxll <- gpde.ll(mle, dat = dat, m = m)
      std.error <-
        sqrt(solve(gpde.infomat(
          par = mle,
          dat = dat,
          method = "exp",
          m = m
        ))[1, 1])
      constr.mle.es <- function(psif) {
        # nloptr::auglag(x0 = 0.1, fn = function(x){-gpde.ll(par = c(psif, x), dat = dat, m = m)},
        # hin = function(x){ c(psif*(1-x)*((m^x-1)/x+1)^(-1),
        # 1+psif*(1-x)*((m^x-1)/x+1)^(-1)/x*xmin, 1+psif*(1-x)*((m^x-1)/x+1)^(-1)/x*xmax)}, lower =
        # -1+1e-10, upper = 1-1e-5)$par}
        optim(
          par = 0.1,
          fn = function(x) {
            -gpde.ll(par = c(psif, x), dat = dat, m = m)
          },
          method = "Brent",
          lower = -1 + 1e-10,
          upper = 1 - 1e-05
        )$par
      }

      # Missing psi vector
      if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
        psirangelow <-
          unique(pmax(
            mean(dat),
            seq(-3, -0.1, length = 15) * std.error + mle[1]
          ))
        lowvals <- sapply(psirangelow, function(par) {
          gpde.ll(c(par, constr.mle.es(par)), m = m, dat = dat)
        }) -
          maxll
        psirangehigh <-
          seq(2.5, 30, length = 12) * std.error + mle[1]
        highvals <- sapply(psirangehigh, function(par) {
          gpde.ll(c(par, constr.mle.es(par)), m = m, dat = dat)
        }) -
          maxll

        lo <- approx(x = lowvals, y = psirangelow, xout = -4)$y
        hi <-
          approx(x = highvals, y = psirangehigh, xout = -4)$y
        lo <-
          ifelse(is.na(lo), lm(psirangelow ~ lowvals)$coef[2] * -4 + mle[1], lo)
        hi <-
          ifelse(
            is.na(hi),
            lm(psirangehigh ~ highvals)$coef[2] * -4.5 + mle[1],
            hi
          )
        psi <-
          c(
            seq(lo, mle[1], length = 15)[-15],
            seq(mle[1], min(mle[1] + 2 * std.error, hi), length = 25)[-1]
          )
        if (mle[1] + 2 * std.error < hi) {
          psi <- c(psi, seq(mle[1] + 2 * std.error, hi, length = 20))
        }
      }
      # Do not remove threshold for expected shortfall

      pars <- cbind(psi, sapply(psi, constr.mle.es))
      # Profile log likelihood values for psi
      profll <- apply(pars, 1, function(par) {
        gpde.ll(par = par, dat = dat, m = m)
      })
      profll[profll == 1e+10] <- NA
      r <- sign(mle[param] - psi) * sqrt(2 * (maxll - profll))
      if ("tem" %in% mod) {
        phi.mle <- gpde.phi(
          par = mle,
          dat = dat,
          m = m,
          V = V
        )
        q2num <- apply(pars, 1, function(par) {
          det(rbind(
            c(
              c(phi.mle) -
                gpde.phi(
                  par = par,
                  dat = dat,
                  V = V,
                  m = m
                )
            ),
            gpde.dphi(
              par = par,
              dat = dat,
              V = V,
              m = m
            )[-1, ]
          ))
        })
        if (isTRUE(any(sign(q2num) * sign(r) < 0, na.rm = TRUE))) {
          warning(
            "Correction factor and likelihood root are of opposite sign - check output"
          )
        }
        logq <- apply(pars, 1, function(par) {
          -0.5 *
            log(gpde.infomat(
              par = par,
              dat = dat,
              method = "obs",
              m = m
            )[-1, -1])
        }) +
          log(abs(q2num))
        qmlecontrib <-
          -log(det(gpde.dphi(
            par = mle,
            dat = dat,
            V = V,
            m = m
          ))) +
          0.5 *
            log(det(gpde.infomat(
              par = mle,
              dat = dat,
              method = "obs",
              m = m
            )))
        logq <- logq + qmlecontrib
        qcor <- sign(q2num) * exp(logq)
        rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r))) / r)

        tem.max.opt <- function(psi, dat = dat) {
          para <- c(psi, constr.mle.es(psi))
          pll <- gpde.ll(par = para, dat = dat, m = m)
          rs <- 2 * (maxll - pll)
          logq <-
            -0.5 *
            log(gpde.infomat(
              par = para,
              dat = dat,
              method = "obs",
              m = m
            )[-1, -1]) +
            qmlecontrib +
            log(abs(det(rbind(
              c(
                c(phi.mle) -
                  gpde.phi(
                    par = para,
                    dat = dat,
                    V = V,
                    m = m
                  )
              ),
              gpde.dphi(
                par = para,
                dat = dat,
                V = V,
                m = m
              )[-1, ]
            ))))
          rs + logq - log(sqrt(abs(rs)))
        }
        tem.max <-
          optim(
            par = mle[1],
            fn = tem.max.opt,
            method = "Brent",
            dat = dat,
            lower = max(1e-05, mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(abstol = 1e-10)
          )$par
      }
      if ("modif" %in% mod) {
        # Tangent exponential model approximation of Fraser and Reid to the profile likelihood
        tem.objfunc.es <- function(par) {
          0.5 *
            log(gpde.infomat(
              par = par,
              dat = dat,
              method = "obs",
              m = m
            )[2, 2]) -
            log(abs(gpde.dphi(
              par = par,
              dat = dat,
              m = m,
              V = V[, 2, drop = FALSE]
            )[2, 1]))
        }
        optim.tem.fn.es <- function(psi) {
          theta.psi.opt <- constr.mle.es(psi)
          param <- c(psi, theta.psi.opt)
          ll <- gpde.ll(param, dat = dat, m = m)
          ll + tem.objfunc.es(param)
        }
        # TEM profile log likelihood values for psi
        proflltem <-
          profll + suppressWarnings(apply(pars, 1, tem.objfunc.es))
        # Maximum objective function for TEM

        tem.mle.opt <-
          optim(
            par = mle[1],
            fn = optim.tem.fn.es,
            method = "Brent",
            lower = max(quantile(dat, 1 - 1 / m), mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(fnscale = -1)
          )
        tem.mle <-
          c(tem.mle.opt$par, constr.mle.es(tem.mle.opt$par))

        # Severini empirical covariance function adjustment to profile likelihood Score function -
        # observation-wise for xi
        gpde.score.f <- function(par, dat, m) {
          xi = par[2]
          es = par[1]
          -((m^xi * log(m) + 1) *
            dat /
            (es * (xi - 1)) -
            dat *
              (m^xi + xi - 1) /
              (es *
                (xi -
                  1)^2)) *
            (1 / xi + 1) /
            (dat * (m^xi + xi - 1) / (es * (xi - 1)) - 1) +
            (m^xi +
              xi -
              1) *
              ((m^xi * log(m) + 1) *
                (xi - 1) *
                xi /
                (m^xi + xi - 1)^2 -
                (xi -
                  1) /
                  (m^xi + xi - 1) -
                xi / (m^xi + xi - 1)) /
              ((xi - 1) * xi) +
            log(pmax(
              0,
              -dat *
                (m^xi + xi - 1) /
                (es * (xi - 1)) +
                1
            )) /
              xi^2
        }
        # Score at MLE (sums to zero)
        score.es.mle <- gpde.score.f(mle, dat, m)
        empcov.objfunc.es <- function(par) {
          0.5 *
            log(gpde.infomat(
              par = par,
              dat = dat,
              method = "obs",
              m = m
            )[2, 2]) -
            log(abs(sum(
              score.es.mle * gpde.score.f(par, dat, m = m)
            )))
        }
        profllempcov <-
          profll + suppressWarnings(apply(pars, 1, empcov.objfunc.es))
        optim.empcov.fn.es <- function(psi) {
          theta.psi.opt <- constr.mle.es(psi)
          param <- c(psi, theta.psi.opt)
          ll <- gpde.ll(param, dat = dat, m = m)
          ll + empcov.objfunc.es(param)
        }
        empcov.mle.opt <-
          optim(
            par = mle[1],
            fn = optim.empcov.fn.es,
            method = "Brent",
            lower = max(quantile(dat, 1 - 1 / m), mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(fnscale = -1)
          )
        empcov.mle <-
          c(empcov.mle.opt$par, constr.mle.es(empcov.mle.opt$par))
      }
    } else if (param == "Nmean") {
      maxll <- gpdN.ll(mle, dat = dat, N = N)
      std.error <-
        sqrt(solve(gpdN.infomat(
          par = mle,
          dat = dat,
          method = "exp",
          N = N
        ))[1])
      constr.mle.Nmean <- function(Nmeant) {
        suppressWarnings(
          alabama::auglag(
            par = 0.01,
            function(lambda, psi, N) {
              -gpdN.ll(par = c(psi, lambda), dat = dat, N = N)
            },
            gr = function(lambda, psi, N) {
              -gpdN.score(par = c(psi, lambda), dat = dat, N = N)[2]
            },
            hin = function(lambda, psi, N) {
              sigma = ifelse(
                abs(lambda > 1e-8),
                psi *
                  lambda /
                  (exp(
                    lgamma(N + 1) + lgamma(-lambda + 1) - lgamma(N - lambda + 1)
                  ) -
                    1),
                psi / (0.57721566490153231044 + psigamma(N + 1))
              )
              if (lambda >= 0) {
                c(1e-8, sigma, 1 - lambda, lambda + 1)
              } else {
                c(-sigma / lambda - xmax, sigma, 1 - lambda, lambda + 1)
              }
            },
            psi = Nmeant,
            N = N,
            control.outer = list(trace = FALSE)
          )$par
        )
      }
      # Missing psi vector
      if (missing(psi) || any(is.null(psi)) || any(is.na(psi))) {
        #compute profile log-lik on a grid left and right of the MLE
        psirangelow <-
          unique(pmax(
            mean(dat),
            seq(-3, -0.25, length = 20) * std.error + mle[1]
          ))
        lowvals <- sapply(psirangelow, function(par) {
          gpdN.ll(c(par, constr.mle.Nmean(par)), dat = dat, N = N)
        }) -
          maxll
        psirangehigh <-
          seq(2.5, 50, length = 20) * std.error + mle[1]
        highvals <- sapply(psirangehigh, function(par) {
          gpdN.ll(c(par, constr.mle.Nmean(par)), dat = dat, N = N)
        }) -
          maxll
        #Try to do linear interpolation - only works if value inside the interval lowvals or highvals
        lo <- approx(x = lowvals, y = psirangelow, xout = -4)$y
        #Else linear interpolation with linear model fitted to lower values
        lo <-
          ifelse(
            is.na(lo),
            spline(x = lowvals, y = psirangelow, xout = -4)$y,
            lo
          )
        psirangelow <- seq(lo, mle[1], length = 20)
        lowvals <- sapply(psirangelow, function(par) {
          gpdN.ll(c(par, constr.mle.Nmean(par)), dat = dat, N = N)
        }) -
          maxll
        #hi <- approx(x = highvals, y = psirangehigh, xout = -4)$y
        #For upper, use spline approx
        hi <-
          spline(x = highvals, y = psirangehigh, xout = -4)$y
        #Recompute the range
        psirangehigh <- seq(psirangehigh[1], hi, length = 30)
        highvals <- sapply(psirangehigh, function(par) {
          gpdN.ll(c(par, constr.mle.Nmean(par)), dat = dat, N = N)
        }) -
          maxll
        #Remove NAs, inf, etc.
        highvals <- highvals[is.finite(highvals)]
        psirangehigh <- psirangehigh[1:length(highvals)]
        #If could not interpolate, use a simple linear model to predict lower value
        #hi <- ifelse(is.na(hi), spline(x = highvals, y = psirangehigh, xout = -4)$y, hi)
        psi <-
          as.vector(c(
            spline(
              x = c(lowvals, 0),
              y = c(psirangelow, mle[1]),
              xout = seq(-4, -0.1, length = 15)
            )$y,
            mle[1],
            spline(
              x = c(0, highvals),
              y = c(mle[1], psirangehigh),
              xout = seq(-0.1, highvals[length(highvals)], length = 20)
            )$y
          ))
      } else {
        psi <- psi - thresh
      }

      if (any(as.vector(psi) < 0)) {
        warning("Negative Nmean values provided.")
        psi <- psi[psi > 0]
        if (length(psi) == 0) {
          psi <- mle[1]
        }
      }

      pars <- cbind(psi, sapply(psi, constr.mle.Nmean))
      # Profile log likelihood values for psi
      profll <- apply(pars, 1, function(par) {
        gpdN.ll(par = par, dat = dat, N = N)
      })
      r <- sign(mle[param] - psi) * sqrt(2 * (maxll - profll))
      if ("tem" %in% mod) {
        phi.mle <- gpdN.phi(
          par = mle,
          dat = dat,
          N = N,
          V = V
        )
        q2num <- apply(pars, 1, function(par) {
          det(rbind(
            c(
              c(phi.mle) -
                gpdN.phi(
                  par = par,
                  dat = dat,
                  V = V,
                  N = N
                )
            ),
            gpdN.dphi(
              par = par,
              dat = dat,
              V = V,
              N = N
            )[-1, ]
          ))
        })
        if (isTRUE(any(sign(q2num) * sign(r) < 0, na.rm = TRUE))) {
          warning(
            "Correction factor and likelihood root are of opposite sign - check output"
          )
        }

        logq <- apply(pars, 1, function(par) {
          -0.5 *
            log(gpdN.infomat(
              par = par,
              dat = dat,
              method = "obs",
              N = N
            )[-1, -1])
        }) +
          log(abs(q2num))
        qmlecontrib <-
          -log(det(gpdN.dphi(
            par = mle,
            dat = dat,
            V = V,
            N = N
          ))) +
          0.5 *
            log(det(gpdN.infomat(
              par = mle,
              dat = dat,
              method = "obs",
              N = N
            )))
        logq <- logq + qmlecontrib
        qcor <- sign(q2num) * exp(logq)
        rstar <- ifelse(r == 0, 0, r + (logq - log(abs(r))) / r)

        tem.max.opt <- function(psi, dat = dat) {
          para <- c(psi, constr.mle.Nmean(psi))
          pll <- gpdN.ll(par = para, dat = dat, N = N)
          rs <- 2 * (maxll - pll)
          logq <-
            -0.5 *
            log(gpdN.infomat(
              par = para,
              dat = dat,
              method = "obs",
              N = N
            )[-1, -1]) +
            qmlecontrib +
            log(abs(det(rbind(
              c(
                c(phi.mle) -
                  gpdN.phi(
                    par = para,
                    dat = dat,
                    V = V,
                    N = N
                  )
              ),
              gpdN.dphi(
                par = para,
                dat = dat,
                V = V,
                N = N
              )[-1, ]
            ))))
          rs + logq - log(sqrt(abs(rs)))
        }
        tem.max <-
          optim(
            par = mle[1],
            fn = tem.max.opt,
            method = "Brent",
            dat = dat,
            lower = max(1e-05, mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(abstol = 1e-10)
          )$par
      }
      if ("modif" %in% mod) {
        # Tangent exponential model approximation of Fraser and Reid to the profile likelihood
        tem.objfunc.Nmean <- function(par) {
          0.5 *
            log(gpdN.infomat(
              par = par,
              dat = dat,
              method = "obs",
              N = N
            )[2, 2]) -
            log(abs(gpdN.dphi(
              par = par,
              dat = dat,
              N = N,
              V = V[, 2, drop = FALSE]
            )[2, 1]))
        }
        optim.tem.fn.Nmean <- function(psi) {
          theta.psi.opt <- constr.mle.Nmean(psi)
          param <- c(psi, theta.psi.opt)
          ll <- gpdN.ll(param, dat = dat, N = N)
          ll + tem.objfunc.Nmean(param)
        }
        # TEM profile log likelihood values for psi
        proflltem <-
          profll + suppressWarnings(apply(pars, 1, tem.objfunc.Nmean))
        # Maximum objective function for TEM
        tem.mle.opt <-
          optim(
            par = mle[1],
            fn = optim.tem.fn.Nmean,
            method = "Brent",
            lower = max(1e-05, mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(fnscale = -1)
          )
        tem.mle <-
          c(tem.mle.opt$par, constr.mle.Nmean(tem.mle.opt$par))

        # Severini empirical covariance function adjustment to profile likelihood Score function -
        # observation-wise for xi
        gpdN.score.f <- function(par, dat, N) {
          z = par[1]
          xi = par[2]
          cst <-
            exp(lgamma(N + 1) + lgamma(1 - xi) - lgamma(N + 1 - xi))
          -(psigamma(N - xi + 1) * cst - psigamma(-xi + 1) * cst) *
            dat *
            (1 /
              xi +
              1) /
            (z *
              (dat * (cst - 1) / z + 1)) +
            ((psigamma(N - xi + 1) *
              cst -
              psigamma(
                -xi +
                  1
              ) *
                cst) *
              xi *
              z /
              (cst - 1)^2 -
              z / (cst - 1)) *
              (cst - 1) /
              (xi * z) +
            log(
              dat *
                (cst - 1) /
                z +
                1
            ) /
              xi^2
        }
        # Score at MLE (sums to zero)
        score.Nmean.mle <- gpdN.score.f(mle, dat, N)
        empcov.objfunc.Nmean <- function(par) {
          0.5 *
            log(gpdN.infomat(
              par = par,
              dat = dat,
              method = "obs",
              N = N
            )[2, 2]) -
            log(abs(sum(
              score.Nmean.mle * gpdN.score.f(par, dat, N = N)
            )))
        }
        profllempcov <-
          profll + suppressWarnings(apply(pars, 1, empcov.objfunc.Nmean))
        optim.empcov.fn.Nmean <- function(psi) {
          theta.psi.opt <- constr.mle.Nmean(psi)
          param <- c(psi, theta.psi.opt)
          ll <- gpdN.ll(param, dat = dat, N = N)
          ll + empcov.objfunc.Nmean(param)
        }
        empcov.mle.opt <-
          optim(
            par = mle[1],
            fn = optim.empcov.fn.Nmean,
            method = "Brent",
            lower = max(1e-05, mle[1] - std.error),
            upper = mle[1] + std.error,
            control = list(fnscale = -1)
          )
        empcov.mle <-
          c(empcov.mle.opt$par, constr.mle.Nmean(empcov.mle.opt$par))
      }
    }
    # Return profile likelihood and quantities of interest (modified likelihoods)
    colnames(pars) <- names(mle)
    ans <-
      list(
        mle = mle,
        pars = pars,
        psi.max = as.vector(mle[param]),
        param = param,
        std.error = std.error,
        psi = psi,
        pll = profll,
        maxpll = maxll,
        r = r
      )
    # Shift by threshold if non-null
    if (shiftres) {
      ans$psi <- ans$psi + thresh
      ans$mle[1] <- ans$psi.max <- ans$mle[1] + thresh
      ans$pars[, 1] <- ans$pars[, 1] + thresh
    }

    if ("tem" %in% mod) {
      ans$q <- qcor
      ans$rstar <- rstar
      ans$tem.psimax <-
        as.vector(tem.max) + ifelse(shiftres, thresh, 0)
      ans$normal <- c(ans$psi.max, ans$std.error)
      if (correction && length(psi) > 10) {
        ans <- spline.corr(ans)
      }
    }
    if ("modif" %in% mod) {
      ans$tem.mle <- ifelse(param == "shape", tem.mle[2], tem.mle[1])
      if (shiftres) {
        ans$tem.mle[1] <- ans$tem.mle[1] + thresh
      }
      ans$tem.pll <- proflltem
      ans$tem.maxpll <- as.vector(tem.mle.opt$value)
      ans$empcov.mle <-
        ifelse(param == "shape", empcov.mle[2], empcov.mle[1])
      if (shiftres) {
        ans$empcov.mle[1] <- ans$empcov.mle[1] + thresh
      }
      ans$empcov.pll <- as.vector(profllempcov)
      ans$empcov.maxpll <- as.vector(empcov.mle.opt$value)
    }
    if ("tem" %in% mod) {
      class(ans) <- c("eprof", "fr")
    } else {
      class(ans) <- "eprof"
    }
    ans$family <- "gpd"
    ans$thresh <- thresh
    ans$param <- param
    if (plot) {
      plot(ans)
    }
    return(invisible(ans))
  }


#' Plot of tangent exponential model profile likelihood
#'
#' This function is adapted from the \code{plot.fr} function from the \code{hoa} package bundle.
#' It differs from the latter mostly in the placement of legends.
#'
#' @param x an object of class \code{fr} returned by \code{\link{gpd.tem}} or \code{\link{gev.tem}}.
#' @param ... further arguments to \code{plot} currently ignored. Providing a numeric vector \code{which} allows for custom selection of the plots. A logical \code{all}. See \strong{Details}.
#' @return graphs depending on argument \code{which}
#' @details Plots produced depend on the integers provided in \code{which}. \code{1} displays the Wald pivot, the likelihood root \code{r}, the modified likelihood root \code{rstar} and the likelihood modification \code{q} as functions of the parameter \code{psi}. \code{2} gives the renormalized profile log likelihood and adjusted form, with the maximum likelihood having ordinate value of zero. \code{3} provides the significance function, a transformation of \code{1}. Lastly, \code{4} plots the correction factor as a function of the likelihood root; it is a diagnostic plot aimed for detecting failure of
#' the asymptotic approximation, often due to poor numerics in a neighborhood of \code{r=0}; the function should be smooth. The function \code{\link{spline.corr}} is designed to handle this by correcting numerically unstable estimates, replacing outliers and missing values with the fitted values from the fit.
#'
#'
#' @references Brazzale, A. R., Davison, A. C. and Reid, N. (2007). \emph{Applied Asymptotics: Case Studies in Small-Sample Statistics}. Cambridge University Press, Cambridge.
#' @export
plot.fr <- function(x, ...) {
  # plot a fraser-reid object
  whichPlot <- c(1:4) #default
  if (length(list(...)) > 0) {
    if ("which" %in% names(list(...))) {
      whichPlot <- list(...)$which
      whichPlot <- (1:4)[c(1:4 %in% whichPlot)]
    } else if ("all" %in% names(list(...))) {
      if (!is.logical(all)) {
        stop("Invalid \"all\" parameter")
      }
      if (list(...)$all) {
        whichPlot <- 1:4
      } else {
        whichPlot <- 1:2
      }
    }
  }
  old.pars <- par(no.readonly = TRUE)
  if (sum(c(1, 2, 3, 4) %in% whichPlot) > 2) {
    par(mfrow = c(2, 2), mar = c(4.8, 4.8, 1, 0.1))
  } else if (sum(c(1, 2, 3, 4) %in% whichPlot) == 2) {
    par(mfrow = c(1, 2))
  }
  fr <- x
  xl <- ifelse(is.null(fr$param), expression(psi), fr$param)

  if (1 %in% whichPlot) {
    plot(
      fr$psi,
      fr$r,
      type = "l",
      xlab = xl,
      ylab = "pivot",
      ylim = c(-4, 4),
      panel.first = abline(
        h = qnorm(c(
          0.005,
          0.025,
          0.05,
          0.5,
          0.95,
          0.975,
          0.995
        )),
        col = "grey",
        lwd = 0.7
      ),
      bty = "l"
    )
    lines(
      fr$psi,
      (fr$normal[1] - fr$psi) / fr$normal[2],
      col = 3,
      lwd = 1.5
    )
    lines(fr$psi, fr$q, col = 2, lwd = 1.5)
    lines(fr$psi, fr$r, lwd = 1.5)
    lines(fr$psi, fr$rstar, col = 4)
    legend(
      x = "topright",
      c("Wald pivot", "lik. root", "modif. root", expression(q(psi))),
      lty = c(1, 1, 1, 1),
      x.intersp = 0.2,
      lwd = 1.5,
      seg.len = 0.5,
      col = c(3, 1, 4, 2),
      bty = "n",
      cex = 0.9,
      xjust = 1
    )
  }
  # top right: log likelihood (and adjusted version, I think?) as a function of psi
  if (2 %in% whichPlot) {
    plot(
      fr$psi,
      -fr$r^2 / 2,
      type = "l",
      xlab = xl,
      ylab = "Profile log likelihood",
      ylim = c(-8, 0),
      panel.first = abline(
        h = -qchisq(c(0.95, 0.99), df = 1) / 2,
        col = "grey",
        lwd = 0.7
      ),
      lwd = 1.5,
      bty = "l"
    )
    lines(fr$psi, -fr$rstar^2 / 2, col = 4, lwd = 1.5)
    legend(
      x = "bottomright",
      c("profile", "tem"),
      lty = c(1, 1),
      x.intersp = 0.2,
      lwd = 1.5,
      seg.len = 0.5,
      col = c("black", 4),
      bty = "n"
    )
    # optional: add diagnostic panels
  }
  if (3 %in% whichPlot) {
    # lower left: plot of Phi(pivot) as a function of psi

    plot(
      fr$psi,
      pnorm(fr$r),
      type = "l",
      xlab = xl,
      ylab = "Significance function",
      ylim = c(0, 1),
      panel.first = abline(
        h = c(0.025, 0.05, 0.5, 0.95, 0.975),
        col = "grey",
        lwd = 0.7
      ),
      lwd = 1.5,
      bty = "l"
    )
    lines(fr$psi, pnorm(fr$q), col = 2, lwd = 1.5)
    lines(fr$psi, pnorm(fr$rstar), col = 4, lwd = 1.5)
    legend(
      x = "topright",
      c("lik. root", "modif. root", expression(q(psi))),
      lty = c(1, 1, 1),
      col = c(1, 4, 2),
      bty = "n",
      x.intersp = 0.2,
      lwd = 1.5,
      seg.len = 0.5,
      cex = 0.9
    )
  }
  # lower right: log(q/r)/r as a function of r (should be smooth)
  if (4 %in% whichPlot) {
    fit.r <-
      stats::smooth.spline(x = na.omit(cbind(fr$r, fr$rstar)), cv = FALSE)
    pr <- predict(fit.r, 0)$y
    plot(
      fr$r,
      fr$rstar - fr$r,
      type = "l",
      xlab = "Likelihood root r",
      ylab = expression(paste("Correction factor log(q/r)/r")),
      panel.first = {
        abline(h = 0, col = "grey", lwd = 0.7)
        abline(v = 0, col = "grey", lwd = 0.7)
        abline(
          v = pr,
          col = "grey",
          lwd = 0.7,
          lty = 2
        )
      },
      lwd = 1.5,
      bty = "l"
    )
  }

  par(old.pars)
}

#' Spline correction for Fraser-Reid approximations
#'
#' The tangent exponential model can be numerically unstable for values close to \eqn{r = 0}.
#' This function corrects these incorrect values, which are interpolated using splines.
#' The function takes as input an object of class \code{fr} and returns the same object with
#' different \code{rstar} values.
#' @section Warning:
#'
#' While penalized (robust) splines often do a good job at capturing and correcting for numerical outliers and \code{NA}, it
#' may also be driven by unusual values lying on the profile log-likelihood the curve or fail to detect outliers (or falsely identifying `correct' values as outliers). The user should always validate by comparing the plots of both the uncorrected (raw) output of the object with that of \code{spline.corr}.
#' @details If available, the function uses \code{cobs} from the eponym package. The latter handles constraints and smoothness penalties, and is more robust than the equivalent \code{\link[stats]{smooth.spline}}.
#'
#' @param fr an object of class \code{fr}, normally the output of \link{gpd.tem} or \link{gev.tem}.
#' @param method string for the method, either \code{cobs} (constrained robust B-spline from eponym package) or \code{smooth.spline}
#' @return an object of class \code{fr}, containing as additional arguments \code{spline} and a modified \code{rstar} argument.
#' @keywords internal
#' @export
spline.corr <- function(fr, method = c("cobs", "smooth.spline")) {
  # Step 1: fit a smoothing spline to rstar If fit failed for some values (for example when
  # shape forced to be < 1) Remove those values
  method <-
    match.arg(method[1], choices = c("cobs", "smooth.spline"))
  if (all(is.nan(fr$q)) || all(is.nan(fr$rstar))) {
    # If could not compute Fraser-Reid correction, abort
    return(fr)
  }
  fitfailed <- which(!is.finite(fr$r))
  if (length(fitfailed) > 0) {
    fr$r <- fr$r[-fitfailed]
    fr$rstar <- fr$rstar[-fitfailed]
    fr$q <- fr$q[-fitfailed]
    fr$psi <- fr$psi[-fitfailed]
  }
  w <- pchisq(fr$r^2, 0.5)
  # If any correction for q failed and returned NA
  corfailed <- which(!is.finite(fr$rstar))
  # If equispaced values for psi between MLE and other, then we have r = 0
  corfailed <- c(corfailed, which(fr$r == 0))
  if (length(corfailed) > 0) {
    resp <- (fr$rstar - fr$r)[-corfailed]
    regr <- fr$r[-corfailed]
    w <- w[-corfailed]
  } else {
    resp <- (fr$rstar - fr$r)
    regr <- fr$r
  }
  if (
    method == "cobs" &&
      requireNamespace("cobs", quietly = TRUE)
  ) {
    spline <-
      cobs::cobs(
        y = resp,
        x = regr,
        w = w,
        constraint = "none",
        lambda = 1,
        nknots = 20,
        print.mesg = FALSE,
        print.warn = FALSE
      )$fitted
  } else {
    spline <-
      rev(
        stats::smooth.spline(
          y = resp,
          x = regr,
          w = w,
          spar = 0.9,
          cv = TRUE
        )$y
      )
  }
  # Compute difference between fitted values and rstar
  departure <- spline - resp
  # Ad-hoc fix of the values close to MLE where the numerical precision causes difficulty
  # Outlier detection via chi-square test From package outliers, (c)Lukasz Komsta
  scores <- function(x, prob = NA) {
    (x - mean(x))^2 / var(x) > qchisq(prob, 1)
  }
  bad <- which(scores(departure, prob = 0.95))

  if (length(bad) > 0) {
    # Exclude those values if they are in the end of the distribution
    bad <-
      bad[intersect(
        which(bad < 0.85 * length(departure)),
        which(bad > 0.15 * length(departure))
      )]
    # Remove outliers and fit again (with less smoothness)
  }
  if (length(bad) > 0) {
    resp[bad] <- NA
    w <- w[-bad]
  }
  if (requireNamespace("cobs", quietly = TRUE)) {
    spline <-
      cobs::cobs(
        y = resp,
        x = regr,
        constraint = "none",
        w = w,
        lambda = -1,
        ic = "SIC",
        knots.add = TRUE,
        repeat.delete.add = TRUE,
        print.mesg = FALSE,
        print.warn = FALSE
      )
    fr$spline <- spline
    fr$rstar <-
      predict(spline, fr$r, interval = "none")[, 2] + fr$r
  } else {
    spline <-
      stats::smooth.spline(
        x = na.omit(cbind(regr, resp)),
        cv = FALSE,
        all.knots = TRUE
      )
    fr$spline <- spline
    fr$rstar <- predict(spline, fr$r)$y + fr$r
  }
  return(fr)
}

#' Bridging the singularity for higher order asymptotics
#'
#' The correction factor \eqn{\log(q/r)/r} for the
#' likelihood root is unbounded in the vincinity of
#' the maximum likelihood estimator. The thesis of
#' Rongcai Li (University of Toronto, 2001)
#' explores different ways of bridging this
#' singularity, notably using asymptotic expansions.
#'
#' The poor man's method used here consists in
#' fitting a robust regression to \eqn{1/q-1/r}
#' as a function of \eqn{r} and using predictions
#' from the model to solve for \eqn{q}. This
#' approach is seemingly superior to that
#' previously used in \link{spline.corr}.
#'
#' @param fr an object of class \code{fr}
#' @param print.warning logical; should warning message be printed? Default to \code{FALSE}
##' @return an object of class \code{fr}, containing as additional arguments \code{spline} and a modified \code{rstar} argument.
#' @keywords internal
#' @export
tem.corr <- function(fr, print.warning = FALSE) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The \"MASS\" package is required for this function to work")
  }
  if (all(is.nan(fr$q)) || all(is.nan(fr$rstar))) {
    # If could not compute Fraser-Reid correction, abort
    return(fr)
  }
  fitfailed <- which(!is.finite(fr$r))
  if (length(fitfailed) > 0) {
    fr$r <- fr$r[-fitfailed]
    fr$rstar <- fr$rstar[-fitfailed]
    fr$q <- fr$q[-fitfailed]
    fr$psi <- fr$psi[-fitfailed]
  }
  if (length(fr$psi) < 25L) {
    if (print.warning) {
      warning(
        "The correction for the tangent exponential model\n is based on less than 25 observations."
      )
    }
  }
  # this is approximately linear for fixed data
  resp <- 1 / fr$q - 1 / fr$r
  # fit a robust regression to discount observations that are outlying
  robust_reg <- MASS::rlm(resp ~ fr$r)
  # Replace q values by predictions
  pred <- predict(robust_reg)
  qhat <- 1 / (pred + 1 / fr$r)
  fr$q_old <- fr$q
  fr$rstar_old <- fr$rstar
  fr$q <- qhat
  fr$rstar <- fr$r + log(fr$q / fr$r) / fr$r
  fr$tem.psimax
  return(fr)
}


#' Tangent exponential model approximation for the GEV distribution
#'
#' The function \code{gev.tem} provides a tangent exponential model (TEM) approximation
#' for higher order likelihood inference for a scalar parameter for the generalized extreme value distribution.
#' Options include location scale and shape parameters as well as value-at-risk (or return levels).
#' The function attempts to find good values for \code{psi} that will
#' cover the range of options, but the fail may fit and return an error.
#'
#'
#' @param param parameter over which to profile
#' @param psi scalar or ordered vector of values for the interest parameter. If \code{NULL} (default), a grid of values centered at the MLE is selected
#' @param N size of block over which to take maxima. Required only for \code{param} \code{Nmean} and \code{Nquant}.
#' @param p tail probability for the (1-p)th quantile (return levels). Required only if \code{param = 'retlev'}
#' @param q probability level of quantile. Required only for \code{param} \code{Nquant}.
#' @param dat sample vector for the GEV distribution
#' @param n.psi number of values of \code{psi} at which the likelihood is computed, if \code{psi} is not supplied (\code{NULL}). Odd values are more prone to give rise to numerical instabilities near the MLE. If \code{psi} is a vector of length 2 and \code{n.psi} is greater than 2, these are taken to be endpoints of the sequence.
#' @param plot logical indicating whether \code{plot.fr} should be called upon exit
#' @param correction logical indicating whether \link{spline.corr} should be called.
#' @author Leo Belzile
#' @return an invisible object of class \code{fr} (see \code{tem} in package \code{hoa}) with elements
#' \itemize{
#' \item \code{normal}: maximum likelihood estimate and standard error of the interest parameter \eqn{\psi}
#' \item \code{par.hat}: maximum likelihood estimates
#' \item \code{par.hat.se}: standard errors of maximum likelihood estimates
#' \item \code{th.rest}: estimated maximum profile likelihood at (\eqn{\psi}, \eqn{\hat{\lambda}})
#' \item \code{r}: values of likelihood root corresponding to \eqn{\psi}
#' \item \code{psi}: vector of interest parameter
#' \item \code{q}: vector of likelihood modifications
#' \item \code{rstar}: modified likelihood root vector
#' \item \code{rstar.old}: uncorrected modified likelihood root vector
#' \item \code{param}: parameter
#' }
#' @export
#' @examples
#' \dontrun{
#' set.seed(1234)
#' dat <- rgev(n = 40, loc = 0, scale = 2, shape = -0.1)
#' gev.tem('shape', dat = dat, plot = TRUE)
#' gev.tem('quant', dat = dat, p = 0.01, plot = TRUE)
#' gev.tem('scale', psi = seq(1, 4, by = 0.1), dat = dat, plot = TRUE)
#' dat <- rgev(n = 40, loc = 0, scale = 2, shape = 0.2)
#' gev.tem('loc', dat = dat, plot = TRUE)
#' gev.tem('Nmean', dat = dat, p = 0.01, N=100, plot = TRUE)
#' gev.tem('Nquant', dat = dat, q = 0.5, N=100, plot = TRUE)
#' }
gev.tem <-
  function(
    param = c(
      "loc",
      "scale",
      "shape",
      "quant",
      "Nmean",
      "Nquant"
    ),
    dat,
    psi = NULL,
    p = NULL,
    q = 0.5,
    N = NULL,
    n.psi = 50,
    plot = TRUE,
    correction = TRUE
  ) {
    if (param %in% c("VaR", "retlev")) {
      param <- "quant"
    } #Compatibility following change of notation 13/07/2017
    if (param == "scale" && !is.null(psi)) {
      if (isTRUE(any(psi < 0))) {
        stop("Invalid argument: scale parameter must be positive")
      }
    }
    tem_out <-
      gev.pll(
        psi = psi,
        param = param,
        mod = "tem",
        dat = dat,
        N = N,
        p = p,
        q = q,
        correction = correction
      )
    if (plot) {
      plot.fr(tem_out)
    }
    return(invisible(tem_out))
  }

#' Tangent exponential model approximation for the GP distribution
#'
#' The function \code{gpd.tem} provides a tangent exponential model (TEM) approximation
#' for higher order likelihood inference for a scalar parameter for the generalized Pareto distribution. Options include
#' scale and shape parameters as well as value-at-risk (also referred to as quantiles, or return levels)
#' and expected shortfall. The function attempts to find good values for \code{psi} that will
#' cover the range of options, but the fit may fail and return an error. In such cases, the user can try to find good
#' grid of starting values and provide them to the routine.
#'
#' As of version 1.11, this function is a wrapper around \code{gpd.pll}.
#'
#' @details The interpretation for \code{m} is as follows: if there are on average \eqn{m_y} observations per year above the threshold, then  \eqn{m = Tm_y} corresponds to \eqn{T}-year return level.
#'
#' @param param parameter over which to profile
#' @param thresh threshold value corresponding to the lower bound of the support or the location parameter of the generalized Pareto distribution.
#' @param psi scalar or ordered vector of values for the interest parameter. If \code{NULL} (default), a grid of values centered at the MLE is selected. If \code{psi} is of length 2 and \code{n.psi}>2, it is assumed to be the minimal and maximal values at which to evaluate the profile log likelihood.
#' @param m number of observations of interest for return levels. See \strong{Details}. Required only for \code{param = 'VaR'} or \code{param = 'ES'}.
#' @param N size of block over which to take maxima. Required only for \code{args} \code{Nmean} and \code{Nquant}.
#' @param p tail probability, equivalent to \eqn{1/m}. Required only for \code{args} \code{quant}.
#' @param q level of quantile for N-block maxima. Required only for \code{args} \code{Nquant}.
#' @param dat sample vector for the GP distribution
#' @param n.psi number of values of \code{psi} at which the likelihood is computed, if \code{psi} is not supplied (\code{NULL}). Odd values are more prone to give rise to numerical instabilities near the MLE
#' @param plot logical indicating whether \code{plot.fr} should be called upon exit
#' @param correction logical indicating whether \link{spline.corr} should be called.
#' @param ... additional arguments, for backward compatibility
#' @author Leo Belzile
#' @return an invisible object of class \code{fr} (see \code{tem} in package \code{hoa}) with elements
#' \itemize{
#' \item \code{normal}: maximum likelihood estimate and standard error of the interest parameter \eqn{\psi}
#' \item \code{par.hat}: maximum likelihood estimates
#' \item \code{par.hat.se}: standard errors of maximum likelihood estimates
#' \item \code{th.rest}: estimated maximum profile likelihood at (\eqn{\psi}, \eqn{\hat{\lambda}})
#' \item \code{r}: values of likelihood root corresponding to \eqn{\psi}
#' \item \code{psi}: vector of interest parameter
#' \item \code{q}: vector of likelihood modifications
#' \item \code{rstar}: modified likelihood root vector
#' \item \code{rstar.old}: uncorrected modified likelihood root vector
#' \item \code{param}: parameter
#' }
#' @export
#' @examples
#' set.seed(123)
#' dat <- rgp(n = 40, scale = 1, shape = -0.1)
#' #with plots
#' m1 <- gpd.tem(param = 'shape', n.psi = 50, dat = dat, plot = TRUE)
#' \dontrun{
#' m2 <- gpd.tem(param = 'scale', n.psi = 50, dat = dat)
#' m3 <- gpd.tem(param = 'VaR', n.psi = 50, dat = dat, m = 100)
#' #Providing psi
#' psi <- c(seq(2, 5, length = 15), seq(5, 35, length = 45))
#' m4 <- gpd.tem(param = 'ES', dat = dat, m = 100, psi = psi, correction = FALSE)
#' mev:::plot.fr(m4, which = c(2, 4))
#' plot(fr4 <- spline.corr(m4))
#' confint(m1)
#' confint(m4, parm = 2, warn = FALSE)
#' m5 <- gpd.tem(param = 'Nmean', dat = dat, N = 100, psi = psi, correction = FALSE)
#' m6 <- gpd.tem(param = 'Nquant', dat = dat, N = 100, q = 0.7, correction = FALSE)
#' }
gpd.tem <-
  function(
    dat,
    param = c(
      "scale",
      "shape",
      "quant",
      "VaR",
      "retlev",
      "ES",
      "Nmean",
      "Nquant"
    ),
    psi = NULL,
    m = NULL,
    thresh = 0,
    n.psi = 50,
    N = NULL,
    p = NULL,
    q = NULL,
    plot = FALSE,
    correction = TRUE,
    ...
  ) {
    args <- list(...)
    if (!is.null(args$threshold)) {
      thresh <- args$threshold
    }
    tem <-
      gpd.pll(
        psi = psi,
        param = param,
        mod = "tem",
        dat = dat,
        N = N,
        m = m,
        mle = NULL,
        q = q,
        p = p,
        thresh = thresh,
        correction = correction
      )
    if (isTRUE(plot)) {
      plot.fr(tem)
    }
    return(invisible(tem))
  }
