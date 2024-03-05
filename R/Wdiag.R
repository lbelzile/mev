#####################################################################
#### Functions to create plots and tests in article
#### Jennifer L. Wadsworth (2016), Technometrics, 'Exploiting structure of maximum likelihood estimators for extreme value threshold selection'
#### Code by J.L. Wadsworth

#' Wadsworth's univariate and bivariate exponential threshold diagnostics
#'
#' Function to produce diagnostic plots and test statistics for the
#' threshold diagnostics exploiting structure of maximum likelihood estimators
#' based on the non-homogeneous Poisson process likelihood
#'
#' @param xdat a numeric vector of data to be fitted.
#' @param model string specifying whether the univariate or bivariate diagnostic should be used. Either \code{nhpp}
#' for the univariate model, \code{exp} (\code{invexp}) for the bivariate exponential model with rate (inverse rate) parametrization. See details.
#' @param u optional; vector of candidate thresholds.
#' @param k number of thresholds to consider (if \code{u} unspecified).
#' @param q1 lowest quantile for the threshold sequence.
#' @param q2 upper quantile limit for the threshold sequence (\code{q2} itself is not used as a threshold,
#'  but rather the uppermost threshold will be at the \eqn{(q_2-1/k): q2-1/k} quantile).
#' @param par parameters of the NHPP likelihood. If \code{missing}, the \code{\link[mev]{fit.pp}} routine will be run to obtain values
#' @param M number of superpositions or 'blocks' / 'years' the process corresponds to (can affect the optimization)
#' @param nbs number of simulations used to assess the null distribution of the LRT, and produce the p-value
#' @param alpha significance level of the LRT
#' @param plots vector of strings indicating which plots to produce; \code{LRT}= likelihood ratio test, \code{WN} = white noise, \code{PS} = parameter stability. Use \code{NULL} if you do not want plots to be produced
#' @param UseQuantiles logical; use quantiles as the thresholds in the plot?
#' @param changepar logical; if \code{TRUE}, the graphical parameters (via a call to \code{par}) are modified.
#' @param ... additional parameters passed to \code{plot}, overriding defaults including
#'
#' @details The function is a wrapper for the univariate (non-homogeneous Poisson process model) and bivariate exponential dependence model.
#' For the latter, the user can select either the rate or inverse rate parameter  (the inverse rate parametrization  works better for uniformity
#' of the p-value distribution under the \code{LR} test.
#'
#' There are two options for the bivariate diagnostic: either provide pairwise minimum of marginally
#' exponentially distributed margins or provide a \code{n} times 2 matrix with the original data, which
#' is transformed to exponential margins using the empirical distribution function.
#'
#' @references Wadsworth, J.L. (2016). Exploiting Structure of Maximum Likelihood Estimators for Extreme Value Threshold Selection, \emph{Technometrics}, \bold{58}(1), 116-126, \code{http://dx.doi.org/10.1080/00401706.2014.998345}.
#'
#' @author Jennifer L. Wadsworth
#' @return plots of the requested diagnostics and an invisible list with components
#' \itemize{
#' \item \code{MLE}: maximum likelihood estimates from all thresholds
#' \item \code{Cov}: joint asymptotic covariance matrix for \eqn{\xi}, \eqn{\eta} or \eqn{1/\eta}.
#' \item \code{WN}: values of the white noise process
#' \item \code{LRT}: values of the likelihood ratio test statistic vs threshold
#' \item \code{pval}: \emph{P}-value of the likelihood ratio test
#' \item \code{k}: final number of thresholds used
#' \item \code{thresh}: threshold selected by the likelihood ratio procedure
#' \item \code{qthresh}: quantile level of threshold selected by the likelihood ratio procedure
#' \item \code{cthresh}: vector of candidate thresholds
#' \item \code{qcthresh}: quantile level of candidate thresholds
#' \item \code{mle.u}: maximum likelihood estimates for the selected threshold
#' \item \code{model}: model fitted
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' # Parameter stability only
#' W.diag(xdat = abs(rnorm(5000)), model = 'nhpp',
#'        k = 30, q1 = 0, plots = "PS")
#' W.diag(rexp(1000), model = 'nhpp', k = 20, q1 = 0)
#' xbvn <- mvrnorm(n = 6000,
#'                 mu = rep(0, 2),
#'                 Sigma = cbind(c(1, 0.7), c(0.7, 1)))
#' # Transform margins to exponential manually
#' xbvn.exp <- -log(1 - pnorm(xbvn))
#' #rate parametrization
#' W.diag(xdat = apply(xbvn.exp, 1, min), model = 'exp',
#'        k = 30, q1 = 0)
#' W.diag(xdat = xbvn, model = 'exp', k = 30, q1 = 0)
#' #inverse rate parametrization
#' W.diag(xdat = apply(xbvn.exp, 1, min), model = 'invexp',
#'        k = 30, q1 = 0)
#' }
#' @export
W.diag <- function(xdat,
                   model = c("nhpp", "exp", "invexp"),
                   u = NULL,
                   k,
                   q1 = 0,
                   q2 = 1,
                   par = NULL,
                   M = NULL,
                   nbs = 1000,
                   alpha = 0.05,
                   plots = c("LRT", "WN", "PS"),
                   UseQuantiles = FALSE,
                   changepar = TRUE,
                   ...) {
    stopifnot(is.logical(UseQuantiles),
              length(UseQuantiles) == 1L,
              is.logical(changepar),
              length(changepar) == 1L)
  model <- match.arg(model)
  if(!is.null(plots)){
    plots <- match.arg(plots,
                       choices = c("LRT", "WN", "PS"),
                       several.ok = TRUE)
  }
  if (ncol(as.matrix(xdat)) == 2 &&
      model %in% c("exp", "invexp")) {
    xdat <- -log(1 - apply(xdat, 2, function(y) {
      rank(y, ties.method = "average") / (length(y) + 1)
    }))
    xdat <- pmin(xdat[, 1], xdat[, 2])
  }
  if (ncol(as.matrix(xdat)) != 1) {
    stop("Invalid input for \"xdat\"")
  }
  switch(
    model,
    nhpp = .NHPP.diag(
      xdat = xdat,
      u = u,
      k = k,
      q1 = q1,
      q2 = q2,
      par = par,
      M = M,
      nbs = nbs,
      alpha = alpha,
      plots = plots,
      UseQuantiles = UseQuantiles,
      changepar = changepar,
      ...
    ),
    exp = .Expl.diag(
      x = xdat,
      u = u,
      k = k,
      q1 = q1,
      q2 = q2,
      nbs = nbs,
      alpha = alpha,
      plots = plots,
      UseQuantiles = UseQuantiles,
      param = "Rate",
      changepar = changepar,
      ...
    ),
    invexp = .Expl.diag(
      x = xdat,
      u = u,
      k = k,
      q1 = q1,
      q2 = q2,
      nbs = nbs,
      alpha = alpha,
      plots = plots,
      UseQuantiles = UseQuantiles,
      param = "InvRate",
      changepar = changepar,
      ...
    )
  )
}


.NHPP.diag <-
  function(xdat,
           u = NULL,
           k,
           q1 = 0,
           q2 = 1,
           par = NULL,
           M = NULL,
           nbs = 1000,
           alpha = 0.05,
           plots = c("LRT", "WN", "PS"),
           UseQuantiles = TRUE,
           changepar = changepar,
           ...) {
    unull <- is.null(u)
    if (unull) {
      thresh <- quantile(xdat, q1)
    } else {
      stopifnot(length(u) > 1)
      thresh <- min(u)
    }
    if (!unull) {
      k <- length(u)
    }
    if (is.null(M)) {
      M <- length(xdat[xdat > thresh])
    }  #why M=nat/3 as default?
    if (is.null(par)) {
      ppf <-
        fit.pp(
          xdat = xdat,
          threshold = quantile(xdat, q1),
          npp = length(xdat) / M,
          show = FALSE
        )
      par <- ppf$estimate
    }
    J1 <-
      .Joint_MLE_NHPP(
        x = xdat,
        u = u,
        k = k,
        q1 = q1,
        q2 = q2,
        par = par,
        M = M
      )
    warn <-
      any(eigen(J1$Cov.xi, only.values = TRUE)$val <= .Machine$double.eps)
    if (!unull && warn) {
      stop("Estimated covariance matrix for xi not positive definite: try different thresholds")
    }

    while (any(eigen(J1$Cov.xi, only.values = TRUE)$val <= .Machine$double.eps)) {
      k <- k - 1
      J1 <-
        .Joint_MLE_NHPP(
          x = xdat,
          k = k,
          q1 = q1,
          q2 = q2,
          par = par,
          M = M
        )
    }
    if (warn) {
      warning(
        paste(
          "Estimated covariance matrix for xi not positive definite for initial k. Final k:",
          k
        )
      )
    }

    if (unull) {
      u <- quantile(xdat, seq(q1, q2, len = k + 1))
    }

    wn <-
      .C1(k) %*% J1$mle[, 3] / sqrt(diag(.C1(k) %*% J1$Cov.xi %*% t(.C1(k))))
    nl <- .norm_LRT(x = wn, u = u[-c(1, k + 1)])

    nlt <- NULL
    for (j in 1:nbs) {
      nlt[j] <- max(.norm_LRT(x = rnorm(k - 1), u[-c(1, k + 1)])[, 2])
    }

    pval <- length(nlt[nlt > max(nl[, 2])]) / nbs

    if (pval < alpha) {
      ustar <- nl[nl[, 2] == max(nl[, 2]), 1]
    } else {
      ustar <- min(u)
    }
    ind <- u[-(k + 1)] == ustar
    theta.hat <- J1$mle[ind,]
    if (isTRUE(unull)) {
      qs <- seq(q1, q2, len = k + 1)[-(k + 1)]
      qlthresh <- mean(xdat <= ustar)
    } else{
      qs <- NULL
      qlthresh <- NULL
    }
  #   if(!is.null(plots)){
  #   # Copy graphical elements from ellipsis
  #   if (changepar) {
  #     old.par <- par(no.readonly = TRUE)
  #     on.exit(par(old.par))
  #     par(mfrow = c(length(plots), 1),
  #         mar = c(4.5, 4.5, 0.1, 0.1))
  #   }
  #   if (is.element("LRT", plots)) {
  #     if (!UseQuantiles) {
  #       plot(
  #         qs,
  #         c(rep(NA, 2), nl[, 2]),
  #         xlab = "quantile",
  #         ylab = "likelihood ratio statistic",
  #         main = paste("p-value:", pval),
  #         ...
  #       )
  #     } else {
  #       plot(
  #         u[-c(k + 1)],
  #         c(rep(NA, 2), nl[, 2]),
  #         bty = "l",
  #         xlab = "threshold",
  #         ylab = "likelihood ratio statistic",
  #         main = paste("p-value:", pval),
  #         ...
  #       )
  #     }
  #   }
  #
  #   if (is.element("WN", plots)) {
  #     if (!UseQuantiles) {
  #       plot(qs,
  #            c(NA, wn),
  #            xlab = "quantile",
  #            ylab = "white noise",
  #            bty = "l",
  #            ...)
  #       abline(h = 0, col = 2)
  #       abline(v = mean(xdat <= ustar), col = 4)
  #     } else {
  #       plot(u[-c(k + 1)],
  #            c(NA, wn),
  #            xlab = "threshold",
  #            ylab = "white noise",
  #            bty = "l",
  #            ...)
  #       abline(h = 0, col = 2)
  #       abline(v = ustar, col = 4)
  #     }
  #   }
  #
  #   if (is.element("PS", plots)) {
  #     TradCI <-
  #       cbind(J1$mle[, 3] - qnorm(0.975) * sqrt(diag(J1$Cov.xi)),
  #             J1$mle[, 3] + qnorm(0.975) * sqrt(diag(J1$Cov.xi)))
  #     if (!UseQuantiles) {
  #       plot(
  #         qs,
  #         J1$mle[, 3],
  #         ylim = c(min(TradCI[, 1]), max(TradCI[, 2])),
  #         xlab = "quantile",
  #         bty = "l",
  #         ylab = "shape",
  #         ...
  #       )
  #       lines(qs, TradCI[, 1], lty = 2)
  #       lines(qs, TradCI[, 2], lty = 2)
  #       abline(v = mean(xdat <= ustar), col = 4)
  #     } else {
  #       plot(
  #         u[-(k + 1)],
  #         J1$mle[, 3],
  #         ylim = c(min(TradCI[, 1]), max(TradCI[, 2])),
  #         bty = "l",
  #         xlab = "threshold",
  #         ylab = "shape",
  #         ...
  #       )
  #       lines(u[-(k + 1)], TradCI[, 1], lty = 2)
  #       lines(u[-(k + 1)], TradCI[, 2], lty = 2)
  #       abline(v = ustar, col = 4)
  #     }
  #   }
  # }
    colnames(J1$mle) <- names(theta.hat) <- c("location", "scale", "shape")
    colnames(J1$Cov.xi) <- rownames(J1$Cov.xi) <- NULL
    ret_list <- list(
      MLE = J1$mle,
      Cov = J1$Cov.xi,
      WN = as.numeric(wn),
      LRT = nl,
      pval = as.numeric(pval),
      k = as.integer(k),
      thresh = as.numeric(ustar),
      qthresh = qlthresh,
      cthresh = as.numeric(u),
      qcthresh = qs,
      mle.u = theta.hat,
      model = "nhpp"
    )
    class(ret_list) <- "mev_thdiag_wadsworth"
    if(!is.null(plots)){
      plot(ret_list,
           plots = plots,
           changepar = changepar,
           UseQuantiles = UseQuantiles,
           ...)
    }
    invisible(ret_list)
  }


#############################################################################################################

.Expl.diag <-
  function(x,
           u = NULL,
           k,
           q1,
           q2 = 1,
           nbs = 1000,
           alpha = 0.05,
           plots = c("LRT", "WN", "PS"),
           UseQuantiles = TRUE,
           param = "InvRate",
           changepar = TRUE,
           ...) {
    unull <- is.null(u)
    if (!unull) {
      k <- length(u)
    }
    J1 <-
      .Joint_MLE_Expl(
        x = x,
        u = u,
        k = k,
        q1 = q1,
        q2 = q2,
        param = param
      )
    warn <- any(eigen(J1$Cov)$val <= .Machine$double.eps)
    if (!unull && warn) {
      stop(
        "Estimated covariance matrix for eta not positive definite: try different thresholds"
      )
    }

    while (any(eigen(J1$Cov)$val <= .Machine$double.eps)) {
      k <- k - 1
      J1 <-
        .Joint_MLE_Expl(
          x = x,
          k = k,
          q1 = q1,
          q2 = q2,
          param = param
        )
    }
    if (warn) {
      warning(
        paste(
          "Estimated covariance matrix for 1/eta not positive definite for initial k. Final k:",
          k
        )
      )
    }
    if (unull & !UseQuantiles) {
      u <- quantile(x, seq(q1, 1, len = k + 1))
    }
    wn <-
      .C1(k) %*% J1$mle / sqrt(diag(.C1(k) %*% J1$Cov %*% t(.C1(k))))
    nl <- .norm_LRT(x = wn, u = u[-c(1, k + 1)])

    nlt <- NULL
    for (j in seq_len(nbs)) {
      nlt[j] <- max(.norm_LRT(x = rnorm(k - 1), u[-c(1, k + 1)])[, 2])
    }

    pval <- length(nlt[nlt > max(nl[, 2])]) / nbs

    if (pval < alpha) {
      ustar <- nl[nl[, 2] == max(nl[, 2]), 1]
    } else {
      ustar <- min(u)
    }
    ind <- u[-(k + 1)] == ustar
    theta.hat <- J1$mle[ind]

    if (isTRUE(unull)) {
      qs <- as.numeric(seq(q1, q2, len = k + 1)[-(k + 1)])
      qlthresh <- mean(x <= ustar)
    } else{
      qs <- NULL
      qlthresh <- NULL
    }
    # if(!is.null(plots)){
    #   if (changepar) {
    #     old.par <- par(no.readonly = TRUE)
    #     on.exit(par(old.par))
    #     par(mfrow = c(length(plots), 1),
    #         mar = c(4.5, 4.5, 0.1, 0.1))
    #   }
    #   if (is.element("LRT", plots)) {
    #     if (unull && UseQuantiles) {
    #       plot(
    #         qs,
    #         c(NA, NA, nl[, 2]),
    #         bty = "l",
    #         xlab = "quantile",
    #         ylab = "likelihood ratio statistic",
    #         main = paste("p-value:", pval),
    #         ...
    #       )
    #     } else {
    #       plot(
    #         u[-c(k + 1)],
    #         c(NA, NA, nl[, 2]),
    #         bty = "l",
    #         xlab = "threshold",
    #         ylab = "likelihood ratio statistic",
    #         main = paste("p-value:", pval),
    #         ...
    #       )
    #     }
    #   }
    #
    #   if (is.element("WN", plots)) {
    #     if (unull && UseQuantiles) {
    #       plot(qs,
    #            c(NA, wn),
    #            xlab = "quantile",
    #            ylab = "white noise",
    #            bty = "l",
    #            ...)
    #       abline(h = 0, col = 2)
    #       abline(v = mean(x <= ustar), col = 4)
    #     } else {
    #       plot(u[-c(k + 1)],
    #            c(NA, wn),
    #            xlab = "threshold",
    #            ylab = "white noise",
    #            bty = "l",
    #            ...)
    #       abline(h = 0, col = 2)
    #       abline(v = ustar, col = 4)
    #     }
    #   }
    #
    #   if (is.element("PS", plots)) {
    #     TradCI <-
    #       cbind(J1$mle - qnorm(0.975) * sqrt(diag(J1$Cov)),
    #             J1$mle + qnorm(0.975) * sqrt(diag(J1$Cov)))
    #     if (UseQuantiles) {
    #       if (param == "InvRate") {
    #         plot(
    #           qs,
    #           J1$mle,
    #           ylim = c(min(TradCI[, 1]), max(TradCI[, 2])),
    #           bty = "l",
    #           xlab = "quantile",
    #           ylab = expression(hat(eta)),
    #           ...
    #         )
    #       } else if (param == "Rate") {
    #         plot(
    #           qs,
    #           J1$mle,
    #           ylim = c(min(TradCI[, 1]), max(TradCI[, 2])),
    #           bty = "l",
    #           xlab = "quantile",
    #           ylab = expression(hat(theta)),
    #           ...
    #         )
    #       }
    #       lines(qs, TradCI[, 1], lty = 2)
    #       lines(qs, TradCI[, 2], lty = 2)
    #       abline(v = mean(x <= ustar), col = 4)
    #     } else {
    #       if (param == "InvRate") {
    #         plot(
    #           u[-(k + 1)],
    #           J1$mle,
    #           bty = "l",
    #           ylim = c(min(TradCI[, 1]), max(TradCI[, 2])),
    #           xlab = "threshold",
    #           ylab = expression(hat(eta)),
    #           ...
    #         )
    #       } else if (param == "InvRate") {
    #         plot(
    #           u[-(k + 1)],
    #           J1$mle,
    #           bty = "l",
    #           ylim = c(min(TradCI[, 1]), max(TradCI[, 2])),
    #           xlab = "threshold",
    #           ylab = expression(hat(theta)),
    #           ...
    #         )
    #       }
    #       lines(u[-(k + 1)], TradCI[, 1], lty = 2)
    #       lines(u[-(k + 1)], TradCI[, 2], lty = 2)
    #       abline(v = ustar, col = 4)
    #     }
    #   }
    # }
    ret_list <- list(
        MLE = J1$mle,
        Cov = J1$Cov,
        WN = as.numeric(wn),
        LRT = nl,
        pval = as.numeric(pval),
        k = as.integer(k),
        thresh = as.numeric(ustar),
        qthresh = qlthresh,
        cthresh = as.numeric(u),
        qcthresh = qs,
        mle.u = as.numeric(theta.hat),
        model = switch(param,
                       "InvRate" = "invexp",
                       "Rate" = "exp")
      )
    class(ret_list) <- "mev_thdiag_wadsworth"
    if(!is.null(plots)){
      plot(ret_list,
           plots = plots,
           changepar = changepar,
           UseQuantiles = UseQuantiles,
           ...)
    }
    invisible(ret_list)
  }







#######################################################################################################

#' Joint maximum likelihood estimation for exponential model
#'
#'
#' Calculates the MLEs of the rate parameter, and joint asymptotic covariance matrix of these MLEs
#' over a range of thresholds as supplied by the user.
#'
#' @param x vector of data
#' @param u vector of thresholds. If not supplied, then \code{k}
#' thresholds between quantiles (\code{q1}, \code{q2}) will be used
#' @param k number of thresholds to consider if u not supplied
#' @param q1 lower quantile to consider for threshold
#' @param q2 upper quantile to consider for threshold
#' @param param character specifying \code{'InvRate'} or \code{'Rate'}
#' for either inverse rate parameter / rate parameter, respectively
#'
#' @author Jennifer L. Wadsworth
#'
#' @return a list with
#' \itemize{
#' \item mle vector of MLEs above the supplied thresholds
#' \item cov joint asymptotic covariance matrix of these MLEs
#' }
#' @keywords internal
.Joint_MLE_Expl <- function(x,
                            u = NULL,
                            k,
                            q1,
                            q2 = 1,
                            param) {
  if (!is.element(param, c("InvRate", "Rate"))) {
    stop("param should be one of InvRate or Rate")
  }
  if (!is.null(u)) {
    k <- length(u)
    x <- x[x > u[1]]
    # add threshold above all data, to 'close' final region
    u <- c(u, max(x) + 1)
  } else {
    u <- quantile(x, seq(q1, q2, len = k + 1))
  }

  I <- n <- m <- thetahat <- NULL

  for (i in 1:k) {
    if (param == "InvRate") {
      thetahat[i] <- mean(x[x >= u[i]] - u[i])
    } else if (param == "Rate") {
      thetahat[i] <- 1 / mean(x[x >= u[i]] - u[i])
    }
    n[i] <- length(x[x >= u[i] & x <= u[i + 1]])
  }
  for (i in 1:k) {
    m[i] <- sum(n[i:k])
    I[i] <- 1 / thetahat[i] ^ 2
  }

  Tcov <- matrix(0, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      Tcov[i, j] <- 1 / (I[min(i, j)] * m[min(i, j)])
    }
  }
  CovT <- Tcov
  return(list(mle = thetahat, Cov = CovT))
}


#####################################################################################


#' Joint maximum likelihood for the non-homogeneous Poisson Process
#'
#' Calculates the MLEs of the parameters (\eqn{\mu}, \eqn{\sigma}, \eqn{\xi}), and joint
#' asymptotic covariance matrix of these MLEs over a range of thresholds as supplied by the user.
#' @param x vector of data
#' @param u optional vector of thresholds. If not supplied, then k thresholds between quantiles (q1, q2) will be used
#' @param k number of thresholds to consider if \code{u} not supplied
#' @param q1 lower  quantile to consider for threshold
#' @param q2 upper quantile to consider for threshold. Default to 1
#' @param par starting values for the optimization
#' @param  M  number of superpositions or 'blocks' / 'years' the process corresponds to.
#' It affects the estimation of \eqn{mu} and \eqn{sigma},
#' but these can be changed post-hoc to correspond to any number)
#'
#' @author Jennifer L. Wadsworth
#' @return a list with components
#' \itemize{
#' \item mle matrix of MLEs above the supplied thresholds; columns are (\eqn{\mu}, \eqn{\sigma}, \eqn{\xi})
#' \item Cov.all joint asymptotic covariance matrix of all MLEs
#' \item Cov.mu joint asymptotic covariance matrix of MLEs for \eqn{\mu}
#' \item Cov.sig joint asymptotic covariance matrix of MLEs for \eqn{\sigma}
#' \item Cov.xi joint asymptotic covariance matrix of MLEs for \eqn{\xi}
#' }
#' @keywords internal
.Joint_MLE_NHPP <- function(x,
                            u = NULL,
                            k,
                            q1,
                            q2 = 1,
                            par,
                            M) {
  if (!is.null(u)) {
    k <- length(u)
    x <- x[x > u[1]]
    # add threshold above all data, to 'close' final region
    u <- c(u, max(x) + 1)
  } else {
    u <- quantile(x, seq(q1, q2, len = k + 1))
  }

  I <- Iinv <- list()
  thetahat <- matrix(NA, ncol = 3, nrow = k)

  for (i in 1:k) {
    opt <- fit.pp(xdat = x,
                  threshold = u[i],
                  np = M)
    thetahat[i,] <- opt$estimate

    ### Deal with xi <- 0.5
    if (thetahat[i, 3] > -0.5) {
      I[[i]] <-
        pp.infomat(
          par = opt$estimate,
          u = u[i],
          np = M,
          method = "exp",
          nobs = 1
        )
      Iinv[[i]] <- solve(I[[i]])
    } else {
      I[[i]] <- Iinv[[i]] <- matrix(0, 3, 3)
    }
  }

  Wcov <- list()
  Wcov1 <- NULL
  for (i in 1:k) {
    Wcov[[i]] <- matrix(0, 3, 3)
    for (j in 1:k) {
      Wcov[[i]] <- cbind(Wcov[[i]], Iinv[[min(i, j)]])
    }
    Wcov1 <- rbind(Wcov1, Wcov[[i]])
  }
  Wcov1 <- Wcov1[, -c(1:3)]

  CovT <- Wcov1
  Cov.mu <- CovT[seq(1, 3 * k, by = 3), seq(1, 3 * k, by = 3)]
  Cov.sig <- CovT[seq(2, 3 * k, by = 3), seq(2, 3 * k, by = 3)]
  Cov.xi <- CovT[seq(3, 3 * k, by = 3), seq(3, 3 * k, by = 3)]

  return(list(
    mle = thetahat,
    Cov.all = CovT,
    Cov.mu = Cov.mu,
    Cov.sig = Cov.sig,
    Cov.xi = Cov.xi
  ))
}


###################################################################################

# norm_LRT

# Details:

# Evaluates the likelihood ratio statistics for testing white noise

# Arguments:

# x - vector of white noise process (WNP, usually normalized estimates of \eqn{xi} or the exponential rate parameter
# \eqn{1/\eta}) u - vector of thresholds that are associated to the WNP


.norm_LRT <- function(x, u) {
  l <- length(u)
  v <-
    u[-c(1)]  # means two or more obs available for std dev calculation
  lr <- NULL
  for (i in 1:length(v)) {
    n1 <- length(x[u <= v[i]])
    num <-
      .nll_norm(theta = c(mean(x[u <= v[i]]), sd(x[u <= v[i]]) * sqrt((n1 - 1) /
                                                                        n1)), x = x[u <= v[i]])
    den <- .nll_norm(theta = c(0, 1), x = x[u <= v[i]])
    lr[i] <- -2 * (num - den)
  }
  return(cbind(v, lr))
}



###################################################################################

# nll_norm - negative log likelihood for the normal distribution

.nll_norm <- function(theta, x) {
  if (theta[2] < 0) {
    return(1e+11)
  } else {
    return(-sum(dnorm(
      x,
      mean = theta[1],
      sd = theta[2],
      log = TRUE
    )))
  }
}



###################################################################################

#' Contrast matrix
#'
#' Produces a contrast matrix with (1,-1) elements running down the two diagonals
#'
#'@param k number of columns (the number of rows is \code{k-1})
#'
#'@return a \code{k-1} x \code{k} contrast matrix
#'@keywords internal
.C1 <- function(k) {
  C <- diag(x = 1, nrow = k - 1, ncol = k)
  C[row(C) + 1 == col(C)] <- -1
  return(C)
}

#' @export
plot.mev_thdiag_wadsworth <-
  function(x,
           plots = c("LRT", "WN", "PS"),
           ...) {
    args <- list(...)
    args$`...` <- NULL #probably not needed anymore
    if(is.null(args$UseQuantiles)){
      UseQuantiles <- FALSE
    } else{
      UseQuantiles <- args$UseQuantiles
      args$UseQuantiles <- NULL
      if(is.null(x$qlevel)){
        UseQuantiles <- FALSE
      }
    }
    model <- match.arg(arg = x$model,
                       choices = c("nhpp","exp", "invexp"),
                       several.ok = FALSE)
    plots <- match.arg(plots,
                       choices = c("LRT", "WN", "PS"),
                       several.ok = TRUE)
    if(length(plots) < 1){
      stop("No choice selected; aborting.")
    }
    # Copy graphical elements from ellipsis
    if (isTRUE(args$changepar)) {
      old.par <- par(no.readonly = TRUE)
      on.exit(par(old.par))
      par(mfrow = c(length(plots), 1),
          mar = c(4.5, 4.5, 0.5, 0.5))
    }
    args$changepar <- NULL
    if(!UseQuantiles){
      xp <- x$cthresh[-c(x$k + 1)]
      xlab <- "threshold"
    } else{
      xp <- x$qthresh
      xlab <- "quantile"
    }
    args$x <- args$y <- args$xlab <- args$ylab <- NULL
    if(is.null(args$bty)){
      args$bty = "l"
    }
    if (is.element("LRT", plots)) {
      do.call(what = plot,
              args = c(list(
        x = xp,
        y = c(rep(NA, 2), x$LRT[, 2]),
        xlab = xlab,
        ylab = "likelihood ratio"
              ), args)
      )
      mtext(cex = 0.8,
            text = paste("p-value:", format.pval(pv = x$pval, eps = 1e-4)),
            side = 3,
            adj = 1)
    }

    if (is.element("WN", plots)) {
      do.call(what = plot,
              args = c(list(
                x = xp,
                 y = c(NA, x$WN),
                 xlab = xlab,
                 ylab = "white noise"),
                args))
      abline(h = 0, col = 2)
      abline(v = ifelse(UseQuantiles,
                        x$qthresh,
                        x$thresh),
             col = 4)
    }
    if (is.element("PS", plots)) {
      col <- switch(model, nhpp = 3, exp = 1, invexp = 1)
      TradCI <-
        cbind(as.matrix(x$MLE)[, col] - qnorm(0.975) * sqrt(diag(x$Cov)),
              as.matrix(x$MLE)[, col] + qnorm(0.975) * sqrt(diag(x$Cov)))
      do.call(plot, args = c(list(
        x = xp,
        y = as.matrix(x$MLE)[, col],
        ylim = c(min(TradCI[, 1]), max(TradCI[, 2])),
        xlab = xlab,
        ylab = switch(model,
                      "nhpp" = "shape",
                      "invexp" = expression(hat(eta)),
                      "exp" = expression(hat(theta)))),
        args))
      lines(xp, TradCI[, 1], lty = 2)
      lines(xp, TradCI[, 2], lty = 2)
      abline(v = ifelse(UseQuantiles,
                        x$qthresh,
                        x$thresh),
             col = 4)
    }
    return(invisible(NULL))
  }


#' @export
print.mev_thdiag_wadsworth <-
  function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("Threshold selection method: Wadsworth's white noise test\n based on sequential Poisson process superposition")
    cat(switch(x$model,
               "nhpp" = "inhomogeneous Poisson process (shape)",
               "invexp" = "coefficient of tail dependence \n(exponential, reciprocal rate)",
               "exp" = "coefficient of tail dependence \n(exponential, rate)"), "\n")
    cat("Selected threshold:", round(x$thresh, digits), "\n")
    return(invisible(NULL))
  }
