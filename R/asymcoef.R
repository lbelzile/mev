#' Coefficient of extremal asymmetry
#'
#' This function implements estimators of the bivariate
#'  coefficient of extremal asymmetry proposed in
#'  Semadeni's (2021) PhD thesis.
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
#' @param  estimator string indicating the estimation method, one of \code{empirical} or empirical likelihood (\code{emplik})
#' @param confint string for the method used to derive confidence intervals, either \code{none} (default) or a nonparametric \code{bootstrap}
#' @param level probability level for confidence intervals, default to 0.95 or bounds for the interval
#' @param B integer; number of bootstrap replicates (if applicable)
#' @param ties.method string; method for handling ties. See the documentation of \link[base]{rank} for available options.
#' @param plot logical; if \code{TRUE}, return a plot.
#' @param ... additional parameters for plots
#' @return an invisible data frame with columns
#' \describe{
#' \item{\code{threshold}}{vector of thresholds on the probability scale}
#' \item{\code{coef}}{extremal asymmetry coefficient estimates}
#' \item{\code{confint}}{either \code{NULL} or a matrix with two columns containing the lower and upper bounds for each threshold}
#' }
#' @references Semadeni, C. (2020). Inference on the Angular Distribution of Extremes, PhD thesis, EPFL, no. 8168.
#' @examples
#' \dontrun{
#' samp <- rmev(n = 1000,
#'              d = 2,
#'              param = 0.2,
#'              model = "log")
#' xasym(samp, confint = "wald")
#' xasym(samp, method = "emplik")
#' }
#' @export
#' @keywords internal
xasym <- function(
  xdat,
  qlev = NULL,
  nq = 40,
  qlim = c(0.8, 0.99),
  estimator = c("emp", "elik"),
  confint = c("none", "wald", "bootstrap"),
  level = 0.95,
  B = 999L,
  # trunc = TRUE,
  ties.method = "random",
  plot = TRUE,
  ...
) {
  args <- list(...)
  if (missing(xdat) & !is.null(args$data)) {
    xdat <- args$data
  } else if (missing(xdat)) {
    stop(
      "Argument \"xdat\" missing: must be a matrix of observations."
    )
  }
  if (!is.null(args$method)) {
    # 2025.10.16 backward compatibility
    estimator <- args$method
  }
  method <- match.arg(
    estimator,
    choices = c("empirical", "emplik", "emp", "elik")
  )
  if (method == "empirical") {
    method <- "emp"
  }
  if (method == "emplik") {
    method = "elik"
  }
  data <- xdat
  if (!is.null(args$u)) {
    qlev <- args$u
  }
  u <- qlev
  if (inherits(data, "data.frame")) {
    data <- as.matrix(data)
  }
  if (!inherits(data, "matrix")) {
    stop("\"data\" must be a matrix.")
  }
  if (ncol(data) != 2L) {
    stop("\"xdat\" must be a matrix with two columns.")
  }
  n <- nrow(data)

  confint <- match.arg(confint)
  if (length(qlim) != 2L) {
    stop("\"qlim\" must be a bivariate numeric vector.")
  }
  if (qlim[1] < 0 || qlim[2] >= 1) {
    stop("\"qlim\" must contain probabilities.")
  }
  ties.method <- match.arg(
    ties.method,
    choices = c("average", "first", "last", "random", "max", "min"),
    several.ok = FALSE
  )
  if (isTRUE(any(c(level < 0, level > 1)))) {
    stop("Invalid \"level\" argument: argument must be in the unit interval.")
  }
  if (length(level) == 2L) {
    qulevels <- sort(level)
  } else if (length(level) == 1L) {
    qulevels <- c((1 - level[1]) / 2, 1 - (1 - level[1]) / 2)
  } else {
    stop("Invalid \"level argument: must be a vector of length 1 or 2.")
  }
  datarank <- apply(X = data, MARGIN = 2, FUN = rank, ties.method = ties.method)
  rowmin <- apply(datarank, 1, min)
  # Sort data by decreasing minimum
  # (this will simplify computational burden)
  # as we only need to keep track
  # of indices of exceedances above the
  # smallest threshold
  od <- order(rowmin, decreasing = TRUE)
  pos_fn <- function(vector, scalar) {
    which.max(vector <= scalar) - 1L
  }
  datarank <- datarank[od, ]
  rowmin <- rowmin[od]
  # Pick quantiles
  eps <- .Machine$double.eps^0.5
  qlim2 <- c(
    min(apply(datarank, 1, max)) / (n + 1) + eps,
    rowmin[1] / (n + 1) - eps
  )
  if (is.null(u)) {
    if (!is.null(qlim)) {
      if (qlim[1] < qlim2[1]) {
        stop("lower quantile limit is too low")
      }
      if (qlim[2] > qlim2[2]) {
        stop("upper quantile limit is too high")
      }
      if (qlim[1] > qlim[2]) {
        stop("lower quantile limit is less than upper quantile limit")
      }
    } else {
      qlim <- qlim2
    }
    u <- seq(qlim[1], qlim[2], length = nq)
  } else {
    u <- sort(u)
  }
  nq <- length(u)
  if (min(u) < qlim2[1] || max(u) > qlim2[2]) {
    warning(
      "Upper quantile limit is too high or lower quantile limit is too low"
    )
  }

  # Compute empirical coefficient
  if (method == "emp") {
    empirical_coef <- function(
      u,
      rkdata,
      rowMin = NULL,
      ordered = TRUE,
      retList = FALSE
    ) {
      asym_coef_v <- vector(mode = "numeric", length = nq)
      nk_v <- psiplus_v <- vector(mode = "integer", length = nq)
      th <- u * (nrow(rkdata) + 1L)
      if (is.null(rowMin)) {
        rowMin <- apply(rkdata, 1, min)
        ordered <- FALSE
      }
      if (!ordered) {
        od <- order(rowMin, decreasing = TRUE)
        rowMin <- rowMin[od]
        rkdata <- rkdata[od, ]
      }
      for (i in seq_along(th)) {
        mind <- pos_fn(rowMin, th[i])
        # Perhaps superfluous, but can have ties
        psiplus <- sum(rkdata[seq_len(mind), 1] < rkdata[seq_len(mind), 2])
        psiminus <- sum(rkdata[seq_len(mind), 2] < rkdata[seq_len(mind), 1])
        nk_v[i] <- sum(psiplus + psiminus)
        psiplus_v[i] <- psiplus
        asym_coef_v[i] <- (psiplus - psiminus) / nk_v[i]
      }
      if (retList) {
        list(m = nk_v, psiplus = psiplus_v, coef = asym_coef_v)
      } else {
        return(asym_coef_v)
      }
    }
    xasym_res <- empirical_coef(
      u = u,
      rkdata = datarank,
      rowMin = rowmin,
      ordered = TRUE,
      retList = TRUE
    )
    est <- xasym_res$coef
    if (confint == "bootstrap") {
      boot_xasym <- matrix(0, ncol = nq, nrow = B + 1L)
      boot_xasym[1, ] <- est
      for (b in seq_len(B)) {
        boot_xasym[b + 1L, ] <-
          empirical_coef(
            u = u,
            rkdata = apply(
              data[sample.int(n, size = n, replace = TRUE), ],
              2,
              rank,
              ties.method = ties.method
            ),
            ordered = FALSE
          )
      }
      conf_int <- t(apply(
        boot_xasym,
        2,
        quantile,
        probs = qulevels,
        na.rm = TRUE
      ))
    } else if (confint == "wald") {
      variance <- with(xasym_res, 4 * psiplus / m * (1 - psiplus / m) / m)
      conf_int <- cbind(
        est + qnorm(qulevels[1]) * sqrt(variance),
        est + qnorm(qulevels[2]) * sqrt(variance)
      )
    }
  } else if (method == "elik") {
    # Function to compute empirical likelihood-based estimator
    xcoef_fun_emplik <- function(angles, weights) {
      (sum(weights * (angles < 0.5)) - 0.5) /
        (0.5 - sum(weights * (1 - 2 * angles) * (angles < 0.5)))
    }

    if (confint == "wald") {
      warning(
        "Method not implemented for\n the empirical likelihood-based estimator."
      )
      confint <- 'none'
    }

    # Execute calculations inside loop
    emplik_coef <- function(u, rkdata, rowMin = NULL, ordered = TRUE) {
      if (is.null(rowMin)) {
        rowMin <- apply(rkdata, 1, min)
        ordered <- FALSE
      }
      if (!ordered) {
        od <- order(rowMin, decreasing = TRUE)
        rowMin <- rowMin[od]
        rkdata <- rkdata[od, ]
      }
      # Use ordering to speed up calculations
      th_ind <- sapply(u * (nrow(rkdata) + 1), function(ui) {
        pos_fn(rowMin, ui)
      })
      # Only transform exceedances above smallest threshold
      frechet <- -1 / (log(rkdata[1:th_ind[1], ]) - log(nrow(rkdata) + 1))
      angs <- frechet[, 1, drop = FALSE] / rowSums(frechet)
      est <- vector("numeric", length = length(u))
      for (i in seq_along(u)) {
        emplik_sol <- try(
          suppressWarnings(
            angmeas(
              xdat = angs[1:th_ind[i], , drop = FALSE],
              thresh = 0,
              wgt = "Empirical",
              is.angle = TRUE
            )
          )
        )
        if (!inherits(emplik_sol, what = "try-error")) {
          est[i] <- xcoef_fun_emplik(
            angles = as.vector(emplik_sol$ang),
            weights = emplik_sol$wts
          )
        } else {
          est[i] <- NA
        }
      }
      return(est)
    }
    est <- emplik_coef(
      u = u,
      rkdata = datarank,
      rowMin = rowmin,
      ordered = TRUE
    )
    if (confint == "bootstrap") {
      boot_xasym <- matrix(0, ncol = nq, nrow = B + 1L)
      boot_xasym[1, ] <- est
      for (b in seq_len(B)) {
        boot_xasym[b + 1L, ] <-
          emplik_coef(
            u = u,
            rkdata = apply(
              data[sample.int(n, size = n, replace = TRUE), ],
              2,
              rank,
              ties.method = ties.method
            ),
            ordered = FALSE
          )
      }
      conf_int <- t(apply(
        boot_xasym,
        2,
        quantile,
        probs = qulevels,
        na.rm = TRUE
      ))
    }
  }
  if (confint == "none") {
    ret <- data.frame(
      threshold = u,
      coef = est
    )
  } else {
    ret <- data.frame(
      qlev = u,
      coef = est,
      lower = pmax(-1, conf_int[, 1]),
      upper = pmin(1, conf_int[, 2])
    )
  }
  attr(ret, "estimator") <- method
  attr(ret, "confint") <- confint
  attr(ret, "measure") <- "xasym"
  class(ret) <- c("mev_xdep", "data.frame")
  if (isTRUE(plot)) {
    plot(ret, ...)
  }
  return(invisible(ret))
}
