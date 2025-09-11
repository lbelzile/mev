#' Ramos and Ledford test of independence
#'
#' The Ramos and Ledford (2005) score test of independence is a modification of tests by Tawn (1988) and Ledford and Tawn (1996) for a logistic model parameter \eqn{\alpha=1}; the latter two have scores with zero expectation, but the variance of the score are infinite, which produces non-regularity and yield test, once suitably normalized, that converge slowly to their asymptotic null distribution. The test, designed for bivariate samples, transforms observations to have unit Frechet margins and considers a bivariate censored likelihood approach for the logistic distribution.
#'
#'
#' @param xdat a \code{n} by 2 matrix of observations
#' @param p probability level for the marginal threshold
#' @param test string; if \code{tawn}, only censor observations in the upper quadrant when both variables are large as in Tawn (1988), otherwise censor marginally for \code{ledford} as in Ledford and Tawn (1996).
#' @return a list with elements
#' \describe{
#' \item{\code{stat}}{value of the score test statistic}
#' \item{\code{pval}}{asymptotic p-value}
#' \item{\code{test}}{\code{test} argument}
#' }
#' @export
#' @examples
#' samp <- rmev(n = 1000, d = 2,
#'     param = 0.99, model = "log")
#' (test.scoreindep(samp, p = 0.9))
test.scoreindep <- function(
  xdat,
  p,
  test = c("ledford", "tawn")
) {
  # Sanity checks:
  test <- match.arg(test)
  xdat <- as.matrix(na.omit(xdat))
  stopifnot(
    ncol(xdat) == 2L,
    is.finite(p),
    isTRUE(all(p > 0, p < 1)),
    length(p) == 1
  )
  n <- nrow(xdat)
  # Transform data to 1/ unit Frechet scale
  X <- -log(rank(xdat[, 1]) / (n + 1))
  Y <- -log(rank(xdat[, 2]) / (n + 1))
  # Threshold on the 1 / Frechet scale
  u <- -log(p)
  if (p < 0.75) {
    warning(
      "Threshold too low: the approximate standard deviation of the score test may be inaccurate."
    )
  }
  # Break down lik by region
  region <- paste0("R", as.integer(X <= u), as.integer(Y <= u))
  nregion <- table(factor(region, levels = c("R00", "R01", "R10", "R11")))
  S11 <- 2 * u * log(2) * exp(-2 * u) / (2 * exp(-u) - exp(-2 * u) - 1)

  if (test == "ledford") {
    Ys <- Y[region == "R01"]
    Xs <- X[region == "R10"]
    S01 <- log(u) * u - (1 - Ys) * log(Ys) + (1 - u - Ys) * log(u + Ys)
    S10 <- log(u) * u - (1 - Xs) * log(Xs) + (1 - u - Xs) * log(u + Xs)
    S00 <- -2 * u * log(2)
    score <- nregion[1] * S00 + sum(S01) + sum(S10) + nregion[4] * S11
    # Standard deviation of score stat
    sigma = 1.107767 + 0.362784 * p - 0.084381 * p^2
  } else if (test == "tawn") {
    Xs <- X[region != "R11"]
    Ys <- Y[region != "R11"]
    # Score statistic for Tawn (1988)
    score <- sum(
      -(1 - Xs) *
        log(Xs) -
        (1 - Ys) * log(Ys) +
        (2 - Xs - Ys) * log(Xs + Ys) -
        1 / (Xs + Ys)
    ) +
      nregion[4] * S11
    # Standard deviation of score
    if (p > 1 - 1e-7) {
      warning(
        "The approximation to the standard deviation of the score statistic is not necessarily good."
      )
    }
    lp <- log(-log(1 - p))
    sigma <- 1.463193 + 0.312435 * lp + 0.132315 * lp^2 + 0.035713 * lp^3
  }
  # Test
  stat <- -as.numeric(score) / sqrt(n) / sigma
  pval <- 2 * pnorm(abs(stat), lower.tail = FALSE)
  res <- list(stat = stat, pval = pval, test = test)
  class(res) <- "mev_test_scoreindep"
  return(invisible(res))
}

#'@export
print.mev_test_scoreindep <- function(
  x,
  digits = min(3, getOption("digits") - 3),
  ...
) {
  cat("Score test of independence of Ramos and Ledford (2005)\n")
  cat(
    "Censoring scheme of",
    switch(x$test, ledford = "Ledford and Tawn (1996)", tawn = "Tawn (1988)"),
    "\n"
  )
  cat("Test statistic:", round(x$stat, digits), "\n")
  cat("P-value:", round(x$pval, digits), "\n")
}


#' Ramos and Ledford test of independence
#'
#' The Ramos and Ledford (2005) score test of independence is a modification of tests by Tawn (1988) and Ledford and Tawn (1996) for a logistic model parameter \eqn{\alpha=1}; the latter two have scores with zero expectation, but the variance of the score are infinite, which produces non-regularity and yield test, once suitably normalized, that converge slowly to their asymptotic null distribution. The test, designed for bivariate samples, transforms observations to have unit Frechet margins and considers a bivariate censored likelihood approach for the logistic distribution.
#'
#'
#' @param xdat a \code{n} by 2 matrix of observations
#' @param p probability level for the marginal threshold
#' @param test string; if \code{tawn}, only censor observations in the upper quadrant when both variables are large as in Tawn (1988), otherwise censor marginally for \code{ledford} as in Ledford and Tawn (1996).
#' @return a list with elements
#' \describe{
#' \item{\code{stat}}{value of the score test statistic}
#' \item{\code{pval}}{asymptotic p-value}
#' \item{\code{test}}{\code{test} argument}
#' }
#' @export
#' @keywords internal
scoreindep <- function(
  xdat,
  p,
  test = c("ledford", "tawn")
) {
  .Deprecated(new = "test.scoreindep", package = "mev")
  test.scoreindep(xdat = xdat, p = p, test = test)
}
