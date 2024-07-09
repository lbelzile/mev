#' Extreme U-statistic Pickands estimator
#'
#' Given a random sample of exceedances, the estimator
#' returns an estimator of the shape parameter or extreme
#' value index using a kernel of order 3, based on
#' \code{m} largest exceedances of \code{xdat}.
#'
#' The calculations are based on the recursions provided in Lemma 4.3 of Oorschot et al.
#' @param xdat vector of observations of length \eqn{n}
#' @param m number of largest order statistics \eqn{3 \leq m \leq n}. Choosing \eqn{m = n} amounts to using only the three largest observations in the sample.
#' @references Oorschot, J, J. Segers and C. Zhou (2023), Tail inference using extreme U-statistics,  Electron. J. Statist. 17(1): 1113-1159. \doi{10.1214/23-EJS2129}
#' @export
#' @examples
#' samp <- rgp(n = 1000, shape = 0.2)
#' PickandsXU(samp, m = 3)
PickandsXU <- function(xdat, m){
  m <- as.integer(m)
  stopifnot(m >= 3,
            is.vector(xdat),
            length(xdat) >= m)
  xdat <- na.omit(as.numeric(xdat))
  n <- length(xdat)
  xdat <- sort(xdat,
               decreasing = TRUE)
  # Initial recursion j=2
  shape <- (2*(n-1)/(m-2)-2)*log(xdat[1] - xdat[2])
  mcst <- 1
  for(j in seq.int(
    from = 3L,
    to = n - m + 3L,
    by = 1L)){
    mcst <- mcst * (n - j - m + 4L) / (n - j + 1L)
     shape <- shape + mcst * (2 * (n - j + 1) / (m - 2) - j) * sum(log(xdat[seq_len(j-1)] - xdat[j]))
  }
   shape <- shape * m * (m - 1) * (m - 2) / (n * (n - 1) * (n - m + 1))
  return(shape)

}
