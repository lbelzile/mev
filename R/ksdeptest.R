#' Test Jeffrey's divergence for equality of distribution
#'
#' Given two multivariate samples
#' \eqn{\boldsymbol{X}_{i,1}, \ldots, \boldsymbol{X}_{n_1,1}} and
#' \eqn{\boldsymbol{X}_{i,2}, \ldots, \boldsymbol{X}_{n_1,2}} from
#' an \eqn{d} dimensional distribution and a risk functional \eqn{r},
#' transform observations to unit Pareto margins above some high threshold
#' and partition the risk region into \eqn{K} distinct sub-regions.
#'
#' The test statistic is the Jeffreys' divergence
#'  (i.e., the symmetrized Kullback-Leibler divergence) between the
#'  two multinomial distributions, defined as
#'  \deqn{\sum_{j=1}^K (\widehat{p}_{1j}-\widehat{p}_{2j})(\log\widehat{p}_{1j}-\log\widehat{p}_{2j})
#' @param xdat a list of matrices with \code{d} columns containing the sample values
#' @param ties.method string indicating the type of method for \code{rank}; see \code{\link[base]{rank}} for a list of options. Default to \code{"random"}
#' @param margtrans marginal transformation; if \code{"none"}, data are assumed to be in uniform margins
#' @param risk string indicating the risk functional, or else a one-homogeneous function
#' @param qlev quantile level for threshold for each dataset, either a scalar or a vector of the same length as \code{xdat}
#' @param region function that returns a table
#' @param B integer; number of bootstrap samples
#' @examples
#' # Example of partition of the risk region
#' region <- function(x){
#'  table(factor(apply(x, 1, which.max), levels = 1:ncol(x)))
#' }
#)
test.jeffdep <- function(
  xdat,
  qlev,
  region,
  nexc = NULL,
  risk = c("sum", "max", "min", "l2"),
  margtrans = c("emp", "none"),
  # ties.method = "random",
  B = 1000
) {
  ties.method <- "random" # Hardcoded for now
  margtrans <- match.arg(margtrans)
  ties.method <- match.arg(
    arg = ties.method,
    choices = c("average", "first", "last", "random", "max", "min")
  )
  boot <- TRUE
  # Data are stored in a list
  stopifnot(
    "\"xdat\" must be a list of matrices" = is.list(xdat),
    "\"xdat\" must contain at least two matrices of samples." = length(xdat) >=
      2L
  )
  xdat <- lapply(xdat, as.matrix)
  D <- length(xdat)
  d <- ncol(xdat[[1]])
  nsub <- sapply(xdat, nrow)
  # Datasets can have different number of rows, but the same columns
  stopifnot(
    "All elements of \"xdat\" must have the same number of columns" = isTRUE(
      length(unique(sapply(xdat, ncol))) == 1
    )
  )
  if (!is.null(nexc)) {
    stopifnot(
      "\"qlev\" must be a vector of length 1 or \"length(xdat)\"." = length(
        qlev
      ) %in%
        c(1L, D),
      "\"qlev\" must contain quantile levels, with entries in [0,1)." = isTRUE(all(
        is.finite(qlev),
        qlev < 1,
        qlev >= 0
      ))
    )
    # Recycle qlev if different
    qlev <- rep(qlev, length = D)
  } else {
    nexc <- as.integer(nexc)
    stopifnot(
      "\"nexc\" must be a scalar or integer vector of the same length as \"xdat\"" = length(
        nexc
      ) %in%
        c(1, D),
      "\"nexc\" must be positive" = isTRUE(all(nexc > 0)),
      "\"nexc\" must be no smaller than the sample size." = isTRUE(all(
        nsub >= nexc
      ))
    )
    qlev <- 1 - nexc / nsub
  }

  if (is.function(risk)) {
    # Check homogeneity of function
    stopifnot(
      "\"risk\" must be an homogeneous function of order 1." = isTRUE(all.equal(
        risk(3 * 1:D),
        3 * risk(1:D)
      ))
    )
  } else {
    risk <- match.arg(risk)
    risk <- switch(
      risk,
      "sum" = sum,
      "max" = max,
      "min" = min,
      "l2" = function(x) {
        sqrt(sum(x^2))
      }
    )
  }

  # Helper functions
  unif2Par <- function(
    xdat,
    margtrans,
    ties.method = ties.method
  ) {
    if (margtrans == "none") {
      for (i in seq_along(xdat)) {
        stopifnot(
          "\"xdat\" must contain uniform values if \"margtrans\" is TRUE." = min(xdat[[
            i
          ]]) >=
            0,
          max(xdat[[i]]) < 1
        )
        # Transform to unit Pareto margins
        xdat[[i]] <- 1 / (1 - xdat[[i]])
      }
    } else {
      for (i in seq_along(xdat)) {
        xdat[[i]] <- 1 /
          (1 -
            apply(xdat[[i]], 2, rank, ties.method = ties.method) /
              (nrow(xdat[[i]]) + 1))
      }
    }
    return(xdat)
  }

  Par2Ang <- function(xdat, risk, qlev) {
    qlev <- rep(qlev, length(xdat))
    xang <- list()
    for (i in seq_along(xdat)) {
      rad <- apply(xdat[[i]], 1, risk)
      xang[[i]] <- (xdat[[i]] / rad)[rad > qlev[i], ]
    }
    return(xang)
  }
  Ang2Classif <- function(xang, region) {
    nexc <- sapply(xang, nrow)
    tabs <- lapply(xang, region)
    stopifnot(
      "\"region\" must classify all observations in a unique bin." = isTRUE(all(
        nexc == sapply(tabs, sum)
      )),
      "\"region\" must return tables of the same size for all subsamples" = length(unique(lapply(
        xang,
        length
      ))) ==
        1L
    )
    D <- length(xang)
    K <- length(tabs[[1]])
    probs <- matrix(nrow = D, ncol = K)
    for (i in seq_along(xang)) {
      probs[i, ] <- tabs[[i]] / sum(tabs[[i]])
    }
    return(probs)
  }
  jeffreysdiv <- function(probs) {
    D <- nrow(probs)
    K <- ncol(probs)
    logprobs <- log(probs)
    combo <- combn(D, 2)
    stat <- 0
    for (i in ncol(combo)) {
      stat <- stat +
        sum(
          (probs[combo[1, i], ] - probs[combo[2, i], ]) *
            (log(probs[combo[1, i], ]) - log(probs[combo[2, i], ]))
        )
    }
    return(stat)
  }
  compute_stat <- function(
    xdat,
    qlev,
    risk,
    region,
    margtrans = "emp",
    ties.method = "random"
  ) {
    jeffreysdiv(
      Ang2Classif(
        xang = Par2Ang(
          xdat = unif2Par(
            xdat,
            margtrans = margtrans,
            ties.method = ties.method
          ),
          risk = risk,
          qlev = qlev
        ),
        region = region
      )
    )
  }
  stat <- compute_stat(
    xdat = xdat,
    qlev = qlev,
    risk = risk,
    region = region,
    margtrans = margtrans,
    ties.method = "random"
  )
  if (isTRUE(boot)) {
    ncumsamp <- cumsum(nsub)
    bootstat <- numeric(B)
    for (b in seq_len(B)) {
      if (margtrans == "emp") {
        # repeat this in case there were ties in the ranks
        xdat2 <- unif2Par(
          xdat = xdat,
          margtrans = "emp",
          ties.method = "random"
        )
      } else {
        xdat2 <- xdat
      }
      pooled <- do.call(rbind, xdat2)
      bootsamp <- pooled[sample.int(nrow(pooled), size = nrow(pooled)), ]
      xboot <- list()
      # Split the data frame into list
      for (i in seq_along(xdat)) {
        xboot[[i]] <- bootsamp[(c(0, nsub)[i] + 1):ncumsamp[i], , drop = FALSE]
      }
      bootstat[b] <- compute_stat(
        xdat = xboot,
        qlev = qlev,
        risk = risk,
        region = region,
        margtrans = margtrans,
        ties.method = "random"
      )
    }
    out <- list(
      stat = stat,
      pval = mean(stat > bootstat)
    )
  }
}
