
#' Intensity function for the extremal Student model
#'
#' The intensity function includes the normalizing constants
#' @param tdat matrix of unit Pareto observations
#' @param df degrees of freedom, must be larger than 1
#' @param Sigma scale matrix
#' @param cholPrecis Cholesky root of the precision matrix \code{solve(Sigma)}. Default to \code{NULL}, meaning the latter is calculated within the function
#' @return log intensity contribution
#' @keywords internal
#' @export
intensXstud <- function(tdat, df, Sigma, cholPrecis = NULL) {
  D <- ncol(Sigma)
  N <- nrow(tdat)
  if (is.null(cholPrecis)) {
    A <- backsolve(chol(Sigma), diag(D))
  } else {
    A <- cholPrecis
  }
  ldet <- -2 * sum(log(diag(A)))

  N * ((1 - D) * log(df) - 0.5 * (D - 1) * log(pi) - 0.5 * ldet - lgamma((df + 1) / 2) + lgamma((df + D) / 2)) +
    (1 / df - 1) * sum(log(abs(tdat))) - 0.5 * (df + D) * sum(log(apply(sign(tdat) * abs(tdat)^(1 / df), 1, function(v) {
      tcrossprod(t(v) %*% A)
    })))
  # faster with sweep?
}


#' Intensity function for the Brown-Resnick model
#'
#' The intensity function includes the normalizing constants
#' @param tdat matrix of unit Pareto observations
#' @param Lambda conditionally negative definite parameter matrix of the Huesler--Reiss model
#' @param cholPrecis Cholesky root of the corresponding precision matrix \code{solve(Sigma)}. Default to \code{NULL}, meaning the latter is calculated within the function
#' @return log intensity contribution
#' @keywords internal
#' @export
intensBR <- function(tdat, Lambda, cholPrecis = NULL) {
  D <- ncol(Lambda)
  N <- nrow(tdat)
  if (is.null(cholPrecis)) {
    A <- backsolve(chol(Lambda2cov(Lambda = Lambda, co = 1, subA = -1)), diag(D - 1))
  } else {
    A <- cholPrecis
  }
  # Om <- t(apply(tdat, 1, function(x){log(x[-1])-log(x[1])+ 2*Lambda[1,-1]}))
  ldet <- -2 * sum(log(diag(A)))
  # mu = -diag(Sigma)/2 == -semivario(di)[-1,1] == -2*Lambda[1,-1]
  N * (-0.5 * (D - 1) * log(pi) - 0.5 * ldet) - sum(log(tdat[, -1])) - 2 * sum(log(tdat[, 1])) - 0.5 * sum(apply(tdat, 1, function(x) {
    tcrossprod(t(log(x[-1]) - log(x[1]) + Lambda[1, -1]) %*% A)
  }))
}


#' Jacobian of the transformation from generalized Pareto to unit Pareto distribution
#'
#' If \code{dat} is a vector, the arguments \code{loc}, \code{scale} and \code{shape} should be numericals of length 1.
#' @param dat vector or matrix of data
#' @param loc vector of location parameters
#' @param scale vector of scale parameters, strictly positive
#' @param shape shape parameter
#' @param lambdau vector of probability of marginal threshold exceedance
#' @param censored a matrix of logical indicating whether the observations are censored
#' @return log-likelihood contribution for the Jacobian
#' @keywords internal
#' @export
jac <- function(dat, loc = 0, scale, shape, lambdau = 1, censored) {
  if (is.vector(dat)) {
    stopifnot(all(is.vector(censored), length(censored) == length(dat), length(scale) == 1, length(shape) == 1))
    return(-(length(dat) - sum(censored)) * log(scale) + (1 / shape - 1) * sum(log(1 + shape / scale * pmax(0, dat[!censored] - loc))))
  } else {
    dat <- as.matrix(dat)
    if (!is.matrix(dat)) {
      stop("\"dat\" must be a matrix")
    }
    N <- nrow(dat)
    if (!all(dim(dat) == dim(censored))) {
      stop("Invalid input in Jacobian contribution.")
    }
    loc <- rep(loc, length.out = ncol(dat))
    scale <- rep(scale, length.out = ncol(dat))
    shape <- rep(shape, length.out = ncol(dat))
    lambdau <- rep(lambdau, length.out = ncol(dat))
    ll <- 0
    for (j in 1:ncol(dat)) {
      ll <- ll - (N - sum(censored[, j])) * (log(scale[j]) + log(lambdau[j])) + (1 / shape[j] - 1) * sum(log(1 + shape[j] / scale[j] * pmax(0, dat[!censored[, j], j] - loc[j])))
    }
    return(ll)
  }
}


#' Transformation from the generalized Pareto to unit Pareto
#'
#' @inheritParams jac
#' @return a vector or matrix of the same dimension as \code{dat} with unit Pareto observations
#' @keywords internal
#' @export
gpdtopar <- function(dat, loc = 0, scale, shape, lambdau = 1) {
  if (is.vector(dat)) {
    return((1 + shape / scale * pmax(dat - loc, 0))^(1 / shape) / lambdau)
  } else {
    loc <- rep(loc, length.out = ncol(dat))
    scale <- rep(scale, length.out = ncol(dat))
    shape <- rep(shape, length.out = ncol(dat))
    sapply(1:ncol(dat), function(j) {
      if (!isTRUE(all.equal(shape[j], 0))) {
        pmax(0, (1 + shape[j] / scale[j] * (dat[, j] - loc[j])))^(1 / shape[j]) / lambdau[j]
      } else {
        exp((dat[, j] - loc[j]) / scale[j]) / lambdau[j]
      }
    })
  }
}


#' Likelihood for multivariate generalized Pareto distribution
#'
#' Likelihood for the Brown--Resnick, extremal Student or logistic vectors over region determined by
#' \deqn{\{y \in F: \max_{j=1}^D \sigma_j \frac{y^\xi_j-1}{\xi_j}+\mu_j  > u\};}
#' where \eqn{\mu} is \code{loc}, \eqn{\sigma} is \code{scale} and \eqn{\xi} is \code{shape}.
#' @param dat matrix of observations
#' @param thresh functional threshold for the maximum
#' @param loc vector of location parameter for the marginal generalized Pareto distribution
#' @param scale vector of scale parameter for the marginal generalized Pareto distribution
#' @param shape vector of shape parameter for the marginal generalized Pareto distribution
#' @param lambdau vector of marginal rate of marginal threshold exceedance.
#' @param likt string indicating the type of likelihood, with an additional contribution for the non-exceeding components: one of  \code{"mgp"}, \code{"binom"} and \code{"pois"}.
#' @param ... additional arguments (see Details)
#' @param par list of parameters: \code{alpha} for the logistic model, \code{Lambda} for the Brown--Resnick model or else \code{Sigma} and \code{df} for the extremal Student.
#' @param model string indicating the model family, one of \code{"log"}, \code{"br"} or \code{"xstud"}
#' @note The location and scale parameters are not identifiable unless one of them is fixed.
#' @details
#' Optional arguments can be passed to the function via \code{...}
#' \itemize{
#' \item \code{cl} cluster instance  created by \code{makeCluster} (default to \code{NULL})
#' \item \code{ncors} number of cores for parallel computing of the likelihood
#' \item \code{mmax} maximum per column
#' \item \code{B1} number of replicates for quasi Monte Carlo integral for the exponent measure
#' \item \code{genvec1} generating vector for the quasi Monte Carlo routine (exponent measure), associated with \code{B1}
#' }
#' @return the value of the log-likelihood with \code{attributes} \code{expme}, giving the exponent measure
#' @export
likmgp <- function(dat, thresh, loc, scale, shape, par, model = c("br", "xstud", "log"),
                   likt = c("mgp", "pois", "binom"), lambdau = 1, ...) {
  # Rename arguments
  tdat <- dat
  N <- nrow(dat)
  D <- ncol(dat)
  A <- rep(scale, length.out = D)
  B <- rep(loc, length.out = D)
  xi <- rep(shape, length.out = D)
  lambdau <- rep(lambdau, length.out = D)
  model <- match.arg(model)
  likt <- match.arg(likt)
  stopifnot(max(lambdau) <= 1, min(lambdau) > 0)
  if (model == "br") {
    Lambda <- par$Lambda
    if (is.null(Lambda)) {
      stop("Invalid \"par\" for \"br\" model.")
    }
  } else if (model == "xstud") {
    Sigma <- par$Sigma
    df <- par$df
    if (any(is.null(df), is.null(Sigma))) {
      stop("Invalid \"par\" for \"xstud\" model.")
    }
  } else if (model == "log") {
    alpha <- par$alpha
    if (is.null(alpha)) {
      stop("Invalid \"par\"")
    }
    alpha <- alpha[1]
    if (alpha > 1) {
      alpha <- 1 / alpha
    }
    if (alpha < 0) {
      stop("Invalid \"par\" for \"log\" model.")
    }
  }
  stopifnot(is.matrix(tdat), ncol(tdat) > 1)
  ellips <- list(...)
  if (likt %in% c("pois", "binom")) {
    ntot <- ellips$ntot
    if (is.null(ntot)) {
      stop("Poisson/binomial likelihood requires the total number of observations above the threshold")
    }
  }


  # Schur complement of submatrix \code{Sigma} excluding indices \code{ind}
  schurcompC <- function(Sigma, ind) {
    stopifnot(c(length(ind) > 0, ncol(Sigma) - length(ind) > 0))
    Sigma[-ind, -ind, drop = FALSE] - Sigma[-ind, ind, drop = FALSE] %*% solve(Sigma[ind, ind, drop = FALSE]) %*% Sigma[ind, -ind, drop = FALSE]
  }

  if (is.null(mmax)) {
    mmax <- apply(tdat, 2, max, na.rm = TRUE)
  }
  if (is.null(mmin)) {
    mmin <- apply(tdat, 2, min, na.rm = TRUE)
  }
  # Check for marginal constraints
  if (!isTRUE(all(ifelse(xi < 0, A + (xi * mmax - B) > 0, A + (xi * mmin - B) > 0)))) {
    return(-1e10)
  }
  # Compute marginal transformation and Jacobian
  tu <- rep(0, D)
  jac <- -N * sum((log(A) + log(lambdau)))
  for (j in 1:D) {
    if (abs(xi[j]) > 1e-5) {
      tdat[, j] <- pmax(0, (1 + xi[j] * (dat[, j] - B[j]) / A[j]))
      jac <- jac + (1 / xi[j] - 1) * sum(log(tdat[, j]))
      tdat[, j] <- tdat[, j]^(1 / xi[j]) / lambdau[j]
      # Map thresholds
      tu[j] <- (1 + xi[j] * (thresh - B[j]) / A[j])^(1 / xi[j]) / lambdau[j]
    } else { # xi is zero
      tdat[, j] <- (dat[, j] - B[j]) / A[j] # this is a transformation onto log scale
      # Jacobian of marginal transformation - for both models
      jac <- jac + sum(tdat[, j])
      tdat[, j] <- exp(tdat[, j] - log(lambdau[j]))
      # Map thresholds
      tu[j] <- exp((thresh - B[j]) / A[j]) / lambdau[j]
    }
  }
  if(model %in% c("br","xstud")){
  if (!requireNamespace("mvPot", quietly = TRUE)) {
    stop(
      "Package \"mvPot\" must be installed to use this function.",
      call. = FALSE
    )
  }
  }

  # Copy from ellipsis
  mmin <- ellips$mmin
  mmax <- ellips$mmax
  genvec1 <- ellips$genvec1
  B1 <- ifelse(is.null(ellips$B1), 1009L, ellips$B1)
  antithetic <- ifelse(is.null(ellips$antithetic), FALSE, ellips$antithetic)
  if (is.null(genvec1)) {
    genvec1 <- mvPot::genVecQMC(B1, ncol(dat) - 1L)$genVec
  }
  M1 <- ifelse(is.null(ellips$M1), 1L, ellips$M1)
  ncores <- ifelse(is.null(ellips$ncores), 1L, ellips$ncores)

  if (model == "br") {
    intens <- intensBR(tdat = tdat, Lambda = Lambda)
    exponentMeasure <- sum(.weightsBR(z = tu, Lambda = Lambda, prime = B1, method = "mvPot", genvec = genvec1, nrep = 1) / tu)
  } else if (model == "xstud") {
    intens <- intensXstud(tdat = tdat, df = df, Sigma = Sigma)
    exponentMeasure <- sum(.weightsXstud(z = tu, Sigma = Sigma, df = df, method = "mvPot", prime = B1, genvec = genvec1, nrep = 1) / tu)
  } else if (model == "log") {
    lVfunlog <- function(x, alpha) {
      if (is.null(dim(x))) {
        alpha * log(sum(x^(-1 / alpha)))
      } else {
        alpha * log(rowSums(x^(-1 / alpha)))
      }
    }
    lVu <- lVfunlog(x = tu, alpha = alpha)
    lVx <- sum(lVfunlog(x = tdat, alpha = alpha))
    lfalfacto1 <- function(x, s) {
      sum(log(abs(seq(x, x - s + 1, by = -1))))
    }
    ldVfunlog <- function(x, alpha, lV) {
      falf <- lfalfacto1(alpha, ncol(x))
      -length(dat) * log(alpha) + nrow(x) * falf -
        (1 / alpha + 1) * sum(log(x)) + (alpha - ncol(x)) * lV / alpha
    }
    intens <- ldVfunlog(x = tdat, alpha = alpha, lV = lVx)
    exponentMeasure <- exp(lVu)
  }
  res <- intens + jac + switch(likt,
    mgp = -N * log(exponentMeasure),
    pois = ntot * exponentMeasure + N * log(ntot) - lgamma(N + 1),
    binom = -(ntot - N) * log(1 - exponentMeasure) + lchoose(ntot, N)
  )
  attributes(res) <- list("expme" = exponentMeasure)
  res
}

#' Censored likelihood for multivariate generalized Pareto distributions
#'
#' Censored likelihood for the logistic distribution and the Brown--Resnick and extremal Student processes.
#'
#' @inheritParams likmgp
#' @param mthresh vector of individuals thresholds under which observations are censored
#' @param ... additional arguments (see Details)
#' @note The location and scale parameters are not identifiable unless one of them is fixed.
#' @details
#' Optional arguments can be passed to the function via \code{...}
#' \itemize{
#' \item \code{censored} matrix of booleans and \code{NA} indicating whether observations \code{dat} fall below the mthreshold \code{mthresh}
#' \item \code{cl} cluster instance  created by \code{makeCluster} (default to \code{NULL})
#' \item \code{ncors} number of cores for parallel computing of the likelihood
#' \item \code{numAbovePerRow} number of observations above mthreshold (non-missing) per row
#' \item \code{numAbovePerCol} number of observations above mthreshold (non-missing) per column
#' \item \code{mmax} maximum per column
#' \item \code{B1} number of replicates for quasi Monte Carlo integral for the exponent measure
#' \item \code{B2} number of replicates for quasi Monte Carlo integral for the censored intensity contribution
#' \item \code{genvec1} generating vector for the quasi Monte Carlo routine (exponent measure), associated with \code{B1}
#' \item \code{genvec2} generating vector for the quasi Monte Carlo routine (individual obs contrib), associated with \code{B2}
#' }
#' @return the value of the log-likelihood with \code{attributes} \code{expme}, giving the exponent measure
#' @export
clikmgp <- function(dat, thresh, mthresh = thresh, loc, scale, shape, par, model = c("br", "xstud", "log"),
                    likt = c("mgp", "pois", "binom"), lambdau = 1, ...) {

  # Rename arguments
  tdat <- dat
  N <- nrow(dat)
  D <- ncol(dat)
  A <- rep(scale, length.out = D)
  B <- rep(loc, length.out = D)
  xi <- rep(shape, length.out = D)
  lambdau <- rep(lambdau, length.out = D)
  model <- match.arg(model)
  likt <- match.arg(likt)
  stopifnot(max(lambdau) <= 1, min(lambdau) > 0)
  if (model == "br") {
    Lambda <- par$Lambda
    if (is.null(Lambda)) {
      stop("Invalid \"par\" for \"br\" model.")
    }
  } else if (model == "xstud") {
    Sigma <- par$Sigma
    df <- par$df
    if (any(is.null(df), is.null(Sigma))) {
      stop("Invalid \"par\" for \"xstud\" model.")
    }
  } else if (model == "log") {
    alpha <- par$alpha
    if (is.null(alpha)) {
      stop("Invalid \"par\"")
    }
    alpha <- alpha[1]
    if (alpha > 1) {
      alpha <- 1 / alpha
    }
    if (alpha < 0) {
      stop("Invalid \"par\" for \"log\" model.")
    }
  }
  stopifnot(is.matrix(tdat), ncol(tdat) > 1)
  if (length(mthresh) < ncol(tdat)) {
    mthresh <- rep(mthresh, length.out = ncol(tdat))
  }
  ellips <- list(...)
  if (likt == "pois") {
    ntot <- ellips$ntot
    if (is.null(ntot)) {
      stop("Poisson likelihood requires the total number of observations above the mthreshold")
    }
  }

  # Copy from ellipsis
  numAbovePerRow <- ellips$numAbovePerRow
  numAbovePerCol <- ellips$numAbovePerCol
  mmax <- ellips$mmax
  genvec1 <- ellips$genvec1
  genvec2 <- ellips$genvec2
  B1 <- ifelse(is.null(ellips$B1), 1009L, ellips$B1)
  B2 <- ifelse(is.null(ellips$B2), 499L, ellips$B2)
  antithetic <- ifelse(is.null(ellips$antithetic), FALSE, ellips$antithetic)
  if (is.null(genvec1)) {
    genvec1 <- mvPot::genVecQMC(B1, ncol(dat) - 1L)$genVec
  }
  if (is.null(genvec2)) {
    genvec2 <- mvPot::genVecQMC(B2, ncol(dat) - 1L)$genVec
  }
  M1 <- ifelse(is.null(ellips$M1), 1L, ellips$M1)
  M2 <- ifelse(is.null(ellips$M2), 1L, ellips$M2)
  ncores <- ifelse(is.null(ellips$ncores), 1L, ellips$ncores)
  censored <- ellips$censored
  if (is.null(censored)) {
    censored <- matrix(FALSE, nrow = nrow(tdat), ncol = ncol(tdat))
    for (j in 1:ncol(tdat)) {
      censored[, j] <- tdat[, j] < mthresh[j]
    }
  }

  # Schur complement of submatrix \code{Sigma} excluding indices \code{ind}
  schurcompC <- function(Sigma, ind) {
    stopifnot(c(length(ind) > 0, ncol(Sigma) - length(ind) > 0))
    Sigma[-ind, -ind, drop = FALSE] - Sigma[-ind, ind, drop = FALSE] %*% solve(Sigma[ind, ind, drop = FALSE]) %*% Sigma[ind, -ind, drop = FALSE]
  }

  if (is.null(numAbovePerRow)) {
    numAbovePerRow <- D - rowSums(censored)
  }
  if (is.null(numAbovePerCol)) {
    numAbovePerCol <- N - colSums(censored)
  }
  if (is.null(mmax)) {
    mmax <- apply(tdat, 2, max, na.rm = TRUE)
  }
  stopifnot(dim(tdat) == dim(censored), length(numAbovePerRow) == N)
  # Check boundary constraints for the support of the GP distribution
  # if(isTRUE(any(ifelse(shape >= 0, margmthresh < (B - shape / A), mmax > (B - shape/A))))){
  #  return(-Inf)
  # }
  # Compute marginal transformation and Jacobian
  tu <- rep(0, D)
  yth <- rep(0, D)
  jac <- -sum(numAbovePerCol * (log(A) + log(lambdau)))
  for (j in 1:D) {
    if (abs(xi[j]) > 1e-5) {
      if (model == "br") {
        tdat[, j] <- log(pmax(0, 1 + xi[j] * (dat[, j] - B[j]) / A[j]))
        jac <- jac + (1 / xi[j] - 1) * sum(tdat[!censored[, j], j])
        tdat[, j] <- (1 / xi[j]) * tdat[, j] - log(lambdau[j])
      } else {
        tdat[, j] <- pmax(0, (1 + xi[j] * (dat[, j] - B[j]) / A[j]))
        jac <- jac + (1 / xi[j] - 1) * sum(log(tdat[!censored[, j], j]))
        tdat[, j] <- tdat[, j]^(1 / xi[j]) / lambdau[j]
      }
      # Map mthresholds
      yth[j] <- (1 + xi[j] * (mthresh[j] - B[j]) / A[j])^(1 / xi[j]) / lambdau[j]
      tu[j] <- (1 + xi[j] * (thresh - B[j]) / A[j])^(1 / xi[j]) / lambdau[j]
    } else { # xi is zero
      tdat[, j] <- (dat[, j] - B[j]) / A[j] # this is a transformation onto log scale
      # Jacobian of marginal transformation - for both models
      jac <- jac + sum(tdat[!censored[, j], j])
      tdat[, j] <- tdat[, j] - log(lambdau[j])
      if (model != "br") {
        tdat[, j] <- exp(tdat[, j])
      }
      # Map mthresholds
      yth[j] <- exp((mthresh[j] - B[j]) / A[j]) / lambdau[j]
      tu[j] <- exp((thresh - B[j]) / A[j]) / lambdau[j]
    }
  }
  # test <- mvPot::censoredLikelihoodBR(obs = split(exp(tdat), row(exp(tdat))), loc = sites, vario = varioloc, thresh = tu, p = B1, vec = genvec1)
  # return(-test + jac)
  #
  # par(mfrow = c(1,1))
  #  wexc <- which(apply(t(t(rain[,stid]) - us), 1, max) > thresh)
  #  plot(log(1/(1-(rank(as.vector(rain[,stid[j]]))/(ellips$ntot+1))))[wexc], tdat[,j], xlab = "theoretical",ylab = "empirical", main = shape[1]);
  #  abline(0,1)
  #  # b <- spunif(x = as.vector(rain[,stid[j]]),mthresh = us[j], scale = scale[j], shape = shape[j])
  #  # plot(log(1/(1-b[wexc])), tdat[,j]); abline(0,1)
  # abline(h = log(yth[j]))
  # Use empirical transformation
  # jac <- 0
  # for(j in 1:D){
  #   tdat[,j] <- log(1/(1-(rank(as.vector(rain[,stid[j]]))/(ellips$ntot+1))))[wexc]
  # }
  if (model != "log") {
    likelihood_xstud <- function(i) {
      if (i < N + 1) {
        k <- numAbovePerRow[i]
        ab <- which(!censored[i, ]) # uncensored
        Zin <- tdat[i, ab]^(1 / df)
        if (k > 1) {
          cholS <- chol(Sigma[ab, ab])
          logdetS <- 2 * sum(log(diag(cholS)))
          Siginv <- chol2inv(cholS)
          schurcomp <- Sigma[-ab, -ab, drop = FALSE] - Sigma[-ab, ab, drop = FALSE] %*% Siginv %*% Sigma[ab, -ab, drop = FALSE]
          vecS <- Siginv %*% Zin
        } else {
          logdetS <- log(Sigma[ab, ab])
          vecS <- Zin / Sigma[ab, ab]
          schurcomp <- Sigma[-ab, -ab, drop = FALSE] - Sigma[-ab, ab, drop = FALSE] %*% Sigma[ab, -ab, drop = FALSE]
          # faster than crossprod
        }
        kst <- c(Zin %*% vecS)
        muC <- c(Sigma[-ab, ab, drop = FALSE] %*% vecS)
        if (D - numAbovePerRow[i] > 0) {
          contribBelow <- mvPot::mvTProbQuasiMonteCarlo(
            p = B2, upperBound = yth[-ab] - muC,
            cov = kst / (length(ab) + df) * schurcomp, nu = df + length(ab),
            genVec = genvec2[1:length(muC)], nrep = M2, antithetic = antithetic
          )[1]
        } else {
          contribBelow <- 1
        } # fixed 19-04-2019 to account for the fact that tdat^(1/df) already
        contribAbove <- -(k + df) / 2 * log(kst) + (1 - df) * sum(log(Zin)) + lgamma((df + k) / 2) -
          lgamma((df + 1) / 2) - 0.5 * logdetS - (k - 1) * log(df) - (k - 1) / 2 * log(pi)
        return(log(contribBelow) + contribAbove)
      } else if (i > N) {
        # Compute exponent measure
        j <- i - N
        mvPot::mvTProbQuasiMonteCarlo(
          p = B1, upperBound = (exp((log(tu[-j]) - log(tu[j])) / df) - Sigma[-j, j]),
          cov = (Sigma[-j, -j] - Sigma[-j, j, drop = FALSE] %*% Sigma[j, -j, drop = FALSE]) / (df + 1),
          nu = df + 1, genVec = genvec1, nrep = M1, antithetic = antithetic
        )[1]
      }
    }
    likelihood_br <- function(i) {
      # Note: the vector tdat is already on the log-scale
      if (i < N + 1) {
        be <- which(censored[i, ])
        ab <- which(!censored[i, ])
        if (length(ab) == 0) {
          warning("Invalid input - non exceedance")
          return(0)
        }
        # remove first non-censored, shift indices by 1
        be2 <- be - I(be > ab[1])
        ab2 <- ab[-1] - I(ab[-1] > ab[1])
        SigmaD <- outer(2 * Lambda[ab[1], -ab[1]], 2 * Lambda[ab[1], -ab[1]], "+") - 2 * Lambda[-ab[1], -ab[1]]
        muD <- -2 * Lambda[ab[1], -ab[1]] + tdat[i, ab[1]]
        if (numAbovePerRow[i] == 1) { # all but one observation fall below mthreshold, so censored
          contribBelow <- mvPot::mvtNormQuasiMonteCarlo(p = B2, upperBound = log(yth[-ab]) - muD, cov = SigmaD, genVec = genvec2[1:length(muD)], nrep = M2, antithetic = antithetic)[1]
          logcontribAbove <- -2 * tdat[i, ab[1]]
        } else if (numAbovePerRow[i] == D) { # all above!
          contribBelow <- 1
          logcontribAbove <- -tdat[i, ab[1]] - sum(tdat[i, ]) + .dmvnorm_arma(x = tdat[i, ab[-1], drop = FALSE], mean = muD, sigma = SigmaD, log = TRUE)
        } else {
          # return(- sum(tdat[i,ab]) - tdat[i,ab[1]] + dcondmvtnorm(x = tdat[i,-ab[1]], ind = be2, ubound = log(yth[-ab]), mu = muD, Sigma = SigmaD, model = "norm", n = 500, log = TRUE))
          logcontribAbove <- -sum(tdat[i, ab]) - tdat[i, ab[1]] +
            .dmvnorm_arma(x = tdat[i, ab[-1], drop = FALSE], mean = as.vector(muD[ab2]), log = TRUE, sigma = as.matrix(SigmaD[ab2, ab2]))
          muC <- c(muD[be2] + SigmaD[be2, ab2] %*% solve(SigmaD[ab2, ab2]) %*% (tdat[i, ab[-1]] - muD[ab2]))
          contribBelow <- mvPot::mvtNormQuasiMonteCarlo(
            p = B2, upperBound = log(yth[-ab]) - muC, cov = schurcompC(SigmaD, ab2),
            genVec = genvec2[1:length(muC)], nrep = M2, antithetic = antithetic
          )[1]
        }
        return(log(contribBelow) + logcontribAbove)
      } else if (i > N) {
        # Compute exponent measure
        j <- i - N
        mvPot::mvtNormQuasiMonteCarlo(p = B1, upperBound = (2 * Lambda[-j, j] + log(tu[-j]) - log(tu[j])), cov = 2 * (outer(Lambda[-j, j], Lambda[j, -j], FUN = "+") - Lambda[-j, -j]), genVec = genvec1[1:(D - 1)], nrep = M1, antithetic = antithetic)[1]
      }
    }
    if (ncores > 1 && !requireNamespace("parallel", quietly = TRUE)) {
      pro <- parallel::mclapply(X = 1:(D + N), FUN = switch(model, br = likelihood_br, xstud = likelihood_xstud))
    } else {
      pro <- lapply(X = 1:(D + N), FUN = switch(model, br = likelihood_br, xstud = likelihood_xstud))
    }
    exponentMeasure <- sum(unlist(pro)[(1 + N):(D + N)] / tu)
    intens <- sum(unlist(pro)[1:N])
  } else { # Model is logistic
    lVfunlog <- function(x, alpha) {
      if (is.null(dim(x))) {
        alpha * log(sum(x^(-1 / alpha)))
      } else {
        alpha * log(rowSums(x^(-1 / alpha)))
      }
    }
    cdat <- t(apply(tdat, 1, function(x) {
      pmax(yth, x)
    }))
    lVu <- lVfunlog(x = tu, alpha = alpha)
    lVx <- lVfunlog(x = cdat, alpha = alpha)
    lfalfacto1 <- function(x, s) {
      sum(log(abs(seq(x, x - s + 1, by = -1))))
    }
    ldVfunlog <- function(x, censored, alpha, numAbovePerRow, lV) {
      falf <- c(log(alpha), sapply(2:D, function(s) {
        lfalfacto1(alpha, s)
      }))
      sum(-numAbovePerRow * log(alpha) + falf[numAbovePerRow]) -
        (1 / alpha + 1) * sum(log(x[!censored])) + sum((alpha - numAbovePerRow) * lV) / alpha
    }
    intens <- ldVfunlog(x = cdat, censored = censored, alpha = alpha, numAbovePerRow = numAbovePerRow, lV = lVx)
    exponentMeasure <- exp(lVu)
  }
  res <- jac + intens + switch(likt,
    mgp = -N * log(exponentMeasure),
    pois = -ntot * exponentMeasure + N * log(ntot) - lgamma(N + 1),
    binom = -(ntot - N) * log(1 - exponentMeasure) + lchoose(ntot, N)
  )
  attributes(res) <- list("expme" = exponentMeasure)
  return(res)
}

#' Exponent measure for multivariate generalized Pareto distributions
#'
#' Integrated intensity over the region defined by \eqn{[0, z]^c} for logistic, Huesler-Reiss, Brown-Resnick and extremal Student processes.
#'
#' @note The list \code{par} must contain different arguments depending on the model. For the Brown--Resnick model, the user must supply the conditionally negative definite matrix \code{Lambda} following the parametrization in Engelke \emph{et al.} (2015) or the covariance matrix \code{Sigma}, following Wadsworth and Tawn (2014). For the Husler--Reiss model, the user provides the mean and covariance matrix, \code{m} and \code{Sigma}. For the extremal student, the covariance matrix \code{Sigma} and the degrees of freedom \code{df}. For the logistic model, the strictly positive dependence parameter \code{alpha}.
#' @param z vector at which to estimate exponent measure
#' @param par list of parameters
#' @param model string indicating the model family
#' @param method string indicating the package from which to extract the numerical integration routine
#' @return numeric giving the measure of the complement of \eqn{[0,z]}.
#' @export
#' @examples
#' \dontrun{
#' # Extremal Student
#' Sigma <- stats::rWishart(n = 1, df = 20, Sigma = diag(10))[, , 1]
#' expme(z = rep(1, ncol(Sigma)), par = list(Sigma = cov2cor(Sigma), df = 3), model = "xstud")
#' # Brown-Resnick model
#' D <- 5L
#' loc <- cbind(runif(D), runif(D))
#' di <- as.matrix(dist(rbind(c(0, ncol(loc)), loc)))
#' semivario <- function(d, alpha = 1.5, lambda = 1) {
#'   (d / lambda)^alpha
#' }
#' Vmat <- semivario(di)
#' Lambda <- Vmat[-1, -1] / 2
#' expme(z = rep(1, ncol(Lambda)), par = list(Lambda = Lambda), model = "br", method = "mvPot")
#' Sigma <- outer(Vmat[-1, 1], Vmat[1, -1], "+") - Vmat[-1, -1]
#' expme(z = rep(1, ncol(Lambda)), par = list(Lambda = Lambda), model = "br", method = "mvPot")
#' }
expme <- function(z, par, model = c("log", "hr", "br", "xstud"),
                  method = c("TruncatedNormal", "mvtnorm", "mvPot")) {
  model <- match.arg(model[1], choices = c("log", "hr", "br", "xstud"))
  if (model != "log") {
    method <- match.arg(method[1], choices = c("mvtnorm", "mvPot", "TruncatedNormal"))
    if (method == "mvtnorm") {
      if (!requireNamespace("mvtnorm", quietly = TRUE)) {
        warning("\"mvtnorm\" package is not installed.")
        method <- "TruncatedNormal"
      }
    } else if (method == "mvPot") {
      if (!requireNamespace("mvPot", quietly = TRUE)) {
        warning("\"mvPot\" package is not installed.")
        method <- "TruncatedNormal"
      }
    } else if(method == "TruncatedNormal"){
      if (!requireNamespace("TruncatedNormal", quietly = TRUE)) {
        stop(
          "Package \"TruncatedNormal\" must be installed to use this function.",
          call. = FALSE
        )
      }
      if(!utils::packageVersion("TruncatedNormal") > "1.1"){
        stop(
          "Please update package \"TruncatedNormal\" to a more recent version.",
          call. = FALSE
        )
      }
    }
  }
  if (model == "log") {
    alpha <- par$alpha
    if (any(c(is.null(alpha), alpha < 0, length(alpha) > 1))) {
      stop("Invalid or missing arguments for the logistic model")
    }
    if (alpha > 1) {
      alpha <- 1 / alpha
    }
    lVfunlog <- function(x, alpha) {
      if (is.null(dim(x))) {
        alpha * log(sum(x^(-1 / alpha)))
      } else {
        alpha * log(rowSums(x^(-1 / alpha)))
      }
    }
    return(lVfunlog(z, alpha))
  } else if (model == "hr") {
    m <- par$m
    Sigma <- par$Sigma
    if (any(c(is.null(m), is.null(Sigma), ncol(Sigma) != nrow(Sigma)))) {
      stop("Invalid or missing arguments for the Huesler-Reiss model")
    }
    D <- ncol(Sigma)
    Sigmainv <- solve(Sigma)
    q <- Sigmainv %*% rep(1, D)
    Q <- (Sigmainv - q %*% t(q) / sum(q))
    l <- c(Sigmainv %*% (((m %*% q - 1) / sum(q))[1] * rep(1, D) - m))
    return(expmeHR(z = z, L = l, Q = Q, method = method))
  } else if (model == "xstud") {
    Sigma <- par$Sigma
    df <- par$df
    if (any(c(is.null(Sigma), is.null(df), ncol(Sigma) != nrow(Sigma)))) {
      stop("Invalid or missing arguments for the extremal Student model")
    }
    return(expmeXS(z = z, Sigma = Sigma, df = df, method = method))
  } else if (model == "br") {
    Lambda <- par$Lambda
    Sigma <- par$Sigma
    if (!is.null(Lambda)) {
      return(expmeBR(z = z, Lambda = Lambda, method = method))
    } else if (!is.null(Sigma)) {
      return(expmeBR_WT(z = z, Sigma = Sigma, method = method))
    } else {
      stop("Invalid or missing arguments for the Brown-Resnick model")
    }
  }
}

expmeBR <- function(z, Lambda, method = c("mvtnorm", "mvPot", "TruncatedNormal")) {
  weights <- .weightsBR(z = z, Lambda = Lambda, method = method)
  sum(weights / z)
}

expmeHR <- function(z, Q, L, method = c("mvtnorm", "mvPot", "TruncatedNormal")) {
  weights <- .weightsHR(z = z, Q = Q, L = L, method = method)
  sum(weights / z)
}

expmeBR_WT <- function(z, Sigma, method = c("mvtnorm", "mvPot", "TruncatedNormal")) {
  weights <- .weightsBR_WT(z = z, Sigma = Sigma, method = method)
  sum(weights / z)
}

expmeXS <- function(z, Sigma, df, method = c("mvtnorm", "mvPot", "TruncatedNormal")) {
  D <- length(z)
  stopifnot(ncol(Sigma) == D | nrow(Sigma) == D | df > 0)
  if (!isTRUE(all.equal(as.vector(diag(Sigma)), rep(1, D)))) {
    warning("Input \"Sigma\" must be a correlation matrix")
    Sigma <- try(cov2cor(Sigma))
    if (inherits(Sigma, what = "try-error")) {
      stop("Could not convert \"Sigma\" to a correlation matrix.")
    }
  }
  weights <- .weightsXstud(z = z, Sigma = Sigma, df = df, method = method)
  sum(weights / z)
}

.weightsHR <- function(z, L, Q, method = c("mvtnorm", "mvPot", "TruncatedNormal")) {
  method <- match.arg(method, choices = c("mvtnorm", "mvPot", "TruncatedNormal"))[1]
  D <- ncol(Q)
  weights <- rep(0, D)
  if (method == "mvPot") {
    genVec <- mvPot::genVecQMC(p = 499, D - 1)
  }
  for (j in 1:D) {
    Qmiinv <- solve(Q[-j, -j])
    weights[j] <- det(Qmiinv)^(0.5) * exp(0.5 * t(L[-j]) %*% Qmiinv %*% L[-j])[1] *
      switch(method,
        mvtnorm = mvtnorm::pmvnorm(upper = c(-Qmiinv %*% L[-j]), sigma = Qmiinv),
        mvPot = mvPot::mvtNormQuasiMonteCarlo(
          p = genVec$primeP,
          upperBound = c(-Qmiinv %*% L[-j]), cov = Qmiinv,
          genVec = genVec$genVec
        )[1],
        TruncatedNormal = TruncatedNormal::mvNqmc(
          l = rep(-Inf, D - 1), n = 1e5,
          u = c(-Qmiinv %*% L[-j]), Sig = Qmiinv
        )$prob
      )
  }
  return(weights)
}


.weightsBR <- function(z, Lambda, method = c("mvtnorm", "mvPot", "TruncatedNormal"), ...) {
  method <- match.arg(method, choices = c("mvtnorm", "mvPot", "TruncatedNormal"))[1]
  ellipsis <- list(...)
  D <- length(z)
  stopifnot(ncol(Lambda) == D | nrow(Lambda) == D)
  weights <- rep(0, D)
  if (method == "mvtnorm") {
    for (j in 1:D) {
      weights[j] <- mvtnorm::pmvnorm(
        lower = rep(-Inf, D - 1), upper = 2 * Lambda[-j, j] + log(z[-j]) - log(z[j]),
        sigma = 2 * (outer(Lambda[-j, j], Lambda[j, -j], FUN = "+") - Lambda[-j, -j])
      )
    }
  } else if (method == "mvPot") {
    if (!is.null(ellipsis$prime)) {
      prime <- ellipsis$prime
      if (!is.null(ellipsis$genvec)) {
        genVec <- ellipsis$genvec
        stopifnot(length(genVec) == D - 1)
      }
    } else {
      prime <- 499L
      genVec <- mvPot::genVecQMC(p = prime, D - 1)$genVec
    }
    for (j in 1:D) {
      weights[j] <- mvPot::mvtNormQuasiMonteCarlo(
        p = prime,
        upperBound = 2 * Lambda[-j, j] + log(z[-j]) - log(z[j]),
        cov = 2 * (outer(Lambda[-j, j], Lambda[j, -j], FUN = "+") - Lambda[-j, -j]),
        genVec = genVec
      )[1]
    }
  } else if (method == "TruncatedNormal") {
    for (j in 1:D) {
      weights[j] <- TruncatedNormal::mvNqmc(
        l = rep(-Inf, D - 1), n = 1e5,
        u = 2 * Lambda[-j, j] + log(z[-j]) - log(z[j]),
        Sig = 2 * (outer(Lambda[-j, j], Lambda[j, -j], FUN = "+") - Lambda[-j, -j])
      )$prob
    }
  }
  return(weights)
}

.weightsBR_WT <- function(z, Sigma, method = c("mvtnorm", "mvPot", "TruncatedNormal"), ...) {
  method <- match.arg(method, choices = c("mvtnorm", "mvPot", "TruncatedNormal"))[1]
  ellipsis <- list(...)
  D <- length(z)
  stopifnot(ncol(Sigma) == D | nrow(Sigma) == D)
  weights <- rep(0, D)
  Ti <- cbind(rep(-1, D - 1), diag(D - 1))
  if (method == "mvPot") {
    genVec <- mvPot::genVecQMC(p = 499, D - 1)
  }
  for (j in 1:D) {
    if (j > 1) {
      Ti[, (j - 1):j] <- Ti[, j:(j - 1)]
    }
    if (method == "mvtnorm") {
      weights[j] <- mvtnorm::pmvnorm(
        lower = rep(-Inf, D - 1),
        upper = log(z[-j] / z[j]) + diag(Sigma)[-j] / 2 + Sigma[j, j] / 2 - Sigma[j, -j],
        sigma = Ti %*% Sigma %*% t(Ti)
      )
    } else if (method == "mvPot") {
      if (!is.null(ellipsis$prime)) {
        prime <- ellipsis$prime
        if (!is.null(ellipsis$genvec)) {
          genVec <- ellipsis$genvec
          stopifnot(length(genVec) == D - 1)
        }
      } else {
        prime <- 499L
        genVec <- mvPot::genVecQMC(p = prime, D - 1)$genVec
      }
      weights[j] <- mvPot::mvtNormQuasiMonteCarlo(
        p = prime,
        upperBound = log(z[-j] / z[j]) + diag(Sigma)[-j] / 2 + Sigma[j, j] / 2 - Sigma[j, -j],
        cov = Ti %*% Sigma %*% t(Ti),
        genVec = genVec
      )[1]
    } else if (method == "TruncatedNormal") {
      weights[j] <- TruncatedNormal::mvNqmc(
        l = rep(-Inf, D - 1), n = 1e5,
        u = log(z[-j] / z[j]) + diag(Sigma)[-j] / 2 + Sigma[j, j] / 2 - Sigma[j, -j],
        Sig = Ti %*% Sigma %*% t(Ti)
      )$prob
    }
  }
  return(weights)
}

.weightsXstud <- function(z, Sigma, df, method = c("mvtnorm", "mvPot", "TruncatedNormal"), ...) {
  method <- match.arg(method, choices = c("mvtnorm", "mvPot", "TruncatedNormal"))[1]
  ellipsis <- list(...)
  D <- nrow(Sigma)
  stopifnot(nrow(Sigma) == length(z))
  weights <- rep(0, D)
  if (method == "mvtnorm") {
    for (j in 1:D) {
      weights[j] <- mvtnorm::pmvt(
        lower = rep(-Inf, D - 1), df = df + 1,
        upper = exp((log(z[-j]) - log(z[j])) / df) - Sigma[-j, j],
        sigma = (Sigma[-j, -j] - Sigma[-j, j, drop = FALSE] %*% Sigma[j, -j, drop = FALSE]) / (df + 1)
      )
    }
  } else if (method == "mvPot") {
    if (!is.null(ellipsis$prime)) {
      prime <- ellipsis$prime
      if (!is.null(ellipsis$genvec)) {
        genVec <- ellipsis$genvec
        stopifnot(length(genVec) == D - 1)
      }
    } else {
      prime <- 499L
      genVec <- mvPot::genVecQMC(p = prime, D - 1)$genVec
    }
    for (j in 1:D) {
      weights[j] <- mvPot::mvTProbQuasiMonteCarlo(
        p = prime,
        upperBound = exp((log(z[-j]) - log(z[j])) / df) - Sigma[-j, j],
        cov = (Sigma[-j, -j] - Sigma[-j, j, drop = FALSE] %*% Sigma[j, -j, drop = FALSE]) / (df + 1),
        nu = df + 1, genVec = genVec
      )[1]
    }
  } else if (method == "TruncatedNormal") {
    for (j in 1:D) {
      weights[j] <- TruncatedNormal::mvTqmc(
        l = rep(-Inf, D - 1), df = df + 1, n = 1e5,
        u = exp((log(z[-j]) - log(z[j])) / df) - Sigma[-j, j],
        Sig = (Sigma[-j, -j] - Sigma[-j, j, drop = FALSE] %*%
          Sigma[j, -j, drop = FALSE]) / (df + 1)
      )$prob
    }
  }
  return(weights)
}
