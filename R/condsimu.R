#' Conditional distribution of Gaussian or Student subvectors
#'
#' The function samples (truncated) Gaussian or Student vectors
#' corresponding to indices \code{ind} given the values at the remaining index.
#' The location vector \code{mu} and the scale matrix \code{Sigma} are those of the \eqn{d+p} vector.
#' The routine relies on the CDF approximation based on minimax exponential tilting implemented in the \code{TruncatedNormal} package.
#'
#' @param ind a \code{d} vector of indices to impute with integer entries in \eqn{\{1, \ldots, d+p\}}
#' @param x a \code{p} vector with the values of process at remaining coordinates
#' @param lbound \code{d} vector of lower bounds
#' @param ubound \code{d} vector of upper bounds for truncated
#' @param mu \code{d+p} vector of location parameters
#' @param Sigma \code{d+p} by \code{d+p} scale matrix
#' @param df degrees of freedom of the \code{d+p} dimensional Student process
#' @param model string indicating family, either \code{norm} for Gaussian or \code{stud} for Student-t
#' @param n sample size for simulations. Default to 500.
#' @return an n by d matrix of conditional simulations
pcondmvtnorm <- function(n = 500, ind, x, lbound = rep(-Inf, length(ind)), ubound = rep(Inf, length(ind)),
                     mu, Sigma, df = NULL, model = c("norm", "stud")){
  if(length(x) + length(ind) != length(mu)){
    stop("Invalid argument")
  }
  stopifnot(length(ind) > 0)
  x <- as.vector(x)
  mu <- as.vector(mu)
  schurcomp <- function(Sigma, ind){
    stopifnot(c(length(ind)>0, ncol(Sigma)-length(ind)>0))
    Sigma[ind, ind, drop = FALSE] - Sigma[ind, -ind, drop = FALSE] %*%
      solve(Sigma[-ind, -ind, drop = FALSE]) %*% Sigma[-ind, ind, drop = FALSE]
  }
  if(all(c(isTRUE(all.equal(lbound, rep(-Inf, length(ind)))),
           isTRUE(all.equal(ubound, rep(Inf, length(ind))))))){
    return(1)
  }
  if(length(x) == 0){ #Unconditional simulation
    switch(model,
             norm = TruncatedNormal::mvNcdf(n = n, l = lbound - mu, u = ubound - mu, Sig = Sigma)$prob,
             stud = TruncatedNormal::mvTcdf(n = n, l = lbound - mu, u = ubound - mu, Sig = Sigma, df = df)$prob
      )
  } else {
      muC <- c(mu[ind] + Sigma[ind, -ind, drop = FALSE] %*% solve(Sigma[-ind, -ind, drop = FALSE]) %*% (x - mu[-ind]))
      switch(model,
             norm = TruncatedNormal::mvNcdf(n = n, l = lbound - muC, u = ubound - muC,
                                               Sig = schurcomp(Sigma, ind))$prob,
             TruncatedNormal::mvTcdf(n = n, l = lbound - muC, u = ubound - muC,
                                        Sig = c(df + t(x- mu[-ind]) %*% solve(Sigma[-ind, -ind, drop = FALSE]) %*% (x- mu[-ind]))/
                                          (df + length(x)) * schurcomp(Sigma, ind),
                                        df = df + length(x))$prob)
  }
}

dcondmvtnorm <- function(x, ind, lbound = rep(-Inf, length(ind)), ubound = rep(Inf, length(ind)),
                         mu, Sigma, df = NULL, model = c("norm", "stud"), n = 500, log = FALSE){
  if(length(x) != length(mu)){
    stop("Invalid argument")
  }
  #stopifnot(length(ind) > 0)
  x <- as.vector(x)
  mu <- as.vector(mu)
  schurcomp <- function(Sigma, ind){
    stopifnot(c(length(ind)>0, ncol(Sigma)-length(ind)>0))
    Sigma[ind, ind, drop = FALSE] - Sigma[ind, -ind, drop = FALSE] %*%
      solve(Sigma[-ind, -ind, drop = FALSE]) %*% Sigma[-ind, ind, drop = FALSE]
  }
  if(length(x[ind]) == 0){ #Unconditional distribution function
    res <- switch(model,
           norm = mvtnorm::dmvnorm(x = x, mean = mu, sigma = as.matrix(Sigma), log = TRUE),
           stud = mvtnorm::dmvt(x = x - mu, sigma = Sigma, df = df, log = TRUE) #TODO check the centering

    )
  } else{
    muC <- c(mu[ind] + Sigma[ind, -ind, drop = FALSE] %*% solve(Sigma[-ind, -ind, drop = FALSE]) %*% (x[-ind] - mu[-ind]))
    sigmaC <- as.matrix(schurcomp(Sigma, ind))
    if(length(ind)==1L){
      res <- switch(model,
                    norm = mvtnorm::dmvnorm(x = x[-ind], mean = mu[-ind], sigma = solve(Sigma[-ind, -ind, drop = FALSE]), log = TRUE) +
                      TruncatedNormal:::lnNpr(a = (lbound - muC)/sqrt(sigmaC), b = (ubound - muC)/sqrt(sigmaC)),
                    stud = mvtnorm::dmvt(x = x[-ind]- mu[-ind], sigma = Sigma[-ind, -ind, drop = FALSE], df = df, log = TRUE) +
                      log(TruncatedNormal::mvTcdf(n = n, l = lbound - muC, u = ubound - muC,
                                                  Sig = c(df + t(x[-ind]- mu[-ind]) %*% solve(Sigma[-ind, -ind, drop = FALSE]) %*% (x[-ind]- mu[-ind]))/
                                                    (df + length(x[-ind])) * sigmaC,
                                                  df = df + length(x[-ind]))$prob))
    } else {
      res <- switch(model,
                    norm = mvtnorm::dmvnorm(x = x[-ind], mean = mu[-ind], sigma = Sigma[-ind,-ind, drop = FALSE], log = TRUE) +
                      log(TruncatedNormal::mvNcdf(n = n, l = lbound - muC, u = ubound - muC, Sig = sigmaC)$prob),
                    stud = mvtnorm::dmvt(x = x[-ind]- mu[-ind], sigma = Sigma[-ind, -ind, drop = FALSE], df = df, log = TRUE) +
                      log(TruncatedNormal::mvTcdf(n = n, l = lbound - muC, u = ubound - muC,
                                                  Sig = c(df + t(x[-ind]- mu[-ind]) %*% solve(Sigma[-ind, -ind, drop = FALSE]) %*% (x[-ind]- mu[-ind]))/
                                                    (df + length(x[-ind])) * sigmaC,
                                                  df = df + length(x[-ind]))$prob))
  }

  }
  return(ifelse(log, res, exp(res)))
}



#' Conditional distribution of Gaussian or Student subvectors
#'
#' The function samples (truncated) Gaussian or Student vectors
#' corresponding to indices \code{ind} given the values at the remaining index.
#' The location vector \code{mu} and the scale matrix \code{Sigma} are those of the \eqn{d+p} vector.
#' The routine relies on the CDF approximation based on minimax exponential tilting implemented in the \code{TruncatedNormal} package.
#'
#' @param ind a \code{d} vector of indices to impute with integer entries in \eqn{\{1, \ldots, d+p\}}
#' @param x a \code{p} vector with the values of process at remaining coordinates
#' @param lbound \code{d} vector of lower bounds
#' @param ubound \code{d} vector of upper bounds for truncated
#' @param mu \code{d+p} vector of location parameters
#' @param Sigma \code{d+p} by \code{d+p} scale matrix
#' @param df degrees of freedom of the \code{d+p} dimensional Student process
#' @param model string indicating family, either \code{norm} for Gaussian or \code{stud} for Student-t
#' @param n sample size for the random vector; default to 1.
#' @return an n by d matrix of conditional simulations
rcondmvtnorm <- function(n = 1L, ind, x, lbound = rep(-Inf, length(ind)), ubound = rep(Inf, length(ind)),
                         mu, Sigma, df = NULL, model = c("norm", "stud"), log = FALSE){
  model <- match.arg(model)
  if(length(x) + length(ind) != length(mu)){
    stop("Invalid argument")
  }
  stopifnot(length(ind) > 0)
  x <- as.vector(x)
  mu <- as.vector(mu)
  schurcomp <- function(Sigma, ind){
    stopifnot(c(length(ind)>0, ncol(Sigma)-length(ind)>0))
    Sigma[ind, ind, drop=FALSE] - Sigma[ind,-ind, drop=FALSE] %*% solve(Sigma[-ind, -ind, drop=FALSE]) %*% Sigma[-ind,ind, drop=FALSE]
  }
 if(length(x) == 0){ #Unconditional simulation
    switch(model,
           norm = TruncatedNormal::mvrandn(n = n, mu = mu, l = lbound, u = ubound, Sig = Sigma),
           stud = TruncatedNormal::mvrandt(n = n, mu = mu, l = lbound, u = ubound, Sig = Sigma, df = df)
    )
  } else {
    muC <- c(mu[ind] + Sigma[ind, -ind, drop=FALSE] %*% solve(Sigma[-ind, -ind, drop=FALSE]) %*% (x - mu[-ind]))
    switch(model,
           norm = TruncatedNormal::mvrandn(n = n, mu = muC, l = lbound, u = ubound,
                                          Sig = schurcomp(Sigma, ind)),
           stud = TruncatedNormal::mvrandt(n = n, l = lbound, mu = muC, u = ubound,
                                   Sig = c(df + t(x- mu[-ind]) %*% solve(Sigma[-ind, -ind, drop=FALSE]) %*% (x- mu[-ind]))/
                                     (df + length(x)) * schurcomp(Sigma, ind),
                                   df = df + length(x))
    )

  }
}


#' Proposals for random walk Metropolis-Hastings
#'
#' This function transforms a vector \code{cur} to an unconstrained scale based on \code{lbound} and \code{ubound},
#' then samples draws from a multivariate Gaussian vector with covariance matrix \code{cov} centered at the current (unconstrained) value.
#' The function then returns a list containing the proposal on the unconstrained scale (\code{trprop}),
#' on the original scale (\code{prop}) along with the log of the jacobian of the transformation (\code{logjac}).
propRWMH<- function(cur, cov, lbound = rep(-Inf, ncol(cov)), ubound = rep(Inf, ncol(cov))){
  d <- length(cur)
  stopifnot(d == dim(cov), length(lbound) == d, length(ubound) == d)
  tr.cur <- rep(0, d);
  prop <- rep(0, d)
  logratio <- rep(0, d)
  expit <- function(x){1/(1+exp(-x))}
  logit <- function(x){log(x)-log(1-x)}
  for(j in 1:d){
    if(lbound[j] == -Inf && ubound[j] == Inf){
      tr.cur[j] <- cur[j]
    } else if(lbound[j] > -Inf && ubound[j] == Inf){
      tr.cur[j] <- log(cur[j] - lbound[j])
    } else if(lbound[j] == -Inf && ubound[j] < Inf){
      tr.cur[j] <- log(ubound[j] - cur[j])
    } else{
      tr.cur[j] <-  logit((cur[j] - lbound[j]) / (ubound[j] - lbound[j]))
    }
  }
  if(d == 1L){
    tr.prop <- rnorm(n = 1, mean = tr.cur, sd =  sqrt(cov))
  } else{
    tr.prop <- mev::mvrnorm(n = 1, mu = tr.cur, Sigma = cov)
  }
  for(j in 1:d){
    if(lbound[j] == -Inf && ubound[j] == Inf){
      prop[j] <- tr.prop[j]
    } else if(lbound[j] > -Inf && ubound[j] == Inf){
      prop[j] <- lbound[j] + exp(tr.prop[j])
      logratio[j] <- tr.prop[j] - tr.cur[j]
    } else if(lbound[j] == -Inf && ubound[j] < Inf){
      prop[j] <- ubound[j] - exp(tr.prop[j])
      logratio[j] <- tr.prop[j] - tr.cur[j]
    } else{
      prop[j] <- lbound[j] + (ubound[j] - lbound[j]) * expit(tr.prop[j])
      logratio[j] <- log(1/(cur[j] - lbound[j]) + 1/(ubound[j] - cur[j])) - log(1/(prop[j] - lbound[j]) + 1/(ubound[j] - prop[j]))
    }
  }
  return(list(prop = prop, trprop = tr.prop, logjac = logratio))
}
