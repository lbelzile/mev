#' Simulation from Pareto processes using composition sampling
#'
#' The algorithm performs forward sampling by simulating first from a
#' mixture, then sample angles conditional on them being less than (max) or greater than (min) one.
#' The resulting sample from the angular distribution is then multiplied by
#' Pareto variates with tail index \code{shape}.
#'
#' For the moment, only exchangeable models and models based n elliptical processes are handled.
#'
#'
#' The parametrization of the Brown--Resnick is in terms of the matrix \code{Lambda}, which
#' is formed by evaluating the semivariogram \eqn{\gamma} at sites \eqn{s_i, s_j}, meaning that
#'  \eqn{\Lambda_{i,j} = \gamma(s_i, s_j)/2}.
#'
#' @author Leo Belzile
#' @param n sample size.
#' @param Lambda parameter matrix for the Brown--Resnick model. See \bold{Details}.
#' @param Sigma correlation matrix if \code{model = 'xstud'}, otherwise
#' the covariance matrix formed from the stationary Brown-Resnick process.
#' @param df degrees of freedom for extremal Student process.
#' @param model string indicating the model family.
#' @param shape tail index of the Pareto variates (reciprocal shape parameter). Must be strictly positive.
#' @param param parameter value for the logistic or negative logistic model
#' @param d dimension of the multivariate model, only needed for logistic or negative logistic models
#' @param risk string indicating the risk functional. Only \code{max} and \code{min} are currently supported.
#' @param ... additional parameters, currently ignored
#' @details The argument \code{Sigma} is ignored for the Brown-Resnick model
#' if \code{Lambda} is provided by the user.
#' @export
#' @seealso \code{\link{rparp}} for general simulation of Pareto processes based on an accept-reject algorithm.
#' @return an \code{n} by \code{d} matrix of samples, where \code{d = ncol(Sigma)}, with \code{attributes} \code{mixt.weights}.
#' @examples
#' \dontrun{
#' #Brown-Resnick, Wadsworth and Tawn (2014) parametrization
#' D <- 20L
#' coord <- cbind(runif(D), runif(D))
#' semivario <- function(d, alpha = 1.5, lambda = 1){0.5 * (d/lambda)^alpha}
#' Lambda <- semivario(as.matrix(dist(coord))) / 2
#' rparpcs(n = 10, Lambda = Lambda, model = 'br', shape = 0.1)
#' #Extremal Student
#' Sigma <- stats::rWishart(n = 1, df = 20, Sigma = diag(10))[,,1]
#' rparpcs(n = 10, Sigma = cov2cor(Sigma), df = 3, model = 'xstud')
#' }
rparpcs <- function(n,
                    model = c("log", "neglog", "br", "xstud"),
                    risk = c("max", "min"),
                    param = NULL,
                    d,
                    Lambda = NULL,
                    Sigma = NULL,
                    df = NULL,
                    shape = 1,
                    ...) {
  args <- list(...)
  if(!is.null(args$riskf)){
    riskf <- args$riskf
  } else{
    riskf <- risk
  }

  riskf <- match.arg(arg = riskf,
                     choices = c("min", "max"),
                     several.ok = TRUE)[1]
  model <- match.arg(model)

  stopifnot(shape > 0)
  if(model %in% c("br", "xstud")){
  if (!requireNamespace("TruncatedNormal", quietly = TRUE)) {
    stop(
      "Package \"TruncatedNormal\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if(!utils::packageVersion("TruncatedNormal") > "1.2"){
    stop(
      "Please update package \"TruncatedNormal\" to a more recent version.",
      call. = FALSE
    )
  }
  }

    if(model %in% c("log", "neglog")){
      if(missing(d)){
        stop("Invalid dimension for the model.")
      }
      if (is.null(param) || param < 0 || d < 1) {
        stop("Invalid parameter value")
      }
      if (length(param) != 1) {
        warning("Only first entry of param vector considered")
        param <- param[1]
      }
      if (model == "log") {
        if (param < 1) {
          param <- 1/param
        }
      }
      weights <- rep(1/d, d)
      out <- matrix(1, nrow = d, ncol = n)
      # Exchangeable model, balanced mixture
      if(model == "log"){
        scale <- 1/gamma(1-1/param)
        F0 <- scale*rgamma(n,
                     shape = 1-1/param)^(-1/param)

        if(riskf == "min"){
          Flow <- rep(
            x = mev::pgev(
              q = F0,
              loc = scale,
              scale = scale/param,
              shape = 1/param),
          each = d-1)
        Fj <- mev::qgev(
          p = Flow + (1-Flow) * runif(n * (d-1)),
              loc = scale,
              scale = scale/param,
              shape = 1/param)
        } else if(riskf == "max"){
          Fj <- mev::qgev(
            p = runif(n * (d-1)) *
              rep(mev::pgev(F0,
                        loc = scale,
                        scale = scale/param,
                        shape = 1/param), d-1),
            loc = scale,
            scale = scale/param,
            shape = 1/param)
        }
        } else if(model == "neglog"){
          scale <- 1/gamma(1+1/param)
          F0 <- scale*rgamma(n,
                             shape = 1+1/param)^(-1/param)

          if(riskf == "min"){
            Flow <- rep(
              x = pweibull(
                q = F0,
                scale = scale,
                shape = param),
              each = d-1)
            Fj <- qweibull(
              p = Flow + (1-Flow) * runif(n * (d-1)),
              scale = scale,
              shape = param)
          } else if(riskf == "max"){
            Fj <- qweibull(
              p = runif(n * (d-1)) *
                rep(pweibull(
                  q = F0,
                  scale = scale,
                  shape = param), d-1),
              scale = scale,
              shape = param)
          }
        }
        # pos <- (sample.int(n = d, size = n, replace = TRUE)-1) * n +  1:n
        pos <- sample.int(n = d, size = n, replace = TRUE) +
          d*(0:(n-1))
        out[-pos] <- Fj/rep(F0, each = d - 1)
        samp <- t(out) * runif(n)^(-shape)
  } else if(model %in% c("br", "xstud")){

    if (model == "xstud") {
        if (is.null(df)) {
            stop("Invalid degree of freedom argument")
        }
        stopifnot(!is.null(Sigma))
        D <- nrow(Sigma)
        if (!isTRUE(all.equal(as.vector(diag(Sigma)), rep(1, D)))) {
            warning("Input \"Sigma\" must be a correlation matrix")
            Sigma <- cov2cor(Sigma)
        }
        weights <- .weightsXstud(z = rep(1, D), Sigma = Sigma, df = df, riskf = riskf)
    } else if (model == "br") {
        if (!is.null(Lambda)) {
            D <- nrow(Lambda)
            weights <- .weightsBR(z = rep(1, D), Lambda = Lambda, riskf = riskf)
        } else if (!is.null(Sigma)) {
            D <- nrow(Sigma)
            T1 <- cbind(rep(-1, D - 1), diag(D - 1))
            weights <- .weightsBR_WT(z = rep(1, D), Sigma = Sigma, riskf = riskf)
        }
    }
    weights <- weights/sum(weights)
    id <- sample.int(n = D, size = n, replace = TRUE, prob = weights)
    tabu <- table(id)
    bins <- as.integer(names(tabu))
    tabu <- as.vector(tabu)
    # Matrix to store samples
    ang <- matrix(1, nrow = D, ncol = n)
    accu <- 0L
    for (i in 1:length(tabu)) {
        j <- bins[i]
        if (riskf == "max") {
            if (model == "xstud") {
                ang[-j, (accu + 1L):(accu + tabu[i])] <- pmax(TruncatedNormal::mvrandt(n = tabu[i], l = rep(-Inf, D - 1), u = rep(1,
                  D - 1), df = df + 1, Sig = (Sigma[-j, -j] - Sigma[-j, j, drop = FALSE] %*% Sigma[j, -j, drop = FALSE])/(df + 1)),
                  0)^df
            } else if (model == "br") {
                if (!is.null(Lambda)) {
                  ang[-j, (accu + 1L):(accu + tabu[i])] <- exp(TruncatedNormal::mvrandn(n = tabu[i], l = rep(-Inf, D - 1), u = 2 *
                    Lambda[j, -j], Sig = 2 * (outer(Lambda[-j, j], Lambda[j, -j], FUN = "+") - Lambda[-j, -j])) - 2 * Lambda[j, -j])
                } else if (!is.null(Sigma)) {
                  me <- diag(Sigma)[-j]/2 + Sigma[j, j]/2 - Sigma[j, -j]
                  Ti <- T1
                  Ti[, c(1L, j)] <- T1[, c(j, 1L)]
                  ang[-j, (accu + 1L):(accu + tabu[i])] <- exp(TruncatedNormal::mvrandn(n = tabu[i], l = rep(-Inf, D - 1), u = -me,
                    Sig = Ti %*% Sigma %*% t(Ti)) + me)
                }
            }
        } else if (riskf == "min") {
            if (model == "xstud") {
                ang[-j, (accu + 1L):(accu + tabu[i])] <- pmax(TruncatedNormal::mvrandt(n = tabu[i], u = rep(Inf, D - 1), l = rep(1,
                  D - 1), df = df + 1, Sig = (Sigma[-j, -j] - Sigma[-j, j, drop = FALSE] %*% Sigma[j, -j, drop = FALSE])/(df + 1)),
                  0)^df
            } else if (model == "br") {
                if (!is.null(Lambda)) {
                  ang[-j, (accu + 1L):(accu + tabu[i])] <- exp(TruncatedNormal::mvrandn(n = tabu[i], u = rep(Inf, D - 1), Sig = 2 *
                    (outer(Lambda[-j, j], Lambda[j, -j], FUN = "+") - Lambda[-j, -j]), l = 2 * Lambda[j, -j]) - 2 * Lambda[j, -j])
                } else if (!is.null(Sigma)) {
                  me <- diag(Sigma)[-j]/2 + Sigma[j, j]/2 - Sigma[j, -j]
                  Ti <- T1
                  Ti[, c(1L, j)] <- T1[, c(j, 1L)]
                  ang[-j, (accu + 1L):(accu + tabu[i])] <- exp(TruncatedNormal::mvrandn(n = tabu[i], u = rep(Inf, D - 1), l = -me,
                    Sig = Ti %*% Sigma %*% t(Ti)) + me)
                }
            }

        }
        accu <- accu + tabu[i]
    }
    samp <- runif(n)^(-shape) * t(ang[, sample.int(n, n, replace = FALSE)])
}
    attr(samp, "mixt.weights") <- weights
    return(samp)
}
#' Simulation of generalized Huesler-Reiss Pareto vectors via composition sampling
#'
#' Sample from the generalized Pareto process associated to Huesler-Reiss spectral profiles.
#' For the Huesler-Reiss Pareto vectors, the matrix \code{Sigma} is utilized to build \eqn{Q} viz.
#' \deqn{Q = \Sigma^{-1} - \frac{\Sigma^{-1}\mathbf{1}_d\mathbf{1}_d^\top\Sigma^{-1}}{\mathbf{1}_d^\top\Sigma^{-1}\mathbf{1}_d}.}
#' The location vector \code{m} and \code{Sigma} are the parameters of the underlying log-Gaussian process.
#'
#' @param n sample size
#' @param u vector of marginal location parameters (must be strictly positive)
#' @param alpha vector of shape parameters (must be strictly positive).
#' @param Sigma covariance matrix of process, used to define \eqn{Q}. See \bold{Details}.
#' @param m location vector of Gaussian distribution.
#' @return \code{n} by d matrix of observations
#' @references Ho, Z. W. O and C. Dombry (2019), Simple models for multivariate regular variations and the
#'   Huesler-Reiss Pareto distribution, Journal of Multivariate Analysis (\bold{173}), p. 525-550, \doi{10.1016/j.jmva.2019.04.008}
#' @export
#' @examples
#' D <- 20L
#' coord <- cbind(runif(D), runif(D))
#' di <- as.matrix(dist(rbind(c(0, ncol(coord)), coord)))
#' semivario <- function(d, alpha = 1.5, lambda = 1){(d/lambda)^alpha}
#' Vmat <- semivario(di)
#' Sigma <- outer(Vmat[-1, 1], Vmat[1, -1], '+') - Vmat[-1, -1]
#' m <- Vmat[-1,1]
#' \dontrun{
#' samp <- rparpcshr(n = 100, u = c(rep(1, 10), rep(2, 10)),
#'           alpha = seq(0.1, 1, length = 20), Sigma = Sigma, m = m)
#' }
rparpcshr <- function(n, u, alpha, Sigma, m) {
  if (!requireNamespace("TruncatedNormal", quietly = TRUE)) {
    stop(
      "Package \"TruncatedNormal\" must be installed to use this function.",
      call. = FALSE
    )
  }
    D <- ncol(Sigma)
    Sigmainv <- solve(Sigma)
    q <- Sigmainv %*% rep(1, D)
    Q <- (Sigmainv - q %*% t(q)/sum(q))
    l <- c(Sigmainv %*% (((m %*% q - 1)/sum(q))[1] * rep(1, D) - m))
    if (any(u <= 0)) {
        stop("Threshold must be strictly positive")
    }
    if (any(alpha <= 0)) {
        stop("Tail indices must all be strictly positive")
    }
    alpha <- rep(alpha, length = D)
    u <- rep(u, length = D)
    Lst <- l - Q %*% diag(alpha) %*% log(u)  #prop 4.2 b/c u \neq 1
    weights <- .weightsHR(z = rep(1, D), L = Lst, Q = Q)
    weights <- weights/sum(weights)
    id <- sample.int(n = D, size = n, replace = TRUE, prob = weights)
    tabu <- table(id)
    bins <- as.integer(names(tabu))
    tabu <- as.vector(tabu)
    # Matrix to store samples
    ang <- matrix(1, nrow = D, ncol = n)
    accu <- 0L
    for (i in 1:length(tabu)) {
        j <- bins[i]
        Qjinv <- solve(Q[-j, -j])
        ang[-j, (accu + 1L):(accu + tabu[i])] <- exp(TruncatedNormal::mvrandn(n = tabu[i], l = rep(-Inf, D - 1), u = c(-Qjinv %*%
            Lst[-j, ]), Sig = Qjinv) + c(Qjinv %*% Lst[-j, ]))
        accu <- accu + tabu[i]
    }
    R <- mev::rgp(n = n, loc = 1, scale = 1, shape = 1)
    samp <- R * t(ang[, sample.int(n, n, replace = FALSE)])
    for (j in 1:ncol(samp)) {
        samp[, j] <- u[j] * (samp[, j]^(alpha[j]))  #coro 4.3
    }
    return(samp)
}

