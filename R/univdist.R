#' Generalized Pareto distribution
#'
#' Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized Pareto distribution
#'
#' @author Leo Belzile
#' @name gpd
#' @param par vector of \code{scale} and \code{shape}
#' @param dat sample vector
#' @param tol numerical tolerance for the exponential model
#' @param method string indicating whether to use the expected  (\code{'exp'}) or the observed (\code{'obs'} - the default) information matrix.
#' @param V vector calculated by \code{gpd.Vfun}
#' @param n sample size
#' @section Usage: \preformatted{gpd.ll(par, dat, tol=1e-5)
#' gpd.ll.optim(par, dat, tol=1e-5)
#' gpd.score(par, dat)
#' gpd.infomat(par, dat, method = c('obs','exp'))
#' gpd.bias(par, n)
#' gpd.Fscore(par, dat, method = c('obs','exp'))
#' gpd.Vfun(par, dat)
#' gpd.phi(par, dat, V)
#' gpd.dphi(par, dat, V)}
#'
#' @section Functions:
#'
#' \itemize{
#' \item{\code{gpd.ll}:} {log likelihood}
#' \item{\code{gpd.ll.optim}:} {negative log likelihood parametrized in terms of \code{log(scale)} and shape
#' in order to perform unconstrained optimization}
#' \item{\code{gpd.score}:} {score vector}
#' \item{\code{gpd.infomat}:} {observed or expected information matrix}
#' \item{\code{gpd.bias}:} {Cox-Snell first order bias}
#' \item{\code{gpd.Fscore}:} {Firth's modified score equation}
#' \item{\code{gpd.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gpd.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gpd.dphi}:} {derivative matrix of the canonical parameter in the local
#' exponential family approximation}
#' }
#' @references Firth, D. (1993). Bias reduction of maximum likelihood estimates, \emph{Biometrika}, \strong{80}(1), 27--38.
#' @references Coles, S. (2001). \emph{An Introduction to Statistical Modeling of Extreme Values}, Springer, 209 p.
#' @references Cox, D. R. and E. J. Snell (1968). A general definition of residuals, \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \strong{30}, 248--275.
#' @references Cordeiro, G. M. and R. Klein (1994). Bias correction in ARMA models, \emph{Statistics and Probability Letters}, \strong{19}(3), 169--176.
#' @references Giles, D. E., Feng, H. and R. T. Godwin (2016).  Bias-corrected maximum likelihood estimation of the  parameters of the generalized Pareto distribution, \emph{Communications in Statistics - Theory and Methods}, \strong{45}(8), 2465--2483.
#'
NULL



#' @title Generalized extreme value distribution
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized extreme value distribution
#'
#' @name gev
#' @param par vector of \code{loc}, \code{scale} and \code{shape}
#' @param dat sample vector
#' @param method string indicating whether to use the expected  (\code{'exp'}) or the observed (\code{'obs'} - the default) information matrix.
#' @param V vector calculated by \code{gev.Vfun}
#' @param n sample size
#' @param p vector of probabilities
#' @section Usage: \preformatted{gev.ll(par, dat)
#' gev.ll.optim(par, dat)
#' gev.score(par, dat)
#' gev.infomat(par, dat, method = c('obs','exp'))
#' gev.retlev(par, p)
#' gev.bias(par, n)
#' gev.Fscore(par, dat, method=c('obs','exp'))
#' gev.Vfun(par, dat)
#' gev.phi(par, dat, V)
#' gev.dphi(par, dat, V)}
#'
#'
#' @section Functions:
#' \itemize{
#' \item{\code{gev.ll}:} {log likelihood}
#' \item{\code{gev.ll.optim}:} {negative log likelihood parametrized in terms of location, \code{log(scale)} and shape
#' in order to perform unconstrained optimization}
#' \item{\code{gev.score}:} {score vector}
#' \item{\code{gev.infomat}:} {observed or expected information matrix}
#' \item{\code{gev.retlev}:} {return level, corresponding to the \eqn{(1-p)}th quantile}
#' \item{\code{gev.bias}:} {Cox-Snell first order bias}
#' \item{\code{gev.Fscore}:} {Firth's modified score equation}
#' \item{\code{gev.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gev.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gev.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
#' @references Firth, D. (1993). Bias reduction of maximum likelihood estimates, \emph{Biometrika}, \strong{80}(1), 27--38.
#' @references Coles, S. (2001). \emph{An Introduction to Statistical Modeling of Extreme Values}, Springer, 209 p.
#' @references Cox, D. R. and E. J. Snell (1968). A general definition of residuals, \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \strong{30}, 248--275.
#' @references Cordeiro, G. M. and R. Klein (1994). Bias correction in ARMA models, \emph{Statistics and Probability Letters}, \strong{19}(3), 169--176.
NULL


#' log likelihood for the generalized Pareto distribution
#'
#' Function returning the density of an \code{n} sample from the GP distribution.
#' \code{gpd.ll.optim} returns the negative log likelihood parametrized in terms of \code{log(scale)} and shape
#' in order to perform unconstrained optimization
#' @seealso \code{\link{gpd}}
#' @inheritParams gpd
#' @export
#' @keywords internal
gpd.ll <- function(par, dat, tol = 1e-05) {
    sigma = par[1]
    if(sigma < 0){
      return(-Inf)
    }
    xi = par[2]
    n <- length(dat)
    if (abs(xi) > tol) {
      if(xi > -1){
        #Most scenarios
        as.vector(-n * log(sigma) - (1 + 1/xi) * sum(log(pmax(1 + xi/sigma * dat, 0))))
      } else if(isTRUE(all.equal(xi, -1, check.attributes = FALSE))){
        -n*log(max(dat))
      } else{
       -Inf
      }
    } else {
      as.vector(-n * log(sigma) - sum(dat)/sigma)
    }
}



#' @rdname gpd.ll
#' @inheritParams gpd
#' @export
gpd.ll.optim <- function(par, dat, tol = 1e-05) {
    sigma = exp(par[1])
    xi = par[2]
    if (abs(xi) > tol) {
        length(dat) * log(sigma) + (1 + 1/xi) * sum(log(pmax(1 + xi/sigma * dat, 0)))
    } else {
        length(dat) * log(sigma) + sum(dat)/sigma
    }
}
#' Score vector for the generalized Pareto distribution
#'
#' @seealso \code{\link{gpd}}
#' @inheritParams gpd
#' @export
#' @keywords internal
gpd.score <- function(par, dat) {
    sigma = par[1]
    xi = as.vector(par[2])
    if (!isTRUE(all.equal(0, xi))) {
        c(sum(dat * (xi + 1)/(sigma^2 * (dat * xi/sigma + 1)) - 1/sigma), sum(-dat * (1/xi + 1)/(sigma *
            (dat * xi/sigma + 1)) + log(pmax(dat * xi/sigma + 1, 0))/xi^2))
    } else {
        c(sum((dat - sigma)/sigma^2), sum(1/2 * (dat - 2 * sigma) * dat/sigma^2))
    }
}

#' Information matrix for the generalized Pareto distribution
#'
#' The function returns the expected or observed information matrix.
#' @seealso \code{\link{gpd}}
#' @inheritParams gpd
#' @param nobs number of observations
#' @export
#' @keywords internal
gpd.infomat <- function(par, dat, method = c("obs", "exp"), nobs = length(dat)) {
    method <- match.arg(method)
    dat <- as.vector(dat)
    sigma <- as.vector(par[1])
    xi <- as.vector(par[2])
    if (method == "obs") {
        if (any((1 + xi * dat/sigma) < 0)) {
            stop("Data outside of range specified by parameter, yielding a zero likelihood")
        }
        c11 <- (length(dat) - (1 + xi) * sum(dat * (2 * sigma + xi * dat)/(sigma + xi * dat)^2))/(sigma^2)
        if (!isTRUE(all.equal(xi, 0))) {
            c12 <- (1 + xi) * sum(dat/(sigma + xi * dat)^2)/xi - sum(dat/(sigma + xi * dat))/(sigma *
                xi)
            c22 <- (2/xi^2) * sum(dat/(sigma + xi * dat)) - 2 * sum(log(pmax(dat * xi/sigma + 1, 0)))/(xi^3) +
                (1 + 1/xi) * sum(dat^2/(sigma + xi * dat)^2)
        } else {
            c12 <- sum(-(dat^2 - dat * sigma)/sigma^3)
            c22 <- sum(-1/3 * (2 * dat - 3 * sigma) * dat^2/sigma^3)
        }
        -matrix(c(c11, c12, c12, c22), nrow = 2, ncol = 2, byrow = TRUE)
    } else if (method == "exp") {
        k22 <- -2/((1 + xi) * (1 + 2 * xi))
        k11 <- -1/(sigma^2 * (1 + 2 * xi))
        k12 <- -1/(sigma * (1 + xi) * (1 + 2 * xi))
        -nobs * cbind(c(k11, k12), c(k12, k22))  #fixed 12-10-2016
    }
}
#' Tangent exponential model statistics for the generalized Pareto distribution
#'
#' Matrix of approximate ancillary statistics, sample space derivative of the
#' log likelihood and mixed derivative for the generalized Pareto distribution.
#' @seealso \code{\link{gpd}}
#' @inheritParams gpd
#' @export
#' @name gpd.temstat
#' @keywords internal
gpd.Vfun <- function(par, dat) {
    sigma = par[1]
    xi = par[2]
    cbind(dat/sigma, sigma * (dat * xi/sigma + 1)^(-1/xi) * (log(dat * xi/sigma + 1)/xi^2 - dat/(sigma *
        (dat * xi/sigma + 1) * xi))/(dat * xi/sigma + 1)^(-1/xi - 1))
}

#' @inheritParams gpd
#' @rdname gpd.temstat
#' @export
gpd.phi <- function(par, dat, V) {
    sigma = par[1]
    xi = par[2]
    rbind(-(xi + 1)/(sigma * (dat * xi/sigma + 1))) %*% V
}

#' @inheritParams gpd
#' @export
#' @rdname gpd.temstat
gpd.dphi <- function(par, dat, V) {
    sigma = par[1]
    xi = as.vector(par[2])
    if (!isTRUE(all.equal(xi, 0))) {
        rbind(xi * (1/xi + 1)/(sigma^2 * (dat * xi/sigma + 1)) - dat * xi^2 * (1/xi + 1)/(sigma^3 * (dat *
            xi/sigma + 1)^2), -(1/xi + 1)/(sigma * (dat * xi/sigma + 1)) + dat * (xi + 1)/(sigma^2 *
            (dat * xi/sigma + 1)^2) + 1/(sigma * (dat * xi/sigma + 1) * xi)) %*% V
    } else {
        rbind(rep(sigma^(-2), length(dat)), (dat - sigma)/sigma^2) %*% V
    }
}


#' log likelihood for the generalized extreme value distribution
#'
#' Function returning the density of an \code{n} sample from the GEV distribution.
#'
#' \code{gev.ll.optim} returns the negative log likelihood parametrized in terms of location, \code{log(scale)} and shape in order to perform unconstrained optimization
#'
#' @inheritParams gev
#' @export
#' @keywords internal
#' @seealso \code{\link{gev}}
gev.ll <- function(par, dat) {
    dat <- as.vector(dat)
    par <- as.vector(par)
    tx <- pmax(1 + par[3]/par[2] * (dat - par[1]), 0)
    if (!isTRUE(all.equal(as.vector(par[3]), 0, tolerance = 1e-07))) {
        sum(-log(par[2]) - (1/par[3] + 1) * log(tx) - tx^(-1/par[3]))
    } else {
        sum(-log(par[2]) - (dat - par[1])/par[2] - exp(-(dat - par[1])/par[2]))
    }
}

#' @rdname gev.ll
#' @inheritParams gev
#' @export
#' @keywords internal
gev.ll.optim <- function(par, dat) {
    tpar = par
    tpar[2] = exp(par[2])
    # parameters with log-scale
    nll = -gev.ll(tpar, dat)
    return(nll)
}

#' Score vector for the generalized extreme value distribution
#'
#' @inheritParams gev
#' @export
#' @keywords internal
gev.score <- function(par, dat) {
    mu = par[1]
    sigma = par[2]
    xi = as.vector(par[3])
    if (!isTRUE(all.equal(xi, 0, tolerance = 1e-10))) {
        c(sum(-(pmax(0, -(mu - dat) * xi/sigma + 1))^(-1/xi - 1)/sigma - xi * (1/xi + 1)/(sigma * ((mu -
            dat) * xi/sigma - 1))),
          sum(-(dat - mu) * ((dat - mu) * xi/sigma + 1)^(-1/xi - 1)/sigma^2 +
            (dat - mu) * (xi + 1)/(sigma^2 * ((dat - mu) * xi/sigma + 1)) - 1/sigma),
          sum(-(mu - dat) * (1/xi + 1)/(sigma * ((mu - dat) * xi/sigma - 1)) -
                (log1p(-(mu - dat) * xi/sigma)/xi^2 - (mu - dat)/(sigma * ((mu - dat) * xi/sigma - 1) * xi))/
                (-(mu - dat) * xi/sigma + 1)^(1/xi) + log1p(-(mu - dat) * xi/sigma)/xi^2))
    } else {
        c(sum(-exp(mu/sigma - dat/sigma)/sigma + 1/sigma),
          sum(mu * exp(mu/sigma - dat/sigma)/sigma^2 -
            dat * exp(mu/sigma - dat/sigma)/sigma^2 - mu/sigma^2 - 1/sigma + dat/sigma^2),
          sum((exp(-(dat/sigma)) *
            (dat - mu) * (exp(mu/sigma) * (mu - dat) + exp(dat/sigma) * (dat - mu - 2 * sigma)))/(2 * sigma^2)))
    }

}

#' Information matrix for the generalized extreme value distribution
#'
#' The function returns the expected or observed information matrix.
#' @inheritParams gev
#' @param nobs number of observations
#' @export
#' @keywords internal
gev.infomat <- function(par, dat, method = c("obs", "exp"), nobs = length(dat)) {
    method <- match.arg(method)
    par <- as.vector(par)
    if (method == "exp") {
        # (Expected) Fisher information does not depend on location parameter
        if (length(par) == 3) {
            sigma = par[2]
            xi = as.vector(par[3])
        } else {
            sigma = par[1]
            xi = as.vector(par[2])
        }
      if(xi < -0.5){
        stop("The Fisher information is not defined if the shape parameter of the GEV is less than -0.5")
      }
        # Limiting case when xi=0 Fix to value at zero
        if (abs(xi) < 0.001) {
            return(nobs * cbind(c(1/sigma^2, -0.422784335098467/sigma^2, 0.41184033042644/sigma), c(-0.422784335098467/sigma^2,
                1.82368066085288/sigma^2, 0.332484907160274/sigma), c(0.41184033042644/sigma, 0.332484907160274/sigma,
                2.42360605517703)))
        }
        p = (1 + xi)^2 * gamma(1 + 2 * xi)
        q = gamma(2 + xi) * (digamma(1 + xi) + (1 + xi)/xi)
        infomat <- nobs * cbind(c(p/sigma^2, -(p - gamma(2 + xi))/(sigma^2 * xi), (p/xi - q)/(sigma *
            xi)), c(-(p - gamma(2 + xi))/(sigma^2 * xi), (1 - 2 * gamma(2 + xi) + p)/(sigma^2 * xi^2),
            -(1 + digamma(1) + (1 - gamma(2 + xi))/xi - q + p/xi)/(sigma * xi^2)), c((p/xi - q)/(sigma *
            xi), -(1 + digamma(1) + (1 - gamma(2 + xi))/xi - q + p/xi)/(sigma * xi^2), (pi^2/6 + (1 +
            digamma(1) + 1/xi)^2 - 2 * q/xi + p/xi^2)/xi^2))
        if (xi < 0.003 & xi > 0) {
            # Linearization because numerically unstable in this region
            y2 <- (pi^2/6 + (1 + digamma(1) + 1/0.003)^2 - 2 * (gamma(2 + 0.003) * (digamma(1 + 0.003) +
                (1 + 0.003)/0.003))/0.003 + (1 + 0.003)^2 * gamma(1 + 2 * 0.003)/0.003^2)/0.003^2
            y1 <- 2.42360605517703
            infomat[3, 3] <- nobs * ((y2 - y1)/0.003 * xi + y1)
        } else if (xi > -0.003 & xi < 0) {
            y2 <- (pi^2/6 + (1 + digamma(1) - 1/0.003)^2 + 2 * (gamma(2 - 0.003) * (digamma(1 - 0.003) -
                (1 - 0.003)/0.003))/0.003 + (1 - 0.003)^2 * gamma(1 - 2 * 0.003)/0.003^2)/0.003^2
            y1 <- 2.42360605517703
            infomat[3, 3] <- nobs * (-(y2 - y1)/0.003 * xi + y1)
        }
        colnames(infomat) <- rownames(infomat) <- c("loc","scale","shape")
        return(infomat)
    } else if (method == "obs") {
        if (length(par) != 3) {
            stop("Invalid parameter vector")
        }
        mu = par[1]
        sigma = par[2]
        xi = as.vector(par[3])
        if(xi < -0.5){
          stop("The observed information matrix is not defined if the shape parameter of the GEV is less than -0.5")
        }
        # Bug fixed 21-10-2016 (parameter were defined after they were used).
        if (any((1 + xi * (dat - mu)/sigma) < 0)) {
            stop("Data outside of range specified by parameter, yielding a zero likelihood")
        }
        infomat <- matrix(0, ncol = 3, nrow = 3)
        if (!isTRUE(all.equal(xi, 0))) {
            infomat[1, 1] <- sum(-((dat - mu) * xi/sigma + 1)^(-1/xi - 2) * (xi + 1)/sigma^2 + xi^2 *
                (1/xi + 1)/(sigma^2 * ((dat - mu) * xi/sigma + 1)^2))
            infomat[1, 2] <- infomat[2, 1] <- sum(-(dat - mu) * ((dat - mu) * xi/sigma + 1)^(-1/xi -
                2) * (xi + 1)/sigma^3 + ((dat - mu) * xi/sigma + 1)^(-1/xi - 1)/sigma^2 - xi * (1/xi +
                1)/(sigma^2 * ((dat - mu) * xi/sigma + 1)) + (dat - mu) * xi^2 * (1/xi + 1)/(sigma^3 *
                ((dat - mu) * xi/sigma + 1)^2))
            infomat[1, 3] <- infomat[3, 1] <- sum((-(mu - dat) * xi/sigma + 1)^(-1/xi - 1) * ((mu - dat) *
                (1/xi + 1)/(sigma * ((mu - dat) * xi/sigma - 1)) - log1p(-(mu - dat) * xi/sigma)/xi^2)/sigma -
                (1/xi + 1)/(sigma * ((mu - dat) * xi/sigma - 1)) + (mu - dat) * (xi + 1)/(sigma^2 * ((mu -
                dat) * xi/sigma - 1)^2) + 1/(sigma * ((mu - dat) * xi/sigma - 1) * xi))
            infomat[2, 3] <- infomat[3, 2] <- sum((mu - dat) * (-(mu - dat) * xi/sigma + 1)^(-1/xi -
                1) * (log1p(-(mu - dat) * xi/sigma)/xi^2 - (mu - dat)/(sigma * ((mu - dat) * xi/sigma -
                1) * xi))/sigma^2 + (mu - dat) * (1/xi + 1)/(sigma^2 * ((mu - dat) * xi/sigma - 1)) -
                (mu - dat)^2 * (xi + 1)/(sigma^3 * ((mu - dat) * xi/sigma - 1)^2) - (mu - dat)/(sigma^2 *
                ((mu - dat) * xi/sigma - 1) * xi) + (mu - dat)^2/(sigma^3 * (-(mu - dat) * xi/sigma +
                1)^(1/xi) * ((mu - dat) * xi/sigma - 1)^2))
            infomat[3, 3] <- sum(-(log1p(-(mu - dat) * xi/sigma)/xi^2 - (mu - dat)/(sigma * ((mu -
                dat) * xi/sigma - 1) * xi))^2/(-(mu - dat) * xi/sigma + 1)^(1/xi) + (2 * log1p(-(mu - dat) *
                xi/sigma)/xi^3 - 2 * (mu - dat)/(sigma * ((mu - dat) * xi/sigma - 1) * xi^2) - (mu -
                dat)^2/(sigma^2 * ((mu - dat) * xi/sigma - 1)^2 * xi))/(-(mu - dat) * xi/sigma + 1)^(1/xi) +
                (mu - dat)^2 * (1/xi + 1)/(sigma^2 * ((mu - dat) * xi/sigma - 1)^2) - 2 * log1p(-(mu -
                dat) * xi/sigma)/xi^3 + 2 * (mu - dat)/(sigma * ((mu - dat) * xi/sigma - 1) * xi^2))
            infomat[2, 2] <- sum(-(mu - dat)^2 * (-(mu - dat) * xi/sigma + 1)^(-1/xi - 2) * (xi + 1)/sigma^4 -
                2 * (mu - dat) * (-(mu - dat) * xi/sigma + 1)^(-1/xi - 1)/sigma^3 - 2 * (mu - dat) *
                (xi + 1)/(sigma^3 * ((mu - dat) * xi/sigma - 1)) + (mu - dat)^2 * xi^2 * (1/xi + 1)/(sigma^4 *
                ((mu - dat) * xi/sigma - 1)^2) + 1/sigma^2)
        } else {
            infomat[1, 1] <- sum(-exp(mu/sigma - dat/sigma)/sigma^2)
            infomat[1, 2] <- infomat[2, 1] <- sum(((mu + sigma) * exp(mu/sigma) - dat * exp(mu/sigma) -
                sigma * exp(dat/sigma)) * exp(-dat/sigma)/sigma^3)
            infomat[2, 2] <- sum(-2 * mu * exp(mu/sigma - dat/sigma)/sigma^3 + 2 * dat * exp(mu/sigma -
                dat/sigma)/sigma^3 + 2 * mu/sigma^3 + 1/sigma^2 - 2 * dat/sigma^3 - (mu^2 * exp(mu/sigma) -
                2 * mu * dat * exp(mu/sigma) + dat^2 * exp(mu/sigma)) * exp(-dat/sigma)/sigma^4)
            infomat[1, 3] <- infomat[3, 1] <- sum(-1/2 * (mu^2 * exp(mu/sigma) + 2 * mu * sigma * exp(mu/sigma) -
                2 * mu * dat * exp(mu/sigma) - 2 * sigma * dat * exp(mu/sigma) + dat^2 * exp(mu/sigma) -
                2 * mu * sigma * exp(dat/sigma) - 2 * sigma^2 * exp(dat/sigma) + 2 * sigma * dat * exp(dat/sigma)) *
                exp(-dat/sigma)/sigma^3)
            infomat[2, 3] <- infomat[3, 2] <- sum(1/2 * (mu^2 * exp(mu/sigma) + 2 * mu * sigma * exp(mu/sigma) -
                2 * mu * dat * exp(mu/sigma) - 2 * sigma * dat * exp(mu/sigma) + dat^2 * exp(mu/sigma) -
                2 * mu * sigma * exp(dat/sigma) - 2 * sigma^2 * exp(dat/sigma) + 2 * sigma * dat * exp(dat/sigma)) *
                (mu - dat) * exp(-dat/sigma)/sigma^4)
            infomat[3, 3] <- sum(exp(-(dat/sigma))*(4*exp(dat/sigma)* sigma*(2*mu + 3*sigma - 2*dat) -
                                 exp(mu/sigma)*(3*mu + 8*sigma - 3*dat)*(mu - dat))*(mu - dat)^2)/(12*sigma^4)
        }
        colnames(infomat) <- rownames(infomat) <- c("loc","scale","shape")
        return(-infomat)
    }
}


#' Tangent exponential model statistics for the generalized extreme value distribution
#'
#' Matrix of approximate ancillary statistics, sample space derivative of the
#' log likelihood and mixed derivative for the generalized extreme value distribution.
#' @seealso \code{\link{gev}}
#' @inheritParams gev
#' @export
#' @name gev.temstat
#' @keywords internal
gev.Vfun <- function(par, dat) {
    cbind(1, (dat - par[1])/par[2], par[2] * (-(par[1] - dat) * par[3]/par[2] + 1)^(-1/par[3]) * (log1p(-(par[1] -
        dat) * par[3]/par[2] )/par[3]^2 - (par[1] - dat)/(par[2] * ((par[1] - dat) * par[3]/par[2] -
        1) * par[3]))/(-(par[1] - dat) * par[3]/par[2] + 1)^(-1/par[3] - 1))
}

#' @rdname gev.temstat
#' @inheritParams gev
#' @export
#' @keywords internal
gev.phi <- function(par, dat, V) {
    mu = par[1]
    sigma = par[2]
    xi = as.vector(par[3])
    if (!isTRUE(all.equal(xi, 0))) {
        t(((dat - mu) * xi/sigma + 1)^(-1/xi - 1)/sigma + xi * (1/xi + 1)/(sigma * ((mu - dat) * xi/sigma -
            1))) %*% V
    } else {
        t(exp((mu - dat)/sigma)/sigma - 1/sigma) %*% V

    }
}


#' @rdname gev.temstat
#' @inheritParams gev
#' @export
#' @keywords internal
gev.dphi <- function(par, dat, V) {
    mu = par[1]
    sigma = par[2]
    xi = as.vector(par[3])
    if (!isTRUE(all.equal(xi, 0))) {
        rbind((-(mu - dat) * xi/sigma + 1)^(-1/xi - 2) * (xi + 1)/sigma^2 - xi^2 * (1/xi + 1)/(sigma^2 *
            ((mu - dat) * xi/sigma - 1)^2), -(mu - dat) * (-(mu - dat) * xi/sigma + 1)^(-1/xi - 2) *
            (xi + 1)/sigma^3 - (-(mu - dat) * xi/sigma + 1)^(-1/xi - 1)/sigma^2 - xi * (1/xi + 1)/(sigma^2 *
            ((mu - dat) * xi/sigma - 1)) + (mu - dat) * xi^2 * (1/xi + 1)/(sigma^3 * ((mu - dat) * xi/sigma -
            1)^2), -(-(mu - dat) * xi/sigma + 1)^(-1/xi - 1) * ((mu - dat) * (1/xi + 1)/(sigma * ((mu -
            dat) * xi/sigma - 1)) - log1p(-(mu - dat) * xi/sigma)/xi^2)/sigma + (1/xi + 1)/(sigma *
            ((mu - dat) * xi/sigma - 1)) - (mu - dat) * (xi + 1)/(sigma^2 * ((mu - dat) * xi/sigma -
            1)^2) - 1/(sigma * ((mu - dat) * xi/sigma - 1) * xi)) %*% V
    } else {
        rbind(exp(-dat/sigma + mu/sigma)/sigma^2, (sigma * exp(dat/sigma) + (dat - mu - sigma) * exp(mu/sigma)) *
            exp(-dat/sigma)/sigma^3, 1/2 * (2 * dat * sigma * exp(dat/sigma) - 2 * mu * sigma * exp(dat/sigma) -
            2 * sigma^2 * exp(dat/sigma) + dat^2 * exp(mu/sigma) - 2 * dat * mu * exp(mu/sigma) + mu^2 *
            exp(mu/sigma) - 2 * dat * sigma * exp(mu/sigma) + 2 * mu * sigma * exp(mu/sigma)) * exp(-dat/sigma)/sigma^3) %*%
            V
    }
}

#' @title Generalized Pareto distribution (expected shortfall parametrization)
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized Pareto distribution parametrized in terms of expected shortfall.
#'
#' The parameter \code{m} corresponds to \eqn{\zeta_u}/(1-\eqn{\alpha}), where \eqn{\zeta_u} is the rate of exceedance over the threshold
#' \code{u} and \eqn{\alpha} is the percentile of the expected shortfall.
#' Note that the actual parametrization is in terms of excess expected shortfall, meaning expected shortfall minus threshold.
#'
#' @details The observed information matrix was calculated from the Hessian using symbolic calculus in Sage.
#'
#' @author Leo Belzile
#' @name gpde
#' @param par vector of length 2 containing \eqn{e_m} and \eqn{\xi}, respectively the expected shortfall at probability 1/(1-\eqn{\alpha}) and the shape parameter.
#' @param dat sample vector
#' @param m number of observations of interest for return levels. See \strong{Details}
#' @param tol numerical tolerance for the exponential model
#' @param method string indicating whether to use the expected  (\code{'exp'}) or the observed (\code{'obs'} - the default) information matrix.
#' @param nobs number of observations
#' @param V vector calculated by \code{gpde.Vfun}
#'
#' @section Usage: \preformatted{gpde.ll(par, dat, m, tol=1e-5)
#' gpde.ll.optim(par, dat, m, tol=1e-5)
#' gpde.score(par, dat, m)
#' gpde.infomat(par, dat, m, method = c('obs', 'exp'), nobs = length(dat))
#' gpde.Vfun(par, dat, m)
#' gpde.phi(par, dat, V, m)
#' gpde.dphi(par, dat, V, m)}
#'
#' @section Functions:
#'
#' \itemize{
#' \item{\code{gpde.ll}:} {log likelihood}
#' \item{\code{gpde.ll.optim}:} {negative log likelihood parametrized in terms of log expected
#' shortfall and shape in order to perform unconstrained optimization}
#' \item{\code{gpde.score}:} {score vector}
#' \item{\code{gpde.infomat}:} {observed information matrix for GPD parametrized in terms of rate of expected shortfall and shape}
#' \item{\code{gpde.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gpde.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gpde.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
NULL

#' @title Generalized Pareto distribution (return level parametrization)
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized Pareto distribution parametrized in terms of return levels.
#'
#' @details The observed information matrix was calculated from the Hessian using symbolic calculus in Sage.
#' @details The interpretation for \code{m} is as follows: if there are on average \eqn{m_y} observations per year above the threshold, then  \eqn{m=Tm_y} corresponds to \eqn{T}-year return level.
#'
#' @author Leo Belzile
#' @name gpdr
#' @param par vector of length 2 containing \eqn{y_m} and \eqn{\xi}, respectively the \eqn{m}-year return level and the shape parameter.
#' @param dat sample vector
#' @param m number of observations of interest for return levels. See \strong{Details}
#' @param tol numerical tolerance for the exponential model
#' @param method string indicating whether to use the expected  (\code{'exp'}) or the observed (\code{'obs'} - the default) information matrix.
#' @param nobs number of observations
#' @param V vector calculated by \code{gpdr.Vfun}
#'
#' @section Usage: \preformatted{gpdr.ll(par, dat, m, tol=1e-5)
#' gpdr.ll.optim(par, dat, m, tol=1e-5)
#' gpdr.score(par, dat, m)
#' gpdr.infomat(par, dat, m, method = c('obs', 'exp'), nobs = length(dat))
#' gpdr.Vfun(par, dat, m)
#' gpdr.phi(par, V, dat, m)
#' gpdr.dphi(par, V, dat, m)}
#' @section Functions:
#'
#' \itemize{
#' \item{\code{gpdr.ll}:} {log likelihood}
#' \item{\code{gpdr.ll.optim}:} {negative log likelihood parametrized in terms of \code{log(scale)} and shape
#' in order to perform unconstrained optimization}
#' \item{\code{gpdr.score}:} {score vector}
#' \item{\code{gpdr.infomat}:} {observed information matrix for GPD parametrized in terms of rate of \eqn{m}-year return level and shape}
#' \item{\code{gpdr.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gpdr.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gpdr.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
NULL

#' @title Generalized extreme value distribution (return level parametrization)
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized extreme value distribution  parametrized in terms of the return level \eqn{z}, scale and shape.
#'
#' @author Leo Belzile
#' @name gevr
#' @param par vector of \code{retlev}, \code{scale} and \code{shape}
#' @param dat sample vector
#' @param p tail probability, corresponding to \eqn{(1-p)}th quantile for \eqn{z}
#' @param method string indicating whether to use the expected  (\code{'exp'}) or the observed (\code{'obs'} - the default) information matrix.
#' @param nobs number of observations
#' @param V vector calculated by \code{gevr.Vfun}
#'
#' @section Usage: \preformatted{gevr.ll(par, dat, p)
#' gevr.ll.optim(par, dat, p)
#' gevr.score(par, dat, p)
#' gevr.infomat(par, dat, p, method = c('obs', 'exp'), nobs = length(dat))
#' gevr.Vfun(par, dat, p)
#' gevr.phi(par, dat, p, V)
#' gevr.dphi(par, dat, p, V)}
#'
#' @section Functions:
#' \itemize{
#' \item{\code{gevr.ll}:} {log likelihood}
#' \item{\code{gevr.ll.optim}:} {negative log likelihood parametrized in terms of return levels, \code{log(scale)} and shape in order to perform unconstrained optimization}
#' \item{\code{gevr.score}:} {score vector}
#' \item{\code{gevr.infomat}:} {observed information matrix}
#' \item{\code{gevr.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gevr.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gevr.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
NULL

#' Negative log likelihood of the generalized Pareto distribution (expected shortfall)
#'
#' @seealso \code{\link{gpde}}
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.ll <- function(par, dat, m) {
    es = par[1]
    xi = as.vector(par[2])
    if (any(xi > 1, es < 0, min(1 + (m^xi - 1 + xi)/((1 - xi) * es) * dat) < 0)) {
        return(-1e+10)
    }
    xizero <- abs(xi) < 1e-9
    sigmae = ifelse(!xizero, es * (1 - xi) * xi/(m^xi - 1 + xi), es/(log(m)+1))
    return(gpd.ll(par = c(sigmae, xi), dat = dat))
}

#' Negative log likelihood of the generalized Pareto distribution (expected shortfall) - optimization
#' The negative log likelihood is parametrized in terms of log expected shortfall and shape in order to perform unconstrained optimization
#' @rdname gpde.ll
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.ll.optim <- function(par, dat, m) {
    es = exp(par[1])
    xi = par[2]
    if (xi > 1) {
        return(1e+10)
    }
    -sum(-log(xi * (1 - xi)/(m^xi - 1 + xi)) - log(es) - (1 + 1/xi) * log1p((m^xi - 1 + xi)/((1 - xi) *
        es) * dat))
}
#' Score vector for the GP distribution (expected shortfall)
#' @seealso \code{\link{gpde}}
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.score <- function(par, dat, m) {
    es = par[1]
    xi = par[2]
    if(missing(m)){
      stop("User must provide `m` parameter in `gpde.score`.")
    }
    xizero <- abs(xi) < 1e-4
    if(!xizero){
    c(sum(-1/es + dat * (m^xi + xi - 1) * (1/xi + 1)/(es^2 * (xi - 1) * (dat * (m^xi + xi - 1)/(es *
        (xi - 1)) - 1))), sum(-((m^xi * log(m) + 1) * dat/(es * (xi - 1)) - dat * (m^xi + xi - 1)/(es *
        (xi - 1)^2)) * (1/xi + 1)/(dat * (m^xi + xi - 1)/(es * (xi - 1)) - 1) + (m^xi + xi - 1) * ((m^xi *
        log(m) + 1) * (xi - 1) * xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi + xi - 1))/((xi -
        1) * xi) + log1p(-dat * (m^xi + xi - 1)/(es * (xi - 1)))/xi^2))
    } else{
      c(sum((dat*log(m) + dat - es)/es^2),
        sum(1/2*((dat^2 - dat*es)*log(m)^3 + (3*dat^2 - 5*dat*es + es^2)*log(m)^2 + dat^2 -
                   4*dat*es + 2*es^2 + (3*dat^2 - 8*dat*es + 2*es^2)*log(m))/(es^2*log(m) + es^2)))
    }
}

#' Observed information matrix for the GP distribution (expected shortfall)
#'
#' The information matrix is parametrized in terms of excess expected shortfall and shape
#' @seealso \code{\link{gpde}}
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.infomat <- function(par, dat, m, method = c("obs", "exp"), nobs = length(dat)) {
    method <- match.arg(method)  #default to observed information
    es = as.vector(par[1])
    xi = as.vector(par[2])
    if(missing(m)){
      stop("User must provide `m` parameter in `gpde.infomat`.")
    }
    if (xi < -0.5) {
        return(matrix(NA, 2, 2))
    }
    xizero <- abs(xi) < 1e-5
    logm <- log(m)
    if (method == "obs") {
      if(!xizero){
        k11 = sum(-1/es^2 + 2 * dat * (m^xi + xi - 1) * (1/xi + 1)/(es^3 * (xi - 1) * (dat * (m^xi +
            xi - 1)/(es * (xi - 1)) - 1)) - dat^2 * (m^xi + xi - 1)^2 * (1/xi + 1)/(es^4 * (xi - 1)^2 *
            (dat * (m^xi + xi - 1)/(es * (xi - 1)) - 1)^2))
        k22 = sum(-((m^xi * logm + 1) * dat/(es * (xi - 1)) - dat * (m^xi + xi - 1)/(es * (xi - 1)^2))^2 *
            (1/xi + 1)/(dat * (m^xi + xi - 1)/(es * (xi - 1)) - 1)^2 + (dat * m^xi * logm^2/(es * (xi -
            1)) - 2 * (m^xi * logm + 1) * dat/(es * (xi - 1)^2) + 2 * dat * (m^xi + xi - 1)/(es * (xi -
            1)^3)) * (1/xi + 1)/(dat * (m^xi + xi - 1)/(es * (xi - 1)) - 1) - (m^xi * (xi - 1) * xi *
            logm^2/(m^xi + xi - 1)^2 - 2 * (m^xi * logm + 1)^2 * (xi - 1) * xi/(m^xi + xi - 1)^3 +
            2 * (m^xi * logm + 1) * (xi - 1)/(m^xi + xi - 1)^2 + 2 * (m^xi * logm + 1) * xi/(m^xi +
            xi - 1)^2 - 2/(m^xi + xi - 1)) * (m^xi + xi - 1)/((xi - 1) * xi) - (m^xi * logm + 1) *
            ((m^xi * logm + 1) * (xi - 1) * xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi +
                xi - 1))/((xi - 1) * xi) + (m^xi + xi - 1) * ((m^xi * logm + 1) * (xi - 1) * xi/(m^xi +
            xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) - xi/(m^xi + xi - 1))/((xi - 1) * xi^2) + (m^xi + xi -
            1) * ((m^xi * logm + 1) * (xi - 1) * xi/(m^xi + xi - 1)^2 - (xi - 1)/(m^xi + xi - 1) -
            xi/(m^xi + xi - 1))/((xi - 1)^2 * xi) - 2 * ((m^xi * logm + 1) * dat/(es * (xi - 1)) -
            dat * (m^xi + xi - 1)/(es * (xi - 1)^2))/(xi^2 * (dat * (m^xi + xi - 1)/(es * (xi - 1)) -
            1)) + 2 * log1p(-dat * (m^xi + xi - 1)/(es * (xi - 1)))/xi^3)
        k12 = sum(-((m^xi * logm + 1) * dat/(es^2 * (xi - 1)) - dat * (m^xi + xi - 1)/(es^2 * (xi -
            1)^2)) * (1/xi + 1)/(dat * (m^xi + xi - 1)/(es * (xi - 1)) - 1) + dat * (m^xi + xi - 1) *
            ((m^xi * logm + 1) * dat/(es * (xi - 1)) - dat * (m^xi + xi - 1)/(es * (xi - 1)^2)) * (1/xi +
            1)/(es^2 * (xi - 1) * (dat * (m^xi + xi - 1)/(es * (xi - 1)) - 1)^2) + dat * (m^xi + xi -
            1)/(es^2 * (xi - 1) * xi^2 * (dat * (m^xi + xi - 1)/(es * (xi - 1)) - 1)))
      } else{
        k11 <- sum(2*dat*es*log(m) + 2*dat*es - es^2)/es^4
        k12 <- 0.5*sum((2*dat*log(m)^2 - es*log(m)^2 + 4*dat*log(m) - 4*es*log(m) + 2*dat - 4*es)*dat)/es^3
        k22 <- sum(1/12*(4*(2*dat^3 - 3*dat^2*es + dat*es^2)*log(m)^5 + (40*dat^3 - 72*dat^2*es + 32*dat*es^2 - es^3)*log(m)^4 +
                           4*(20*dat^3 - 45*dat^2*es + 25*dat*es^2 - es^3)*log(m)^3 + 8*dat^3 - 36*dat^2*es + 48*dat*es^2 -
                           12*es^3 + 4*(20*dat^3 - 57*dat^2*es + 42*dat*es^2 - 3*es^3)*log(m)^2 +
                           8*(5*dat^3 - 18*dat^2*es + 18*dat*es^2 - 3*es^3)*log(m))/(es^3*log(m)^2 + 2*es^3*log(m) + es^3))
      }
        return(cbind(c(k11, k12), c(k12, k22)))
    } else if (method == "exp") {
        sigmae = ifelse(!xizero, es * (1 - xi) * xi/(m^xi - 1 + xi), es/(log(m)+1))
        if(!xizero){
          Jac <- rbind(c((1 - xi) * xi/(m^xi - 1 + xi), (m^xi * logm + 1) * es * (xi - 1) * xi/(m^xi +
            xi - 1)^2 - es * (xi - 1)/(m^xi + xi - 1) - es * xi/(m^xi + xi - 1)), c(0, 1))
        } else{
         Jac <-  rbind(c(1/(logm+1), -1/2*(log(m)^2 + 2*log(m) + 2)*es/(log(m)^2 + 2*log(m) + 1)), c(0, 1))
        }
        return(t(Jac) %*% gpd.infomat(par = c(sigmae, xi), dat = dat, method = "exp", nobs = nobs) %*%
            Jac)
    }
}
#' Tangent exponential model statistics for the generalized Pareto distribution (expected shortfall)
#'
#' Vector implementing conditioning on approximate ancillary statistics for the TEM
#' @seealso \code{\link{gpde}}
#' @name gpde.temstat
#' @inheritParams gpde
#' @keywords internal
#' @export
gpde.Vfun <- function(par, dat, m) {
    es = par[1]
    xi = par[2]
    cbind(dat/es, es * (xi - 1) * xi * (-dat * (m^xi + xi - 1)/(es * (xi - 1)) + 1)^(-1/xi) * (((m^xi *
        log(m) + 1) * dat/(es * (xi - 1)) - dat * (m^xi + xi - 1)/(es * (xi - 1)^2))/(xi * (dat * (m^xi +
        xi - 1)/(es * (xi - 1)) - 1)) - log1p(-dat * (m^xi + xi - 1)/(es * (xi - 1)))/xi^2)/((m^xi +
        xi - 1) * (-dat * (m^xi + xi - 1)/(es * (xi - 1)) + 1)^(-1/xi - 1)))
}

#' Canonical parameter in the local exponential family approximation
#'
#' @inheritParams gpde
#' @rdname gpde.temstat
#' @keywords internal
#' @export
gpde.phi <- function(par, dat, V, m) {
    es = par[1]
    xi = par[2]
    rbind(-(m^xi + xi - 1) * (1/xi + 1)/(es * (xi - 1) * (dat * (m^xi + xi - 1)/(es * (xi - 1)) - 1))) %*%
        V
}

## Derivative matrix of the canonical parameter in the local exponential family approximation
#' @inheritParams gpde
#' @rdname gpde.temstat
#' @keywords internal
#' @export
gpde.dphi <- function(par, dat, V, m) {
    es = par[1]
    xi = par[2]
    rbind((m^xi + xi - 1) * (1/xi + 1)/(es^2 * (xi - 1) * (dat * (m^xi + xi - 1)/(es * (xi - 1)) - 1)) -
        dat * (m^xi + xi - 1)^2 * (1/xi + 1)/(es^3 * (xi - 1)^2 * (dat * (m^xi + xi - 1)/(es * (xi -
            1)) - 1)^2), (m^xi + xi - 1) * ((m^xi * log(m) + 1) * dat/(es * (xi - 1)) - dat * (m^xi +
        xi - 1)/(es * (xi - 1)^2)) * (1/xi + 1)/(es * (xi - 1) * (dat * (m^xi + xi - 1)/(es * (xi - 1)) -
        1)^2) - (m^xi * log(m) + 1) * (1/xi + 1)/(es * (xi - 1) * (dat * (m^xi + xi - 1)/(es * (xi -
        1)) - 1)) + (m^xi + xi - 1) * (1/xi + 1)/(es * (xi - 1)^2 * (dat * (m^xi + xi - 1)/(es * (xi -
        1)) - 1)) + (m^xi + xi - 1)/(es * (xi - 1) * xi^2 * (dat * (m^xi + xi - 1)/(es * (xi - 1)) -
        1))) %*% V
}

#' Negative log likelihood of the generalized Pareto distribution (return levels)
#'
#' @seealso \code{\link{gpdr}}
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.ll <- function(par, dat, m) {
  if(missing(m)){
    stop("Need to specify the reciprocal tail probability `m`")
  }
    ym = par[1]
    xi = par[2]
    if (par[1] < 0) {
        return(-1e+15)
    }
    nn = length(dat)
    pr <- m^xi - 1
    sigmar <- ifelse(abs(xi) > 1e-7, ym * xi/(m^xi - 1), ym/log(m))
    return(gpd.ll(par = c(sigmar, xi), dat = dat))
    # if (abs(xi) > 1e-7) {
    #   -nn * log(xi/pr) - nn * log(ym) - (1 + 1/xi) * sum(log(pmax(1 + exp(log(dat) - log(ym)) * pr,0)))
    # } else {
    #   nn * log(log(m)) - nn * log(ym) + sum(-dat * log(m)/ym)
    # }
}

#' Negative log likelihood parametrized in terms of log return level and shape in order to perform unconstrained optimization
#' @inheritParams gpdr
#' @keywords internal
#' @rdname gpdr.ll
#' @export
gpdr.ll.optim <- function(par, dat, m) {
    tol <- 1e-05
    nn = length(dat)
    ym = exp(par[1])
    xi = par[2]
    p <- m^xi - 1
    if (abs(xi) > tol) {
        -(-nn * log(xi/p) - nn * log(ym) - (1 + 1/xi) * sum(log(pmax(1 + exp(log(dat) - log(ym)) * p,
            0))))
    } else {
        -(nn * log(log(m)) - nn * log(ym) + sum(-dat * log(m)/ym))
    }
}

#' Score of the profile log likelihood for the GP distribution (return levels parametrization)
#' @seealso \code{\link{gpdr}}
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.score <- function(par, dat, m) {
  if(missing(m)){
    stop("Need to specify the reciprocal tail probability `m`")
  }
    nn = length(dat)
    xi = par[2]
    ym = par[1]
    p <- (m)^xi - 1
    xizero <- abs(xi) < 1e-6
    if(!xizero){
    as.vector(c(-nn/ym + (1 + 1/xi) * sum(dat * p/(ym + dat * p))/ym, -nn/xi + exp(log(nn) + xi * log(m) + log(log(m)))/p +
        sum(log1p(dat/ym * p))/xi^2 - (1 + 1/xi) * sum(exp(log(dat) + log(log(m)) + xi * log(m) - log(ym +
        dat * p)))))
    } else{
      as.vector(c(sum((dat*log(m) - ym))/ym^2,
        0.5*sum((dat^2*log(m)^2 + ym^2*log(m) - (dat*log(m)^2 + 2*dat*log(m))*ym))/ym^2))
    }
}

#' Observed information matrix for GP distribution (return levels)
#'
#'The information matrix is parametrized in terms of rate of \code{m}-year return level and shape
#' @seealso \code{\link{gpdr}}
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.infomat <- function(par, dat, m, method = c("obs", "exp"), nobs = length(dat)) {
  if(missing(m)){
    stop("Need to specify the reciprocal tail probability `m`")
  }
    xi = as.vector(par[2])
    r = as.vector(par[1])
    method <- method[1]
    if (xi < -0.5) {
        return(matrix(NA, 2, 2))
    }
    xizero <- abs(xi) < 1e-5
    logm <- log(m)
    if (method == "obs") {
      if(!xizero){
        info <- matrix(ncol = 2, nrow = 2)
        info[1, 1] <- sum(dat^2 * (m^xi - 1)^2 * (1/xi + 1)/((dat * (m^xi - 1)/r + 1)^2 * r^4) - 2 *
            dat * (m^xi - 1) * (1/xi + 1)/((dat * (m^xi - 1)/r + 1) * r^3) + 1/r^2)
        info[2, 2] <- sum(dat^2 * (m^xi)^2 * (1/xi + 1) * logm^2/((dat * (m^xi - 1)/r + 1)^2 * r^2) -
            dat * m^xi * (1/xi + 1) * logm^2/((dat * (m^xi - 1)/r + 1) * r) + m^xi * (m^xi * xi * logm/(m^xi -
            1)^2 - 1/(m^xi - 1)) * logm/xi + (m^xi * xi * logm^2/(m^xi - 1)^2 - 2 * (m^xi)^2 * xi *
            logm^2/(m^xi - 1)^3 + 2 * m^xi * logm/(m^xi - 1)^2) * (m^xi - 1)/xi - (m^xi - 1) * (m^xi *
            xi * logm/(m^xi - 1)^2 - 1/(m^xi - 1))/xi^2 + 2 * dat * m^xi * logm/((dat * (m^xi - 1)/r +
            1) * r * xi^2) - 2 * log1p(dat * (m^xi - 1)/r)/xi^3)
        info[2, 1] <- info[1, 2] <- sum(-dat^2 * (m^xi - 1) * m^xi * (1/xi + 1) * logm/((dat * (m^xi -
            1)/r + 1)^2 * r^3) + dat * m^xi * (1/xi + 1) * logm/((dat * (m^xi - 1)/r + 1) * r^2) -
            dat * (m^xi - 1)/((dat * (m^xi - 1)/r + 1) * r^2 * xi^2))
        return(-info)
      } else{
        k11 <- sum((2*dat*r*log(m) - r^2))/r^4
        k12 <- 0.5*sum((2*dat*log(m) - r*log(m) - 2*r)*dat)*log(m)/r^3
        k22 <- sum((8*dat^3*log(m)^3 - r^3*log(m)^2 + 4*(dat*log(m)^3 + 3*dat*log(m)^2)*r^2 - 12*(dat^2*log(m)^3 + dat^2*log(m)^2)*r)/(12*r^3))
        return(matrix(c(k11, k12, k12, k22), byrow = TRUE, ncol = 2, nrow = 2))
      }
    } else if (method == "exp") {
        sigmar <- ifelse(!xizero, r * xi/(m^xi - 1), r/logm)
        if(!xizero){
          Jac <- rbind(c(xi/(m^xi - 1), -m^xi * r * xi * logm/(m^xi - 1)^2 + r/(m^xi - 1)), c(0, 1))
        } else{
          Jac <-  rbind(c(1/logm, -0.5*r), c(0, 1))
        }
        return(t(Jac) %*% gpd.infomat(par = c(sigmar, xi), method = "exp", dat = dat, nobs = length(dat)) %*% Jac)
    }

}

#' Tangent exponential model statistics for the generalized Pareto distribution (return level)
#'
#' Vector implementing conditioning on approximate ancillary statistics for the TEM
#' @seealso \code{\link{gpdr}}
#' @name gpdr.temstat
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.Vfun <- function(par, dat, m) {
    xi = par[2]
    ym = par[1]
    p <- m^xi - 1
    cbind(dat/ym, -(ym + dat * p)/p * ((m)^xi * log(m)/(p + ym/dat) - log1p(dat/ym * p)/xi))
}

#' Canonical parameter in the local exponential family approximation
#' @seealso \code{\link{gpdr}}
#' @rdname gpdr.temstat
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.phi <- function(par, dat, V, m) {
    xi = par[2]
    r = par[1]
    matrix(-(1 + 1/xi) * (m^xi - 1)/(r + dat * (m^xi - 1)), nrow = 1) %*% V
}

#' Derivative of the canonical parameter \eqn{\phi(\theta)} in the local exponential family approximation
#' @rdname gpdr.temstat
#' @inheritParams gpdr
#' @keywords internal
#' @export
gpdr.dphi <- function(par, dat, V, m) {
    xi = par[2]
    r = par[1]
    p = m^xi - 1
    rbind((1 + 1/xi) * p/(r + dat * p)^2, p/(xi^2 * (r + dat * p)) - (1 + 1/xi) * log(m) * m^xi * r/(r +
        dat * p)^2) %*% V

}
####################################### GEV in terms of return levels ###


#' Return level for the generalized extreme value distribution
#'
#' This function returns the \eqn{1-p}th quantile of the GEV.
#' @inheritParams gev
#' @seealso \code{\link{gev}}
#' @keywords internal
#' @export
gev.retlev <- function(par, p) {
    mu = par[1]
    sigma = par[2]
    xi = as.vector(par[3])
    if (!isTRUE(all.equal(xi, 0))) {
        mu - sigma/xi * (1 - (-log(1 - p))^(-xi))
    } else {
        mu - sigma * log(-log(-p + 1))
    }
}

#' This function returns the mean of N observations from the GEV.
#' @inheritParams gevN
#' @seealso \code{\link{gevN}}
#' @keywords internal
#' @export
gevN.mean <- function(par, N) {
    mu = par[1]
    sigma = par[2]
    xi = as.vector(par[3])
    if (!isTRUE(all.equal(xi, 0))) {
        mu - sigma/xi * (1 - N^xi * gamma(1 - xi))
    } else {
        mu + sigma * (log(N) - psigamma(1))
    }
}

gevN.quant <- function(par, N, q) {
    mu = par[1]
    sigma = par[2]
    xi = as.vector(par[3])
    if (!isTRUE(all.equal(xi, 0))) {
        mu - sigma/xi * (1 - (N/log(1/q))^xi)
    } else {
        mu + sigma * (log(N) - log(log(1/q)))
    }
}

#' This function returns the mean of N maxima from the GP.
#' @inheritParams gpdN
#' @seealso \code{\link{gpdN}}
#' @keywords internal
#' @export
gpdN.mean <- function(par, N) {
    sigma = par[1]
    xi = as.vector(par[2])
    if (!isTRUE(all.equal(xi, 0))) {
        (exp(lgamma(N + 1) + lgamma(1 - xi) - lgamma(N + 1 - xi)) - 1) * sigma/xi
    } else {
        (-psigamma(1) + psigamma(N + 1)) * sigma
    }
}

#' This function returns the qth percentile of N maxima from the GP.
#' @inheritParams gpdN
#' @seealso \code{\link{gpdN}}
#' @keywords internal
#' @export
gpdN.quant <- function(par, q, N) {
    sigma = par[1]
    xi = as.vector(par[2])
    if (!isTRUE(all.equal(xi, 0))) {
        sigma/xi * ((1 - q^(1/N))^(-xi) - 1)
    } else {
        sigma * (-log(-q^(1/N) + 1))
    }
}


#' Negative log likelihood of the generalized extreme value distribution (return levels)
#'
#' @seealso \code{\link{gevr}}
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.ll <- function(par, dat, p) {
    z = par[1]
    sigma = par[2]
    xi = as.vector(par[3])
    muf = ifelse(!isTRUE(all.equal(xi, 0)),
                 z + sigma/xi * (1 - (-log(1 - p))^(-xi)),
                 sigma * log(-log(1 - p)) + z)
    if (!isTRUE(all.equal(xi, 0, tolerance = 1e-08))) {
        sum(-log(sigma) - (1/xi + 1) * log(pmax(1 + xi * (dat - muf)/sigma, 0)) - pmax(1 + xi * (dat -
            muf)/sigma, 0)^(-1/xi))
    } else {
        sum((muf - dat)/sigma - exp((muf - dat)/sigma) - log(sigma))
    }

}

#' Negative log likelihood parametrized in terms of location, log return level and shape in order to perform unconstrained optimization
#' @rdname gevr.ll
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.ll.optim <- function(par, dat, p) {
    tpar = par
    tpar[2] = exp(par[2])
    # if(tpar[3] <= -1){return(1e20)}
    -gevr.ll(tpar, dat, p)
}

#' Score of the log likelihood for the GEV distribution (return levels)
#' @seealso \code{\link{gevr}}
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.score <- function(par, dat, p) {
    z = par[1]
    sigma = par[2]
    xi = par[3]
    log1mp <- log(1 - p)
    if(abs(xi) > 1e-4){
    c(sum((((dat - z) * xi/sigma + 1/(-log1mp)^xi)^(-1/xi - 2) * (xi + 1) * exp(-1/((dat - z) *
        xi/sigma + 1/(-log1mp)^xi)^(1/xi))/sigma^2 - ((dat - z) * xi/sigma + 1/(-log1mp)^xi)^(-2/xi -
        2) * exp(-1/((dat - z) * xi/sigma + 1/(-log1mp)^xi)^(1/xi))/sigma^2) * sigma * ((dat - z) *
        xi/sigma + 1/(-log1mp)^xi)^(1/xi + 1) * exp(1/(((dat - z) * xi/sigma + 1/(-log1mp)^xi)^(1/xi)))),
        sum(-(dat * (-log1mp)^xi - z * (-log1mp)^xi - (dat * (-log1mp)^xi - z * (-log1mp)^xi -
            sigma) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma *
            (-log1mp)^xi))^(1/xi))/((sigma * dat * xi * (-log1mp)^xi - sigma * xi * z * (-log1mp)^xi
                                     + sigma^2) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma *
            (-log1mp)^xi))^(1/xi))), sum(-(xi * z * (-log1mp)^xi - (dat * (-log1mp)^xi -
            sigma * log(-log1mp)) * xi + ((dat * (-log1mp)^xi - sigma * log(-log1mp)) *
            xi^2 + (dat * (-log1mp)^xi - sigma * log(-log1mp)) * xi - (xi^2 * (-log1mp)^xi +
            xi * (-log1mp)^xi) * z) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi +
            sigma)/(sigma * (-log1mp)^xi))^(1/xi) + (dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi
                                                     - (dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma) * ((dat * xi *
            (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi) +
            sigma) * log((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma *
            (-log1mp)^xi)))/((dat * xi^3 * (-log1mp)^xi - xi^3 * z * (-log1mp)^xi + sigma *
            xi^2) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi))))
    } else{
     as.vector(c(sum(
       (1 +  exp((-dat + z)/sigma)*log1mp)/sigma),
       sum((-sigma + dat - z + exp((-dat + z)/sigma)*(dat - z)*log1mp))/sigma^2,
       sum((1/(2*sigma^2))*(( exp(z/sigma)*(-dat + z)*log1mp* (-dat + z + 2*sigma*log(-log1mp)) +  exp(dat/sigma)*((dat - z)*(-2*sigma + dat - z) + 2*sigma*(sigma - dat + z)*log(-log1mp)))/ exp(dat/sigma)))
     ))
    }
}

#' Observed information matrix for GEV distribution (return levels)
#'
#'The information matrix is parametrized in terms of return level (\eqn{(1-p)}th quantile), scale and shape.
#' @seealso \code{\link{gevr}}
#' @inheritParams gevr
#' @param nobs integer number of observations
#' @keywords internal
#' @export
gevr.infomat <- function(par, dat, method = c("obs", "exp"), p, nobs = length(dat)) {
    method <- match.arg(method)
    z = par[1]
    sigma = par[2]
    xi = par[3]
    log1mp <- log(1 - p)
    logmlog1mp <- log(-log1mp)
    if (method == "obs") {
        infomat <- matrix(0, ncol = 3, nrow = 3)
        if(abs(xi) > 1e-4){
        infomat[1, 1] <- sum(-(xi * (-log1mp)^(2 * xi) - (xi^2 * (-log1mp)^(2 * xi) + xi *
            (-log1mp)^(2 * xi)) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi +
            sigma)/(sigma * (-log1mp)^xi))^(1/xi) + (-log1mp)^(2 * xi))/((dat^2 * xi^2 * (-log1mp)^(2 * xi) +
             xi^2 * z^2 * (-log1mp)^(2 * xi) + 2 * sigma * dat * xi * (-log1mp)^xi +
            sigma^2 - 2 * (dat * xi^2 * (-log1mp)^(2 * xi) + sigma * xi * (-log1mp)^xi) * z) *
            ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi)))

        infomat[1, 2] <- infomat[2, 1] <- sum(-(dat * (-log1mp)^(2 * xi) - z * (-log1mp)^(2 *
            xi) - sigma * (-log1mp)^xi + (sigma * xi * (-log1mp)^xi + sigma * (-log1mp)^xi) *
            ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi))/
              ((sigma * dat^2 * xi^2 * (-log1mp)^(2 * xi) + sigma * xi^2 * z^2 *
            (-log1mp)^(2 * xi) + 2 * sigma^2 * dat * xi * (-log1mp)^xi + sigma^3 - 2 * (sigma *
            dat * xi^2 * (-log1mp)^(2 * xi) + sigma^2 * xi * (-log1mp)^xi) * z) * ((dat * xi *
            (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi)))
        infomat[1, 3] <- infomat[3, 1] <- sum(-((sigma * (-log1mp)^xi * logmlog1mp - dat *
            (-log1mp)^(2 * xi)) * xi^2 + (sigma * (-log1mp)^xi * logmlog1mp - dat *
            (-log1mp)^(2 * xi)) * xi + (xi^2 * (-log1mp)^(2 * xi) + xi * (-log1mp)^(2 *
            xi)) * z - (sigma * xi^3 * (-log1mp)^xi * logmlog1mp + xi^2 * z *
            (-log1mp)^(2 * xi) + (sigma * (-log1mp)^xi * (logmlog1mp + 1) - dat * (-log1mp)^(2 *
            xi)) * xi^2) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma *
            (-log1mp)^xi))^(1/xi) + (dat * xi * (-log1mp)^(2 * xi) - xi * z * (-log1mp)^(2 *
            xi) + sigma * (-log1mp)^xi) * log((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/
                                                (sigma * (-log1mp)^xi)))/((dat^2 * xi^4 * (-log1mp)^(2 * xi) +
            xi^4 * z^2 * (-log1mp)^(2 * xi) + 2 * sigma * dat * xi^3 * (-log1mp)^xi + sigma^2 *
            xi^2 - 2 * (dat * xi^4 * (-log1mp)^(2 * xi) + sigma * xi^3 * (-log1mp)^xi) * z) *
            ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi)))
        infomat[2, 2] <- sum((dat^2 * xi * (-log1mp)^(2 * xi) + (xi * (-log1mp)^(2 * xi) -
            (-log1mp)^(2 * xi)) * z^2 - dat^2 * (-log1mp)^(2 * xi) + 2 * sigma * dat * (-log1mp)^xi -
              2 * (dat * xi * (-log1mp)^(2 * xi) - dat * (-log1mp)^(2 * xi) + sigma *
            (-log1mp)^xi) * z - (dat^2 * xi * (-log1mp)^(2 * xi) + xi * z^2 * (-log1mp)^(2 *
            xi) + 2 * sigma * dat * (-log1mp)^xi - sigma^2 - 2 * (dat * xi * (-log1mp)^(2 *
            xi) + sigma * (-log1mp)^xi) * z) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/
                                                  (sigma * (-log1mp)^xi))^(1/xi))/((sigma^2 * dat^2 * xi^2 *
             (-log1mp)^(2 * xi) + sigma^2 * xi^2 * z^2 * (-log1mp)^(2 * xi) + 2 * sigma^3 * dat * xi *
            (-log1mp)^xi + sigma^4 - 2 * (sigma^2 * dat * xi^2 * (-log1mp)^(2 * xi) + sigma^3 *
            xi * (-log1mp)^xi) * z) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi +
            sigma)/(sigma * (-log1mp)^xi))^(1/xi)))
        infomat[3, 2] <- infomat[2, 3] <- sum(-((sigma * dat * (-log1mp)^xi * logmlog1mp -
            dat^2 * (-log1mp)^(2 * xi)) * xi^2 - (xi^2 * (-log1mp)^(2 * xi) + xi *
            (-log1mp)^(2 * xi)) * z^2 + (sigma * dat * (-log1mp)^xi * logmlog1mp - dat^2 *
            (-log1mp)^(2 * xi)) * xi - ((sigma * (-log1mp)^xi * logmlog1mp - 2 * dat * (-log(-p +
            1))^(2 * xi)) * xi^2 + (sigma * (-log1mp)^xi * logmlog1mp - 2 * dat *
            (-log1mp)^(2 * xi)) * xi) * z - (sigma * dat * xi^3 * (-log1mp)^xi * logmlog1mp - xi^2 *
            z^2 * (-log1mp)^(2 * xi) + (sigma * dat * (-log1mp)^xi * (logmlog1mp + 1) -
            dat^2 * (-log1mp)^(2 * xi)) * xi^2 - (sigma * xi^3 * (-log1mp)^xi * logmlog1mp +
            (sigma * (-log1mp)^xi * (logmlog1mp + 1) - 2 * dat * (-log1mp)^(2 *
            xi)) * xi^2) * z) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma *
            (-log1mp)^xi))^(1/xi) + (dat^2 * xi * (-log1mp)^(2 * xi) + xi * z^2 *
             (-log1mp)^(2 * xi) + sigma * dat * (-log1mp)^xi - (2 * dat * xi * (-log1mp)^(2 * xi) +
            sigma * (-log1mp)^xi) * z) * log((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/
                                               (sigma * (-log1mp)^xi)))/((sigma * dat^2 * xi^4 * (-log1mp)^(2 *
            xi) + sigma * xi^4 * z^2 * (-log1mp)^(2 * xi) + 2 * sigma^2 * dat * xi^3 *
              (-log1mp)^xi + sigma^3 * xi^2 - 2 * (sigma * dat * xi^4 * (-log1mp)^(2 * xi) + sigma^2 * xi^3 *
            (-log1mp)^xi) * z) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma *
            (-log1mp)^xi))^(1/xi)))
        infomat[3, 3] <- sum((sigma * dat * xi^4 * (-log1mp)^xi * logmlog1mp^2 + (4 * sigma *
            dat * (-log1mp)^xi * logmlog1mp - 3 * dat^2 * (-log1mp)^(2 * xi)) * xi^3 +
            (2 * sigma * dat * (-log1mp)^xi * (logmlog1mp - 1) - (logmlog1mp^2 - 2 *
                logmlog1mp) * sigma^2 - dat^2 * (-log1mp)^(2 * xi)) * xi^2 - (3 * xi^3 *
            (-log1mp)^(2 * xi) + xi^2 * (-log1mp)^(2 * xi)) * z^2 - (dat^2 * xi^2 *
            (-log1mp)^(2 * xi) + xi^2 * z^2 * (-log1mp)^(2 * xi) + 2 * sigma * dat * xi * (-log1mp)^xi +
            sigma^2 - 2 * (dat * xi^2 * (-log1mp)^(2 * xi) + sigma * xi * (-log1mp)^xi) * z) *
            log((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma *
            (-log1mp)^xi))^2 - (sigma * xi^4 * (-log1mp)^xi * logmlog1mp^2 + 2 * (2 * sigma *
            (-log1mp)^xi * logmlog1mp - 3 * dat * (-log1mp)^(2 * xi)) * xi^3 + 2 * (sigma *
            (-log1mp)^xi * (logmlog1mp - 1) - dat * (-log1mp)^(2 * xi)) * xi^2) * z -
            (sigma * dat * xi^5 * (-log1mp)^xi * logmlog1mp^2 + ((logmlog1mp^2 + 2 *
                logmlog1mp) * sigma * dat * (-log1mp)^xi - dat^2 * (-log1mp)^(2 * xi)) *
                xi^4 + (4 * sigma * dat * (-log1mp)^xi * logmlog1mp - 3 * dat^2 *
                          (-log1mp)^(2 * xi)) * xi^3 - 2 * (sigma * dat * (-log1mp)^xi - sigma^2 * logmlog1mp) *
               xi^2 - (xi^4 * (-log1mp)^(2 * xi) + 3 * xi^3 * (-log1mp)^(2 * xi)) *
                z^2 - (sigma * xi^5 * (-log1mp)^xi * logmlog1mp^2 + ((logmlog1mp^2 +
                2 * logmlog1mp) * sigma * (-log1mp)^xi - 2 * dat * (-log1mp)^(2 * xi)) *
                xi^4 + 2 * (2 * sigma * (-log1mp)^xi * logmlog1mp - 3 * dat * (-log1mp)^(2 * xi)) *
                  xi^3 - 2 * sigma * xi^2 * (-log1mp)^xi) * z) * ((dat * xi * (-log1mp)^xi -
                  xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi) + 2 *
            (dat^2 * xi^3 * (-log1mp)^(2 * xi) - (sigma * dat * (-log1mp)^xi * (log(-log(-p +
                1)) - 2) - dat^2 * (-log1mp)^(2 * xi)) * xi^2 + (xi^3 * (-log1mp)^(2 * xi) +
                xi^2 * (-log1mp)^(2 * xi)) * z^2 + (sigma * dat * (-log1mp)^xi - sigma^2 *
                (logmlog1mp - 1)) * xi - (2 * dat * xi^3 * (-log1mp)^(2 * xi) - (sigma *
                (-log1mp)^xi * (logmlog1mp - 2) - 2 * dat * (-log1mp)^(2 * xi)) * xi^2 +
                sigma * xi * (-log1mp)^xi) * z - (dat^2 * xi^3 * (-log1mp)^(2 * xi) + xi^3 *
                z^2 * (-log1mp)^(2 * xi) + 2 * sigma * dat * xi^2 * (-log1mp)^xi + sigma^2 *
                xi - 2 * (dat * xi^3 * (-log1mp)^(2 * xi) + sigma * xi^2 * (-log1mp)^xi) *
                z) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma *
                (-log1mp)^xi))^(1/xi)) * log((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma *
            (-log1mp)^xi)))/((dat^2 * xi^6 * (-log1mp)^(2 * xi) + xi^6 * z^2 * (-log1mp)^(2 *
            xi) + 2 * sigma * dat * xi^5 * (-log1mp)^xi + sigma^2 * xi^4 - 2 * (dat * xi^6 *
           (-log1mp)^(2 * xi) + sigma * xi^5 * (-log1mp)^xi) * z) * ((dat * xi * (-log1mp)^xi -
            xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi)))
        } else{ #xi numerically zero
          infomat[1, 1] <- sum(( exp((-dat + z)/sigma)*log1mp))/sigma^2
          infomat[1, 2] <- sum(-((sigma + exp((-dat + z)/sigma)*(sigma - dat + z)*log1mp)))/sigma^3
          infomat[2, 1] <- infomat[1, 2]
          infomat[1, 3] <- sum((1/(2*sigma^3))*((2* exp(dat/sigma)* sigma*(sigma - dat + z + sigma*logmlog1mp) +  exp(z/sigma)* log1mp*((dat - z)*(-2*sigma + dat - z) + 2*sigma*(sigma - dat + z)* logmlog1mp))/ exp(dat/sigma)))
          infomat[3, 1] <- infomat[1, 3]
          infomat[2, 2] <- sum(sigma*(sigma - 2*dat + 2*z) + exp((-dat + z)/sigma)*(dat - z)*(-2*sigma + dat - z) * log1mp)/sigma^4
          infomat[2, 3] <- sum((1/(2*sigma^4))*(((dat - z)*(2*exp(dat/sigma)*sigma*(sigma - dat + z + sigma*logmlog1mp) +   exp(z/sigma)*  log1mp*((dat - z)*(-2*sigma + dat - z) +   2*sigma*(sigma - dat + z)* logmlog1mp)))/exp(dat/sigma)))
          infomat[3, 2] <- infomat[2, 3]
          infomat[3, 3] <- sum((1/(12*sigma^4))*((dat - z)*(4*sigma*((dat - z)*(3*sigma - 2*dat + 2*z) - 3*sigma*logmlog1mp*(2*(sigma - dat + z) + sigma*logmlog1mp)) +  exp((-dat + z)/sigma)*log1mp*((-8*sigma + 3*dat - 3*z)*(dat - z)^2 + 12*sigma*logmlog1mp*((dat - z)*(2*sigma - dat + z) - sigma*(sigma - dat + z)*logmlog1mp)))))
        }
        return(-infomat)
        # else expected information matrix
    } else {
        muf = z + sigma/xi * (1 - (-log(1 - p))^(-xi))
        if (!isTRUE(all.equal(xi, 0, tolerance = 1e-07))) {
            Jac <- rbind(c(1, -((-log1mp)^(-xi) - 1)/xi, sigma * (-log1mp)^(-xi) * log(-log1mp)/xi +
                             sigma * ((-log1mp)^(-xi) - 1)/xi^2), c(0, 1, 0), c(0, 0, 1))
        } else {
            Jac <- rbind(c(1, log(-log1mp), -sigma * log(-log1mp)^2 + 1/2 * sigma * (-log1mp)^(-xi)
                           * log(-log1mp)^2), c(0, 1, 0), c(0, 0, 1))
        }
        return(t(Jac) %*% gev.infomat(dat = dat, method = "exp", par = c(muf, sigma, xi)) %*% Jac)

    }
}


#' Tangent exponential model statistics for the GEV distribution (return level)
#'
#' Vector implementing conditioning on approximate ancillary statistics for the TEM
#' @seealso \code{\link{gevr}}
#' @name gevr.temstat
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.Vfun <- function(par, dat, p) {
    z = par[1]
    sigma = par[2]
    xi = par[3]
    log1mp <- log(1-p)
    cbind(1, (dat - z)/sigma, sigma * ((dat - z) * xi/sigma + (-log1mp)^(-xi))^(-1/xi) * (((-log1mp)^(-xi) * log(-log1mp) - (dat - z)/sigma)/(((dat - z) * xi/sigma + (-log1mp)^(-xi)) *
        xi) + log((dat - z) * xi/sigma + (-log1mp)^(-xi))/xi^2)/((dat - z) * xi/sigma + (-log1mp)^(-xi))^(-1/xi - 1))
}

#' Canonical parameter in the local exponential family approximation
#' @rdname gevr.temstat
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.phi <- function(par, dat, p, V) {
    z = par[1]
    sigma = par[2]
    xi = par[3]
    log1mp <- log(1-p)
    t(-((xi * (-log1mp)^xi + (-log1mp)^xi) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi) - (-log1mp)^xi)/((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi))) %*% V
}

#' Derivative of the canonical parameter \eqn{\phi(\theta)} in the local exponential family approximation
#' @rdname gevr.temstat
#' @inheritParams gevr
#' @keywords internal
#' @export
gevr.dphi <- function(par, dat, p, V) {
    z = par[1]
    sigma = par[2]
    xi = par[3]
    log1mp <- log(1 - p)
    rbind((xi * (-log1mp)^(2 * xi) - (xi^2 * (-log1mp)^(2 * xi) + xi * (-log1mp)^(2 *
        xi)) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi) + (-log1mp)^(2 * xi))/((dat^2 * xi^2 * (-log1mp)^(2 * xi) + xi^2 *
        z^2 * (-log1mp)^(2 * xi) + 2 * sigma * dat * xi * (-log1mp)^xi + sigma^2 - 2 * (dat *
        xi^2 * (-log1mp)^(2 * xi) + sigma * xi * (-log1mp)^xi) * z) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi)), (dat * (-log1mp)^(2 * xi) - z * (-log1mp)^(2 * xi) - sigma * (-log1mp)^xi + (sigma * xi * (-log1mp)^xi + sigma * (-log1mp)^xi) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi +
        sigma)/(sigma * (-log1mp)^xi))^(1/xi))/((sigma * dat^2 * xi^2 * (-log1mp)^(2 * xi) +
        sigma * xi^2 * z^2 * (-log1mp)^(2 * xi) + 2 * sigma^2 * dat * xi * (-log1mp)^xi + sigma^3 -
        2 * (sigma * dat * xi^2 * (-log1mp)^(2 * xi) + sigma^2 * xi * (-log1mp)^xi) * z) *
        ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi)),
        ((sigma * (-log1mp)^xi * log(-log1mp) - dat * (-log1mp)^(2 * xi)) * xi^2 + (sigma *
            (-log1mp)^xi * log(-log1mp) - dat * (-log1mp)^(2 * xi)) * xi + (xi^2 * (-log1mp)^(2 * xi) + xi *
        (-log1mp)^(2 * xi)) * z - (sigma * xi^3 * (-log1mp)^xi * log(-log1mp) + xi^2 * z * (-log1mp)^(2 * xi) + (sigma * (-log1mp)^xi * (log(-log1mp) +
            1) - dat * (-log1mp)^(2 * xi)) * xi^2) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi))^(1/xi) + (dat * xi * (-log1mp)^(2 * xi) -
            xi * z * (-log1mp)^(2 * xi) + sigma * (-log1mp)^xi) * log((dat * xi * (-log1mp)^xi - xi * z *
        (-log1mp)^xi + sigma)/(sigma * (-log1mp)^xi)))/((dat^2 * xi^4 *
            (-log1mp)^(2 * xi) + xi^4 * z^2 * (-log1mp)^(2 * xi) + 2 * sigma * dat * xi^3 *
            (-log1mp)^xi + sigma^2 * xi^2 - 2 * (dat * xi^4 * (-log1mp)^(2 * xi) + sigma *
            xi^3 * (-log1mp)^xi) * z) * ((dat * xi * (-log1mp)^xi - xi * z * (-log1mp)^xi +
            sigma)/(sigma * (-log1mp)^xi))^(1/xi))) %*% V
}


#' @title Generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized Pareto distribution parametrized in terms of average maximum of \code{N} exceedances.
#'
#' The parameter \code{N} corresponds to the number of threshold exceedances of interest over which the maxima is taken.
#' \eqn{z} is the corresponding expected value of this block maxima.
#' Note that the actual parametrization is in terms of excess expected mean, meaning expected mean minus threshold.
#'
#' @details The observed information matrix was calculated from the Hessian using symbolic calculus in Sage.
#'
#' @author Leo Belzile
#' @name gpdN
#' @param par vector of length 2 containing \eqn{z} and \eqn{\xi}, respectively the mean excess of the maxima of N exceedances above the threshold and the shape parameter.
#' @param dat sample vector
#' @param N block size for threshold exceedances.
#' @param tol numerical tolerance for the exponential model
#' @param V vector calculated by \code{gpdN.Vfun}
#' @section Usage: \preformatted{gpdN.ll(par, dat, N, tol=1e-5)
#' gpdN.score(par, dat, N)
#' gpdN.infomat(par, dat, N, method = c('obs', 'exp'), nobs = length(dat))
#' gpdN.Vfun(par, dat, N)
#' gpdN.phi(par, dat, N, V)
#' gpdN.dphi(par, dat, N, V)}
#'
#' @section Functions:
#'
#' \itemize{
#' \item{\code{gpdN.ll}:} {log likelihood}
#' \item{\code{gpdN.score}:} {score vector}
#' \item{\code{gpdN.infomat}:} {observed information matrix for GP parametrized in terms of mean of the maximum of \code{N} exceedances and shape}
#' \item{\code{gpdN.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gpdN.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gpdN.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
NULL

#' Negative log likelihood of the generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#'
#' @seealso \code{\link{gpdN}}
#' @inheritParams gpdN
#' @keywords internal
#' @export
gpdN.ll <- function(par, dat, N) {
    xi = par[2]
    z = par[1]
    euler_gamma = 0.57721566490153231044
    sigma = ifelse(abs(xi) > 1e-8,
                   z * xi/(exp(lgamma(N + 1) + lgamma(-xi + 1) - lgamma(N - xi + 1)) - 1),
                   z / (euler_gamma + psigamma(N + 1)))
    gpd.ll(par = c(sigma, xi), dat = dat)
}

#' Score vector of the generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#' @seealso \code{\link{gpdN}}
#' @inheritParams gpdN
#' @keywords internal
#' @export
gpdN.score <- function(par, dat, N) {
    z = par[1]
    xi = par[2]
    cst <- exp(lgamma(N + 1) + lgamma(1 - xi) - lgamma(N + 1 - xi))
    xizero <- abs(xi) < 1e-6
    if(!xizero){
      as.vector(c(sum(dat * (cst - 1) * (1/xi + 1)/(z^2 * (dat * (cst - 1)/z + 1)) - 1/z), sum(-(psigamma(N - xi +
          1) * cst - psigamma(-xi + 1) * cst) * dat * (1/xi + 1)/(z * (dat * (cst - 1)/z + 1)) + ((psigamma(N -
          xi + 1) * cst - psigamma(-xi + 1) * cst) * xi * z/(cst - 1)^2 - z/(cst - 1)) * (cst - 1)/(xi *z) +
            log(dat * (cst - 1)/z + 1)/xi^2)))
    } else{
      euler_gamma <- 0.57721566490153231044
      psi1pN <- psigamma(1 + N)
      psip1pN <- psigamma(1 + N, deriv = 1)
      as.vector(c( sum(dat*euler_gamma - z + dat*psi1pN)/z^2,
       sum((6*dat^2*euler_gamma^3 - 12*dat*euler_gamma^2*z - 6*dat*euler_gamma^3*z -
          dat*euler_gamma*pi^2*z + 6*euler_gamma^2*z^2 + pi^2*z^2 +
          6*(3*dat^2*euler_gamma - dat*(2 + 3*euler_gamma)*z + z^2)*
          psigamma(1 + N)^2 +
          6*dat*(dat - z)*psigamma(1 + N)^3 + 6*(dat*euler_gamma - z)*z*
          psigamma(deriv = 1, 1+N) + psigamma(1 + N)*(18*dat^2*euler_gamma^2 -
          dat*(24*euler_gamma + 18*euler_gamma^2 + pi^2)*z + 12*euler_gamma*z^2 +
          6*dat*z*psigamma(deriv = 1, 1+N)))/(12*z^2*(euler_gamma + psigamma(1 + N))))))
    }
}

#' Information matrix of the generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#' @seealso \code{\link{gpdN}}
#' @inheritParams gpdN
#' @keywords internal
#' @export
gpdN.infomat <- function(par, dat, N, method = c("obs", "exp"), nobs = length(dat)) {
    z = as.vector(par[1])
    xi = as.vector(par[2])
    if (xi < -0.5) {
        return(matrix(NA, 2, 2))
    }
    method <- method[1]  #default corresponds to observed information rather than Fisher information
    xizero <- abs(xi) < 1e-4
    # Fisher information matrix
    cst <- exp(lgamma(N + 1) + lgamma(1 - xi) - lgamma(N + 1 - xi))
    if (method == "exp") {
      if(!xizero){
        sigmaf <- z * xi/(cst - 1)
        Jac <- rbind(c(xi/(cst - 1),
                       -(digamma(N - xi + 1) * cst - digamma(-xi + 1) * cst) * xi * z/(cst - 1)^2 + z/(cst - 1)), c(0, 1))
      } else{
        euler_gamma <- 0.57721566490153231044
        psi1pN <- psigamma(1 + N)
        sigmaf <- z/(euler_gamma + psi1pN)
        Jac <- rbind(c(1/(euler_gamma + psi1pN),
                       -(z*(6*euler_gamma^2 + pi^2 + 12*euler_gamma*psigamma(1 + N) +
                               6*psigamma(1 + N)^2 - 6*psigamma(1 + N, deriv = 1)))/
                           (12*(euler_gamma + psigamma(1 + N))^2)), c(0,1))
      }
        return(t(Jac) %*% gpd.infomat(par = c(sigmaf, xi), dat = dat, method = "exp") %*% Jac)

    } else if (method == "obs") {
        # Observed information
        if(!xizero){
        k11 <- sum(2 * dat * (cst - 1) * (1/xi + 1)/(z^3 * (dat * (cst - 1)/z + 1)) - dat^2 * (cst -
            1)^2 * (1/xi + 1)/(z^4 * (dat * (cst - 1)/z + 1)^2) - 1/z^2)
        k12 <- sum(-(digamma(N - xi + 1) * cst - digamma(-xi + 1) * cst) * dat * (1/xi + 1)/(z^2 * (dat *
            (cst - 1)/z + 1)) + (digamma(N - xi + 1) * cst - digamma(-xi + 1) * cst) * dat^2 * (cst -
            1) * (1/xi + 1)/(z^3 * (dat * (cst - 1)/z + 1)^2) + dat * (cst - 1)/(xi^2 * z^2 * (dat *
            (cst - 1)/z + 1)))
        k22 <- sum(-(digamma(N - xi + 1) * cst - digamma(-xi + 1) * cst)^2 * dat^2 * (1/xi + 1)/(z^2 *
            (dat * (cst - 1)/z + 1)^2) + (digamma(N - xi + 1)^2 * cst - 2 * digamma(N - xi + 1) * digamma(-xi +
            1) * cst + digamma(-xi + 1)^2 * cst - psigamma(deriv = 1, N - xi + 1) * cst + psigamma(deriv = 1,
            -xi + 1) * cst) * dat * (1/xi + 1)/(z * (dat * (cst - 1)/z + 1)) - (digamma(N - xi + 1) *
            cst - digamma(-xi + 1) * cst) * ((digamma(N - xi + 1) * cst - digamma(-xi + 1) * cst) * xi *
            z/(cst - 1)^2 - z/(cst - 1))/(xi * z) + (2 * (digamma(N - xi + 1) * cst - digamma(-xi + 1) *
            cst)^2 * xi * z/(cst - 1)^3 - (digamma(N - xi + 1)^2 * cst - 2 * digamma(N - xi + 1) * digamma(-xi +
            1) * cst + digamma(-xi + 1)^2 * cst - psigamma(deriv = 1, N - xi + 1) * cst + psigamma(deriv = 1,
            -xi + 1) * cst) * xi * z/(cst - 1)^2 - 2 * (digamma(N - xi + 1) * cst - digamma(-xi + 1) *
            cst) * z/(cst - 1)^2) * (cst - 1)/(xi * z) + ((digamma(N - xi + 1) * cst - digamma(-xi +
            1) * cst) * xi * z/(cst - 1)^2 - z/(cst - 1)) * (cst - 1)/(xi^2 * z) - 2 * (digamma(N - xi +
            1) * cst - digamma(-xi + 1) * cst) * dat/(xi^2 * z * (dat * (cst - 1)/z + 1)) + 2 * log(dat *
            (cst - 1)/z + 1)/xi^3)
        } else{
          euler_gamma <- 0.57721566490153231044
          psigamma1pN <- psigamma(1 + N)
          zeta3 <- 1.2020569031595944587
          psigammap1pN <- psigamma(1 + N, deriv = 1)
          psigammapp1pN <-  psigamma(1 + N, deriv = 2)
          k11 <- -sum(-2*dat*euler_gamma + z - 2*dat*psigamma1pN)/z^3
          k12 <- -sum((1/(12*z^3))*(dat*(-12*dat*euler_gamma^2 + (6*euler_gamma*(2 + euler_gamma) + pi^2)*z +
                     12*(-2*dat*euler_gamma + z + euler_gamma*z)*psigamma1pN +
                     6*(-2*dat + z)*psigamma1pN^2 -  6*z*psigammap1pN)))
          k22 <- -sum((-96*dat^3*euler_gamma^5 + 24*dat^2*euler_gamma^3*(6*euler_gamma*(1 + euler_gamma) + pi^2)*z -
                   24*dat*euler_gamma^2*(2*euler_gamma^2*(3 + euler_gamma) + (1 + euler_gamma)*pi^2)*
                   z^2 + (12*euler_gamma^4 + 12*euler_gamma^2*pi^2 - pi^4)*z^3 -
                   12*((40*dat^3*euler_gamma + 4*dat*(3 + 5*euler_gamma)*z^2 - z^3 -
                   12*dat^2*(z + 5*euler_gamma*z))*psigamma1pN^4 +
                   4*dat*(2*dat^2 - 3*dat*z + z^2)*psigamma1pN^5 +
                   2*psigamma1pN^3*(40*dat^3*euler_gamma^2 -  dat^2*(12*euler_gamma*(2 + 5*euler_gamma) + pi^2)*z +
                   dat*(4*euler_gamma*(6 + 5*euler_gamma) + pi^2)*z^2 -2*euler_gamma*z^3 + 6*dat*(dat - z)*z*psigammap1pN) +
                   z*((12*dat^2*euler_gamma^3 - 12*dat*euler_gamma^2*(1 + euler_gamma)*z + (6*euler_gamma^2 - pi^2)*z^2)*
                   psigammap1pN + 3*z^2*psigammap1pN^2 + 4*euler_gamma*(dat*euler_gamma - z)*z*
                   (psigammapp1pN + 2*zeta3)) + psigamma1pN^2*(80*dat^3*euler_gamma^3 -
                    6*dat^2*euler_gamma*(4*euler_gamma*(3 + 5*euler_gamma) + pi^2)*z + 2*dat*(4*euler_gamma^2*(9 + 5*euler_gamma) +
                     (1 + 3*euler_gamma)*pi^2)*z^2 - (6*euler_gamma^2 + pi^2)*z^3 +6*z*(6*dat^2*euler_gamma - 2*dat*(1 + 3*euler_gamma)*z + z^2)*
                     psigammap1pN + 4*dat*z^2*(psigammapp1pN + 2*zeta3)) +2*psigamma1pN*(euler_gamma*(20*dat^3*euler_gamma^3 -
                     3*dat^2*euler_gamma*(2*euler_gamma*(4 + 5*euler_gamma) + pi^2)*z + dat*(2*euler_gamma^2*(12 + 5*euler_gamma) +
                    (2 + 3*euler_gamma)*pi^2)*z^2 - (2*euler_gamma^2 + pi^2)*z^3) +6*euler_gamma*z*(3*dat^2*euler_gamma - dat*(2 + 3*euler_gamma)*z + z^2)*
                    psigammap1pN - 2*z^2*(-2*dat*euler_gamma + z)*(psigammapp1pN + 2*zeta3))))/(144*z^3*(psigamma1pN + euler_gamma)^2))
        }
        return(cbind(c(k11, k12), c(k12, k22)))
    }
}

#' Tangent exponential model statistics for the generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#'
#' Vector implementing conditioning on approximate ancillary statistics for the TEM
#' @seealso \code{\link{gpdN}}
#' @name gpdN.temstat
#' @inheritParams gpdN
#' @keywords internal
#' @export
gpdN.Vfun <- function(par, dat, N) {
    xi = par[2]
    z = par[1]
    cst <- exp(lgamma(N + 1) + lgamma(-xi + 1) - lgamma(N - xi + 1))
    cbind(dat/z, -xi * z * (dat * (cst - 1)/z + 1)^(-1/xi) * ((psigamma(N - xi + 1) * cst - psigamma(-xi +
        1) * cst) * dat/(xi * z * (dat * (cst - 1)/z + 1)) - log(dat * (cst - 1)/z + 1)/xi^2)/((dat *
        (cst - 1)/z + 1)^(-1/xi - 1) * (cst - 1)))
}

#' Canonical parameter in the local exponential family approximation
#'
#' @inheritParams gpdN
#' @rdname gpdN.temstat
#' @keywords internal
#' @export
gpdN.phi <- function(par, dat, N, V) {
    xi = par[2]
    z = par[1]
    cst <- exp(lgamma(N + 1) + lgamma(-xi + 1) - lgamma(N - xi + 1))
    matrix(-(cst - 1) * (1/xi + 1)/(z * (dat * (cst - 1)/z + 1)), nrow = 1) %*% V
}

## Derivative matrix of the canonical parameter in the local exponential family approximation
#' @inheritParams gpdN
#' @rdname gpdN.temstat
#' @keywords internal
#' @export
gpdN.dphi <- function(par, dat, N, V) {
    xi = par[2]
    z = par[1]
    cst <- exp(lgamma(N + 1) + lgamma(-xi + 1) - lgamma(N - xi + 1))
    rbind((cst - 1) * (1/xi + 1)/(z^2 * (dat * (cst - 1)/z + 1)) - dat * (cst - 1)^2 * (1/xi + 1)/(z^3 *
        (dat * (cst - 1)/z + 1)^2), -(psigamma(N - xi + 1) * cst - psigamma(-xi + 1) * cst) * (1/xi +
        1)/(z * (dat * (cst - 1)/z + 1)) + (psigamma(N - xi + 1) * cst - psigamma(-xi + 1) * cst) * dat *
        (cst - 1) * (1/xi + 1)/(z^2 * (dat * (cst - 1)/z + 1)^2) + (cst - 1)/(xi^2 * z * (dat * (cst -
        1)/z + 1))) %*% V
}



############################################################################################

############################################################################################

#' @title Generalized extreme value distribution (quantile/mean of N-block maxima parametrization)
#'
#' @description Likelihood, score function and information matrix, bias,
#' approximate ancillary statistics and sample space derivative
#' for the generalized extreme value distribution  parametrized in terms of the
#' quantiles/mean of N-block maxima parametrization \eqn{z}, scale and shape.
#'
#' @author Leo Belzile
#' @name gevN
#' @param par vector of \code{loc}, quantile/mean of N-block maximum and \code{shape}
#' @param dat sample vector
#' @param V vector calculated by \code{gevN.Vfun}
#' @param q probability, corresponding to \eqn{q}th quantile of the \code{N}-block maximum
#' @param qty string indicating whether to calculate the \code{q} quantile or the mean
#' @section Usage: \preformatted{gevN.ll(par, dat, N, q, qty = c('mean', 'quantile'))
#' gevN.ll.optim(par, dat, N, q = 0.5, qty = c('mean', 'quantile'))
#' gevN.score(par, dat, N, q = 0.5, qty = c('mean', 'quantile'))
#' gevN.infomat(par, dat, qty = c('mean', 'quantile'), method = c('obs', 'exp'), N, q = 0.5, nobs = length(dat))
#' gevN.Vfun(par, dat, N, q = 0.5, qty = c('mean', 'quantile'))
#' gevN.phi(par, dat, N, q = 0.5, qty = c('mean', 'quantile'), V)
#' gevN.dphi(par, dat, N, q = 0.5, qty = c('mean', 'quantile'), V)}
#'
#' @section Functions:
#' \itemize{
#' \item{\code{gevN.ll}:} {log likelihood}
#' \item{\code{gevN.score}:} {score vector}
#' \item{\code{gevN.infomat}:} {expected and observed information matrix}
#' \item{\code{gevN.Vfun}:} {vector implementing conditioning on approximate ancillary statistics for the TEM}
#' \item{\code{gevN.phi}:} {canonical parameter in the local exponential family approximation}
#' \item{\code{gevN.dphi}:} {derivative matrix of the canonical parameter in the local exponential family approximation}
#' }
NULL

#' Negative log likelihood of the generalized extreme value distribution (quantile/mean of N-block maxima parametrization)
#'
#' @seealso \code{\link{gevN}}
#' @inheritParams gevN
#' @keywords internal
#' @export
gevN.ll <- function(par, dat, N, q = 0.5, qty = c("mean", "quantile")) {
    qty <- match.arg(qty)
    mu = par[1]
    z = par[2]
    xi = as.vector(par[3])
    if (!isTRUE(all.equal(xi, 0))) {
        sigma <- switch(qty,
                        quantile = (z - mu) * xi/(N^xi * (log(1/q))^(-xi) - 1),
                        mean = (z - mu) * xi/(N^xi * gamma(1 - xi) - 1))
    } else {
        sigma <- switch(qty,
                        quantile = (z - mu)/(log(N) - log(-log(q))),
                        mean = (z - mu)/(-psigamma(1) + log(N)))
    }
    gev.ll(par = c(mu, sigma, xi), dat = dat)
}

#' Score of the generalized extreme value distribution (quantile/mean of N-block maxima parametrization)
#'
#' @seealso \code{\link{gevN}}
#' @inheritParams gevN
#' @keywords internal
#' @export
gevN.score <- function(par, dat, N, q = 0.5, qty = c("mean", "quantile")) {
    qty <- match.arg(qty)
    mu = par[1]
    z = par[2]
    xi = par[3]
    log1q = log(1/q)
    if (qty == "quantile") {
        # quantiles at prob. q
        if(abs(xi) > 1e-4){
        as.vector(c(sum((-(N^xi * log1q^(-xi) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi - 1) * ((N^xi * log1q^(-xi) -
            1) * (dat - mu)/(mu - z)^2 + (N^xi * log1q^(-xi) - 1)/(mu - z))/xi + ((N^xi * log1q^(-xi) -
            1) * (dat - mu)/(mu - z)^2 + (N^xi * log1q^(-xi) - 1)/(mu - z)) * (1/xi + 1)/((N^xi *
            log1q^(-xi) - 1) * (dat - mu)/(mu - z) - 1) - 1/(mu - z)), sum(-(N^xi * log1q^(-xi) -
            1) * (dat - mu) * (-(N^xi * log1q^(-xi) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu -
            z)^2 * xi) - (N^xi * log1q^(-xi) - 1) * (dat - mu) * (1/xi + 1)/(((N^xi * log1q^(-xi) -
            1) * (dat - mu)/(mu - z) - 1) * (mu - z)^2) + 1/(mu - z)), sum((-(N^xi * log1q^(-xi) -
            1) * (dat - mu)/(mu - z) + 1)^(-1/xi) * ((N^xi * log1q^(-xi) * log(N) - N^xi * log1q^(-xi) *
            log(log1q)) * (dat - mu)/(((N^xi * log1q^(-xi) - 1) * (dat - mu)/(mu - z) - 1) * (mu -
            z) * xi) - log(-(N^xi * log1q^(-xi) - 1) * (dat - mu)/(mu - z) + 1)/xi^2) - (N^xi * log1q^(-xi) *
            log(N) - N^xi * log1q^(-xi) * log(log1q)) * (dat - mu) * (1/xi + 1)/(((N^xi * log1q^(-xi) -
            1) * (dat - mu)/(mu - z) - 1) * (mu - z)) + (N^xi * log1q^(-xi) - 1) * ((N^xi * log1q^(-xi) *
            log(N) - N^xi * log1q^(-xi) * log(log1q)) * (mu - z) * xi/(N^xi * log1q^(-xi) -
            1)^2 - (mu - z)/(N^xi * log1q^(-xi) - 1))/((mu - z) * xi) + log1p(-(N^xi * log1q^(-xi) -
            1) * (dat - mu)/(mu - z))/xi^2)))
    } else{
      logN <- log(N)
      loglog1q <- log(log1q)
      as.vector(c(
        sum((-mu + z + (-1 + exp(((mu - dat)*(-logN + loglog1q))/(mu - z)))*(dat - z)*(logN - loglog1q))/(mu - z)^2),
        sum((mu - z + (mu - dat)*(-1 +  N^((-mu + dat)/(mu - z))*log1q^((mu - dat)/(mu - z)))*(logN -  loglog1q))/
              (mu - z)^2),
        sum((1/(2*(mu - z)^2))*(((-(mu - z))*(mu - 2*dat + z) + (mu - dat)*(dat - z)* (-1 +  N^((-mu + dat)/(mu - z))*
            log1q^((mu - dat)/(mu - z)))*(log( N) - loglog1q))*(logN - loglog1q)))
      ))
      }
    } else {
        # Mean
      if(abs(xi) > 1e-4){
        as.vector(c(sum((-(N^xi * gamma(-xi + 1) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi - 1) * ((N^xi * gamma(-xi +
            1) - 1) * (dat - mu)/(mu - z)^2 + (N^xi * gamma(-xi + 1) - 1)/(mu - z))/xi + ((N^xi * gamma(-xi +
            1) - 1) * (dat - mu)/(mu - z)^2 + (N^xi * gamma(-xi + 1) - 1)/(mu - z)) * (1/xi + 1)/((N^xi *
            gamma(-xi + 1) - 1) * (dat - mu)/(mu - z) - 1) - 1/(mu - z)), sum(-(N^xi * gamma(-xi + 1) -
            1) * (dat - mu) * (-(N^xi * gamma(-xi + 1) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu -
            z)^2 * xi) - (N^xi * gamma(-xi + 1) - 1) * (dat - mu) * (1/xi + 1)/(((N^xi * gamma(-xi +
            1) - 1) * (dat - mu)/(mu - z) - 1) * (mu - z)^2) + 1/(mu - z)), sum((-(N^xi * gamma(-xi +
            1) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi) * ((N^xi * log(N) * gamma(-xi + 1) - N^xi * psigamma(-xi +
            1) * gamma(-xi + 1)) * (dat - mu)/(((N^xi * gamma(-xi + 1) - 1) * (dat - mu)/(mu - z) - 1) *
            (mu - z) * xi) - log1p(-(N^xi * gamma(-xi + 1) - 1) * (dat - mu)/(mu - z))/xi^2) - (N^xi *
            log(N) * gamma(-xi + 1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1)) * (dat - mu) * (1/xi +
            1)/(((N^xi * gamma(-xi + 1) - 1) * (dat - mu)/(mu - z) - 1) * (mu - z)) + (N^xi * gamma(-xi +
            1) - 1) * ((N^xi * log(N) * gamma(-xi + 1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1)) *
            (mu - z) * xi/(N^xi * gamma(-xi + 1) - 1)^2 - (mu - z)/(N^xi * gamma(-xi + 1) - 1))/((mu -
            z) * xi) + log1p(-(N^xi * gamma(-xi + 1) - 1) * (dat - mu)/(mu - z))/xi^2)))
      } else{
        logN <- log(N)
        euler_gamma <- 0.57721566490153231044
        as.vector(c(
         sum((1/(mu - z)^2)*(-mu - euler_gamma*dat + z + euler_gamma*z - dat*log(N) +
                               z*log(N) +
                               exp(((mu - dat)*(euler_gamma + log(N)))/(-mu + z))*(dat -
                                                                                  z)*(euler_gamma + log(N)))),
         sum((1/(mu - z)^2)*(mu - euler_gamma*mu + euler_gamma*dat - z - mu*log(N) +
                               dat*log(N) +
                               exp(((mu - dat)*(euler_gamma + log(N)))/(-mu + z))*(mu -
                                                                                  dat)*(euler_gamma + log(N)))),
         sum((exp((mu*(euler_gamma + log(N)))/(-mu + z))*
                (exp((dat*(euler_gamma + log(N)))/(mu - z))*(mu - dat)*(euler_gamma + log(N))*
                   (pi^2*(mu - z) + 6*euler_gamma^2*(dat - z) + 6*(dat - z)*log(N)*
                      (2*euler_gamma + log(N))) + exp((mu*(euler_gamma + log(N)))/(mu - z))*
                   ((-euler_gamma)*pi^2*(mu - dat)*(mu - z) + pi^2*(mu - z)^2 - 6*euler_gamma^3*(mu - dat)*(dat - z) -
                      6*euler_gamma^2*(mu - z)*(mu - 2*dat + z) + log(N)*((-pi^2)*(mu - dat)*(mu - z) -
                      18*euler_gamma^2*(mu - dat)*(dat - z) - 12*euler_gamma*(mu - z)*(mu - 2*dat + z) +
               6*log(N)*(-mu^2 + 3*euler_gamma*dat*(dat - z) + z*(-2*dat + z) +
               mu*((2 - 3*euler_gamma)*dat + 3*euler_gamma*z) - (mu - dat)*(dat - z)*log(N))))))/
               (12*(mu - z)^2*(euler_gamma + log(N))))
        ))
      }
    }
}

#' Information matrix of the generalized extreme value distribution (quantile/mean of N-block maxima parametrization)
#' @seealso \code{\link{gevN}}
#' @inheritParams gevN
#' @keywords internal
#' @export
gevN.infomat <- function(par, dat, method = c("obs", "exp"), qty = c("mean", "quantile"), N, q = 0.5,  nobs = length(dat)) {
    method <- match.arg(method)
    par <- as.vector(par)
    mu = par[1]
    z = par[2]
    xi = par[3]
    logN = log(N)
    log1q = log(1/q)
    loglog1q = log(log1q)
    qty <- match.arg(qty)
    xizero <- isTRUE(all.equal(xi, 0, tolerance = 1e-05))
    euler_gamma <- 0.57721566490153231044
    zeta3 <- 1.2020569031595944587
    if (qty == "mean") {
        if (!xizero) {
            # z = mu + sigma/xi*(N^xi*gamma(1-xi)-1)
            sigmaq <- (z - mu) * xi/(N^xi * gamma(1 - xi) - 1)
        } else {
            # z = mu + sigma*(log(N) + euler_gamma)
            sigmaq <- (z - mu)/(logN + euler_gamma)
        }

        if (method == "obs") {
          if(!xizero){
            k11 <- sum(((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 2) * ((N^xi *
                gamma(-xi + 1) - 1) * (mu - dat)/(mu - z)^2 - (N^xi * gamma(-xi + 1) - 1)/(mu - z))^2 *
                (1/xi + 1)/xi - 2 * ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi -
                1) * ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z)^3 - (N^xi * gamma(-xi + 1) -
                1)/(mu - z)^2)/xi - ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z)^2 - (N^xi * gamma(-xi +
                1) - 1)/(mu - z))^2 * (1/xi + 1)/((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) +
                1)^2 + 2 * ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z)^3 - (N^xi * gamma(-xi +
                1) - 1)/(mu - z)^2) * (1/xi + 1)/((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) +
                1) - 1/(mu - z)^2)
            k12 <- sum(-(N^xi * gamma(-xi + 1) - 1) * ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu -
                z) + 1)^(-1/xi - 2) * (mu - dat) * ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z)^2 -
                (N^xi * gamma(-xi + 1) - 1)/(mu - z)) * (1/xi + 1)/((mu - z)^2 * xi) + ((N^xi * gamma(-xi +
                1) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 1) * (2 * (N^xi * gamma(-xi + 1) - 1) * (mu -
                dat)/(mu - z)^3 - (N^xi * gamma(-xi + 1) - 1)/(mu - z)^2)/xi - (2 * (N^xi * gamma(-xi +
                1) - 1) * (mu - dat)/(mu - z)^3 - (N^xi * gamma(-xi + 1) - 1)/(mu - z)^2) * (1/xi + 1)/((N^xi *
                gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1) + (N^xi * gamma(-xi + 1) - 1) * (mu -
                dat) * ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z)^2 - (N^xi * gamma(-xi + 1) -
                1)/(mu - z)) * (1/xi + 1)/(((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^2 *
                (mu - z)^2) + 1/(mu - z)^2)
            k13 <- sum(-((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 1) * ((N^xi *
                logN * gamma(-xi + 1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1)) * (mu - dat) * (1/xi +
                1)/(((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1) * (mu - z)) - log1p((N^xi *
                gamma(-xi + 1) - 1) * (mu - dat)/(mu - z))/xi^2) * ((N^xi * gamma(-xi + 1) - 1) *
                (mu - dat)/(mu - z)^2 - (N^xi * gamma(-xi + 1) - 1)/(mu - z))/xi + ((N^xi * gamma(-xi +
                1) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 1) * ((N^xi * logN * gamma(-xi + 1) - N^xi *
                psigamma(-xi + 1) * gamma(-xi + 1)) * (mu - dat)/(mu - z)^2 - (N^xi * logN * gamma(-xi +
                1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1))/(mu - z))/xi - ((N^xi * logN * gamma(-xi +
                1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1)) * (mu - dat)/(mu - z)^2 - (N^xi * logN *
                gamma(-xi + 1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1))/(mu - z)) * (1/xi + 1)/((N^xi *
                gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1) + (N^xi * logN * gamma(-xi + 1) - N^xi *
                psigamma(-xi + 1) * gamma(-xi + 1)) * (mu - dat) * ((N^xi * gamma(-xi + 1) - 1) * (mu -
                dat)/(mu - z)^2 - (N^xi * gamma(-xi + 1) - 1)/(mu - z)) * (1/xi + 1)/(((N^xi * gamma(-xi +
                1) - 1) * (mu - dat)/(mu - z) + 1)^2 * (mu - z)) - ((N^xi * gamma(-xi + 1) - 1) * (mu -
                dat)/(mu - z) + 1)^(-1/xi - 1) * ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z)^2 -
                (N^xi * gamma(-xi + 1) - 1)/(mu - z))/xi^2 + ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu -
                z)^2 - (N^xi * gamma(-xi + 1) - 1)/(mu - z))/(((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu -
                z) + 1) * xi^2))
            k22 <- sum((N^xi * gamma(-xi + 1) - 1)^2 * ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu -
                z) + 1)^(-1/xi - 2) * (mu - dat)^2 * (1/xi + 1)/((mu - z)^4 * xi) - 2 * (N^xi * gamma(-xi +
                1) - 1) * ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 1) * (mu -
                dat)/((mu - z)^3 * xi) - (N^xi * gamma(-xi + 1) - 1)^2 * (mu - dat)^2 * (1/xi + 1)/(((N^xi *
                gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^2 * (mu - z)^4) + 2 * (N^xi * gamma(-xi +
                1) - 1) * (mu - dat) * (1/xi + 1)/(((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) +
                1) * (mu - z)^3) - 1/(mu - z)^2)
            k23 <- sum((N^xi * gamma(-xi + 1) - 1) * ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu -
                z) + 1)^(-1/xi - 1) * (mu - dat) * ((N^xi * logN * gamma(-xi + 1) - N^xi * psigamma(-xi +
                1) * gamma(-xi + 1)) * (mu - dat) * (1/xi + 1)/(((N^xi * gamma(-xi + 1) - 1) * (mu -
                dat)/(mu - z) + 1) * (mu - z)) - log1p((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z))/
                  xi^2)/((mu - z)^2 * xi) - (N^xi * logN * gamma(-xi + 1) - N^xi * psigamma(-xi +
                1) * gamma(-xi + 1)) * ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi -
                1) * (mu - dat)/((mu - z)^2 * xi) - (N^xi * logN * gamma(-xi + 1) - N^xi * psigamma(-xi +
                1) * gamma(-xi + 1)) * (N^xi * gamma(-xi + 1) - 1) * (mu - dat)^2 * (1/xi + 1)/(((N^xi *
                gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^2 * (mu - z)^3) + (N^xi * logN * gamma(-xi +
                1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1)) * (mu - dat) * (1/xi + 1)/(((N^xi * gamma(-xi +
                1) - 1) * (mu - dat)/(mu - z) + 1) * (mu - z)^2) + (N^xi * gamma(-xi + 1) - 1) * ((N^xi *
                gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 1) * (mu - dat)/((mu - z)^2 *
                xi^2) - (N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(((N^xi * gamma(-xi + 1) - 1) * (mu -
                dat)/(mu - z) + 1) * (mu - z)^2 * xi^2))
            k33 <- sum(((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi) * ((N^xi * logN *
                gamma(-xi + 1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1)) * (mu - dat)/(((N^xi * gamma(-xi +
                1) - 1) * (mu - dat)/(mu - z) + 1) * (mu - z) * xi) - log1p((N^xi * gamma(-xi + 1) - 1) *
                (mu - dat)/(mu - z))/xi^2)^2 + ((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) +
                1)^(-1/xi) * ((N^xi * logN * gamma(-xi + 1) - N^xi * psigamma(-xi + 1) * gamma(-xi +
                1))^2 * (mu - dat)^2/(((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^2 * (mu -
                z)^2 * xi) - (N^xi * logN^2 * gamma(-xi + 1) - 2 * N^xi * logN * psigamma(-xi + 1) *
                gamma(-xi + 1) + N^xi * psigamma(-xi + 1)^2 * gamma(-xi + 1) + N^xi * psigamma(1, -xi +
                1) * gamma(-xi + 1)) * (mu - dat)/(((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) +
                1) * (mu - z) * xi) + 2 * (N^xi * logN * gamma(-xi + 1) - N^xi * psigamma(-xi + 1) *
                gamma(-xi + 1)) * (mu - dat)/(((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1) *
                (mu - z) * xi^2) - 2 * log1p((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z))/xi^3) -
                (N^xi * logN * gamma(-xi + 1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1))^2 * (mu -
                  dat)^2 * (1/xi + 1)/(((N^xi * gamma(-xi + 1) - 1) * (mu - dat)/(mu - z) + 1)^2 * (mu -
                  z)^2) + (N^xi * logN^2 * gamma(-xi + 1) - 2 * N^xi * logN * psigamma(-xi + 1) *
                gamma(-xi + 1) + N^xi * psigamma(-xi + 1)^2 * gamma(-xi + 1) + N^xi * psigamma(1, -xi +
                1) * gamma(-xi + 1)) * (mu - dat) * (1/xi + 1)/(((N^xi * gamma(-xi + 1) - 1) * (mu -
                dat)/(mu - z) + 1) * (mu - z)) + (N^xi * gamma(-xi + 1) - 1) * (2 * (N^xi * logN *
                gamma(-xi + 1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1))^2 * (mu - z) * xi/(N^xi *
                gamma(-xi + 1) - 1)^3 - (N^xi * logN^2 * gamma(-xi + 1) - 2 * N^xi * logN * psigamma(-xi +
                1) * gamma(-xi + 1) + N^xi * psigamma(-xi + 1)^2 * gamma(-xi + 1) + N^xi * psigamma(1,
                -xi + 1) * gamma(-xi + 1)) * (mu - z) * xi/(N^xi * gamma(-xi + 1) - 1)^2 - 2 * (N^xi *
                logN * gamma(-xi + 1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1)) * (mu - z)/(N^xi *
                gamma(-xi + 1) - 1)^2)/((mu - z) * xi) - (N^xi * logN * gamma(-xi + 1) - N^xi * psigamma(-xi +
                1) * gamma(-xi + 1)) * ((N^xi * logN * gamma(-xi + 1) - N^xi * psigamma(-xi + 1) *
                gamma(-xi + 1)) * (mu - z) * xi/(N^xi * gamma(-xi + 1) - 1)^2 - (mu - z)/(N^xi * gamma(-xi +
                1) - 1))/((mu - z) * xi) + (N^xi * gamma(-xi + 1) - 1) * ((N^xi * logN * gamma(-xi +
                1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1)) * (mu - z) * xi/(N^xi * gamma(-xi + 1) -
                1)^2 - (mu - z)/(N^xi * gamma(-xi + 1) - 1))/((mu - z) * xi^2) - 2 * (N^xi * logN *
                gamma(-xi + 1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1)) * (mu - dat)/(((N^xi * gamma(-xi +
                1) - 1) * (mu - dat)/(mu - z) + 1) * (mu - z) * xi^2) + 2 * log1p((N^xi * gamma(-xi + 1) -
                1) * (mu - dat)/(mu - z))/xi^3)
          } else{
            k11 <- -sum((1/(mu - z)^4)*((-exp(((mu - dat)*(euler_gamma + logN))/(-mu +  z)))*(dat - z)*(euler_gamma + logN)*
                                     (2*mu + euler_gamma*dat - (2 + euler_gamma)*z + (dat - z)*logN) +
                                     (mu - z)*(mu + 2*euler_gamma*dat - z - 2*euler_gamma*z +  2*(dat - z)*logN)))
            k22 <- -sum(-((1/(mu - z)^4)*(exp(((mu - dat)*(euler_gamma + logN))/(-mu + z))*(mu - dat)*(euler_gamma + logN)*
                                        ((-2 + euler_gamma)*mu - euler_gamma*dat + 2*z + (mu - dat)*logN) +
                                        (mu - z)*((-1 + 2*euler_gamma)*mu - 2*euler_gamma*dat + z + 2*(mu - dat)*logN))))
            k12 <- -sum((1/(mu - z)^4)*((-exp(((mu - dat)*(euler_gamma + logN))/(-mu + z)))*(euler_gamma + logN)*
                                         (mu^2 + (-2 + euler_gamma)*mu*dat - euler_gamma*dat^2 - euler_gamma*mu*z + (2 + euler_gamma)*dat*z - z^2 +
                                            (mu - dat)*(dat - z)*logN) + (mu - z)*((-1 + euler_gamma)*mu - 2*euler_gamma*dat + z + euler_gamma*z +
                                        (mu - 2*dat + z)*logN)))
            k13 <- -sum((1/(12*(mu - z)^4))*(exp(((2*mu - dat)*(euler_gamma + logN))/(-mu + z))*(dat - z)*
                    (exp(((2*mu - dat)*(euler_gamma + logN))/(mu - z))*
                       (mu - z)*(12*euler_gamma*(-mu + z) + pi^2*(-mu + z) + 6*euler_gamma^2*(mu - 2*dat + z)) +
                       exp((mu*(euler_gamma + logN))/(mu - z))*
                       ((-euler_gamma)*pi^2*(mu - dat)*(mu - z) + pi^2*(mu - z)^2 -
                          6*euler_gamma^3*(mu - dat)*(dat - z) -
                          6*euler_gamma^2*(mu - z)*(mu - 2*dat + z)) + logN*(12*exp(((2*mu - dat)*
                       (euler_gamma + logN))/(mu - z))*(mu - z)*((-1 + euler_gamma)*mu + z + euler_gamma*(-2*dat + z)) +
                       exp((mu*(euler_gamma + logN))/(mu - z))*((-pi^2)*(mu - dat)*(mu - z) -
                       18*euler_gamma^2*(mu - dat)*(dat - z) - 12*euler_gamma*(mu - z)*(mu - 2*dat + z)) +
                       6*logN*(exp(((2*mu - dat)*(euler_gamma + logN))/(mu -z))*(mu - z)*(mu - 2*dat + z) +
                       exp((mu*(euler_gamma + logN))/(mu - z))*(-mu^2 + 3*euler_gamma*dat*(dat - z) + z*(-2*dat + z) +
                        mu*((2 - 3*euler_gamma)*dat + 3*euler_gamma*z) - (mu - dat)*(dat - z)*logN))))))
           k23 <- -sum((1/(12*(mu - z)^4))*(exp(((2*mu - dat)*(euler_gamma + logN))/(-mu +  z))*(mu - dat)*
              (exp(((2*mu - dat)*(euler_gamma + logN))/(mu - z))*(mu - z)*(12*euler_gamma*(-mu + z) +
                  pi^2*(-mu + z) + 6*euler_gamma^2*(mu - 2*dat + z)) + exp((mu*(euler_gamma + logN))/(mu - z))*
                 ((-euler_gamma)*pi^2*(mu - dat)*(mu - z) + pi^2*(mu - z)^2 -   6*euler_gamma^3*(mu - dat)*(dat - z) -  6*euler_gamma^2*(mu - z)*(mu - 2*dat + z)) +
                 logN*(12*exp(((2*mu - dat)*(euler_gamma + logN))/(mu - z))*(mu - z)*((-1 + euler_gamma)*mu + z + euler_gamma*(-2*dat + z)) +
                           exp((mu*(euler_gamma + logN))/(mu - z))* ((-pi^2)*(mu - dat)*(mu - z) -
                        18*euler_gamma^2*(mu - dat)*(dat - z) - 12*euler_gamma*(mu - z)* (mu - 2*dat + z)) +
                           6*logN*(exp(((2*mu - dat)*(euler_gamma + logN))/(mu - z))*(mu - z)*(mu - 2*dat + z) +  exp((mu*(euler_gamma + logN))/(mu - z))*
                                       (-mu^2 + 3*euler_gamma*dat*(dat - z) +   z*(-2*dat + z) + mu*((2 - 3*euler_gamma)*dat + 3*euler_gamma*z) - (mu - dat)*(dat - z)*logN))))))
           k33 <- -sum((exp(((mu - dat)*(euler_gamma + logN))/(-mu + z))* (24*((-(mu - dat))*(dat - z)*(2*euler_gamma*pi^2*(mu - dat)*(mu - z) -
                      pi^2*(mu - z)^2 + 30*euler_gamma^3*(mu - dat)*(dat - z) + 20*euler_gamma^2*(mu - z)*(mu - 2*dat + z)) +  exp(((mu - dat)*(euler_gamma + logN))/
                      (mu - z))*(mu -  z)*((-pi^2)*(mu - dat)*(mu - z)*(dat - z) + 20*euler_gamma^2*(mu - dat)*(dat - z)*(mu - 2*dat + z) +
                      2*euler_gamma*(mu - z)*  (mu^2 - 12*mu*dat + 12*dat^2 + 10*mu*z - 12*dat*z + z^2)))*logN^3 +  12*((-(mu - dat))*(dat - z)*(pi^2*(mu - dat)*(mu - z) +
                      45*euler_gamma^2*(mu - dat)*(dat - z) + 20*euler_gamma*(mu - z)*(mu - 2*dat + z)) +
                        exp(((mu - dat)*(euler_gamma + logN))/(mu - z))*(mu - z)* (mu^3 + 40*euler_gamma*dat^3 - 12*(1 + 5*euler_gamma)*dat^2*z + 4*(3 + 5*euler_gamma)*dat*z^2 -
                        z^3 + mu^2*(4*(-3 + 5*euler_gamma)*dat + (9 - 20*euler_gamma)*z) + mu*((12 - 60*euler_gamma)*dat^2 +
                         80*euler_gamma*dat*z - (9 + 20*euler_gamma)*z^2)))* logN^4 +  24*(mu - dat)*(dat - z)*(-2*mu^2 + 4*mu*dat - 9*euler_gamma*mu*dat +
                      9*euler_gamma*dat^2 + 9*euler_gamma*mu*z - 4*dat*z -9*euler_gamma*dat*z + 2*z^2 +
                        2*exp(((mu - dat)*(euler_gamma + logN))/(mu - z))*(mu - z)* (mu - 2*dat + z))*logN^5 - 36*(mu - dat)^2*(dat - z)^2*logN^6 -
                        euler_gamma^2*(mu - dat)*(mu^3*pi^4 + 48*euler_gamma^3*mu^2*(dat - z) + 12*euler_gamma^2*mu^2*pi^2*(dat - z) +
                        36*euler_gamma^4*mu*(dat - z)^2 - 36*euler_gamma^4*dat*(dat - z)^2 +48*euler_gamma*mu*pi^2*(dat - z)*z + 12*euler_gamma^2*pi^2*dat*(dat - z)*z +
                          48*euler_gamma^3*(dat - z)*(2*dat - z)*z - pi^4*dat*z^2 + 24*euler_gamma*mu^2*pi^2*(-dat + z) +
                        96*euler_gamma^3*mu*dat*(-dat + z) + 24*euler_gamma*pi^2*z^2*(-dat + z) - 12*euler_gamma^2*mu*pi^2*(dat - z)*(dat + z) + mu*pi^4*z*(2*dat + z) -
                          mu^2*pi^4*(dat + 2*z) - 96*(mu - z)^3*zeta3) + exp(((mu - dat)*(euler_gamma + logN))/(mu - z))*(mu - z)* ((-pi^4)*(mu - z)^3 -
                        24*euler_gamma^3*pi^2*(mu - dat)*(mu - z)*(dat - z) + 48*euler_gamma^5*(mu - dat)*(dat - z)*(mu - 2*dat + z) +
                        12*euler_gamma^4*(mu - z)* (mu^2 + 12*dat^2 - 12*dat*z + z^2 + 2*mu*(-6*dat + 5*z)) + 96*euler_gamma*(mu - z)^3*zeta3 -
                          12*euler_gamma^2*(mu - z)^2*(pi^2*(mu - 2*dat + z) + 8*(mu - dat)*zeta3)) +  logN^2*((-(mu - dat))*(mu^3*pi^4 +
                      480*euler_gamma^3*mu^2*(dat - z) + 72*euler_gamma^2*mu^2*pi^2*  (dat - z) + 540*euler_gamma^4*mu*(dat - z)^2 -
                        540*euler_gamma^4*dat*(dat - z)^2 + 144*euler_gamma*mu*pi^2*(dat - z)*z +  72*euler_gamma^2*pi^2*dat*(dat - z)*z + 480*euler_gamma^3*(dat - z)*(2*dat - z)*z -
                        pi^4*dat*z^2 +72*euler_gamma*mu^2*pi^2*(-dat + z) + 960*euler_gamma^3*mu*dat*(-dat + z) +  72*euler_gamma*pi^2*z^2*(-dat + z) -
                        72*euler_gamma^2*mu*pi^2*(dat - z)*(dat + z) +  mu*pi^4*z*(2*dat + z) - mu^2*pi^4*(dat + 2*z) - 96*(mu - z)^3*zeta3) +
                        12*exp(((mu - dat)*(euler_gamma + logN))/(mu - z))*(mu - z)* (-6*euler_gamma*pi^2*(mu - dat)*(mu - z)*(dat - z) +                                                                       40*euler_gamma^3*(mu - dat)*(dat - z)*  (mu - 2*dat + z) +  6*euler_gamma^2*(mu - z)*(mu^2 + 12*dat^2 - 12*dat*z + z^2 +
                        2*mu*(-6*dat + 5*z)) - (mu - z)^2*(pi^2*(mu - 2*dat + z) + 8*(mu - dat)*zeta3))) +
                        2* logN*((-euler_gamma)*(mu - dat)*(mu^3*pi^4 +  120*euler_gamma^3*mu^2*(dat - z) + 24*euler_gamma^2*mu^2*pi^2*(dat - z) +
                      108*euler_gamma^4*mu*(dat - z)^2 - 108*euler_gamma^4*dat*(dat - z)^2 +  72*euler_gamma*mu*pi^2*(dat - z)*z +  24*euler_gamma^2*pi^2*dat*(dat - z)*  z +
                        120*euler_gamma^3*(dat - z)*(2*dat - z)*z -  pi^4*dat*z^2 + 36*euler_gamma*mu^2*pi^2*(-dat + z) + 240*euler_gamma^3*mu*dat*(-dat + z) +
                        36*euler_gamma*pi^2*z^2*(-dat + z) - 24*euler_gamma^2*mu*pi^2*(dat - z)*(dat + z) +  mu*pi^4*z*(2*dat + z) -
                        mu^2*pi^4*(dat + 2*z) - 96*(mu - z)^3*zeta3) + 12*exp(((mu - dat)*(euler_gamma + logN))/(mu - z))*(mu - z)*
                        (-3*euler_gamma^2*pi^2*(mu - dat)*(mu - z)*(dat - z) +
                        10*euler_gamma^4*(mu - dat)*(dat - z)*  (mu - 2*dat + z) +  2*euler_gamma^3*(mu - z)*(mu^2 + 12*dat^2 - 12*dat*z + z^2 + 2*mu*(-6*dat + 5*z)) +
                          4*(mu - z)^3*zeta3 - euler_gamma*(mu - z)^2*(pi^2*(mu - 2*dat + z) + 8*(mu - dat)*zeta3)))))/  (144*(mu - z)^4*(euler_gamma + logN)^2))
          }
            return(cbind(c(k11, k12, k13), c(k12, k22, k23), c(k13, k23, k33)))
        } else if (method == "exp") {
            if (!xizero) {
                Jac <- rbind(c(1, 0, 0), c(-xi/(N^xi * gamma(-xi + 1) - 1), xi/(N^xi * gamma(-xi + 1) -
                  1), (N^xi * logN * gamma(-xi + 1) - N^xi * psigamma(-xi + 1) * gamma(-xi + 1)) *
                  (mu - z) * xi/(N^xi * gamma(-xi + 1) - 1)^2 - (mu - z)/(N^xi * gamma(-xi + 1) - 1)),
                  c(0, 0, 1))
                return(t(Jac) %*% gev.infomat(par = c(mu, sigmaq, xi), method = "exp", nobs = nobs) %*%
                  Jac)
            } else {
                Jac <- rbind(c(1, 0, 0), c(-1/(euler_gamma + logN), 1/(euler_gamma + logN), 1/12 *
                  (6 * euler_gamma^2 + pi^2 + 12 * euler_gamma * logN + 6 * logN^2) * (mu - z)/(euler_gamma +
                  logN)^2), c(0, 0, 1))

                return(t(Jac) %*% gev.infomat(par = c(mu, sigmaq, 0), method = "exp", nobs = nobs) %*%
                  Jac)

            }
        }
    } else if (qty == "quantile") {
        if (!xizero) {
            # z = mu + sigma/xi*(N^xi*log(1/q)^(-xi)-1)
            sigmaq <- (z - mu) * xi/(N^xi * (log1q)^(-xi) - 1)
        } else {
            # z = mu + sigma*log(N/log1q)
            sigmaq <- (z - mu)/log(N/log1q)
        }
        # Quantiles, observed information
        if (method == "obs") {
          if(!xizero){
            k11 <- sum(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 2) * ((N^xi *
                log1q^(-xi) - 1) * (mu - dat)/(mu - z)^2 - (N^xi * log1q^(-xi) - 1)/(mu - z))^2 *
                (1/xi + 1)/xi - 2 * ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi -
                1) * ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z)^3 - (N^xi * log1q^(-xi) -
                1)/(mu - z)^2)/xi - ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z)^2 - (N^xi * log1q^(-xi) -
                1)/(mu - z))^2 * (1/xi + 1)/((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1)^2 +
                2 * ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z)^3 - (N^xi * log1q^(-xi) - 1)/(mu -
                  z)^2) * (1/xi + 1)/((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1) - 1/(mu -
                z)^2)
            k12 <- sum(-(N^xi * log1q^(-xi) - 1) * ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu -
                z) + 1)^(-1/xi - 2) * (mu - dat) * ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z)^2 -
                (N^xi * log1q^(-xi) - 1)/(mu - z)) * (1/xi + 1)/((mu - z)^2 * xi) + ((N^xi * log1q^(-xi) -
                1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 1) * (2 * (N^xi * log1q^(-xi) - 1) * (mu -
                dat)/(mu - z)^3 - (N^xi * log1q^(-xi) - 1)/(mu - z)^2)/xi - (2 * (N^xi * log1q^(-xi) -
                1) * (mu - dat)/(mu - z)^3 - (N^xi * log1q^(-xi) - 1)/(mu - z)^2) * (1/xi + 1)/((N^xi *
                log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1) + (N^xi * log1q^(-xi) - 1) * (mu -
                dat) * ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z)^2 - (N^xi * log1q^(-xi) -
                1)/(mu - z)) * (1/xi + 1)/(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1)^2 *
                (mu - z)^2) + 1/(mu - z)^2)
            k13 <- sum(-((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 1) * ((N^xi *
                log1q^(-xi) * logN - N^xi * log1q^(-xi) * loglog1q) * (mu - dat) * (1/xi +
                1)/(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1) * (mu - z)) - log1p((N^xi *
                log1q^(-xi) - 1) * (mu - dat)/(mu - z))/xi^2) * ((N^xi * log1q^(-xi) - 1) *
                (mu - dat)/(mu - z)^2 - (N^xi * log1q^(-xi) - 1)/(mu - z))/xi + ((N^xi * log1q^(-xi) -
                1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 1) * ((N^xi * log1q^(-xi) * logN - N^xi *
                log1q^(-xi) * loglog1q) * (mu - dat)/(mu - z)^2 - (N^xi * log1q^(-xi) * logN -
                N^xi * log1q^(-xi) * loglog1q)/(mu - z))/xi - ((N^xi * log1q^(-xi) * logN -
                N^xi * log1q^(-xi) * loglog1q) * (mu - dat)/(mu - z)^2 - (N^xi * log1q^(-xi) *
                logN - N^xi * log1q^(-xi) * loglog1q)/(mu - z)) * (1/xi + 1)/((N^xi * log1q^(-xi) -
                1) * (mu - dat)/(mu - z) + 1) + (N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) *
                loglog1q) * (mu - dat) * ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z)^2 -
                (N^xi * log1q^(-xi) - 1)/(mu - z)) * (1/xi + 1)/(((N^xi * log1q^(-xi) - 1) * (mu -
                dat)/(mu - z) + 1)^2 * (mu - z)) - ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) +
                1)^(-1/xi - 1) * ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z)^2 - (N^xi * log1q^(-xi) -
                1)/(mu - z))/xi^2 + ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z)^2 - (N^xi * log1q^(-xi) -
                1)/(mu - z))/(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1) * xi^2))
            k22 <- sum((N^xi * log1q^(-xi) - 1)^2 * ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu -
                z) + 1)^(-1/xi - 2) * (mu - dat)^2 * (1/xi + 1)/((mu - z)^4 * xi) - 2 * (N^xi * log1q^(-xi) -
                1) * ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 1) * (mu - dat)/((mu -
                z)^3 * xi) - (N^xi * log1q^(-xi) - 1)^2 * (mu - dat)^2 * (1/xi + 1)/(((N^xi * log1q^(-xi) -
                1) * (mu - dat)/(mu - z) + 1)^2 * (mu - z)^4) + 2 * (N^xi * log1q^(-xi) - 1) * (mu -
                dat) * (1/xi + 1)/(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1) * (mu - z)^3) -
                1/(mu - z)^2)
            k23 <- sum((N^xi * log1q^(-xi) - 1) * ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu -
                z) + 1)^(-1/xi - 1) * (mu - dat) * ((N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) *
                loglog1q) * (mu - dat) * (1/xi + 1)/(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu -
                z) + 1) * (mu - z)) - log1p((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) )/xi^2)/((mu -
                z)^2 * xi) - (N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) * loglog1q) *
                ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi - 1) * (mu - dat)/((mu -
                z)^2 * xi) - (N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) * loglog1q) *
                (N^xi * log1q^(-xi) - 1) * (mu - dat)^2 * (1/xi + 1)/(((N^xi * log1q^(-xi) - 1) *
                (mu - dat)/(mu - z) + 1)^2 * (mu - z)^3) + (N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) *
                loglog1q) * (mu - dat) * (1/xi + 1)/(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu -
                z) + 1) * (mu - z)^2) + (N^xi * log1q^(-xi) - 1) * ((N^xi * log1q^(-xi) - 1) *
                (mu - dat)/(mu - z) + 1)^(-1/xi - 1) * (mu - dat)/((mu - z)^2 * xi^2) - (N^xi * log1q^(-xi) -
                1) * (mu - dat)/(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1) * (mu - z)^2 *
                xi^2))
            k33 <- sum(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi) * ((N^xi * log1q^(-xi) *
                logN - N^xi * log1q^(-xi) * loglog1q) * (mu - dat)/(((N^xi * log1q^(-xi) -
                1) * (mu - dat)/(mu - z) + 1) * (mu - z) * xi) - log1p((N^xi * log1q^(-xi) - 1) * (mu -
                dat)/(mu - z))/xi^2)^2 + ((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1)^(-1/xi) *
                ((N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) * loglog1q)^2 * (mu - dat)^2/(((N^xi *
                  log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1)^2 * (mu - z)^2 * xi) - (N^xi * log1q^(-xi) *
                  logN^2 - 2 * N^xi * log1q^(-xi) * logN * loglog1q + N^xi * log1q^(-xi) *
                  loglog1q^2) * (mu - dat)/(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) +
                  1) * (mu - z) * xi) + 2 * (N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) *
                  loglog1q) * (mu - dat)/(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1) *
                  (mu - z) * xi^2) - 2 * log1p((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z))/xi^3) -
                (N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) * loglog1q)^2 * (mu - dat)^2 *
                  (1/xi + 1)/(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1)^2 * (mu - z)^2) +
                (N^xi * log1q^(-xi) * logN^2 - 2 * N^xi * log1q^(-xi) * logN * loglog1q +
                  N^xi * log1q^(-xi) * loglog1q^2) * (mu - dat) * (1/xi + 1)/(((N^xi * log1q^(-xi) -
                  1) * (mu - dat)/(mu - z) + 1) * (mu - z)) + (N^xi * log1q^(-xi) - 1) * (2 * (N^xi *
                log1q^(-xi) * logN - N^xi * log1q^(-xi) * loglog1q)^2 * (mu - z) * xi/(N^xi *
                log1q^(-xi) - 1)^3 - (N^xi * log1q^(-xi) * logN^2 - 2 * N^xi * log1q^(-xi) *
                logN * loglog1q + N^xi * log1q^(-xi) * loglog1q^2) * (mu - z) * xi/(N^xi *
                log1q^(-xi) - 1)^2 - 2 * (N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) *
                loglog1q) * (mu - z)/(N^xi * log1q^(-xi) - 1)^2)/((mu - z) * xi) - (N^xi * log1q^(-xi) *
                logN - N^xi * log1q^(-xi) * loglog1q) * ((N^xi * log1q^(-xi) * logN -
                N^xi * log1q^(-xi) * loglog1q) * (mu - z) * xi/(N^xi * log1q^(-xi) - 1)^2 -
                (mu - z)/(N^xi * log1q^(-xi) - 1))/((mu - z) * xi) + (N^xi * log1q^(-xi) - 1) *
                ((N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) * loglog1q) * (mu - z) *
                  xi/(N^xi * log1q^(-xi) - 1)^2 - (mu - z)/(N^xi * log1q^(-xi) - 1))/((mu - z) *
                xi^2) - 2 * (N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) * loglog1q) *
                (mu - dat)/(((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z) + 1) * (mu - z) * xi^2) +
                2 * log1p((N^xi * log1q^(-xi) - 1) * (mu - dat)/(mu - z))/xi^3)
          } else{
           k11 <- -sum((1/(mu - z)^4)*((mu - z)^2 + (-dat + z)*(-2*mu + 2*z + N^((-mu + dat)/(mu - z))*
                  log1q^((mu - dat)/(mu - z))*(2*(mu - z) + (dat - z)*(logN - loglog1q)))*
                  (logN - loglog1q)))
           k12 <- -sum(-((1/(mu - z)^4)*((mu - z)*(mu -z - (mu - 2*dat + z)*(logN - loglog1q)) +  N^((-mu + dat)/(mu - z))*
                       log1q^((mu - dat)/(mu - z))*((mu - z)*(mu - 2*dat + z) + (mu - dat)*(dat - z)*(logN - loglog1q))*(logN -loglog1q))))
           k13 <- -sum((1/(2*(mu - z)^4))*(2*(mu - z)^2*(-mu + z - (dat - z)*(logN - loglog1q))*(logN - loglog1q) +
                     2*(mu - dat)*(mu - z)*(mu - z + (dat - z)*(logN - loglog1q))*(logN - loglog1q) +
                     2*(mu - dat)*(mu - z)*(dat - z)*(logN - loglog1q)^2 + (mu - z)*(mu - 2*dat + z)*(-dat + z)*(logN - loglog1q)^2 +
                     exp(((mu - dat)*(-logN + loglog1q))/(mu - z))*(mu - z)*(mu - 2*dat + z)*(-dat + z)*
                     (logN - loglog1q)^2 - N^((-mu + dat)/(mu - z))*(mu - dat)*(dat - z)^2*log1q^((mu - dat)/(mu - z))*
                     (logN - loglog1q)^3))
           k22 <- -sum((1/(mu - z)^4)*((mu - z)^2 + (mu - dat)*(-2*mu + 2*z + N^((-mu + dat)/(mu - z))*
                      log1q^((mu - dat)/(mu - z))*(2*(mu - z) - (mu - dat)*(logN - loglog1q)))*(logN - loglog1q)))
           k23 <- -sum((1/(2*(mu - z)^4))*((mu - dat)*((mu - z)*(-2*mu +  2*z + (mu - 2*dat + z)*(logN - loglog1q)) +
                      N^((-mu + dat)/(mu - z))*log1q^((mu - dat)/(mu - z))*((-(mu - z))*(mu - 2*dat + z) - (mu - dat)*(dat - z)*(logN - loglog1q))*(logN -
                      loglog1q))*(logN - loglog1q)))
           k33 <- -sum((1/12)*((logN - loglog1q)^2 +  (12*(mu - dat)*(dat - z)*(-mu +
                   z + (mu - 2*dat + z)*(logN - loglog1q))*(logN -  loglog1q)^2)/
                  (mu - z)^3 - (4*exp(((mu - dat)*(-logN + loglog1q))/(mu - z))*(mu - dat)*(dat - z)*(mu - 2*dat + z)*
                                  (logN - loglog1q)^3)/(mu - z)^3 - (3*    N^((-mu + dat)/(mu - z))*(mu - dat)^2*(dat - z)^2*
                   log1q^((mu - dat)/(mu - z))*(logN - loglog1q)^4)/(mu - z)^4 +
                                 (8*(mu - dat)*(dat - z)*(mu - 2*dat + z)*(-logN + loglog1q)^3)/(mu - z)^3))

          }
            return(cbind(c(k11, k12, k13), c(k12, k22, k23), c(k13, k23, k33)))

        } else if (method == "exp") {
            if (!xizero) {
                Jac <- rbind(c(1, 0, 0), c(-xi/(N^xi * log1q^(-xi) - 1), xi/(N^xi * log1q^(-xi) -
                  1), (N^xi * log1q^(-xi) * logN - N^xi * log1q^(-xi) * loglog1q) * (mu -
                  z) * xi/(N^xi * log1q^(-xi) - 1)^2 - (mu - z)/(N^xi * log1q^(-xi) - 1)), c(0,
                  0, 1))
                return(t(Jac) %*% gev.infomat(par = c(mu, sigmaq, xi), method = "exp", nobs = nobs) %*%
                  Jac)
            } else {
                Jac <- rbind(c(1, 0, 0), c(-1/(logN - log(-log(q))), 1/(logN - log(-log(q))), 1/2 *
                  (mu - z)), c(0, 0, 1))
                return(t(Jac) %*% gev.infomat(par = c(mu, sigmaq, 0), method = "exp", nobs = nobs) %*%
                  Jac)
            }
        }
    }

}



#' Tangent exponential model statistics for the generalized Pareto distribution (mean of maximum of N exceedances parametrization)
#'
#' Vector implementing conditioning on approximate ancillary statistics for the TEM
#' @seealso \code{\link{gevN}}
#' @name gevN.temstat
#' @inheritParams gevN
#' @keywords internal
#' @export
gevN.Vfun <- function(par, dat, N, q = 0.5, qty = c("mean", "quantile")) {
    # Quantiles, profiling by substituting sigma by zq
    qty <- match.arg(qty)
    mu = par[1]
    z = par[2]
    xi = par[3]
    log1q = log(1/q)
    if (qty == "quantile") {
        cbind((mu - z) * ((N^xi * log1q^(-xi) - 1) * (dat - mu)/(mu - z)^2 + (N^xi * log1q^(-xi) -
            1)/(mu - z))/(N^xi * log1q^(-xi) - 1), -(dat - mu)/(mu - z), (-(N^xi * log1q^(-xi) -
            1) * (dat - mu)/(mu - z) + 1)^(-1/xi) * (mu - z) * xi * ((N^xi * log1q^(-xi) * log(N) -
            N^xi * log1q^(-xi) * log(log1q)) * (dat - mu)/(((N^xi * log1q^(-xi) - 1) * (dat -
            mu)/(mu - z) - 1) * (mu - z) * xi) - log(-(N^xi * log1q^(-xi) - 1) * (dat - mu)/(mu -
            z) + 1)/xi^2)/((N^xi * log1q^(-xi) - 1) * (-(N^xi * log1q^(-xi) - 1) * (dat - mu)/(mu -
            z) + 1)^(-1/xi - 1)))
    } else {
        cbind((mu - z) * ((N^xi * gamma(-xi + 1) - 1) * (dat - mu)/(mu - z)^2 + (N^xi * gamma(-xi + 1) -
            1)/(mu - z))/(N^xi * gamma(-xi + 1) - 1), -(dat - mu)/(mu - z), (-(N^xi * gamma(-xi + 1) -
            1) * (dat - mu)/(mu - z) + 1)^(-1/xi) * (mu - z) * xi * ((N^xi * log(N) * gamma(-xi + 1) -
            N^xi * psigamma(-xi + 1) * gamma(-xi + 1)) * (dat - mu)/(((N^xi * gamma(-xi + 1) - 1) * (dat -
            mu)/(mu - z) - 1) * (mu - z) * xi) - log(-(N^xi * gamma(-xi + 1) - 1) * (dat - mu)/(mu -
            z) + 1)/xi^2)/((N^xi * gamma(-xi + 1) - 1) * (-(N^xi * gamma(-xi + 1) - 1) * (dat - mu)/(mu -
            z) + 1)^(-1/xi - 1)))
    }
}

#' Canonical parameter in the local exponential family approximation
#'
#' @inheritParams gevN
#' @rdname gevN.temstat
#' @keywords internal
#' @export
gevN.phi <- function(par, dat, N, q = 0.5, qty = c("mean", "quantile"), V) {
    qty <- match.arg(qty)
    mu = par[1]
    z = par[2]
    xi = par[3]
    if (qty == "mean") {
        matrix(-(N^xi * gamma(-xi + 1) - 1) * (-(N^xi * gamma(-xi + 1) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi -
            1)/((mu - z) * xi) - (N^xi * gamma(-xi + 1) - 1) * (1/xi + 1)/(((N^xi * gamma(-xi + 1) -
            1) * (dat - mu)/(mu - z) - 1) * (mu - z)), nrow = 1) %*% V
    } else {
        matrix(-(N^xi * log(1/q)^(-xi) - 1) * (-(N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi -
            1)/((mu - z) * xi) - (N^xi * log(1/q)^(-xi) - 1) * (1/xi + 1)/(((N^xi * log(1/q)^(-xi) -
            1) * (dat - mu)/(mu - z) - 1) * (mu - z)), nrow = 1) %*% V
    }
}

## Derivative matrix of the canonical parameter in the local exponential family approximation
#' @inheritParams gevN
#' @rdname gevN.temstat
#' @keywords internal
#' @export
gevN.dphi <- function(par, dat, N, q = 0.5, qty = c("mean", "quantile"), V) {
    qty <- match.arg(qty)
    mu = par[1]
    z = par[2]
    xi = par[3]
    if (qty == "mean") {
        rbind(-(N^xi * mu * xi^2 * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat -
            z)/(mu - z))^(1/xi) * gamma(-xi + 1) - N^xi * xi^2 * z * ((N^xi * mu * gamma(-xi + 1) - (N^xi *
            gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi) * gamma(-xi + 1) + N^xi * mu * xi * ((N^xi *
            mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi) * gamma(-xi +
            1) - N^xi * xi * z * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu -
            z))^(1/xi) * gamma(-xi + 1) - N^xi * mu * xi * gamma(-xi + 1) + N^xi * xi * z * gamma(-xi +
            1) - N^xi * dat * gamma(-xi + 1) + N^xi * z * gamma(-xi + 1) + dat - z) * (N^xi * gamma(-xi +
            1) - 1)/((N^xi * dat * gamma(-xi + 1) - N^xi * mu * gamma(-xi + 1) - dat + z)^2 * (mu - z) *
            xi^2 * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi)),
            (mu * xi^2 * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu -
                z))^(1/xi) - xi^2 * z * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) *
                dat - z)/(mu - z))^(1/xi) + mu * xi * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi +
                1) - 1) * dat - z)/(mu - z))^(1/xi) - xi * z * ((N^xi * mu * gamma(-xi + 1) - (N^xi *
                gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi) - N^xi * dat * gamma(-xi + 1) + N^xi *
                mu * gamma(-xi + 1) - mu * xi + xi * z + dat - mu) * (N^xi * gamma(-xi + 1) - 1)/((N^xi *
                dat * gamma(-xi + 1) - N^xi * mu * gamma(-xi + 1) - dat + z)^2 * (mu - z) * xi^2 * ((N^xi *
                mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi)), (N^xi *
                mu * xi^3 * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu -
                z))^(1/xi) * log(N) * gamma(-xi + 1) - N^xi * xi^3 * z * ((N^xi * mu * gamma(-xi + 1) -
                (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi) * log(N) * gamma(-xi + 1) - N^xi *
                mu * xi^3 * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu -
                z))^(1/xi) * psigamma(-xi + 1) * gamma(-xi + 1) + N^xi * xi^3 * z * ((N^xi * mu * gamma(-xi +
                1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi) * psigamma(-xi + 1) * gamma(-xi +
                1) + N^xi * mu * xi^2 * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) *
                dat - z)/(mu - z))^(1/xi) * log(N) * gamma(-xi + 1) - N^xi * xi^2 * z * ((N^xi * mu *
                gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi) * log(N) * gamma(-xi +
                1) - N^xi * mu * xi^2 * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) *
                dat - z)/(mu - z))^(1/xi) * psigamma(-xi + 1) * gamma(-xi + 1) + N^xi * xi^2 * z * ((N^xi *
                mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi) * psigamma(-xi +
                1) * gamma(-xi + 1) - N^xi * mu * xi^2 * log(N) * gamma(-xi + 1) + N^xi * xi^2 * z *
                log(N) * gamma(-xi + 1) + N^xi * mu * xi^2 * psigamma(-xi + 1) * gamma(-xi + 1) - N^xi *
                xi^2 * z * psigamma(-xi + 1) * gamma(-xi + 1) + N^(2 * xi) * dat * xi * ((N^xi * mu *
                gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi) * gamma(-xi +
                1)^2 - N^(2 * xi) * mu * xi * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) -
                1) * dat - z)/(mu - z))^(1/xi) * gamma(-xi + 1)^2 - N^(2 * xi) * dat * xi * log(N) *
                gamma(-xi + 1)^2 + N^(2 * xi) * mu * xi * log(N) * gamma(-xi + 1)^2 + N^(2 * xi) * dat *
                xi * psigamma(-xi + 1) * gamma(-xi + 1)^2 - N^(2 * xi) * mu * xi * psigamma(-xi + 1) *
                gamma(-xi + 1)^2 - 2 * N^xi * dat * xi * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi +
                1) - 1) * dat - z)/(mu - z))^(1/xi) * gamma(-xi + 1) + N^xi * mu * xi * ((N^xi * mu *
                gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi) * gamma(-xi +
                1) + N^xi * xi * z * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat -
                z)/(mu - z))^(1/xi) * gamma(-xi + 1) + N^xi * dat * xi * log(N) * gamma(-xi + 1) - N^xi *
                mu * xi * log(N) * gamma(-xi + 1) - N^xi * dat * xi * psigamma(-xi + 1) * gamma(-xi +
                1) + N^xi * mu * xi * psigamma(-xi + 1) * gamma(-xi + 1) - N^(2 * xi) * dat * xi * gamma(-xi +
                1)^2 + N^(2 * xi) * mu * xi * gamma(-xi + 1)^2 + N^(2 * xi) * dat * log((N^xi * mu *
                gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z)) * gamma(-xi + 1)^2 -
                N^(2 * xi) * mu * log((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat -
                  z)/(mu - z)) * gamma(-xi + 1)^2 + 2 * N^xi * dat * xi * gamma(-xi + 1) - N^xi * mu *
                xi * gamma(-xi + 1) - N^xi * xi * z * gamma(-xi + 1) - 2 * N^xi * dat * log((N^xi * mu *
                gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z)) * gamma(-xi + 1) +
                N^xi * mu * log((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu -
                  z)) * gamma(-xi + 1) + N^xi * z * log((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi +
                1) - 1) * dat - z)/(mu - z)) * gamma(-xi + 1) + dat * xi * ((N^xi * mu * gamma(-xi +
                1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi) - xi * z * ((N^xi * mu *
                gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi) - dat * xi +
                xi * z + dat * log((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat -
                z)/(mu - z)) - z * log((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat -
                z)/(mu - z)))/((N^xi * dat * gamma(-xi + 1) - N^xi * mu * gamma(-xi + 1) - dat + z)^2 *
                xi^3 * ((N^xi * mu * gamma(-xi + 1) - (N^xi * gamma(-xi + 1) - 1) * dat - z)/(mu - z))^(1/xi))) %*%
            V
    } else {
        rbind((N^xi * log(1/q)^(-xi) - 1) * (-(N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi -
            2) * ((N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu - z)^2 + (N^xi * log(1/q)^(-xi) - 1)/(mu -
            z)) * (1/xi + 1)/((mu - z) * xi) - (N^xi * log(1/q)^(-xi) - 1) * ((N^xi * log(1/q)^(-xi) -
            1) * (dat - mu)/(mu - z)^2 + (N^xi * log(1/q)^(-xi) - 1)/(mu - z)) * (1/xi + 1)/(((N^xi *
            log(1/q)^(-xi) - 1) * (dat - mu)/(mu - z) - 1)^2 * (mu - z)) + (N^xi * log(1/q)^(-xi) - 1) *
            (-(N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)^2 * xi) +
            (N^xi * log(1/q)^(-xi) - 1) * (1/xi + 1)/(((N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu -
                z) - 1) * (mu - z)^2), -(N^xi * log(1/q)^(-xi) - 1)^2 * (dat - mu) * (-(N^xi * log(1/q)^(-xi) -
            1) * (dat - mu)/(mu - z) + 1)^(-1/xi - 2) * (1/xi + 1)/((mu - z)^3 * xi) - (N^xi * log(1/q)^(-xi) -
            1) * (-(N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu - z)^2 * xi) +
            (N^xi * log(1/q)^(-xi) - 1)^2 * (dat - mu) * (1/xi + 1)/(((N^xi * log(1/q)^(-xi) - 1) * (dat -
                mu)/(mu - z) - 1)^2 * (mu - z)^3) - (N^xi * log(1/q)^(-xi) - 1) * (1/xi + 1)/(((N^xi *
            log(1/q)^(-xi) - 1) * (dat - mu)/(mu - z) - 1) * (mu - z)^2), (N^xi * log(1/q)^(-xi) - 1) *
            (-(N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi - 1) * ((N^xi * log(1/q)^(-xi) *
            log(N) - N^xi * log(1/q)^(-xi) * log(log(1/q))) * (dat - mu) * (1/xi + 1)/(((N^xi * log(1/q)^(-xi) -
            1) * (dat - mu)/(mu - z) - 1) * (mu - z)) - log(-(N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu -
            z) + 1)/xi^2)/((mu - z) * xi) - (N^xi * log(1/q)^(-xi) * log(N) - N^xi * log(1/q)^(-xi) *
            log(log(1/q))) * (-(N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu -
            z) * xi) + (N^xi * log(1/q)^(-xi) * log(N) - N^xi * log(1/q)^(-xi) * log(log(1/q))) * (N^xi *
            log(1/q)^(-xi) - 1) * (dat - mu) * (1/xi + 1)/(((N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu -
            z) - 1)^2 * (mu - z)^2) - (N^xi * log(1/q)^(-xi) * log(N) - N^xi * log(1/q)^(-xi) * log(log(1/q))) *
            (1/xi + 1)/(((N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu - z) - 1) * (mu - z)) + (N^xi *
            log(1/q)^(-xi) - 1) * (-(N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu - z) + 1)^(-1/xi - 1)/((mu -
            z) * xi^2) + (N^xi * log(1/q)^(-xi) - 1)/(((N^xi * log(1/q)^(-xi) - 1) * (dat - mu)/(mu -
            z) - 1) * (mu - z) * xi^2)) %*% V
    }
}

#' Simulate r-largest observations from point process of extremes
#'
#' Simulate the \code{r}-largest observations from a Poisson point process with intensity
#' \deqn{\Lambda(x) = (1+\xi(x-\mu)/\sigma)^{-1/\xi}}.
#'
#' @param n sample size
#' @param r number of observations per block
#' @param loc location parameter
#' @param scale scale parameter
#' @param shape shape parameter
#' @return an \code{n} by \code{r} matrix of samples from the point process, ordered from largest to smallest in each row.
#' @export
rrlarg <- function(n, r, loc, scale, shape){
  U <- matrix(rexp(n*r), nrow = n, ncol = r)
  U <- t(apply(U, 1, cumsum))
  if(r == 1){
   U <- matrix(U, ncol = 1, nrow = n)
  }
  if(!isTRUE(all.equal(shape, 0))){
   loc + scale*(U^(-shape)-1)/shape
  } else{
   loc - log(U)*scale
  }
}

#' Distribution of the r-largest observations
#'
#' Likelihood, score function and information matrix for the r-largest observations likelihood.
#'
#' @author Leo Belzile
#' @name rlarg
#' @param par vector of \code{loc}, \code{scale} and \code{shape}
#' @param dat an \code{n} by \code{r} sample matrix, ordered from largest to smallest in each row
#' @param method string indicating whether to use the expected  (\code{'exp'}) or the observed (\code{'obs'} - the default) information matrix.
#' @param nobs number of observations for the expected information matrix. Default to \code{length(dat)} if \code{dat} is provided.
#' @param r number of order statistics kept. Default to \code{ncol(dat)}
#' @section Usage:
#' \preformatted{
#' rlarg.ll(par, dat, u, np)
#' rlarg.score(par, dat)
#' rlarg.infomat(par, dat, method = c('obs', 'exp'), nobs = nrow(dat), r = ncol(dat))}
#'
#' @section Functions:
#'
#' \itemize{
#' \item{\code{rlarg.ll}:} {log likelihood}
#' \item{\code{rlarg.score}:} {score vector}
#' \item{\code{rlarg.infomat}:} {observed or expected information matrix}
#' }
#' @references Coles, S. (2001). \emph{An Introduction to Statistical Modeling of Extreme Values}, Springer, 209 p.
#' @references Smith, R.L. (1986).  Extreme value theory based on the r largest annual events, \emph{Journal of Hydrology}, \bold{86}(1-2), 27--43, \code{http://dx.doi.org/10.1016/0022-1694(86)90004-1}.
NULL


#' Log-likelihood of the point process of r-largest observations
#'
#' @param par a vector of location, scale and shape parameter
#' @param dat an \code{n} by \code{r} matrix of observations, ordered from largest to smallest in each row.
#' @export
#' @keywords internal
rlarg.ll <- function(par, dat){
 #Maximum is in first column, ordered
  dat <- as.matrix(dat)
  if(!isTRUE(all(dat[1,1] >= dat[1,]))){
    stop("Observations in `dat` must be ordered from largest to smallest in each row.")
  }
  r <- ncol(dat)
  n <- nrow(dat)
  if(abs(par[3]) > 1e-7){
  - sum((1+par[3]*(dat[,r]-par[1])/par[2])^(-1/par[3])) - n*r*log(par[2]) -
      (1/par[3]+1)*sum(log1p(par[3]*(dat-par[1])/par[2]))
  } else{
    -sum(exp((par[1]-dat[,r])/par[2])) - n*r*log(par[2]) - sum((dat-par[1])/par[2])
  }
}
#' Score of the r-largest observations
#'
#' The score is computed via linear interpolation for the shape parameter in a neighborhood of zero
#' @inheritParams rlarg.ll
#' @return a vector of size 3
#' @export
#' @keywords internal
rlarg.score <- function(par, dat){
  par <- as.vector(par)
  mu <- par[1]; sigma <- par[2]; xi <- par[3]
  dat <- as.matrix(dat)
  if(!isTRUE(all(dat[1,1] >= dat[1,]))){
    stop("Observations in `dat` must be ordered from largest to smallest in each row.")
  }
  r <- ncol(dat)
  n <- nrow(dat)
  # xi \neq 0
  if(abs(xi) > 1e-7){
  score <- cbind(
    -(-(mu - dat[,r])*xi/sigma + 1)^(-1/xi - 1)/sigma -
      rowSums(xi*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1))),
    -r/sigma + (mu - dat[,r])*(-(mu - dat[,r])*xi/sigma + 1)^(-1/xi - 1)/sigma^2 +
      rowSums((mu - dat)*xi*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1))),
    -(-(mu - dat[,r])*xi/sigma + 1)^(-1/xi)*(log1p(-(mu - dat[,r])*xi/sigma)/xi^2 -
                                                (mu - dat[,r])/(sigma*((mu - dat[,r])*xi/sigma - 1)*xi)) +
      rowSums( - (mu - dat)*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) +
                 log1p(-(mu - dat)*xi/sigma)/xi^2))
 return(colSums(score))
 } else{
   #Linearly interpolate for values between 1e-7*c(-1,1)
   yxi <- sapply(xivals <- c(-1e-6,-1e-7,1e-7,1e-6),  function(xi){
       sum(-(-(mu - dat[,r])*xi/sigma + 1)^(-1/xi)*(log1p(-(mu - dat[,r])*xi/sigma)/xi^2 -
                                                      (mu - dat[,r])/(sigma*((mu - dat[,r])*xi/sigma - 1)*xi))) +
         sum( - (mu - dat)*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1)) +
                log1p(-(mu - dat)*xi/sigma)/xi^2)})
  scorexi <- approx(x = xivals, y = yxi, xout = xi)$y
 score <- c(sum(-exp(-dat[,r]/sigma + mu/sigma)/sigma) + length(dat)/sigma,
            -n*r/sigma + sum((mu-dat[,r])/sigma^2*exp((mu-dat[,r])/sigma)) + sum((dat-mu)/sigma^2),
            scorexi)
   return(score)
  }

}


#' Information matrix for the r-largest observations.
#'
#' The function returns the expected or observed information matrix.
#' @seealso \code{\link{rlarg}}
#' @inheritParams rlarg
#' @rdname rlarg
#' @export
#' @keywords internal
rlarg.infomat <- function(par, dat, method = c("obs", "exp"), nobs = nrow(dat), r = ncol(dat)){
  eta <- par[1]
  tau <- par[2]
  xi <- par[3]
  r <- as.integer(r)
  xizero <- abs(xi) < 1e-5
  if(method == "exp"){
    m11a <- -(((r + xi*(2 + xi))*gamma(r + 2*xi))/(tau^2*gamma(r)))
    if(!xizero){
      m12a <- (-gamma(1 + r + xi) + (r + xi*(2 + xi))*gamma(r + 2*xi))/(tau^2*xi*gamma(r))
      m13a <- ((-r)*xi*gamma(r + xi) - (r + xi*(2 + xi))*gamma(r + 2*xi) +
                 gamma(1 + r + xi)*(1 + xi + xi*digamma(1 + r + xi)))/(tau*xi^2*gamma(r))
      m22a <- -((gamma(1 + r) - 2*gamma(1 + r + xi) + (r + xi*(2 + xi))*gamma(r + 2*xi))/(tau^2* xi^2*gamma(r)))
      m23a <- (1/(tau*xi^3*gamma(r)))*((r + xi*(2 + xi))*gamma(r + 2*xi) +
                                         gamma(1 + r)*(1 + xi*digamma(1 + r)) - gamma(r + xi)*(2*r + xi*(2 + xi) +
                                                                                                 xi*(r + xi)*digamma(1 + r + xi)))


    } else{
      m12a <- (1 + r*digamma(r))/tau^2
      m13a <- -((digamma(r)*(2 + r*digamma(r)) + r*psigamma(deriv = 1, r))/(2*tau))
      m22a <- -((1 + digamma(r)*(2 + r*digamma(r)) + r*psigamma(deriv = 1, r))/tau^2)
      m23a <-  (1/(2*tau))*(3*psigamma(deriv = 1, r) +  digamma(r)*(2 + digamma(r)*(3 + r*digamma(r)) +
                                                                      3*r*psigamma(deriv = 1, r)) + r*psigamma(deriv = 2, r))
    }
    if(abs(xi) > 1e-3){
      m33a <- -((1/(xi^4*gamma(r)))*((r + xi*(2 + xi))*gamma(r + 2*xi) -
                                       2*gamma(r + xi)*(r + xi + xi^2 + xi*(r + xi)*digamma(1 + r + xi)) +
                                       gamma(r)*(r + xi*(2 + xi) + r*xi*(2*digamma(r) + xi*(digamma(1 + r)^2 + psigamma(deriv = 1, 1 + r))))))
    } else{
      m33a <- -(4*psigamma(deriv=0, r)^3 + r*psigamma(deriv=0, r)^4 +
                  psigamma(deriv=1, r)*(4 + 3*r*psigamma(deriv=1, r)) + psigamma(deriv=0, r)^2*
                  (4 + 6*r*psigamma(deriv=1, r)) + 4*psigamma(deriv=2, r) +
                  4*psigamma(deriv=0, r)*(3*psigamma(deriv=1, r) + r*psigamma(deriv=2, r)) +
                  r*psigamma(deriv=3, r))/4
    }
    ma <- matrix(c(m11a, m12a, m13a, m12a, m22a, m23a, m13a, m23a, m33a), nrow = 3, ncol = 3)
    if(r == 1L){
      return(nobs*ma)
    } else if(r > 1L){
      if(!xizero){
        m11b <- -((xi^2*gamma(r + 2*xi))/(tau^2*(1 + 2*xi)*gamma(r)))
        m12b <- (xi*gamma(r + 2*xi))/(tau^2*(1 + 2*xi)*gamma(r))
        m13b <- (gamma(r + xi)/(1 + xi) - gamma(r + 2*xi)/(1 + 2*xi))/(tau*gamma(r))
        m22b <- -(gamma(r + 2*xi)/(tau^2*(1 + 2*xi)*gamma(r)))
        m23b <- (-(gamma(r + xi)/(1 + xi)) + gamma(r + 2*xi)/(1 + 2*xi))/(tau*xi*gamma(r))
        m33b <- (-1 + (2*gamma(r + xi))/(gamma(r) + xi*gamma(r)) -
                   gamma(r + 2*xi)/(gamma(r) + 2*xi*gamma(r)))/xi^2
      } else{
        m11b <- 0
        m12b <- 0
        m13b <- 0
        m22b <- -(1/tau^2)
        m23b <- (-1 + digamma(r))/tau
        m33b <- -2 + 2*digamma(r) - digamma(r)^2 - psigamma(deriv = 1, r)
      }
      mb <- matrix(c(m11b, m12b, m13b, m12b, m22b, m23b, m13b, m23b, m33b), nrow = 3, ncol = 3)
      return(-nobs*(ma + (r-1)*mb))
    }

  } else if(method == "obs"){
    #Partial check for ordering in first column
    if(!isTRUE(all(dat[1,1] >= dat[1,]))){
      stop("Observations in `dat` must be ordered from largest to smallest in each row.")
    }
    #Observed information matrix
    dat <- as.matrix(dat)
    r <- ncol(dat)
    yr <- dat[,r]

    if(!xizero){
      d11a <- sum(xi^2*(r/xi + 1)/(tau^2*((eta - yr)*xi/tau - 1)^2) + (-(eta - yr)*xi/tau + 1)^(1/xi - 2)*xi*(1/xi - 1)/(tau^2*(-(eta - yr)*xi/tau + 1)^(2/xi)) - 2*(-(eta - yr)*xi/tau + 1)^(2/xi - 2)/(tau^2*(-(eta - yr)*xi/tau + 1)^(3/xi)))
      d12a <- sum(xi*(r/xi + 1)/(tau^2*((eta - yr)*xi/tau - 1)) - (eta - yr)*xi^2*(r/xi + 1)/(tau^3*((eta - yr)*xi/tau - 1)^2) - (eta - yr)*(-(eta - yr)*xi/tau + 1)^(1/xi - 2)*xi*(1/xi - 1)/(tau^3*(-(eta - yr)*xi/tau + 1)^(2/xi)) + (-(eta - yr)*xi/tau + 1)^(1/xi - 1)/(tau^2*(-(eta - yr)*xi/tau + 1)^(2/xi)) + 2*(eta - yr)*(-(eta - yr)*xi/tau + 1)^(2/xi - 2)/(tau^3*(-(eta - yr)*xi/tau + 1)^(3/xi)))
      d13a <- sum(-(r/xi + 1)/(tau*((eta - yr)*xi/tau - 1)) + (eta - yr)*xi*(r/xi + 1)/(tau^2*((eta - yr)*xi/tau - 1)^2) - (-(eta - yr)*xi/tau + 1)^(1/xi - 1)*((eta - yr)*(1/xi - 1)/(tau*((eta - yr)*xi/tau - 1)) - log(-(eta - yr)*xi/tau + 1)/xi^2)/(tau*(-(eta - yr)*xi/tau + 1)^(2/xi)) - 2*(-(eta - yr)*xi/tau + 1)^(1/xi - 1)*(log(-(eta - yr)*xi/tau + 1)/xi^2 - (eta - yr)/(tau*((eta - yr)*xi/tau - 1)*xi))/(tau*(-(eta - yr)*xi/tau + 1)^(2/xi)) + r/(tau*((eta - yr)*xi/tau - 1)*xi))
      d22a <- sum(-2*(eta - yr)*xi*(r/xi + 1)/(tau^3*((eta - yr)*xi/tau - 1)) + (eta - yr)^2*xi^2*(r/xi + 1)/(tau^4*((eta - yr)*xi/tau - 1)^2) + (eta - yr)^2*(-(eta - yr)*xi/tau + 1)^(1/xi - 2)*xi*(1/xi - 1)/(tau^4*(-(eta - yr)*xi/tau + 1)^(2/xi)) + 1/tau^2 - 2*(eta - yr)*(-(eta - yr)*xi/tau + 1)^(1/xi - 1)/(tau^3*(-(eta - yr)*xi/tau + 1)^(2/xi)) - 2*(eta - yr)^2*(-(eta - yr)*xi/tau + 1)^(2/xi - 2)/(tau^4*(-(eta - yr)*xi/tau + 1)^(3/xi)))
      d23a <- sum((eta - yr)*(r/xi + 1)/(tau^2*((eta - yr)*xi/tau - 1)) - (eta - yr)^2*xi*(r/xi + 1)/(tau^3*((eta - yr)*xi/tau - 1)^2) + (eta - yr)*(-(eta - yr)*xi/tau + 1)^(1/xi - 1)*(log(-(eta - yr)*xi/tau + 1)/xi^2 - (eta - yr)/(tau*((eta - yr)*xi/tau - 1)*xi))/(tau^2*(-(eta - yr)*xi/tau + 1)^(2/xi)) - (eta - yr)*r/(tau^2*((eta - yr)*xi/tau - 1)*xi) + (eta - yr)^2/(tau^3*(-(eta - yr)*xi/tau + 1)^(1/xi)*((eta - yr)*xi/tau - 1)^2))
      d33a <- sum(-(log(-(eta - yr)*xi/tau + 1)/xi^2 - (eta - yr)/(tau*((eta - yr)*xi/tau - 1)*xi))^2/(-(eta - yr)*xi/tau + 1)^(1/xi) + (2*log(-(eta - yr)*xi/tau + 1)/xi^3 - 2*(eta - yr)/(tau*((eta - yr)*xi/tau - 1)*xi^2) - (eta - yr)^2/(tau^2*((eta - yr)*xi/tau - 1)^2*xi))/(-(eta - yr)*xi/tau + 1)^(1/xi) + (eta - yr)^2*(r/xi + 1)/(tau^2*((eta - yr)*xi/tau - 1)^2) - 2*r*log(-(eta - yr)*xi/tau + 1)/xi^3 + 2*(eta - yr)*r/(tau*((eta - yr)*xi/tau - 1)*xi^2))
      da <- - matrix(c(d11a, d12a, d13a, d12a, d22a, d23a, d13a, d23a, d33a), nrow = 3, ncol = 3)
    } else {
      d110a <- sum(-exp(eta/tau - yr/tau)/tau^2)
      d120a <- sum(-(r*tau*exp(yr/tau) - (eta + tau)*exp(eta/tau) + yr*exp(eta/tau))*exp(-yr/tau)/tau^3)
      d130a <- sum(1/2*(2*(eta + tau)*yr*exp(eta/tau) - yr^2*exp(eta/tau) - (eta^2 + 2*eta*tau)*exp(eta/tau) + 2*(eta*r*tau - r*tau*yr + tau^2)*exp(yr/tau))*exp(-yr/tau)/tau^3)
      d220a <- sum((2*(eta + tau)*yr*exp(eta/tau) - yr^2*exp(eta/tau) - (eta^2 + 2*eta*tau)*exp(eta/tau) + (2*eta*r*tau - 2*r*tau*yr + tau^2)*exp(yr/tau))*exp(-yr/tau)/tau^4)
      d230a <- sum(1/2*((3*eta + 2*tau)*yr^2*exp(eta/tau) - yr^3*exp(eta/tau) - (3*eta^2 + 4*eta*tau)*yr*exp(eta/tau) + (eta^3 + 2*eta^2*tau)*exp(eta/tau) - 2*(eta^2*r*tau + r*tau*yr^2 + eta*tau^2 - (2*eta*r*tau + tau^2)*yr)*exp(yr/tau))*exp(-yr/tau)/tau^4)
      d330a <- sum(1/12*(4*(3*eta + 2*tau)*yr^3*exp(eta/tau) - 3*yr^4*exp(eta/tau) - 6*(3*eta^2 + 4*eta*tau)*yr^2*exp(eta/tau) + 12*(eta^3 + 2*eta^2*tau)*yr*exp(eta/tau) - (3*eta^4 + 8*eta^3*tau)*exp(eta/tau) + 4*(2*eta^3*r*tau - 2*r*tau*yr^3 + 3*eta^2*tau^2 + 3*(2*eta*r*tau + tau^2)*yr^2 - 6*(eta^2*r*tau + eta*tau^2)*yr)*exp(yr/tau))*exp(-yr/tau)/tau^4)
      da <- - matrix(c(d110a, d120a, d130a, d120a, d220a, d230a, d130a, d230a, d330a), nrow = 3, ncol = 3)
    }
    if(r == 1){
      return(da)
    } else if(r > 1){
      y <- dat[,-r]
      if(!xizero){
        d11b <- sum(-2*xi^3*(y - yr)*(1/xi + 1)/(((eta - yr)*xi - tau)^3*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)) + xi^4*(y - yr)^2*(1/xi + 1)/(((eta - yr)*xi - tau)^4*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)^2) + xi^2/((eta - yr)*xi - tau)^2)
        d12b <- sum(2*xi^2*(y - yr)*(1/xi + 1)/(((eta - yr)*xi - tau)^3*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)) - xi^3*(y - yr)^2*(1/xi + 1)/(((eta - yr)*xi - tau)^4*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)^2) - xi/((eta - yr)*xi - tau)^2)
        d22b <- sum(-2*xi*(y - yr)*(1/xi + 1)/(((eta - yr)*xi - tau)^3*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)) + xi^2*(y - yr)^2*(1/xi + 1)/(((eta - yr)*xi - tau)^4*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)^2) + 1/((eta - yr)*xi - tau)^2)
        d13b <- sum(-2*(eta - yr)*xi^2*(y - yr)*(1/xi + 1)/(((eta - yr)*xi - tau)^3*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)) + xi^2*((eta - yr)*xi*(y - yr)/((eta - yr)*xi - tau)^2 - (y - yr)/((eta - yr)*xi - tau))*(y - yr)*(1/xi + 1)/(((eta - yr)*xi - tau)^2*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)^2) + (eta - yr)*xi/((eta - yr)*xi - tau)^2 + 2*xi*(y - yr)*(1/xi + 1)/(((eta - yr)*xi - tau)^2*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)) - 1/((eta - yr)*xi - tau) - (y - yr)/(((eta - yr)*xi - tau)^2*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)))
        d23b <- sum(2*(eta - yr)*xi*(y - yr)*(1/xi + 1)/(((eta - yr)*xi - tau)^3*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)) - xi*((eta - yr)*xi*(y - yr)/((eta - yr)*xi - tau)^2 - (y - yr)/((eta - yr)*xi - tau))*(y - yr)*(1/xi + 1)/(((eta - yr)*xi - tau)^2*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)^2) - (eta - yr)/((eta - yr)*xi - tau)^2 - (y - yr)*(1/xi + 1)/(((eta - yr)*xi - tau)^2*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)) + (y - yr)/(((eta - yr)*xi - tau)^2*xi*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)))
        d33b <- sum(((eta - yr)*xi*(y - yr)/((eta - yr)*xi - tau)^2 - (y - yr)/((eta - yr)*xi - tau))^2*(1/xi + 1)/(xi*(y - yr)/((eta - yr)*xi - tau) - 1)^2 - 2*((eta - yr)^2*xi*(y - yr)/((eta - yr)*xi - tau)^3 - (eta - yr)*(y - yr)/((eta - yr)*xi - tau)^2)*(1/xi + 1)/(xi*(y - yr)/((eta - yr)*xi - tau) - 1) + (eta - yr)^2/((eta - yr)*xi - tau)^2 - 2*((eta - yr)*xi*(y - yr)/((eta - yr)*xi - tau)^2 - (y - yr)/((eta - yr)*xi - tau))/(xi^2*(xi*(y - yr)/((eta - yr)*xi - tau) - 1)) - 2*log(-xi*(y - yr)/((eta - yr)*xi - tau) + 1)/xi^3)
      } else{
        d11b <- 0
        d12b <- 0
        d22b <- sum((tau - 2*y + 2*yr)/tau^3)
        d13b <- sum((tau - y + yr)/tau^2)
        d23b <- sum(-(eta*tau - (2*eta + tau)*y + y^2 + 2*eta*yr - yr^2)/tau^3)
        d33b <- sum(1/3*(3*eta^2*tau + 3*(2*eta + tau)*y^2 - 2*y^3 + 6*eta^2*yr - 6*eta*yr^2 + 2*yr^3 - 6*(eta^2 + eta*tau)*y)/tau^3)
      }
      db <-  - matrix(c(d11b, d12b, d13b, d12b, d22b, d23b, d13b, d23b, d33b), nrow = 3, ncol = 3)
      return(da + db)
    }
  }
}

#' Poisson process of extremes.
#'
#' Likelihood, score function and information matrix for the Poisson process likelihood.
#'
#' @author Leo Belzile
#' @name pp
#' @param par vector of \code{loc}, \code{scale} and \code{shape}
#' @param dat sample vector
#' @param u threshold
#' @param method string indicating whether to use the expected  (\code{'exp'}) or the observed (\code{'obs'} - the default) information matrix.
#' @param np number of periods of observations. This is a \emph{post hoc} adjustment for the intensity so that the parameters of the model coincide
#' with those of a generalized extreme value distribution with block size \code{length(dat)/np}.
#' @param nobs number of observations for the expected information matrix. Default to \code{length(dat)} if \code{dat} is provided.
#' @section Usage:
#' \preformatted{pp.ll(par, dat)
#' pp.ll(par, dat, u, np)
#' pp.score(par, dat)
#' pp.infomat(par, dat, method = c('obs', 'exp'))}
#'
#' @section Functions:
#'\itemize{
#' \item{\code{pp.ll}:} {log likelihood}
#' \item{\code{pp.score}:} {score vector}
#' \item{\code{pp.infomat}:} {observed or expected information matrix}
#' }
#' @references Coles, S. (2001). \emph{An Introduction to Statistical Modeling of Extreme Values}, Springer, 209 p.
#' @references Wadsworth, J.L. (2016). Exploiting Structure of Maximum Likelihood Estimators for Extreme Value Threshold Selection, \emph{Technometrics}, \bold{58}(1), 116-126, \code{http://dx.doi.org/10.1080/00401706.2014.998345}.
#'
NULL



#Poisson likelihood

#' Poisson process log-likelihood
#' @export
#' @rdname pp
pp.ll <- function(par, dat, u, np = 1){
  par <- as.vector(par)
  if(abs(par[3]) > 1e-6){
  -np[1]*(1+par[3]*(u[1]-par[1])/par[2])^(-1/par[3]) -
    length(dat)*log(par[2]) - (1+1/par[3]) * sum(log1p(par[3]*(dat-par[1])/par[2]))
  } else{
    -np[1]*exp((par[1]-u)/par[2]) - length(dat)*log(par[2]) + sum(par[1] - dat)/par[2]

  }
}

#' Score vector for the Poisson process above u
#'
#' @seealso \code{\link{pp}}
#' @inheritParams pp
#' @export
#' @rdname pp
#' @keywords internal
pp.score <- function(par, dat, u, np = 1){
  mu <- par[1]; sigma <- par[2]; xi <- par[3]; nu <- length(dat)
  c(  -sum(xi*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1))) -  np*(-(mu - u)*xi/sigma + 1)^(-1/xi - 1)/sigma,
       sum((mu - dat)*xi*(1/xi + 1)/(sigma^2*((mu - dat)*xi/sigma - 1))) - nu/sigma + (mu - u)*np*((u - mu)*xi/sigma + 1)^(1/xi - 1)/(sigma^2*((u - mu)*xi/sigma + 1)^(2/xi)),
       - sum((mu - dat)*(1/xi + 1)/(sigma*((mu - dat)*xi/sigma - 1))) + sum(log1p(-(mu - dat)*xi/sigma)/xi^2) -
         np*(log1p((u - mu)*xi/sigma)/xi^2 - (mu - u)/(sigma*((mu - u)*xi/sigma - 1)*xi))/(-(mu - u)*xi/sigma + 1)^(1/xi))
}



#' Information matrix for the Poisson process likelihood
#'
#' The function returns the expected or observed information matrix.
#'
#' @note For the expected information matrix, the number of points above the threshold is random, but should correspond to
#' \code{np}\eqn{\Lambda}. In this sense, unless \code{np} is linear in the sample size, the
#'
#' @seealso \code{\link{pp}}
#' @inheritParams pp
#' @export
#' @rdname pp
#' @keywords internal
pp.infomat <- function(par, dat, method = c("obs", "exp"), u, np = 1, nobs = length(dat)) {
  method <- match.arg(method)
  if(method == "obs"){
  dat <- as.vector(dat)
  if(sum(dat >u) != length(dat)){
   warning("`dat` contains non-exceedances")
    dat <- dat[dat >u]
  }
  }
  mu <- par[1]
  sigma <- par[2]
  xi <- par[3]
  r1 <- nobs;
  m <- np

  vm <- (u - mu)/sigma
  xizero <- abs(xi) < 1e-5

   if(method == "obs"){
     zm <- (dat - mu)/sigma
     if(!xizero){
       d11 <- -sum(xi^2*(1/xi + 1)/((xi*zm + 1)^2*sigma^2)) - (vm*xi + 1)^(1/xi - 2)*m*xi*(1/xi - 1)/((vm*xi + 1)^(2/xi)*sigma^2) + 2*(vm*xi + 1)^(2/xi - 2)*m/((vm*xi + 1)^(3/xi)*sigma^2)
       d12 <- (vm*xi + 1)^(-1/xi - 2)*m*vm*xi*(1/xi + 1)/sigma^2 - xi^2*sum(zm*(1/xi + 1)/((xi*zm + 1)^2*sigma^2)) -
         (vm*xi + 1)^(-1/xi - 1)*m/sigma^2 + xi*(1/xi + 1)*sum(1/(xi*zm + 1))/sigma^2
       d13 <- -(vm*xi + 1)^(-1/xi - 1)*m*(vm/((vm*xi + 1)*xi) - log(vm*xi + 1)/xi^2)/sigma + sum(xi*zm*(1/xi + 1)/((xi*zm + 1)^2*sigma)) -
         (vm*xi + 1)^(-1/xi)*m*vm/((vm*xi + 1)^2*sigma) - sum((1/xi + 1)/((xi*zm + 1)*sigma)) + sum(1/((xi*zm + 1)*sigma*xi))
       d22 <- (vm*xi + 1)^(-1/xi - 2)*m*vm^2*xi*(1/xi + 1)/sigma^2 - sum(xi^2*zm^2*(1/xi + 1)/((xi*zm + 1)^2*sigma^2)) - 2*(vm*xi + 1)^(-1/xi - 1)*m*vm/sigma^2 + 2*xi*sum(zm*(1/xi + 1)/((xi*zm + 1)*sigma^2)) - r1/sigma^2
       d23 <- -(vm*xi + 1)^(-1/xi - 1)*m*vm*(vm/((vm*xi + 1)*xi) - log(vm*xi + 1)/xi^2)/sigma + sum(xi*zm^2*(1/xi + 1)/((xi*zm + 1)^2*sigma)) - (vm*xi + 1)^(-1/xi)*m*vm^2/((vm*xi + 1)^2*sigma) - sum(zm*(1/xi + 1)/((xi*zm + 1)*sigma)) + sum(zm/((xi*zm + 1)*sigma*xi))
       d33 <- (vm*xi + 1)^(-1/xi)*m*(vm/((vm*xi + 1)*xi) - log(vm*xi + 1)/xi^2)^2 + (vm*xi + 1)^(-1/xi)*m*(vm^2/((vm*xi + 1)^2*xi) + 2*vm/((vm*xi + 1)*xi^2) - 2*log(vm*xi + 1)/xi^3) - sum(zm^2*(1/xi + 1)/(xi*zm + 1)^2) - 2*sum(zm/((xi*zm + 1)*xi^2)) + 2*sum(log(xi*zm + 1))/xi^3
     } else{ #Observed information, xi = 0
       d11 <- m*exp(-vm)/sigma^2
       d12 <- m*vm*exp(-vm)/sigma^2 - m*exp(-vm)/sigma^2 + r1/sigma^2
       d22 <- m*vm^2*exp(-vm)/sigma^2 - 2*m*vm*exp(-vm)/sigma^2 + 2*sum(zm)/sigma^2 - r1/sigma^2
       d13 <- 0.5*(m*vm^2 + 2*sum(zm)*exp(vm) - 2*m*vm - 2*r1*exp(vm))*exp(-vm)/sigma
       d23 <- 0.5*(m*vm^3 + 2*sum(zm^2)*exp(vm) - 2*m*vm^2 - 2*sum(zm)*exp(vm))*exp(-vm)/sigma
       d33 <- ((1/12)*(m*vm^3*(-8 + 3*vm) + 4*exp(vm)*sum(zm^2*(-3 + 2*zm))))/exp(vm)
     }
     return(matrix(c(d11,d12,d13,d12,d22,d23,d13,d23,d33), nrow = 3, ncol = 3, byrow = TRUE))
} else{
  r2 <- r1 <- 1
  if(!xizero){
    e11 <- (vm*xi + 1)^(-1/xi)*m*r2*(xi + 1)*xi/((sigma*vm*xi + sigma)^2*(2*xi + 1)) - (vm*xi + 1)^(-1/xi - 2)*m*xi*(1/xi + 1)/sigma^2
    e12 <- -(vm*xi + 1)^(-1/xi - 2)*m*vm*xi*(1/xi + 1)/sigma^2 - (vm*xi + 1)^(-1/xi)*m*r2*(xi + 1)/((sigma*vm*xi + sigma)^2*(2*xi + 1)) +
          (vm*xi + 1)^(-1/xi - 1)*m/sigma^2
    e13 <- (vm*xi + 1)^(-1/xi - 1)*m*(vm*(1/xi + 1)/(vm*xi + 1) - log(vm*xi + 1)/xi^2)/sigma -
      (2*vm*xi + vm - xi)*(vm*xi + 1)^(-1/xi)*m*r2/((vm*xi + 1)^2*(2*xi^2 + 3*xi + 1)*sigma)
    e22 <- -(vm*xi + 1)^(-1/xi - 2)*m*vm^2*xi*(1/xi + 1)/sigma^2 + 2*(vm*xi + 1)^(-1/xi - 1)*m*vm/sigma^2 +
      ((vm*xi + 1)^2*r1*(2*xi + 1) - ((vm*xi + 2)*vm*(2*xi + 1) + 2)*r2*(xi + 1))*(vm*xi + 1)^(-1/xi)*m/((sigma*vm*xi + sigma)^2*(2*xi + 1))
    e23 <- (vm*xi + 1)^(-1/xi - 1)*m*vm*(vm*(1/xi + 1)/(vm*xi + 1) - log(vm*xi + 1)/xi^2)/sigma -
      ((vm*xi + vm + 1)*vm*(2*xi + 1) + 1)*(vm*xi + 1)^(-1/xi)*m*r2/((vm*xi + 1)^2*sigma*(2*xi + 1)*(xi + 1))
    e33 <- -(vm*xi + 1)^(-1/xi)*m*(vm/((vm*xi + 1)*xi) - log(vm*xi + 1)/xi^2)^2 -
      (vm*xi + 1)^(-1/xi)*m*(vm^2/((vm*xi + 1)^2*xi) + 2*vm/((vm*xi + 1)*xi^2) - 2*log(vm*xi + 1)/xi^3) -
      (2*(vm*xi + 1)^2*(2*xi + 1)*(xi + 1)*log(vm*xi + 1) +
         (vm^2*(2*xi + 1)*(xi + 1)*(xi - 3)*xi + 2*((2*xi^2 - xi - 3)*xi - 1)*vm + 2*xi^2)*xi)*
      (vm*xi + 1)^(-1/xi)*m*r2/((vm*xi + 1)^2*(2*xi + 1)*(xi + 1)*xi^3)
  } else{
    e11 <- -m*exp(-vm)/sigma^2
    e12 <- -m*r2*exp(-vm)/sigma^2 - m*vm*exp(-vm)/sigma^2 + m*exp(-vm)/sigma^2
    e13 <- -1/2*(2*r2*vm + vm^2 - 2*vm)*m*exp(-vm)/sigma
    e22 <- -2*m*r2*vm*exp(-vm)/sigma^2 - m*vm^2*exp(-vm)/sigma^2 + m*r1*exp(-vm)/sigma^2 - 2*m*r2*exp(-vm)/sigma^2 + 2*m*vm*exp(-vm)/sigma^2
    e23 <- -1/2*(2*r2*vm^2 + vm^3 + 2*r2*vm - 2*vm^2 + 2*r2)*m*exp(-vm)/sigma
    e33 <- -1/12*(8*r2*vm^3 + 3*vm^4 + 12*r2*vm^2 - 8*vm^3 + 24*r2*vm + 24*r2)*m*exp(-vm)
}
  return(-matrix(c(e11,e12,e13,e12,e22,e23,e13,e23,e33), nrow = 3, ncol = 3, byrow = TRUE))
}

}
