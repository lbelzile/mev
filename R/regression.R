#' General linear models with generalized extreme value distribution
#'
#' This function returns an object of class \code{mev_gev}, with default methods for printing and quantile-quantile plots.
#' The default starting values are the solution of the probability weighted moments.
#' @param data a data frame with response and observations
#' @param start named list of starting values
#' @param method string indicating the outer optimization routine for the augmented Lagrangian. One of \code{nlminb} or \code{BFGS}.
#' @return a list containing the following components:
#' \itemize{
#' \item \code{estimate} a vector containing the maximum likelihood estimates.
#' \item \code{std.err} a vector containing the standard errors.
#' \item \code{vcov} the variance covariance matrix, obtained as the numerical inverse of the observed information matrix.
#' \item \code{method} the method used to fit the parameter.
#' \item \code{nllh} the negative log-likelihood evaluated at the parameter \code{estimate}.
#' \item \code{convergence} components taken from the list returned by \code{\link[alabama]{auglag}}.
#' Values other than \code{0} indicate that the algorithm likely did not converge.
#' \item \code{counts} components taken from the list returned by \code{\link[alabama]{auglag}}.
#' \item \code{xdat} vector of data
#' }
#' @examples
#' data(venice, package = "mev")
#' proto <- lm.gev(
#' formula = list(
#'   loc = ~year, scale = ~1, shape = ~1),
#'   response = "r1",
#'   data = venice)
lm.gev <- function(
  formula = list(loc = ~1, scale = ~1, shape = ~1),
  response,
  data = NULL,
  link = list(loc = identity, scale = log, shape = identity),
  invlink = list(loc = identity, scale = exp, shape = identity),
  start = NULL,
  scale = TRUE,
  method = c("nlminb", "BFGS"),
  show = FALSE,
  ...
) {
  param_names <- c("loc", "scale", "shape")
  # if (!isTRUE(all(sapply(link, is.function)))) {
  #   stop("Argument \"link\" should be a named list of functions.")
  # }
  if (!isTRUE(all(sort(names(link)) == param_names))) {
    stop(
      "Argument \"link\" should have elements with names \"loc\", \"scale\" and \"shape\"."
    )
  }
  if (!isTRUE(all(sapply(invlink, is.function)))) {
    stop("Argument \"invlink\" should be a named list of functions.")
  }
  if (!isTRUE(all(sort(names(invlink)) == param_names))) {
    stop(
      "Argument \"invlink\" should have elements with names \"loc\", \"scale\" and \"shape\"."
    )
  }
  # Check that link functions match their inverses
  x <- runif(1)
  matches <- c(
    link$loc(invlink$loc(x)),
    invlink$loc(link$loc(x)),
    link$scale(invlink$scale(x)),
    invlink$scale(link$scale(x)),
    link$shape(invlink$shape(x)),
    invlink$shape(link$shape(x))
  )
  if (
    !isTRUE(all(sapply(matches, function(xv) {
      isTRUE(all.equal(xv, x, check.attributes = FALSE))
    })))
  ) {
    stop("Some link function do not match their inverse.")
  }
  if (
    !isTRUE(all(sapply(formula, FUN = function(x) {
      inherits(x, "formula")
    })))
  ) {
    stop("Argument \"formula\" should be a named list of formulas.")
  }
  if (!isTRUE(all(sort(names(formula)) == param_names))) {
    stop(
      "Argument \formula\" should have elements with names \"loc\", \"scale\" and \"shape\"."
    )
  }
  y <- NULL
  # Cast matrix to data frame
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  if (is.data.frame(data)) {
    # By default, 'data' is a function, here NULL argument
    y <- data[[response]]
    mloc <- model.matrix(formula$loc, data = data)
    mscale <- model.matrix(formula$scale, data = data)
    mshape <- model.matrix(formula$shape, data = data)
  } else {
    # Call without the 'data' argument
    mloc <- model.matrix(formula$loc)
    mscale <- model.matrix(formula$scale)
    mshape <- model.matrix(formula$shape)
  }
  dimn <- c(
    dimnames(mloc)[[2]][1],
    dimnames(mscale)[[2]][1],
    dimnames(mshape)[[2]][1]
  )
  # Some arguments may be NULL if the model matrix is empty
  # then, concatenating vectors yields length < 3
  if (!isTRUE(all(dimn == "(Intercept)", length(dimn) == 3L))) {
    stop("Invalid formula: all formulas should minimally include an intercept.")
  }
  # Get response variable from either data frame or environment
  if (is.null(y)) {
    y <- try(get(response), silent = TRUE)
    if (inherits(y, "try-error")) {
      stop(
        "Could not find \"response\" in either current environment or the data frame \"data\"."
      )
    }
  }
  y <- as.numeric(y)
  nobs <- length(y)
  loc_loc <- loc_scale <- loc_shape <- NULL
  sc_loc <- sc_scale <- sc_shape <- NULL
  # Standardize matrices, or make into scalar
  if (isTRUE(dim(mloc)[2] == 1L)) {
    msc_loc <- 1
    n_loc <- 1L
  } else if (isTRUE(dim(mloc)[2] > 1)) {
    if (isTRUE(scale)) {
      msc_loc <- scale(mloc[, -1])
      loc_loc <- attr(msc_loc, "scaled:center")
      sc_loc <- attr(msc_loc, "scaled:scale")
      msc_loc <- cbind(1, msc_loc)
    } else {
      msc_loc <- mloc
    }
    stopifnot(nrow(msc_loc) %in% c(1L, nobs))
    n_loc <- ncol(msc_loc)
  }
  if (isTRUE(dim(mscale)[2] == 1L)) {
    msc_scale <- 1
    n_scale <- 1L
  } else if (isTRUE(dim(mscale)[2] > 1)) {
    if (isTRUE(scale)) {
      msc_scale <- scale(mscale[, -1])
      loc_scale <- attr(msc_scale, "scaled:center")
      sc_scale <- attr(msc_scale, "scaled:scale")
      msc_scale <- cbind(1, msc_scale)
    } else {
      msc_scale <- mscale
    }
    stopifnot(nrow(msc_scale) %in% c(1L, nobs))
    n_scale <- ncol(msc_scale)
  }
  if (isTRUE(dim(mshape)[2] == 1L)) {
    msc_shape <- 1
    n_shape <- 1L
  } else if (isTRUE(dim(mshape)[2] > 1)) {
    if (isTRUE(scale)) {
      msc_shape <- scale(mshape[, -1])
      loc_shape <- attr(msc_shape, "scaled:center")
      sc_shape <- attr(msc_shape, "scaled:scale")
      msc_shape <- cbind(1, msc_shape)
    } else {
      msc_shape <- mshape
    }
    stopifnot(nrow(msc_shape) %in% c(1L, nobs))
    n_shape <- 1L
  }
  # Scale response (GEV is a location-scale family)
  # BUT back-transformation is complicated if
  # non-identity link for say scale
  # Should it affect only the intercept?
  # y <- scale(y)
  # resp_loc <- attr(y, "scaled:center")
  # resp_scale <- attr(y, "scaled:scale")
  par_cst <- 1L ==
    c(
      length(msc_loc),
      length(msc_scale),
      length(msc_shape)
    )
  if (isTRUE(all(par_cst))) {
    warning(
      "No covariate is specified; use \"fit.gev\" to estimate the model instead."
    )
  }
  fitted <- list() # container
  method <- match.arg(method)
  y <- as.double(y[is.finite(y)])
  if (length(y) != nobs) {
    stop("Invalid arguments: response contains missing or infinite values.")
  }
  npars <- c(n_loc, n_scale, n_shape)
  ntpars <- sum(npars)
  nps <- c(1L, cumsum(npars)[-3] + 1L)
  # TODO figure out link function and it's inverse
  if (is.null(start)) {
    start <- fit.gev(xdat = y)
    init <- coef(start)
    spars <- c(
      link$loc(init[1]),
      rep(0, n_loc - 1L),
      link$scale(init[2]),
      rep(0, n_scale - 1L),
      link$shape(init[3]),
      rep(0, n_shape - 1L)
    )
  } else {
    stopifnot(
      is.list(start),
      isTRUE(all.equal(sort(names(start)), param_names))
    )
    init_loc <- start$loc
    stopifnot(length(init_loc) == n_loc)
    init_scale <- start$scale
    stopifnot(length(init_scale) == n_scale)
    init_shape <- start$shape
    stopifnot(length(init_shape) == nshape)
    spars <- c(init_loc, init_scale, init_shape)
  }
  m2p <- function(spars, invlink, msc_loc, msc_scale, msc_shape...) {
    if (par_cst[1]) {
      mu <- spars[nps[1]]
    } else {
      mu <- c(
        msc_loc %*%
          spars[seq(nps[1], length.out = npars[1], by = 1L)]
      )
    }
    mu <- invlink$loc(mu)
    if (par_cst[2]) {
      sigma <- spars[nps[2]]
    } else {
      sigma <- c(
        msc_scale %*%
          spars[seq(nps[2], length.out = npars[2], by = 1L)]
      )
    }
    sigma <- invlink$scale(sigma)
    if (par_cst[3]) {
      xi <- spars[nps[3]]
    } else {
      xi <- c(
        msc_shape %*%
          spars[seq(nps[3], length.out = npars[3], by = 1L)]
      )
    }
    xi <- invlink$shape(xi)
    return(list(mu = mu, sigma = sigma, xi = xi))
  }

  # Check boundary constraints
  mle <- try(suppressWarnings(
    alabama::auglag(
      par = spars,
      fn = function(spars, invlink, ...) {
        pars <- m2p(
          spars = spars,
          invlink = invlink,
          msc_loc = msc_loc,
          msc_scale = msc_scale,
          msc_shape = msc_shape
        )
        mu <- pars$mu
        sigma <- pars$sigma
        xi <- pars$xi
        if (isTRUE(any(sigma <= 0) | any(xi < -1))) {
          return(1e10)
        }
        nll <- -sum(revdbayes::dgev(
          x = y,
          loc = mu,
          scale = sigma,
          shape = xi,
          log = TRUE
        ))
        ifelse(is.finite(nll), nll, 1e10)
      },
      hin = function(spars, invlink, ...) {
        pars <- m2p(
          spars = spars,
          invlink = invlink,
          msc_loc = msc_loc,
          msc_scale = msc_scale,
          msc_shape = msc_shape
        )
        mu <- pars$mu
        sigma <- pars$sigma
        xi <- pars$xi
        c(
          min(sigma + xi * (y - mu)),
          min(sigma),
          min(xi) + 1
        )
      },
      y = y,
      invlink = invlink,
      msc_loc = msc_loc,
      msc_scale = msc_scale,
      msc_shape = msc_shape,
      control.outer = list(method = method, trace = FALSE),
      control.optim = switch(
        method,
        nlminb = list(
          iter.max = 500L,
          rel.tol = 1e-10,
          step.min = 1e-10
        ),
        list(maxit = 1000L, reltol = 1e-10)
      )
    )
  ))
  # Special case of MLE on the boundary xi = -1
  if (inherits(mle, what = "try-error")) {
    stop("Optimization routine for the GEV did not converge.")
  }
  fitted$nllh <- mle$value
  fitted$estimate <- mle$par
  fitted$hessian <- mle$hessian
  fitted$parameters <- m2p(
    spars = fitted$estimate,
    invlink = invlink,
    msc_loc = msc_loc,
    msc_scale = msc_scale,
    msc_shape = msc_shape
  )
  beta_mu <- fitted$estimate[seq(from = nps[1], length.out = npars[1], by = 1)]
  beta_sigma <- fitted$estimate[seq(
    from = nps[2],
    length.out = npars[2],
    by = 1
  )]
  beta_xi <- fitted$estimate[seq(from = nps[3], length.out = npars[3], by = 1)]
  jac <- diag(sum(npars))
  if (scale) {
    # browser()
    # Back-transform covariates and coefficients
    if (n_loc > 1L) {
      beta_mu[1] <- beta_mu[1] + sum(loc_loc / sc_loc * beta_mu[-1])
      beta_mu[-1] <- beta_mu[-1] / sc_loc
      jac[nps[1], 1:npars[1]] <- c(1, loc_loc / sc_loc)
      jac[2:npars[1], 2:npars[1]] <- diag(
        x = 1 / sc_loc,
        nrow = npars[1] - 1,
        ncol = npars[1] - 1
      )
    }
    if (n_scale > 1L) {
      beta_sigma[1] <- beta_sigma[1] -
        sum(loc_scale / sc_scale * beta_sigma[-1])
      beta_sigma[-1] <- beta_sigma[-1] / sc_scale
      jac[nps[2], seq(nps[2], length.out = npars[2], by = 1L)] <- c(
        1,
        -loc_scale / sc_scale
      )
      jac[nps[2] + 2:npars[2], nps[2] + 2:npars[2]] <- diag(
        x = 1 / sc_scale,
        nrow = npars[2] - 1,
        ncol = npars[2] - 1
      )
    }
    if (n_shape > 1L) {
      beta_xi[1] <- beta_xi[1] - sum(loc_shape / sc_shape * beta_xi[-1])
      beta_xi[-1] <- beta_xi[-1] / sc_shape
      jac[nps[3], seq(nps[3], length.out = npars[3], by = 1L)] <- c(
        1,
        -loc_shape / sc_shape
      )
      jac[nps[3] + 2:npars[3], nps[3] + 2:npars[3]] <- diag(
        x = 1 / sc_shape,
        nrow = npars[3] - 1,
        ncol = npars[3] - 1
      )
    }
    fitted$estimate <- c(beta_mu, beta_sigma, beta_xi)
    optfn <- function(spars, invlink, ...) {
      pars <- m2p(
        spars = spars,
        invlink = invlink,
        msc_loc = mloc,
        msc_scale = mscale,
        msc_shape = mshape
      )
      mu <- pars$mu
      sigma <- pars$sigma
      xi <- pars$xi
      if (isTRUE(any(sigma <= 0) | any(xi < -1))) {
        return(1e10)
      }
      nll <- -sum(revdbayes::dgev(
        x = y,
        loc = mu,
        scale = sigma,
        shape = xi,
        log = TRUE
      ))
      ifelse(is.finite(nll), nll, 1e10)
    }
    # browser()
    # hess <- numDeriv::(
    #   func = optfn,
    #   x = fitted$estimate,
    #   invlink = invlink
    # )

    mle2 <- try(suppressWarnings(
      alabama::auglag(
        par = fitted$estimate,
        fn = function(spars, invlink, ...) {
          pars <- m2p(
            spars = spars,
            invlink = invlink,
            msc_loc = mloc,
            msc_scale = mscale,
            msc_shape = mshape
          )
          mu <- pars$mu
          sigma <- pars$sigma
          xi <- pars$xi
          if (isTRUE(any(sigma <= 0) | any(xi < -1))) {
            return(1e10)
          }
          nll <- -sum(revdbayes::dgev(
            x = y,
            loc = mu,
            scale = sigma,
            shape = xi,
            log = TRUE
          ))
          ifelse(is.finite(nll), nll, 1e10)
        },
        hin = function(spars, invlink, ...) {
          pars <- m2p(
            spars = spars,
            invlink = invlink,
            msc_loc = mloc,
            msc_scale = mscale,
            msc_shape = mshape
          )
          mu <- pars$mu
          sigma <- pars$sigma
          xi <- pars$xi
          c(
            min(sigma + xi * (y - mu)),
            min(sigma),
            min(xi) + 1
          )
        },
        y = y,
        invlink = invlink,
        msc_loc = mloc,
        msc_scale = mscale,
        msc_shape = mshape,
        control.outer = list(method = "BFGS", trace = FALSE),
        control.optim = list(
          maxit = 1000L,
          parscale = fitted$estimate,
          reltol = 1e-10
        )
      )
    ))
    fitted$hessian <- jac %*% mle$hessian %*% t(jac)
  }
  fitted$link <- link
  fitted$invlink <- invlink
  #Observed information matrix and standard errors

  fitted$vcov <- matrix(NA, ncol = length(mle$par), nrow = length(mle$par))
  fitted$std.err <- rep(NA, length(mle$par))
  if (fitted$param[3] > -0.5) {
    vcovt <- try(solve(fitted$hessian))
    if (!inherits(vcovt, what = "try-error")) {
      fitted$vcov <- vcovt
      fitted$std.err <- try(sqrt(diag(fitted$vcov)), silent = TRUE)
    }
  }
  names(fitted$std.err) <- names(fitted$estimate) <- c(
    paste0("loc ", colnames(mloc)),
    paste0("scale ", colnames(mscale)),
    paste0("shape ", colnames(mshape))
  )
  fitted$method <- "auglag"
  fitted$nobs <- length(y)
  if (isTRUE(all(mle$kkt1, mle$kkt2))) {
    fitted$convergence <- "successful"
  } else {
    fitted$convergence <- "dubious convergence"
  }
  fitted$counts <- mle$counts
  fitted$xdat <- y
  fitted$start <- spars #if start not provided, this is PWM
  class(fitted) <- "mev_gev_reg"
  if (show) {
    print(fitted)
  }
  invisible(fitted)
}
