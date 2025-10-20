# Map data to standard Laplace scale
qlaplace <- function(x) {
  ifelse(x < 0.5, log(2 * x), -log(2 * (1 - x)))
}

# Compute constants for nonlinear inequality constraints
condex.dataconstr <- function(
  par,
  data,
  thresh = 0.95,
  index = 1:ncol(data),
  ...
) {
  zpos <- c()
  zneg <- c()
  z <- c()
  v <- c()
  alpha <- par[1]
  beta <- par[2]
  u <- qexp(2 * thresh - 1)
  for (ii in seq_along(index)) {
    i <- index[ii]
    X <- data[data[, i] > u, -i, drop = FALSE]
    X0 <- data[data[, i] > u, i]
    Z <- (X - X0 * alpha) / (X0^beta)
    v <- max(c(v, X0))
    zpos <- range(c(
      zpos,
      min(apply(X, 1, min) - X0),
      max(apply(X, 1, max) - X0)
    ))
    z <- range(c(z, Z))
    zneg <- range(c(
      zneg,
      min(apply(X, 1, min) + X0),
      max(apply(X, 1, max) + X0)
    ))
  }
  list(
    z = z,
    zpos = zpos,
    zneg = zneg,
    v = v + 0.1
  )
}

# Adapted from function texmex:::ConstraintsAreSatisfied
# by Harry Southworth and Janet Heffernan
# Adaptive barrier; returns negative value when violated
optim_constr_ht <- function(par, constants, ...) {
  a <- par[1]
  b <- par[2]
  z <- constants$z
  zpos <- constants$zpos
  zneg <- constants$zneg
  v <- constants$v
  C1e <- min(
    min(
      1,
      1 - b * min(z) * v^(b - 1),
      1 -
        v^(b - 1) * min(z) +
        min(zpos) /
          v
    ) -
      a,
    min(
      1,
      1 - b * max(z) * v^(b - 1),
      1 -
        v^(b - 1) * max(z) +
        max(zpos) /
          v
    ) -
      a
  )
  # In 'texmex', the constraint is violated as soon as we have a FALSE,
  # thanks to lazy evaluation - this means that NAs are ignored if the first arguments
  # are FALSE...
  C1o <- min(
    a - (1 - b * min(z) * v^(b - 1)),
    a - (1 - b * max(z) * v^(b - 1)),
    (1 - 1 / b) *
      (b * min(z))^(1 / (1 - b)) *
      (1 - a)^(-b / (1 - b)) +
      min(zpos),
    (1 - 1 / b) *
      (b * max(z))^(1 / (1 - b)) *
      (1 - a)^(-b / (1 - b)) +
      max(zpos)
  )
  C2e <- min(
    1 + b * v^(b - 1) * min(z),
    1 + v^(b - 1) * min(z) - min(zneg) / v,
    1 + b * v^(b - 1) * max(z),
    1 + v^(b - 1) * max(z) - max(zneg) / v
  ) +
    a
  C2o <- min(
    -a - (1 + b * v^(b - 1) * min(z)),
    -a - (1 + b * v^(b - 1) * max(z)),
    (1 - 1 / b) *
      (-b * min(z))^(1 / (1 - b)) *
      (1 + a)^(-b / (1 - b)) -
      min(zneg),
    (1 - 1 / b) *
      (-b * max(z))^(1 / (1 - b)) *
      (1 + a)^(-b / (1 - b)) -
      max(zneg)
  )
  C1 <- suppressWarnings(max(c(C1e, C1o), na.rm = TRUE))
  C2 <- suppressWarnings(max(c(C2e, C2o), na.rm = TRUE))
  min(C1, C2)
}


#' Conditional extremes pseudo log likelihood
#'
#' @param xdat matrix or data frame of observations on the Laplace scale
#' @param par vector of parameters, \eqn{\alpha}, \eqn{\beta}
#' @param index vector of indices over which to condition, by default the first component
#' @param thresh quantile level (uniform scale)
#' @param type estimating equation, either \code{norm} for normal or \code{skewnorm} for skew-normal distribution
#' @param constrain logical; if \code{TRUE}, apply the constraints from Keef, Papastathopoulos and Tawn (2013)
#' @param ... additional arguments, currently ignored
#' @export
#' @return a numeric value with the pseudo log likelihood, with additional arguments \code{"loc"} and \code{"scale"}
#' @keywords internal
condex.pseudoll <- function(
  par,
  xdat,
  index = 1,
  thresh = 0.95,
  constrain = FALSE,
  ...
) {
  args <- list(...)
  stopifnot(is.logical(constrain), length(constrain) == 1L)
  stopifnot(isTRUE(all(index %in% seq_len(ncol(data)))))
  if (!is.null(thresh)) {
    stopifnot(length(thresh) == 1L, thresh >= 0, thresh < 1)
  } else {
    stopifnot(length(index) == 1L)
  }
  if (abs(par[1]) > 1 | par[2] > 1) {
    return(-1e16)
  }
  data <- as.matrix(xdat)
  # Define parameters
  alpha <- par[1]
  beta <- par[2]
  # Vectors for constraints
  zpos <- c()
  zneg <- c()
  z <- c()
  v <- c()
  # For-loop
  u <- qexp(2 * thresh - 1)
  obj <- rep(0, length(index))

  loc <- scale <- matrix(nrow = length(index), ncol = ncol(data) - 1L)
  for (ii in seq_along(index)) {
    i <- index[ii]
    exc <- data[, i] > u & !is.na(data[, i])
    X <- data[exc, -i, drop = FALSE]
    X0 <- data[exc, i]
    n <- length(X0)
    stopifnot(length(par) == 2L)
    Z <- (X - X0 * alpha) / (X0^beta)
    for (j in seq_len(ncol(X))) {
      loc[ii, j] <- mean(Z[, j], na.rm = TRUE)
      scale[ii, j] <- sd(Z, na.rm = TRUE) * sqrt((n - 1) / n)
      scale_mod <- scale[ii, j] * (X0^beta)
      loc_mod <- alpha * X0 + (X0^beta) * loc[ii, j]
      obj[ii] <- obj[ii] +
        sum(dnorm(
          na.omit(X[, j]),
          mean = loc_mod,
          sd = scale_mod,
          log = TRUE
        ))
    }
    if (constrain) {
      Z <- (X - X0 * alpha) / (X0^beta)
      zpos <- range(
        c(
          zpos,
          min(apply(X, 1, min) - X0),
          max(apply(X, 1, max) - X0)
        ),
        na.rm = TRUE
      )
      z <- range(c(z, Z), na.rm = TRUE)
      zneg <- range(
        c(
          zneg,
          min(apply(X, 1, min) + X0),
          max(apply(X, 1, max) + X0)
        ),
        na.rm = TRUE
      )
      v <- max(c(v, X0))
    }
  }
  if (constrain) {
    const <- condex.dataconstr(
      a = alpha,
      b = beta,
      z = z,
      zpos = zpos,
      zneg = zneg,
      v = v + 0.1
    ) >=
      0
    if (!isTRUE(const)) {
      obj <- -1e20
    }
  }
  obj <- sum(obj)
  attr(obj, "loc") <- loc
  attr(obj, "scale") <- scale
  return(obj)
}


#' Heffernan-Tawn model profile pseudo log likelihood
#'
#' @param alpha_grid grid of values for \eqn{\alpha}
#' @param beta_grid grid of values for \eqn{\beta}
#' @export
profile.condex <- function(
  fitted,
  which = c("alpha", "beta"),
  data,
  maxsteps = 101L,
  index = 1,
  thresh = 0.95
) {
  which <- match.arg(
    arg = which,
    choices = c("alpha", "beta"),
    several.ok = TRUE
  )
  # TODO does not use fitted
  # calls undefined function 'texmex_ConstraintsAreSatisfied'
  # Profile is by default bivariate - want to get univariate
  # for profile-based confidence intervals
  #
  stopifnot(length(index) == 1L)
  alpha_grid = seq(-1, 1, length.out = maxsteps)
  beta_grid = seq(-3, 1, length.out = maxsteps)
  type <- match.arg(type)
  alpha_grid <- alpha_grid[alpha_grid >= -1 & alpha_grid <= 1]
  stopifnot(length(alpha_grid) >= 1)
  beta_grid <- beta_grid[beta_grid <= 1]
  stopifnot(length(beta_grid) >= 1)
  stopifnot(isTRUE(all(
    index %in%
      seq_len(ncol(
        data
      ))
  )))
  stopifnot(length(thresh) == 1L, thresh > 0.5, thresh < 1)
  data <- as.matrix(data)
  u <- qexp(2 * thresh - 1)
  profloglik <- matrix(
    data = NA,
    nrow = length(alpha_grid),
    ncol = length(beta_grid)
  )
  D <- ncol(data)
  for (i in seq_along(alpha_grid)) {
    for (j in seq_along(beta_grid)) {
      alpha <- alpha_grid[i]
      beta <- beta_grid[j]
      Zs <- c()
      Xmin <- c()
      Xmax <- c()
      X0 <- c()
      for (ii in seq_along(index)) {
        k <- index[ii]
        X <- data[data[, k] > u, -k, drop = FALSE]
        X0s <- data[data[, k] > u, k]
        X0 <- c(X0, X0s)
        for (l in seq_len(ncol(X))) {
          Zs <- c(Zs, (X[, l] - X0s * alpha) / (X0s^beta))
        }
        Xmin <- c(Xmin, apply(X, 1, min))
        Xmax <- c(Xmax, apply(X, 1, max))
      }
      zpos <- c(min(Xmin - X0), max(Xmax - X0))
      z <- range(Zs)
      zneg <- c(min(Xmin + X0), max(Xmax + X0))
      constraint <- optim_constr_ht(
        a = alpha,
        b = beta,
        z = z,
        zpos = zpos,
        zneg = zneg,
        v = max(X0) + 0.1
      )
      if (isTRUE(any(constraint < 0))) {
        profloglik[i, j] <- -beta *
          (D - 1) *
          sum(log(X0)) +
          sum(dnorm(
            Zs,
            mean = mean(Zs),
            sd = sd(Zs),
            log = TRUE
          ))
      }
    }
  }
  return(list(
    alpha = alpha_grid,
    beta = beta_grid,
    pll = profloglik
  ))
}

#' Residuals from conditional extremes model
#'
#' @param alpha scale parameter
#' @param beta power parameter
#' @param data matrix of observations
#' @param thresh vector of thresholds on the uniform scale
#' @param testIndep logical; if \code{TRUE}, compute a test of independence using energy statistics for each conditioning variable in turn
#' @param group integer vector for the group; only observations in group 1 are considered considered to build residuals.
residuals.condex <- function(
  fitted,
  data,
  thresh
) {
  p <- ncol(data)
  stopifnot(length(group) == p, isTRUE(all(group %in% 1:2)))
  ng <- sum(group == 1)
  # Reorder so that observations in group 1 are first
  od <- c(which(group == 1), which(group != 1))
  data <- data[, od]
  thresh <- rep(thresh, ng)
  stopifnot(isTRUE(all(thresh > 0, thresh < 1, is.finite(thresh))))
  qlaplace <- function(x) {
    ifelse(x < 0.5, log(2 * x), -log(2 * (1 - x)))
  }
  u <- qlaplace(thresh)
  alpha <- rep(alpha, p)
  beta <- rep(beta, p)

  # Step 1: Identify exceedances and create vectors of residuals
  # Step 2: Compute maximum of each component (see constraint)
  # Step 3: Check constraint
  sdata <- data[
    apply(data[, 1:ng], 1, function(x) {
      isTRUE(any(x > u))
    }),
  ]
  nexc <- rowSums(t(sdata[, 1:ng]) >= u)
  res <- matrix(0, nrow = sum(nexc), ncol = p - 1)
  ind <- 0
  pvalindep <- rep(0, ng)
  for (i in seq_len(ng)) {
    exc <- which(sdata[, i] > u[i])
    Y0 <- sdata[exc, i]
    Z <- (sdata[exc, -i] - alpha[i] * Y0) / (Y0^beta[i])
    if (testIndep) {
      pvalindep[i] <- energy::indep.test(Z, Y0, R = 999)$p.value
    }
    res[(ind + 1):(ind + nexc[i]), ] <- Z
    ind <- ind + nexc[i]
  }
  if (!testIndep) {
    pvalindep <- NULL
  }
  return(list(
    pvalindep = pvalindep,
    ng = ng,
    res = res,
    nexc = nexc
  ))
}


#' Prediction from Heffernan-Tawn model
#'
#' @param alpha vector of \eqn{\alpha} parameters, recycled if necessary
#' @param beta vector of \eqn{\beta} parameters, recycled if necessary
#' @param data matrix of data
#' @param thresh vector of uniform thresholds
#' @param region risk region on uniform scale
#' @param B number of Monte Carlo replications
#' @param constraint string indicating whether margins are supposed equiprobable
#' or the relative frequency with which a variable is largest is estimated empirically (multinomial)
#' @param type string; should
predict.condex <- function(
  fitted,
  B = 1e5,
  nm = 1e4
) {
  p <- ncol(data)
  constraint <- match.arg(constraint)
  thresh <- rep(thresh, p)
  stopifnot(isTRUE(all(thresh > 0, thresh < 1, is.finite(thresh))))
  qlaplace <- function(x) {
    ifelse(x < 0.5, log(2 * x), -log(2 * (1 - x)))
  }
  u <- qlaplace(thresh)
  alpha <- rep(alpha, p)
  beta <- rep(beta, p)
  region <- rep(region, p)
  # Step 1: Identify exceedances and create vectors of residuals
  # Step 2: Compute maximum of each component (see constraint)
  # Step 3: Check constraint
  sdata <- data[
    apply(data, 1, function(x) {
      isTRUE(any(x > u))
    }),
  ]
  wmax <- table(apply(sdata, 1, which.max))
  nexc <- rowSums(t(sdata) >= u)
  res <- matrix(0, nrow = sum(nexc), ncol = p - 1)
  ind <- 0
  for (i in seq_len(p)) {
    exc <- which(sdata[, i] > u[i])
    res[(ind + 1):(ind + nexc[i]), ] <- (sdata[exc, -i] -
      alpha[i] * sdata[exc, i]) /
      (sdata[exc, i]^beta[i])
    ind <- ind + nexc[i]
  }
  prob <- rep(0, p)
  for (i in seq_len(p)) {
    nsim <- 0
    while (nsim < B) {
      # Generate exceedances above threshold
      y0 <- rexp(n = nm) + region[i]
      #
      newY <- alpha[i] *
        y0 +
        res[
          sample.int(n = nrow(res), size = nm, replace = TRUE),
          ,
          drop = FALSE
        ] *
          y0^beta[i]
      # Check if the conditioning variable is the largest
      condIsMax <- which(y0 > apply(newY, 1, max))
      # Increment the number of simulations completed
      nsim <- nsim + length(condIsMax)
      # Increment probability by the sum of points
      #  falling in the risk region
      prob[i] <- prob[i] +
        sum(apply(newY[condIsMax, ], 1, function(x) {
          isTRUE(all(x > region[-i]))
        }))
    }
    # Compute the Monte Carlo average of conditional probability
    # multiplied by margtransf probability of exceedance (Y0)
    prob[i] <- prob[i] / nsim * (1 - thresh[i])
  }
  pred <- ifelse(
    constraint == "equiprob",
    mean(prob),
    sum(wmax / sum(wmax) * prob)
  )
  return(list(
    pred = pred,
    prob = prob,
    res = res,
    nexc = nexc
  ))
}


#' Pseudo maximum likelihood estimation for the conditional extremes model
#'
#' @description Numerical optimization of the Gaussian pseudo likelihood for the conditional extremes model of Heffernan and Tawn, with constraints.
#'
#' @export
fit.condex <- function(
  xdat,
  thresh,
  mthresh,
  index = 1L,
  ties.method = c(
    "average",
    "first",
    "last",
    "random",
    "max",
    "min"
  ),
  margtransf = c(
    "empirical",
    "none",
    "semipar"
  ),
  use = c("complete", "pairwise"),
  constraint = TRUE,
  ...
) {
  constraint <- isTRUE(constraint)
  margtransf <- match.arg(margtransf)
  ties.method <- match.arg(ties.method)
  use <- match.arg(use)
  if (is.data.frame(xdat)) {
    xdat <- as.matrix(xdat)
  }
  stopifnot(is.matrix(xdat))
  D <- ncol(xdat)
  index <- unique(as.integer(index))
  stopifnot(is.finite(index), isTRUE(all(index %in% 1:D)))
  if (use == "complete") {
    xdat <- na.omit(xdat)
  }
  if (margtransf == "empirical") {
    sdata <- qlaplace(
      apply(xdat, 2, function(x) {
        rank(x, ties.method = "random", na.last = "keep") /
          (sum(!is.na(x)) + 1L)
      })
    )
  } else if (margtransf == "semipar") {
    stopifnot(
      is.vector(mthresh),
      length(mthresh) == ncol(xdat)
    )
    sdata <- qlaplace(
      mev::spunif(xdat, thresh = rep(mthresh, ncol(xdat)))
    )
  } else {
    sdata <- xdat
  }
  # TODO the index above can be multiple indices
  # TODO figure out constraints are enforced or not, based on argument
  fn <- function(par, ...) {
    -condex.pseudoll(
      par = par,
      index = index,
      thresh = thresh,
      data = sdata,
      constrain = FALSE,
      ...
    )
  }
  opt <- Rsolnp::solnp(
    pars = c(0.2, 0.1),
    fun = fn,
    ineqfun = function(par, ...) {
      optim_constr_ht(
        par = par,
        constants = condex.dataconstr(
          par = par,
          thresh = thresh,
          data = sdata
        )
      )
    },
    ineqLB = 0,
    ineqUB = Inf,
    LB = c(-1 + 1e-10, -1e3),
    UB = c(rep(1 - 1e-10, 2)),
    control = list(trace = FALSE)
  )
  f <- fn(opt$par, savePars = TRUE)
  hessian <- try(numDeriv::hessian(func = fn, x = opt$par))
  vcov <- solve(hessian)
  std.err <- sqrt(diag(vcov))
  loc <- attr(f, "loc")
  scale <- attr(f, "scale")
  rownames(loc) <- rownames(scale) <- index
  list(
    estimate = opt$par,
    std.err = std.err,
    value = as.numeric(f),
    vcov = vcov,
    pars = opt$par,
    loc = loc,
    scale = scale
  )
}

# TODO simulate from conditional extremes model above


# predict method with importance sampling and algorithm from paper
