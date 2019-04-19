
#' Simulation from R-Pareto processes
#'
#' @details For \code{riskf=max} and \code{riskf=min}, the procedure uses rejection sampling based on Pareto variates
#' sampled from \code{sum} and may be slow if \code{d} is large.
#'
#' @inheritParams rmev
#' @param shape shape tail index of Pareto variable
#' @param riskf string indicating risk functional.
#' @param siteindex integer between 1 and d specifying the index of the site or variable
#' @return an \code{n} by \code{d} sample from the R-Pareto process, with \code{attributes}
#' \code{accept.rate} if the procedure uses rejection sampling.
#' @export
#' @examples
#' rparp(n=10, riskf = 'site', siteindex=2, d=3, param=2.5, model='log')
#' rparp(n=10, riskf = 'min', d=3, param=2.5, model='neglog')
#' rparp(n=10, riskf = 'max', d=4, param=c(0.2,0.1,0.9,0.5), model='bilog')
#' rparp(n=10, riskf = 'sum', d=3, param=c(0.8,1.2,0.6, -0.5), model='sdir')
#' vario <- function(x, scale=0.5, alpha=0.8){ scale*x^alpha }
#' grid.loc <- as.matrix(expand.grid(runif(4), runif(4)))
#' rparp(n=10, riskf = 'max', vario=vario,loc=grid.loc, model='br')
rparp <- function(n, shape = 1, riskf = c("sum", "site", "max", "min", "l2"), siteindex = NULL, d, param, sigma, model = c("log", "neglog", "bilog", "negbilog", "hr", "br", "xstud", "smith", "schlather", "ct", "sdir", "dirmix"), weights, vario, loc, ...) {
  stopifnot(shape > 0)
  riskf <- match.arg(riskf)
  if (is.null(siteindex) && riskf == "site") {
    stop("For exceedances of site, the user needs to provide an index between 1 and d")
  }
  # Body of rmevspec
  models <- c("log", "neglog", "bilog", "negbilog", "hr", "br", "xstud", "smith", "schlather", "ct", "sdir", "dirmix", "negdir",
              "dir")
  model <- match.arg(model, models)[1]
  if (model == "schlather") {
    if (!missing(param))
      warning("Parameter value (degrees of freedom) set to one for Schlather model")
    param <- 1
    model <- "xstud"
  }
  # Define model families
  m1 <- c("log", "neglog")
  m2 <- c("bilog", "negbilog")
  m3 <- c("br", "xstud", "smith", "isbr")
  m4 <- c("ct", "dir", "negdir", "sdir")

  # Sanity checks
  if (model %in% c(m1, m2, m4) && (!missing(param) && mode(param) != "numeric")) {
    stop("Invalid parameter")
  }

  if (model %in% m1) {
    d <- as.integer(d)
    sigma = cbind(0)
    if (missing(param) || param < 0 || d < 1) {
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
  } else if (model %in% m2) {
    d <- as.integer(d)
    sigma = cbind(0)
    if (missing(param) || length(param) != d)
      stop("Invalid parameter value")
    # Check whether arguments are valid
    if (model == "bilog" && all(param >= 1)) {
      param <- 1/param
    }
    if (model == "negbilog" && all(param >= 0)) {
      param <- -param
    }
    if (any(param > 1))
      stop("Invalid param vector for bilogistic or negative bilogistic")
    if (any(param < 0) && model == "bilog")
      warning("Negative parameter values in bilogistic")
    if (any(param > 0) && model == "negbilog")
      warning("Positive parameter values in negative bilogistic")
  } else if (model %in% m4) {
    sigma = cbind(0)
    if (missing(param)) {
      stop("Invalid parameter value")
    }
    if (model == "ct") {
      if (length(param) != d) {
        if (length(param) == (d + 1)) {
          warning("Use `sdir' model for the scaled extremal Dirichlet model.")
          model = "sdir"
        } else {
          stop("Invalid arguments for the Coles and Tawn (extremal Dirichlet) model.")
        }
      }
      if (isTRUE(any(param < 0))) {
        stop("Invalid arguments for the Coles and Tawn (extremal Dirichlet) model.")
      }
    }
    if (model != "ct") {
      if (length(param) != (d + 1)) {
        stop("Invalid arguments for the Coles and Tawn (extremal Dirichlet) model.")
      }
      if (model == "negdir" && param[d + 1] > 0) {
        param[d + 1] <- -param[d + 1]
      }
      if (param[d + 1] < 0 && param[d + 1] <= -min(param[-(d + 1)])) {
        stop("Invalid parameters for the scaled Dirichlet. rho must be greater than -min(alpha)")
      }
      if (isTRUE(any(param[-(d + 1)] < 0))) {
        stop("Invalid arguments for the scaled Dirichlet model - alpha must be positive.")
      }
    }
    model = "sdir"
  } else if (model %in% m3) {
    # Smith, Brown-Resnick, extremal student
    if (model == "br") {
      if (missing(sigma) && !missing(vario) && !missing(loc)) {
        if (is.vector(loc))
          loc <- matrix(loc, ncol = 1)  #1 dimensional process
        stopifnot(is.function(vario))
        if (model == "br") {
          model = "isbr"
          m3 <- c(m3, model)
          if (vario(0, ...) > 1e-15) {
            stop("Cannot have a nugget term in the variogram for the Brown-Resnick process")
          }
          semivario2mat <- function(loc, semivario, ...) {
            di <- as.matrix(dist(loc))  #fields::rdist(loc) is faster...
            covmat <- matrix(0, nrow = nrow(di), ncol = ncol(di))
            covmat[lower.tri(covmat)] <- semivario(di[lower.tri(di)], ...)
            covmat[upper.tri(covmat)] <- t(covmat)[upper.tri(covmat)]
            return(covmat)
          }
          sigma <- semivario2mat(loc, vario, ...)/2
          # changed 14-05-2018 Matrix is half of Semivariogram, quarter of variogram
        }
      }
    }
    if (model != "isbr") {
      if (missing(sigma) || ncol(sigma) != nrow(sigma))
        stop("Invalid covariance matrix")
      if (any(diag(sigma) <= 0))
        stop("Degenerate covariance matrix; negative or zero entries found")
    }
    if (model == "xstud" && any(diag(sigma) != 1)) {
      warning("Extremal student requires correlation matrix")
      sigma <- cov2cor(sigma)
    }
    if (model == "xstud" && (missing(param) || length(param) != 1)) {
      stop("Degrees of freedom argument missing or invalid")
    }
    if (model == "smith" && missing(loc))
      stop("Location should be provided for the Smith model")
    if (model == "smith" && ncol(as.matrix(loc)) != ncol(sigma)) {
      stop("Covariance matrix of the Smith model should be
           of the same dimension as dimension of location vector")
    }
    d <- switch(model, xstud = ncol(sigma), br = ncol(sigma), smith = nrow(loc), isbr = ncol(sigma))
    if (model %in% c("smith", "br", "isbr")) {
      param <- 0
    }
    } else if (model == "dirmix") {
      if (any(missing(param), length(weights) != ncol(param) && ncol(param) != 1, any(param < 0))) {
        stop("Invalid arguments for the Dirichlet mixture")
      }
      if (!missing(weights)) {
        if (any(weights < 0))
          stop("Negative weights provided")
        if (sum(weights) != 1)
          warning("weights do not sum to one")
        weights <- weights/sum(weights)
      }
      if (missing(d)) {
        d <- nrow(param)
      } else if (d != nrow(param)) {
        stop("Dimension of d and provided param do not match")
      }
      # Checking for the mean constraints
      mar_mean <- colSums(t(param)/colSums(param) * weights)
      if (!isTRUE(all.equal(mar_mean, rep(1/d, d), tolerance = .Machine$double.eps^0.5))) {
        stop("Invalid mixture components")
      }
      # Switching parameters around to pass them to Rcpp function
      sigma <- param
      param <- weights
    } else if (model == "hr") {
      param = 0
      d <- ncol(sigma)
    }
  if (!model == "smith") {
    loc <- cbind(0)
  }
  # Model
  mod <- switch(model, log = 1, neglog = 2, dirmix = 3, bilog = 4, negbilog = 4, xstud = 5, br = 6, sdir = 7, smith = 8, hr = 9,
                isbr = 9)
  if (riskf == "sum") {
    # Generate from spectral measure
    return(evd::rgpd(n = n, loc = 1, scale = 1, shape = shape) * .rmevspec_cpp(n = n, d = d, para = param, model = mod, Sigma = sigma,
                                                                               loc = loc))
  } else if (riskf == "site") {
    # Check now that siteindex corresponds to a particular site
    # Dimension d could have been modified earlier for spatial models
    siteindex <- as.integer(siteindex)
    if (siteindex < 1 || siteindex > d) {
      stop("Invalid site index")
    }
    return(evd::rgpd(n = n, loc = 1, scale = 1, shape = shape) * .rPsite(n = n, j = siteindex, d = d, para = param, model = mod,
                                                                         Sigma = sigma, loc = loc))
  } else if (riskf %in% c("max", "min", "l2")) {
    ustar <- switch(riskf, max = 1, min = d, l2 = 1)
    ind <- 0L
    ntotsim <- 0L
    ntotacc <- 0L
    nsim <- ceiling(ifelse(n < 10, 4 * n, n))
    samp <- matrix(0, nrow = n, ncol = d)
    while (ind < n) {
      candidate <- evd::rgpd(n = nsim, loc = 1, scale = 1, shape = shape) * .rmevspec_cpp(n = nsim, d = d, para = param, model = mod,
                                                                                          Sigma = sigma, loc = loc)/ustar
      accept <- switch(riskf,
                       max = apply(candidate, 1, function(x) { max(x) > 1 }),
                       min = apply(candidate, 1, function(x) { min(x) > 1 }),
                       l2 = apply(candidate, 1, function(x) { sum(x^2) > 1 }))
      sum_accept <- sum(accept)
      ntotacc <- ntotacc + sum_accept
      ntotsim <- ntotsim + nsim
      if (sum_accept > 0) {
        if (sum_accept < (n - ind)) {
          samp[(ind + 1L):(ind + sum_accept), ] <- candidate[accept, ]
          ind <- ind + sum_accept
          nsim <- min(1e+06, ceiling(1.25 * (nsim/sum_accept) * (n - ind)))
        } else {
          samp[(ind + 1L):n, ] <- as.matrix(candidate[accept, ])[1:(n - ind), ]
          ind <- n
        }
      } else {
        nsim <- min(1e+06, ceiling(1.25 * nsim))
      }
    }
    attr(samp, "accept.rate") <- ntotacc/ntotsim
    return(samp)
  } else {
    stop("Model not implemented")
  }
}



#' Simulation from generalized R-Pareto processes
#'
#' The generalized R-Pareto process is supported on \code{(B - A / shape, Inf)} if \code{shape > 0},
#' or \code{(-Inf, B - A / shape)} for negative shape parameters, conditional on \eqn{(X-r(B))/r(A)>0}.
#' The standard Pareto process corresponds to \code{A = B = rep(1, d)}.
#'
#'
#' @inheritParams rmev
#' @param shape shape parameter of the generalized Pareto variable
#' @param riskf string indicating the risk functional.
#' @param thresh univariate threshold for the exceedances of risk functional
#' @param siteindex integer between 1 and d specifying the index of the site or variable
#' @param A scale vector
#' @param B location vector
#' @return an \code{n} by \code{d} sample from the generalized R-Pareto process, with \code{attributes}
#' \code{accept.rate} if the procedure uses rejection sampling.
#' @export
#' @examples
#' rgparp(n = 10, riskf = 'site', siteindex = 2, d = 3, param = 2.5,
#'    model = 'log', A = c(1, 2, 3), B = c(2, 3, 4))
#' rgparp(n = 10, riskf = 'max', d = 4, param = c(0.2, 0.1, 0.9, 0.5),
#'    A = 1:4, B = 1:4, model = 'bilog')
#' rgparp(n = 10, riskf = 'sum', d = 3, param = c(0.8, 1.2, 0.6, -0.5),
#'    A = 1:3, B = 1:3, model = 'sdir')
#' vario <- function(x, scale = 0.5, alpha = 0.8){ scale*x^alpha }
#' grid.loc <- as.matrix(expand.grid(runif(4), runif(4)))
#' rgparp(n = 10, riskf = 'max', vario = vario,loc = grid.loc,
#'    model = 'br', A = runif(16), B = rnorm(16))
rgparp <- function(n, shape = 1, thresh = 1, riskf = c("mean", "sum", "site", "max", "min", "l2"), siteindex = NULL, d, A, B, param, sigma, model = c("log", "neglog","bilog", "negbilog", "hr", "br", "xstud", "smith", "schlather", "ct", "sdir", "dirmix"), weights, vario, loc, ...) {
  riskf <- match.arg(riskf)
  #shape <- as.vector(shape[1])
  #stopifnot(is.numeric(shape))
  if (is.null(siteindex) && riskf == "site") {
    stop("For exceedances of site, the user needs to provide an index between 1 and d")
  }
  # Body of rmevspec
  models <- c("log", "neglog", "bilog", "negbilog", "hr", "br", "xstud", "smith", "schlather", "ct", "sdir", "dirmix", "negdir",
              "dir")
  model <- match.arg(model, models)[1]
  if (model == "schlather") {
    if (!missing(param))
      warning("Parameter value (degrees of freedom) set to one for Schlather model")
    param <- 1
    model <- "xstud"
  }
  # Define model families
  m1 <- c("log", "neglog")
  m2 <- c("bilog", "negbilog")
  m3 <- c("br", "xstud", "smith", "isbr")
  m4 <- c("ct", "dir", "negdir", "sdir")

  # Sanity checks
  if (model %in% c(m1, m2, m4) && (!missing(param) && mode(param) != "numeric")) {
    stop("Invalid parameter")
  }

  if (model %in% m1) {
    d <- as.integer(d)
    sigma = cbind(0)
    if (missing(param) || param < 0 || d < 1) {
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
  } else if (model %in% m2) {
    d <- as.integer(d)
    sigma = cbind(0)
    if (missing(param) || length(param) != d)
      stop("Invalid parameter value")
    # Check whether arguments are valid
    if (model == "bilog" && all(param >= 1)) {
      param <- 1/param
    }
    if (model == "negbilog" && all(param >= 0)) {
      param <- -param
    }
    if (any(param > 1))
      stop("Invalid param vector for bilogistic or negative bilogistic")
    if (any(param < 0) && model == "bilog")
      warning("Negative parameter values in bilogistic")
    if (any(param > 0) && model == "negbilog")
      warning("Positive parameter values in negative bilogistic")
  } else if (model %in% m4) {
    sigma = cbind(0)
    if (missing(param)) {
      stop("Invalid parameter value")
    }
    if (model == "ct") {
      if (length(param) != d) {
        if (length(param) == (d + 1)) {
          warning("Use `sdir' model for the scaled extremal Dirichlet model.")
          model = "sdir"
        } else {
          stop("Invalid arguments for the Coles and Tawn (extremal Dirichlet) model.")
        }
      }
      if (isTRUE(any(param < 0))) {
        stop("Invalid arguments for the Coles and Tawn (extremal Dirichlet) model.")
      }
    }
    if (model != "ct") {
      if (length(param) != (d + 1)) {
        stop("Invalid arguments for the Coles and Tawn (extremal Dirichlet) model.")
      }
      if (model == "negdir" && param[d + 1] > 0) {
        param[d + 1] <- -param[d + 1]
      }
      if (param[d + 1] < 0 && param[d + 1] <= -min(param[-(d + 1)])) {
        stop("Invalid parameters for the scaled Dirichlet. rho must be greater than -min(alpha)")
      }
      if (isTRUE(any(param[-(d + 1)] < 0))) {
        stop("Invalid arguments for the scaled Dirichlet model - alpha must be positive.")
      }
    }
    model = "sdir"
  } else if (model %in% m3) {
    # Smith, Brown-Resnick, extremal student
    if (model == "br") {
      if (missing(sigma) && !missing(vario) && !missing(loc)) {
        if (is.vector(loc))
          loc <- matrix(loc, ncol = 1)  #1 dimensional process
        stopifnot(is.function(vario))
        if (model == "br") {
          model = "isbr"
          m3 <- c(m3, model)
          if (vario(0, ...) > 1e-15) {
            stop("Cannot have a nugget term in the variogram for the Brown-Resnick process")
          }
          semivario2mat <- function(loc, semivario, ...) {
            di <- as.matrix(dist(loc))  #fields::rdist(loc) is faster...
            covmat <- matrix(0, nrow = nrow(di), ncol = ncol(di))
            covmat[lower.tri(covmat)] <- semivario(di[lower.tri(di)], ...)
            covmat[upper.tri(covmat)] <- t(covmat)[upper.tri(covmat)]
            return(covmat)
          }
          sigma <- semivario2mat(loc, vario, ...)/2
          # changed 14-05-2018 Matrix is half of Semivariogram, quarter of variogram
        }
      }
    }
    if (model != "isbr") {
      if (missing(sigma) || ncol(sigma) != nrow(sigma))
        stop("Invalid covariance matrix")
      if (any(diag(sigma) <= 0))
        stop("Degenerate covariance matrix; negative or zero entries found")
    }
    if (model == "xstud" && any(diag(sigma) != 1)) {
      warning("Extremal student requires correlation matrix")
      sigma <- cov2cor(sigma)
    }
    if (model == "xstud" && (missing(param) || length(param) != 1)) {
      stop("Degrees of freedom argument missing or invalid")
    }
    if (model == "smith" && missing(loc))
      stop("Location should be provided for the Smith model")
    if (model == "smith" && ncol(as.matrix(loc)) != ncol(sigma)) {
      stop("Covariance matrix of the Smith model should be
           of the same dimension as dimension of location vector")
    }
    d <- switch(model, xstud = ncol(sigma), br = ncol(sigma), smith = nrow(loc), isbr = ncol(sigma))
    if (model %in% c("smith", "br", "isbr")) {
      param <- 0
    }
    } else if (model == "dirmix") {
      if (any(missing(param), length(weights) != ncol(param) && ncol(param) != 1, any(param < 0))) {
        stop("Invalid arguments for the Dirichlet mixture")
      }
      if (!missing(weights)) {
        if (any(weights < 0))
          stop("Negative weights provided")
        if (sum(weights) != 1)
          warning("weights do not sum to one")
        weights <- weights/sum(weights)
      }
      if (missing(d)) {
        d <- nrow(param)
      } else if (d != nrow(param)) {
        stop("Dimension of d and provided param do not match")
      }
      # Checking for the mean constraints
      mar_mean <- colSums(t(param)/colSums(param) * weights)
      if (!isTRUE(all.equal(mar_mean, rep(1/d, d), tolerance = .Machine$double.eps^0.5))) {
        stop("Invalid mixture components")
      }
      # Switching parameters around to pass them to Rcpp function
      sigma <- param
      param <- weights
    } else if (model == "hr") {
      param = 0
      d <- ncol(sigma)
    }
  if (!model == "smith") {
    loc <- cbind(0)
  }
  # Model
  mod <- switch(model, log = 1, neglog = 2, dirmix = 3, bilog = 4, negbilog = 4, xstud = 5, br = 6, sdir = 7, smith = 8, hr = 9,
                isbr = 9)
  # Additional checks and arguments for accept-reject algorithm for generalized R-Pareto process
  # Scale vector
  if(missing(A)){
    stop("Missing scale function")
  } else {
    if(length(A) == 1){
      A <- rep(A, length.out = d)
    }
    stopifnot(length(A) == d, all(A > 0))
  }
  if(missing(B)){
    stop("Missing location function")
  } else{
    if(length(B) == 1){
      B <- rep(B, length.out = d)
    } else{
   stopifnot(length(B) == d)
    }
  }

  if (riskf == "site") {
    # Check that siteindex corresponds to a particular site
    # Dimension d could have been modified earlier for spatial models
    siteindex <- as.integer(siteindex)
    if (siteindex < 1 || siteindex > d) {
      stop("Invalid site index")
    }
  }
  #rB <- switch(riskf, sum = sum(B), max = max(B), min = min(B), site = B[siteindex], l2 = sqrt(sum(B^2)))
  #rA <- switch(riskf, sum = sum(A), max = max(A), min = min(A), site = A[siteindex], l2 = sqrt(sum(A^2)))
  #B <- B - rB #identifiability constraint r(B) = 0, r(B) is the threshold
  #A <- A / rA #identifiability constraint r(A) = 1, r(A) is scale of GP
  # Process becomes after simulation with A, B return rA*X + rB
  #Compute threshold for the l1 norm
  shape <- rep(shape, length.out = d)
  us <- thresh
  stopifnot(length(us) == 1L)
  if(riskf %in% c("max", "l2", "sum", "mean")){
    if(isTRUE(all(ifelse(shape < 0, -A/shape+B < us, FALSE)))){
     stop("Invalid input: the threshold selected is above the upper endpoint of all of the marginal distributions.")
    }
  } else if(riskf == "min"){
    if(isTRUE(any(ifelse(shape < 0, -A/shape+B < us, FALSE)))){
      stop("Invalid input: the threshold selected is above the upper endpoint of some the marginal distributions.")
    }
  } else if(riskf == "site"){
   if(shape[siteindex] < 0 && -A[siteindex]/shape[siteindex]+B[siteindex] < us){
     stop("Invalid input: the threshold selected is above the upper endpoint of the marginal distribution at the selected site.")
   }
  }
  if(riskf %in% c("max", "l2")){
    ustar <- min(ifelse(sapply(shape, function(xi){isTRUE(all.equal(xi, 0))}),
                    exp(( us - B) / A),
                    (1 + shape * (us - B) / A)^(1/shape)), na.rm = TRUE)
  } else if(riskf %in% c("sum", "mean", "l2")){ #no bound for l2, contained in l1
    if(riskf == "mean"){
      riskf <- "sum"
      us <- d*us
    }
    ustar <- 0
    zeroshape <- sapply(shape, function(xi){isTRUE(all.equal(xi, 0))})
    if(us - sum(B) < 0){
      ustar <- d
    } else if(all(shape <0)){
    inter <- 1-(us-sum(B))/sum(A)*min(abs(shape))
      if(inter > 0){
      ustar <- d*inter^(-1/max(abs(shape)))
      }
    } else if(all(shape > 0)){
      ustar <- (min(shape)*(us-sum(B))/sum(A)+1)^(1/max(shape))
    } else if (all(zeroshape)){
      ustar <- D*(exp((us-sum(B))/sum(A))^(1/d))
    } else {
      xiprime <- ifelse(zeroshape, 1, abs(shape))
      ustar <- (min(xiprime)*(us-sum(B))/sum(A)+1)^(1/max(xiprime))/2
    }
    ustar <- max(d, ustar)
    } else if(riskf == "min"){
    ustar <- sum(ifelse(sapply(shape, function(xi){isTRUE(all.equal(xi, 0))}),
                    exp(( us - B) / A),
                    (1 + shape * (us - B) / A)^(1/shape)))
  } else if(riskf == "site"){
    #for this, simulate directly from angular measure P0
    ustar <- (1+shape*(us-B[siteindex])/A[siteindex])^(1/shape)
  } else if(riskf == "mean"){ # difference vs max is 1/d factor
    ustar <- min(ifelse(sapply(shape, function(xi){isTRUE(all.equal(xi, 0))}),
                        exp(( us - B) / A),
                        (1 + shape * (us - B) / A)^(1/shape)), na.rm = TRUE)/d
  }
  #Algorithm 1
  # Nonlinear risk functionals
  if (riskf %in% c("sum", "mean", "max", "min", "l2")) {
    ind <- 0L
    ntotsim <- 0L
    ntotacc <- 0L
    nsim <- ceiling(ifelse(n < 10, 10 * n, n))
    samp <- matrix(0, nrow = n, ncol = d)
    while (ind < n) {
      candidate <- ustar / runif(nsim) * .rmevspec_cpp(n = nsim, d = d, para = param, model = mod, Sigma = sigma, loc = loc)
        for(j in 1:d){
          if(!isTRUE(all.equal(shape[j], 0))){
          candidate[,j] <- (candidate[,j]^shape[j] - 1) / shape[j] * A[j] + B[j]
        } else{
          candidate[,j]  <- A[j] * log(candidate[,j]) + B[j]
        }
      }
      accept <- switch(riskf,
                       site = apply(candidate, 1, function(x){ x[siteindex] > us}),
                       max = apply(candidate, 1, function(x) { max(x) > us}),
                       min = apply(candidate, 1, function(x) { min(x) > us}),
                       l2 =  apply(candidate, 1, function(x) { sum(x^2) > us}),
                       sum = apply(candidate, 1, function(x) { sum(x) > us}),
                       mean = apply(candidate, 1, function(x) { mean(x) > us}))
      sum_accept <- sum(accept)
      ntotacc <- ntotacc + sum_accept
      ntotsim <- ntotsim + nsim
      if (sum_accept > 0) {
        if (sum_accept < (n - ind)) {
          samp[(ind + 1L):(ind + sum_accept), ] <- candidate[accept, , drop = FALSE]
          ind <- ind + sum_accept
          nsim <- min(1e+06, ceiling(1.25 * (nsim/sum_accept) * (n - ind)))
        } else {
          samp[(ind + 1L):n, ] <- candidate[accept, ,drop = FALSE][1:(n - ind), , drop = FALSE]
          ind <- n
        }
      } else {
        nsim <- min(1e+06, ceiling(1.25 * nsim))
      }
    }
    #samp <- rA * samp + rB
    attr(samp, "accept.rate") <- ntotacc/ntotsim
    return(samp)
  } else {
     #Acceptance rate is 1
     samp <- ustar / runif(n) * .rPsite(n = n, j = siteindex, d = d, para = param, model = mod, Sigma = sigma, loc = loc)
     for(j in 1:d){
     if(!isTRUE(all.equal(shape[j], 0))){
           samp[,j] <- (samp[,j]^shape[j] - 1) / shape[j] * A[j] +  B[j]
     } else{
        samp[,j]  <- A[j] * log(samp[,j]) + B[j]
       }
     }
    #samp <- rA * samp + rB
    attr(samp, "accept.rate") <- 1
    return(samp)
  }
}
