#' Exact simulations of multivariate extreme value distributions
#'
#' Implementation of the random number generators for multivariate extreme-value distributions
#' and max-stable processes based on the two algorithms described in Dombry, Engelke and Oesting (2016).
#'
#' @param n number of observations
#' @param d dimension of sample
#' @param param parameter vector for the logistic, bilogistic, negative bilogistic and extremal Dirichlet (Coles and Tawn) model.
#' Parameter matrix for the Dirichlet mixture. Degree of freedoms for extremal student model. See \bold{Details}.
#' @param sigma covariance matrix for Brown-Resnick and extremal Student-t distributions. Symmetric matrix of squared  coefficients \eqn{\lambda^2} for the Husler-Reiss model, with zero diagonal elements.
#' @param asy list of asymmetry parameters, as in function \code{rmvevd} from package \code{evd}, of \eqn{2^d-1} vectors of size corresponding to the power set of \code{d}, with sum to one constraints.
#' @param alg algorithm, either simulation via extremal function (\code{'ef'}) or via the spectral measure (\code{'sm'}). Default to \code{ef}.
#' @param model for multivariate extreme value distributions, users can choose between 1-parameter logistic and negative logistic, asymmetric logistic and negative logistic, bilogistic, Husler-Reiss, extremal Dirichlet model (Coles and Tawn) or the Dirichlet mixture. Spatial models include
#' the Brown-Resnick, Smith, Schlather and extremal Student max-stable processes. Max linear models are also supported
#' @param vario semivariogram function whose first argument must be distance. Used only if provided in conjunction with \code{coord} and if \code{sigma} is missing
#' @param coord \code{d} by \code{k} matrix of coordinates, used as input in the variogram \code{vario} or as parameter for the Smith model. If \code{grid} is \code{TRUE}, unique entries should be supplied.
#' @param weights vector of length \code{m} for the \code{m} mixture components that sum to one. For the \code{"maxlin"} model, weights should be a matrix with \code{d} columns that represent the weight of the components and whose column sum to one (if provided, this argument overrides \code{asy}).
#' @param grid Logical. \code{TRUE} if the coordinates are two-dimensional grid points (spatial models).
#' @param dist symmetric matrix of pairwise distances. Default to \code{NULL}.
#' @param ... additional arguments for the \code{vario} function
#' @author Leo Belzile
#' @details The vector param differs depending on the model
#' \itemize{
#'  \item \code{log}: one dimensional parameter greater than 1
#'  \item \code{alog}: \eqn{2^d-d-1} dimensional parameter for \code{dep}. Values are recycled if needed.
#'  \item \code{neglog}: one dimensional positive parameter
#'  \item \code{aneglog}: \eqn{2^d-d-1} dimensional parameter for \code{dep}. Values are recycled if needed.
#'  \item \code{bilog}: \code{d}-dimensional vector of parameters in \eqn{[0,1]}
#'  \item \code{negbilog}: \code{d}-dimensional vector of negative parameters
#'  \item \code{ct, dir, negdir, sdir}: \code{d}-dimensional vector of positive (a)symmetry parameters. For \code{dir} and \code{negdir}, a \eqn{d+1}
#'  vector consisting of the \code{d} Dirichlet parameters and the last entry is an index of regular variation in \eqn{(-\min(\alpha_1, \ldots, \alpha_d), 1]} treated as shape parameter
#'  \item \code{xstud}: one dimensional parameter corresponding to degrees of freedom \code{alpha}
#'  \item \code{dirmix}: \code{d} by \code{m}-dimensional matrix of positive (a)symmetry parameters
#'  \item \code{pairbeta, pairexp}: \code{d(d-1)/2+1} vector of parameters, containing the concentration parameter and the coefficients of the pairwise beta, in lexicographical order e.g., \eqn{\beta_{12}, \beta_{13}, \ldots}
#'  \item \code{wdirbs, wexpbs}: \code{2d} vector of \code{d} concentration parameters followed by the \code{d} Dirichlet parameters
#' }
#'
#' Stephenson points out that the multivariate asymmetric negative logistic model given in e.g. Coles and Tawn (1991) is not a valid distribution function in dimension \eqn{d>3} unless additional constraints are imposed on the parameter values.
#' The implementation in \code{mev} uses the same construction as the asymmetric logistic distribution (see the vignette). As such it does not match the bivariate implementation of \link[evd]{rbvevd}.
#'
#' The dependence parameter of the \code{evd} package for the Husler-Reiss distribution can be recovered taking
#' for the Brown--Resnick model  \eqn{2/r=\sqrt(2\gamma(h))} where \eqn{h} is the lag vector between sites and \eqn{r=1/\lambda} for the Husler--Reiss.
#'
#' @section Warning:
#'As of version 1.8 (August 16, 2016), there is a distinction between models \code{hr} and \code{br}. The latter is meant to be used in conjunction with variograms. The parametrization differs between the two models.
#'
#'
#'The family of scaled Dirichlet is now parametrized by a parameter in \eqn{-\min(\alpha)}{-min(\alpha)} appended to the the \code{d} vector \code{param} containing the parameter \code{alpha}
#'of the Dirichlet model. Arguments \code{model='dir'} and \code{model='negdir'} are still supported internally, but not listed in the options.
#'
#' @return an \code{n} by \code{d} exact sample from the corresponding multivariate extreme value model
#'
#' @export
#' @references
#'Dombry, Engelke and Oesting (2016). Exact simulation of max-stable processes, \emph{Biometrika}, \bold{103}(2), 303--317.
#' @seealso \link{rmevspec}, \link[evd]{rmvevd}, \link[evd]{rbvevd}
#' @examples
#'set.seed(1)
#'rmev(n=100, d=3, param=2.5, model='log', alg='ef')
#'rmev(n=100, d=4, param=c(0.2,0.1,0.9,0.5), model='bilog', alg='sm')
#'## Spatial example using power variogram
#'#NEW: Semi-variogram must take distance as argument
#'semivario <- function(x, scale, alpha){ scale*x^alpha }
#'#grid specification
#'grid.coord <- as.matrix(expand.grid(runif(4), runif(4)))
#'rmev(n=100, vario=semivario, coord=grid.coord, model='br', scale = 0.5, alpha = 1)
#'#using the Brown-Resnick model with a covariance matrix
#' vario2cov <- function(coord, semivario,...){
#'  sapply(1:nrow(coord), function(i) sapply(1:nrow(coord), function(j)
#'   semivario(sqrt(sum((coord[i,])^2)), ...) +
#'   semivario(sqrt(sum((coord[j,])^2)), ...) -
#'   semivario(sqrt(sum((coord[i,]-coord[j,])^2)), ...)))
#' }
#' rmev(n=100, sigma=vario2cov(grid.coord, semivario = semivario, scale = 0.5, alpha = 1), model='br')
#' # asymmetric logistic model - see function 'rmvevd' from package 'evd '
#' asy <- list(0, 0, 0, 0, c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0),
#'   c(.2,.1,.2), c(.1,.1,.2), c(.3,.4,.1), c(.2,.2,.2), c(.4,.6,.2,.5))
#' rmev(n=1, d=4, param=0.3, asy=asy, model="alog")
#'#Example with a grid (generating an array)
#'rmev(n=10, sigma=cbind(c(2,1), c(1,3)), coord=cbind(runif(4), runif(4)), model='smith', grid=TRUE)
#'## Example with Dirichlet mixture
#'alpha.mat <- cbind(c(2,1,1),c(1,2,1),c(1,1,2))
#'rmev(n=100, param=alpha.mat, weights=rep(1/3,3), model='dirmix')
#'rmev(n=10, param=c(0.1,1,2,3), d=3, model='pairbeta')
rmev <- function(
    n,
    d,
    param,
    asy,
    sigma,
    model = c("log",
              "alog",
              "neglog",
              "aneglog",
              "bilog",
              "negbilog",
              "hr",
              "br",
              "xstud",
              "smith",
              "schlather",
              "ct",
              "sdir",
              "dirmix",
              "pairbeta",
              "pairexp",
              "wdirbs",
              "wexpbs",
              "maxlin"),
    alg = c("ef", "sm"),
    weights = NULL,
    vario = NULL,
    coord = NULL,
    grid = FALSE,
    dist = NULL,
    ...) {
 model <- match.arg(model)
  # Create gridded values if specification is for random field discretization
  ellips <- list(...)
  if(is.null(coord) && !is.null(ellips$loc)){
    coord <- ellips$loc
  }
  if(missing(d)){ d <- NULL }
  if(missing(param)){ param <- NULL }
  if(missing(asy)){ asy <- NULL }
  if(missing(sigma)){ sigma <- NULL }
  out <- .rmev_checks(
    n = n,
    d = d,
    param = param,
    asy = asy,
    sigma = sigma,
    model =  model,
    alg = alg,
    weights = weights,
    vario = vario,
    coord = coord,
    grid = grid,
    dist = dist,
    ...)
    mod <- switch(out$model,
                  log = 1,
                  neglog = 2,
                  dirmix = 3,
                  bilog = 4,
                  negbilog = 4,
                  xstud = 5,
                  br = 6,
                  sdir = 7,
                  smith = 8,
                  hr = 9,
                  isbr = 9,
                  pairbeta = 10,
                  pairexp = 11,
                  wdirbs = 12,
                  wexpbs = 13)
    if (out$model %in% c("alog", "aneglog","maxlin")) {
    mod <- switch(out$model,
                  alog = 1,
                  aneglog = 2,
                  maxlin = 14)
    .rmevasy(n = n, d = out$d, par = out$param,
             asym = out$asym, ncompo = out$ncompo,
             Sigma = out$sigma, model = mod)
    } else {
        if (model != "smith") {
            coordat <- cbind(0)
        } else {
            coordat <- out$coord
        }
        if (out$model %in% c("br", "xstud", "smith", "isbr") && grid == TRUE) {
            npdim <- out$d^(1/ncol(out$coord))
            if (!all.equal(npdim, as.integer(npdim))) {
                stop("The dimension of the input grid does not match (square) covariance matrix")
            }
            array(t(switch(out$alg,
                           ef = .rmevA2(n = n, d = out$d, par = out$param, model = mod, Sigma = out$sigma, coordat),
                           sm = .rmevA1(n = n, d = out$d, par = out$param, model = mod, Sigma = out$sigma, coordat))), dim = c(rep(npdim, ncol(coord)), n))
        } else {
            switch(out$alg,
                   ef = .rmevA2(n = n, d = out$d, par = out$param, model = mod, Sigma = out$sigma, coordat),
                   sm = .rmevA1(n = n, d = out$d, par = out$param, model = mod, Sigma = out$sigma, coordat))
        }
    }
}




#' Random samples from spectral distributions of multivariate extreme value models.
#'
#' Generate from \eqn{Q_i}{Qi}, the spectral measure of a given multivariate extreme value model based on the L1 norm.
#'
#' @section Note:
#'  This functionality can be useful to generate for example Pareto processes with marginal exceedances.
#'
#' @inheritParams rmev
#'
#' @author Leo Belzile
#' @details The vector param differs depending on the model
#' \itemize{
#'  \item \code{log}: one dimensional parameter greater than 1
#'  \item \code{neglog}: one dimensional positive parameter
#'  \item \code{bilog}: \code{d}-dimensional vector of parameters in \eqn{[0,1]}
#'  \item \code{negbilog}: \code{d}-dimensional vector of negative parameters
#'  \item \code{ct}, \code{dir}, \code{negdir}: \code{d}-dimensional vector of positive (a)symmetry parameters. Alternatively, a \eqn{d+1}
#'  vector consisting of the \code{d} Dirichlet parameters and the last entry is an index of regular variation in \eqn{(0, 1]} treated as scale
#'  \item \code{xstud}: one dimensional parameter corresponding to degrees of freedom \code{alpha}
#'  \item \code{dirmix}: \code{d} by \code{m}-dimensional matrix of positive (a)symmetry parameters
#'  \item \code{pairbeta, pairexp}: \code{d(d-1)/2+1} vector of parameters, containing the concentration parameter and the coefficients of the pairwise beta, in lexicographical order e.g., \eqn{\beta_{1,2}, \beta_{1,3}, \ldots}
#'  \item \code{wdirbs, wexpbs}: \code{2d} vector of \code{d} concentration parameters followed by the \code{d} Dirichlet parameters
#' }
#' @return an \code{n} by \code{d} exact sample from the corresponding multivariate extreme value model
#'
#' @references
#'Dombry, Engelke and Oesting (2016). Exact simulation of max-stable processes, \emph{Biometrika}, \bold{103}(2), 303--317.
#' @references Boldi (2009). A note on the representation of parametric models for multivariate extremes.
#' \emph{Extremes} \bold{12}, 211--218.
#'
#' @examples
#' set.seed(1)
#' rmevspec(n=100, d=3, param=2.5, model='log')
#' rmevspec(n=100, d=3, param=2.5, model='neglog')
#' rmevspec(n=100, d=4, param=c(0.2,0.1,0.9,0.5), model='bilog')
#' rmevspec(n=100, d=2, param=c(0.8,1.2), model='ct') #Dirichlet model
#' rmevspec(n=100, d=2, param=c(0.8,1.2,0.5), model='sdir') #with additional scale parameter
#'#Variogram gamma(h) = scale*||h||^alpha
#'#NEW: Variogram must take distance as argument
#'vario <- function(x, scale=0.5, alpha=0.8){ scale*x^alpha }
#' #grid specification
#' grid.coord <- as.matrix(expand.grid(runif(4), runif(4)))
#' rmevspec(n=100, vario=vario,coord=grid.coord, model='br')
#' ## Example with Dirichlet mixture
#' alpha.mat <- cbind(c(2,1,1),c(1,2,1),c(1,1,2))
#' rmevspec(n=100, param=alpha.mat, weights=rep(1/3,3), model='dirmix')
#' @export
rmevspec <- function(n, d, param, sigma, model = c("log", "neglog", "bilog", "negbilog", "hr", "br", "xstud", "smith", "schlather",
    "ct", "sdir", "dirmix","pairbeta","pairexp","wdirbs","wexpbs","maxlin"),
    weights = NULL, vario = NULL, coord = NULL, grid = FALSE, dist = NULL, ...) {
  # Dummy algorithm argument for internal checks
  alg <- "sm"
   model <- match.arg(model)
  if(model %in% c("alog","aneglog", "maxlik")){
    stop("Invalid model: cannot simulate from angular distribution of asymmetric models.")
  }
  #TODO add angular measure of the max linear model

  asy <- NULL
  if(missing(d)){ d <- NULL }
  if(missing(param)){ param <- NULL }
  if(missing(sigma)){ sigma <- NULL }
    out <- .rmev_checks(n = n, d = d, param = param, asy = asy, sigma = sigma,
                 model = model, alg = alg, weights = weights,
                 vario = vario, coord = coord, grid = grid, dist = dist,
                 ...)
    mod <- switch(out$model, log = 1, neglog = 2, dirmix = 3,
                  bilog = 4, negbilog = 4, xstud = 5, br = 6,
                  sdir = 7, smith = 8, hr = 9, isbr = 9,
                  pairbeta = 10, pairexp = 11, wdirbs = 12, wexpbs = 13)
    # Generate from spectral measure
    .rmevspec_cpp(n = n,
                  d = out$d,
                  par = out$param,
                  model = mod,
                  Sigma = out$sigma,
                  loc = out$coord)
}


#' Internal function
#'
#' Takes a list of asymmetry parameters with an associated dependence vector and returns
#' the corresponding asymmetry matrix for the asymmetric logistic and asymmetric negative logistic models
#'
#' This function is extracted from the evd package and modified
#' (C) Alec Stephenson
#'
#' @param asy a list of \eqn{2^d-1} asymmetry components, as in Stephenson bvevd functions
#' @param dep vector of \eqn{2^d-d-1} values for the dependence parameter
#' @param d dimension of the model
#' @param model, either \code{alog} for the asymmetric logistic or \code{aneglog}
#' for the asymmetric negative logistic
#'
#' @return a matrix of asymmetry components, enumerating all possible \eqn{2^d-1} subsets of
#' the power set
#' @keywords internal
.mvasym.check <- function(asy, dep, d, model = c("alog", "aneglog","maxlin")) {
  if(is.null(d)){
    stop("Argument \"d\" must be provided.")
  }
    # Function subset is an internal function from the evd package
    subsets <- function(d) {
        x <- 1:d
        k <- NULL
        for (m in x) k <- rbind(cbind(TRUE, k), cbind(FALSE, k))
        pset <- apply(k, 1, function(z) x[z])
        pset[sort.list(unlist(lapply(pset, length)))[-1]]
    }
    if (model == "alog") {
        if (mode(dep) != "numeric" || any(dep <= 0) || any(dep > 1))
            stop("invalid argument for \"dep\"")
    } else if(model == "aneglog"){
        if (mode(dep) != "numeric" || any(dep <= 0))
            stop("invalid argument for \"dep\"")
    }
    nb <- 2^d - 1
    if (mode(asy) != "list" || length(asy) != nb)
        stop(paste("\"asy\" should be a list of length", nb))
    tasy <- function(theta, b) {
        trans <- matrix(0, nrow = nb, ncol = d)
        for (i in 1:nb) trans[i, (1:d %in% b[[i]])] <- theta[[i]]
        trans
    }
    b <- subsets(d)
    if (any(sapply(asy, length) != sapply(b, length)))
        stop("\"asy\" is not of the correct form")
    asy <- tasy(asy, b)
    if (!is.matrix(asy) || mode(asy) != "numeric")
        stop("\"asy\" is not of the correct form")
    if (min(asy) < 0 || max(asy) > 1)
        stop("\"asy\" must contain parameters in [0,1]")
    if (any(apply(asy, 2, sum) != 1) || any(asy[c(rep(FALSE, d), dep == 1), ] != 0) || any(apply(asy[-(1:d), , drop = FALSE], 1, function(x) sum(x !=
        0)) == 1))
        stop("\"asy\" does not satisfy the appropriate constraints")
    asy
}

#' Is the matrix conditionally negative semi-definite?
#' Function adapted from 'is.CNSD' in the CEGO package, v 2.1.0
#' @author Martin Zaefferer
#' @param X a symmetric matrix
#' @param tol tolerance value; eigenvalues between \code{-tol} and \code{tol} are assumed to be zero.
#' @keywords internal
.is.CNSD <- function(X, tol = 1e-08) {
    n <- nrow(X)
    P <- diag(n)
    if (n > 2) {
        diag(P[, -1]) <- -1
    } else if (n == 2) {
        # error with one dimensional case...
        P[1, -1] <- -1
    }
    Xhat <- P %*% X %*% t(P)
    eigs <- eigen(Xhat[-n, -n], TRUE, TRUE)$values
    !eigs[1] > tol
}

#' @keywords internal
.rmev_checks <- function(n = n, d = d, param = param, asy = asy, sigma = sigma,
                   model =  model, alg = alg, weights = weights, vario = vario,
                   coord = coord, grid = grid, dist = dist, ...){
  models <- c("log", "alog", "neglog", "aneglog",
              "bilog", "negbilog", "hr", "br",
              "isbr", "xstud", "smith", "schlather",
              "ct", "sdir", "dirmix", "negdir", "dir",
              "pairbeta","pairexp","wdirbs","wexpbs","maxlin")
  alg <- match.arg(alg)
  asym <- NULL
  ncompo <- NULL
  model <- match.arg(model, models)[1]
  if (!model %in% c("alog", "aneglog","maxlin")) {
    if (!is.null(asy)) {
      warning("Asymmetry parameter ignored")
    } else {
      asym <- matrix(TRUE, ncol = 1, nrow = 1)
    }
  }
  if (model == "schlather") {
    if (!is.null(param) && !isTRUE(all.equal(param, 1, check.attributes = FALSE))){
      warning("Parameter value (degrees of freedom) set to one for Schlather model")
    }
    param <- 1
    model <- "xstud"
  }
  # Define model families
  m1 <- c("log", "neglog")
  m2 <- c("bilog", "negbilog")
  m3 <- c("br", "xstud", "smith", "isbr")
  m4 <- c("ct", "dir", "negdir", "sdir")
  m5 <- c("pairbeta","pairexp","wdirbs","wexpbs")
  # Check that parameters are provided
  if (model %in% c(m1, m2, m4, m5) && (!is.null(param) && mode(param) != "numeric")) {
    stop("Invalid parameter")
  }
  # Check spatial requirements
  if ((!is.null(coord) || !is.null(vario) || !is.null(dist)) && !model %in% m3) {
    if (model == "hr") {
      warning("Obsolete. Please use \"model=br\" instead of \"hr\" for spatial models.")
      model <- "br"
    } else {
      warning("Unused arguments \"vario\", \"dist\" or \"coord\"; only implemented for extremal Student-t or Brown-Resnick processes.")
    }
  }
  if(model %in% m3){
    if(!is.logical(grid) || length(grid) != 1){
      stop("Argument \"grid\" must be a logical, either TRUE or FALSE.")
    }
    if(!is.null(dist)){
      if(!is.matrix(dist) || !isSymmetric(dist) || any(dist[upper.tri(dist, diag = FALSE)] <= 0) || !isTRUE(all(diag(dist) == 0))){
        stop("Invalid pairwise distance matrix \"dist\": distances must be nonnegative, locations unique and the matrix must be symmetric.")
      }
      if(!is.null(coord) || grid){
        warning("Arguments \"coord\" and \"grid\" are ignored if pairwise distance matrix \"dist\" is provided.")
        coord <- NULL
        grid <- FALSE
      }
    }
    if (!is.null(coord)) {
      coord <- as.matrix(coord)
      if (ncol(coord) == 1)
        grid <- FALSE
      if (grid) {
        if (all(sapply(1:ncol(coord), function(i) {
          length(unique(coord[, i])) == nrow(coord)
        }))) {
          coord <- matrix(unlist(expand.grid(apply(coord, 2, as.list))), ncol = ncol(coord))
        } else {
          stop("Duplicate values in \"coord\" using \"grid=TRUE\" not allowed")
        }
      }
    }
  }
  # One-parameter families
  if (model %in% m1) {
    d <- as.integer(d)
    if(is.null(d)){
      stop("No dimension \"d\" provided. Aborting.")
    }
    sigma = cbind(0)
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
  } else if (model %in% m2) {
    d <- as.integer(d)
    sigma = cbind(0)
    # if(model %in% c('bilog','negbilog')){
    if (is.null(param) || length(param) != d)
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
    # Scaled dirichlet family
  } else if (model %in% m4) {
    sigma = cbind(0)
    if (is.null(param)) {
      stop("Invalid parameter value")
    }
    if (model == "ct") {
      if (length(param) != d) {
        if (length(param) == (d + 1)) {
          warning("Use \"sdir\" model for the scaled extremal Dirichlet model.")
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
        stop("Invalid parameters for the scaled Dirichlet. rho must be greater than min(alpha)")
      }
      if (isTRUE(any(param[-(d + 1)] < 0))) {
        stop("Invalid arguments for the scaled Dirichlet model - alpha must be positive.")
      }
    }
    model = "sdir"
  } else if (model %in% m3) {
    if(!is.null(dist)){

    }
    # Smith, Brown-Resnick, extremal student
    if (is.null(sigma) && !is.null(vario) && (!is.null(coord) || !is.null(dist))) {
      if(!is.null(coord) && is.vector(coord)){
        coord <- matrix(coord, ncol = 1)  #1 dimensional process
      }
      stopifnot(is.function(vario))
      if (model == "br") {
        model = "isbr"
        m3 <- c(m3, model)
        if (vario(0, ...) > 1e-15) {
          stop("Cannot have a nugget term in the variogram for the Brown-Resnick process")
        }
        semivario2mat <- function(coord, semivario, dist, ...) {
          if(is.null(dist)){
            di <- as.matrix(dist(coord))
          } else{
            di <- dist
          } #fields::rdist(coord) is faster...
          covmat <- matrix(0, nrow = nrow(di), ncol = ncol(di))
          covmat[lower.tri(covmat)] <- semivario(di[lower.tri(di)], ...)
          covmat[upper.tri(covmat)] <- t(covmat)[upper.tri(covmat)]
          return(covmat)
        }
        sigma <- semivario2mat(coord, vario, dist, ...)/2
        # changed 08-03-2018 Matrix is half of Semivariogram, quarter of variogram
      }
    }
    if (model == "xstud") {
      sigma <- cov2cor(sigma)
    }
    if (model != "isbr") {
      if (is.null(sigma) || ncol(sigma) != nrow(sigma))
        stop("Invalid covariance matrix")
      if (any(diag(sigma) <= 0))
        stop("Degenerate covariance matrix; negative or zero entries found")
    }
    if (model == "xstud" && any(diag(sigma) != 1))
      warning("Extremal student requires correlation matrix")
    if (model == "xstud" && (is.null(param) || length(param) != 1)) {
      stop("Degrees of freedom argument missing or invalid")
    }
    if (model == "smith" && is.null(coord))
      stop("Location should be provided for the Smith model")
    if (model == "smith" && ncol(as.matrix(coord)) != ncol(sigma)) {
      stop("Covariance matrix of the Smith model should be
           of the same dimension as the number of columns of the coordinates matrix.")
    }
    d <- switch(model, xstud = ncol(sigma), br = ncol(sigma), smith = nrow(coord), isbr = ncol(sigma))
    if (model %in% c("smith", "br", "isbr")) {
      param <- 0
    }
  } else if (model %in% c("alog", "aneglog","maxlin")) {
    # Sigma will be index of logistic sub-mixtures param is vector of dependence parameters
    if(model %in% c("alog","aneglog")){
    if (any(param < 0))
      stop("Parameter vector must be positive")
    param <- rep(param, length.out = 2^d - 1 - d)  #if vector too short, recycle dep arguments
    if (model == "alog") {
      if (isTRUE(all.equal(param >= 1, rep(TRUE, length(param))))) {
        param <- 1/param  #sampler works with
      }
      sigma <- .mvasym.check(asy, param, d = d, model = "alog")  #check parameters constraints
      param <- 1/param
    } else {
      sigma <- .mvasym.check(asy, param, d = d, model = "aneglog")  #check parameters constraints
    }
    } else{
      # Prefered argument for maxlin is weights
      if(!is.null(weights) & !is.null(d)){
        weights <- as.matrix(weights)
        if(ncol(weights) != d){
          stop("Invalid \"weights\" argument: must be a matrix with \"d\" columns.")
        }
        weights <- apply(weights, 2, function(x){x/sum(x)})
        if(!isTRUE(all(is.finite(weights)))){
          stop("Invalid arguments in \"weights\": either columns do not normalize to one or non-finite values.")
        }
      sigma <- weights
     } else{
      sigma <- .mvasym.check(asy, param,
                             d = d,
                             model = "maxlin")
     }
      # Shed matrix to remove zero lines
      sigma <- sigma[rowSums(sigma) > 0, ]
      param <- 1
     }
    # Transform list to matrix, with correct correspondance
    asym <- sigma > 0
    if(model != "maxlin"){
    # Shed output to remove zero weight combinations
    if (d == 2) {
      param <- c(sigma[1, 1] != 0, sigma[2, 2] != 0, param)  # not a matrix for d=2
    } else {
      zero_line <- which(rowSums(sigma[-(1:d), ]) == 0) + d  #Possibly empty,
      param <- c(!1:d %in% zero_line, param)  #set dummies for point masses on edges
      if (length(zero_line) > 0) {
        sigma <- sigma[-zero_line, , drop = FALSE]
        asym <- asym[-zero_line, , drop = FALSE]
        param <- param[-zero_line]
      }
    }
    stopifnot(isTRUE(all.equal(nrow(sigma), nrow(asym), length(param))), isTRUE(all.equal(colSums(sigma), rep(1, ncol(sigma)))))
}
    ncompo <- rowSums(asym)
  } else if (model == "dirmix") {
    if (any(is.null(param), length(weights) != ncol(param) && ncol(param) != 1, any(param < 0))) {
      stop("Invalid arguments for the Dirichlet mixture")
    }
    if (!is.null(weights)) {
      if (any(weights < 0))
        stop("Negative weights provided")
      if (sum(weights) != 1)
        warning("Weights do not sum to one")
      weights <- weights/sum(weights)
    }
    if (is.null(d)) {
      d <- nrow(param)
    } else if (d != nrow(param)) {
      stop("\"d\" and dimension of provided \"param\" vector do not match.")
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
    d <- nrow(sigma)
    param <- 0
    if (any(c(sigma < 0, diag(sigma) != rep(0, ncol(sigma)), ncol(sigma) != nrow(sigma), ncol(sigma) != d))) {
      stop("Invalid parameter matrix for the Husler-Reiss model.
           Note: for Brown-Resnick model, please use model=\"br\" instead.")
    }
    if (!.is.CNSD(sigma)) {
      stop("Parameter matrix for the Huesler-Reiss model is not conditionally negative definite")
    }
  } else if(model %in% m5){
    sigma = matrix(0,1,1)
    alg <- "sm"
    # pairwise beta of Cooley et al. (2010), JMVA
    # pairwise exponential, weighted Dirichlet
    # and weighted exponential models of Ballani and Schlather (2011), Biometrika
    if(model %in% c("pairbeta","pairexp")){
      if(length(param) != d*(d-1)/2 + 1L){
        stop(paste("Invalid dimension for the parameter vector of the", model, "model."))
      }
      if((model=="pairexp")&(param[1] <= 0)){
        stop(paste("The parameter alpha of the pairwise exponential model must be positive."))
      }
      if((model == "pairbeta")&(isTRUE(any(param <= 0)))){
        stop(paste("The parameters of the pairwise beta model must be positive."))
      }
    } else if(model %in% c("wdirbs","wexpbs")){
      if(length(param) != 2*d){
        stop(paste("Invalid dimension for the parameter vector of the", model, "model."))
      }
      if(isTRUE(any(param[1:d] <= 0))){
        stop(paste("The parameters alpha of the", model, "model must be positive."))
      }
      if((model == "wdirbs")&(any(param[(d+1):(2*d)] <= 0))){
        stop(paste("The parameters beta of the", model, "model must be positive."))
      }
    }
  } else{
    stop("Model currently not implemented.")
  }
  if(is.null(coord)){ coord <- matrix(1,0,0)}
  return(list(n = n, d = d, param = param, asy = asy, sigma = sigma,
              model =  model, alg = alg, weights = weights,
              vario = vario, coord = coord, asym = asym, ncompo = ncompo))
}
