#' Estimators of the tail coefficient
#'
#' Estimators proposed by Krupskii and Joe under second order expansion
#' for the coefficient of tail dependence \eqn{\eta} and the
#' joint tail orthant probability
#' @examples
#' rho <- runif(1, -1, 1)
#' d <- 2
#' data <- mev::mvrnorm(
#'    n = 1e5,
#'    mu = rep(0, d),
#'   Sigma = matrix(rho, d, d) + diag(1 - rho, d))
#' q <- seq(0.95, 0.995, by = 0.005)
#' taildep <- kjtail(data = data, q = q)
#' with(taildep,
#'  plot(x = 1-pat,
#'       y = eta,
#'       ylim = c(0,1),
#'       panel.first = {abline(h = (1+rho)/2)}))
kjtail <- function(
    data,
    q,
    ptail = NULL,
    mqu,
    type = 1,
    ties.method = eval(formals(rank)$ties.method),
    ...){
  data <- na.omit(as.matrix(data))
  d <- ncol(data)
  n <- nrow(data)
  if(isTRUE(any(q <= 0, q > 1))){
    stop("\"q\" must be a vector of probabilities.")
  }
  if(isTRUE(any(q < 0.5))){
    warning("Small quantile selected.")
  }
  # Transform to tail quantiles
  q <- 1 - q
  if(missing(mqu)){
    unif <- apply(data, 2, rank, ties.method = ties.method)/(n + 1L)
  } else{
  if(!length(mqu) %in% c(1L, d)){
    stop("Invalid marginal threshold vector \"mqu\"")
  }
  if(isTRUE(any(mqu <= 0, mqu > 1))){
      stop("\"mqu\" must be a vector of probabilities.")
  }
  mqu <- rep(mqu, length.out = d)
  qulvl <- numeric(d)
  for(i in seq_len(d)){
    qulvl[i] <- quantile(data[,i], probs = mqu[i])
  }
  unif <- mev::spunif(x = data, thresh = qulvl)
  }
  # Consider minimum of tail
  Ud <- apply(unif, 1, function(x){max(1-x)})
  q <- sort(q, decreasing = TRUE)
  etahat <- etasd <- k1hat <- nabove <- jtprob <- plev <- rep(NA_real_, length(q))
  neglik_fn <- function(par, xdat, Tp){
    k1 <- par[1]
    eta <- par[2]
    - sum(log(pmax(0, k1+(1-k1*Tp)/eta*exp((1/eta-1)*log(xdat) - log(Tp)/eta))))
  }
  for(j in seq_along(q)){
    plev[j] <- Tp <- quantile(Ud, q[j])
    sUd <- Ud[Ud < Tp]
    nabove[j] <- length(sUd)
    if(j == 1L){
      start <- c(1/Tp, 1/(d-0.5))
    } else{
      start <- est$pars
    }
  est <- Rsolnp::solnp(
     pars = c(1/Tp, 1/(d-0.5)),
     fun = neglik_fn,
     ineqfun = function(par, xdat, Tp){
       # Probability must be positive
       # k1 must be positive?
      par[1]*range(xdat)+(1-par[1]*Tp)*exp((log(range(xdat)) - log(Tp)/par[2]))
    },
    ineqLB = rep(0, 2),
    ineqUB = rep(1, 2),
    LB = rep(0, 2),
    UB = c(Inf, 1),
    Tp = Tp,
    xdat = sUd,
    control = list(trace = 0))
  est$hessian <- numDeriv::hessian(func = neglik_fn, x = est$pars, Tp = Tp, xdat = sUd)
    if(isTRUE(est$convergence == 0)){
     k1hat[j] <- est$pars[1]
     etahat[j] <- est$pars[2]
     etasd_j <- try(silent = TRUE, suppressWarnings(sqrt(solve(est$hessian[1:2,1:2])[2,2])))
     if(!inherits(etasd_j, "try-error")){
       etasd[j] <- etasd_j
     }
    }
  }
  lambdahat <- nabove/n*k1hat
  ret_list <- list(
    p = 1-as.numeric(plev),
    eta = as.numeric(etahat),
    eta_sd = etasd,
    pat = nabove/n,
    lambda = as.numeric(lambdahat))
  if(!is.null(ptail)){
    if(isTRUE(all(length(ptail) == 1L, ptail > 0, ptail < min(plev)))){
    ret_list$tailprob <- as.numeric(nabove/n*(k1hat*ptail + (1-k1hat*plev)*(ptail/plev)^(1/etahat)))
    } else{
      warning("Invalid \"ptail\" argument. Must be a scalar smaller than the probability selected.")
    }
  }
  return(ret_list)
}
