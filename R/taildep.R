#' Coefficient of tail correlation and tail dependence
#'
#' For data with unit Pareto margins, the coefficient of tail dependence \eqn{\eta} is defined  via \deqn{\Pr(\min(X) > x) = L(x)x^{-1/\eta},}
#' where \eqn{L(x)} is a slowly varying function. Ignoring the latter, several estimators of \eqn{\eta} can be defined. In unit Pareto margins, \eqn{\eta} is a nonnegative shape parameter that can be estimated by fitting a generalized Pareto distribution above a high threshold. In exponential margins, \eqn{\eta} is a scale parameter and the maximum likelihood estimator of the latter is the Hill estimator. Both methods are based on peaks-over-threshold and the user can choose between pointwise confidence confint obtained through a likelihood ratio test statistic (\code{"lrt"}) or the Wald statistic (\code{"wald"}).
#'
#' The most common approach for estimation is the empirical survival copula, by evaluating the proportion of sample minima with uniform margins that exceed a given \eqn{x}. An alternative estimator uses a smoothed estimator of the survival copula using Bernstein polynomial, resulting in the so-called \code{betacop} estimator. Approximate pointwise confidence confint for the latter are obtained by assuming the proportion of points is binomial.
#'
#' The coefficient of tail correlation \eqn{\chi} is
#' \deqn{\chi = \lim_{u \to 1} \frac{\Pr(F_1(X_1)>u, \ldots, F_D(X_D)>u)}{1-u}.}
#' Asymptotically independent vectors have \eqn{\chi = 0}. The estimator uses an estimator of the survival copula
#' @note As of version 1.15, the percentiles used are from the minimum variable. This ensures that, regardless of the number of variables,
#' there is no error message returned because the quantile levels are too low for there to be observations
#' @export
#' @param data an \eqn{n} by \eqn{d} matrix of multivariate observations
#' @param u vector of percentiles between 0 and 1
#' @param nq number of quantiles of the structural variable at which to form a grid; only used if \code{u = NULL}.
#' @param qlim limits for the sequence \code{u} of the structural variable
#' @param depmeas dependence measure, either of \code{"eta"} or \code{"chi"}
#' @param confint string indicating the type of confidence interval for \eqn{\eta}, one of \code{"wald"} or \code{"lrt"}
#' @param level the confidence level required (default to 0.95).
#' @param trunc logical indicating whether the estimates and confidence intervals should be truncated in \eqn{[0,1]}
#' @param ties.method string indicating the type of method for \code{rank}; see \code{\link[base]{rank}} for a list of options. Default to \code{"random"}
#' @param empirical.transformation logical indicating whether observations should be transformed to pseudo-uniform scale (default to \code{TRUE}); otherwise, they are assumed to be uniform
#' @param method named list giving the estimation method for \code{eta} and \code{chi}. Default to \code{"emp"} for both.
#' @param plot logical; should graphs be plotted?
#' @param ... additional arguments passed to \code{plot}; current support for \code{main}, \code{xlab}, \code{ylab}, \code{add} and further \code{pch}, \code{lty}, \code{type}, \code{col} for points; additional arguments for confidence intervals are handled via \code{cipch}, \code{cilty}, \code{citype}, \code{cicol}.
#' @return a named list with elements
#' \itemize{
#' \item \code{u}: a \code{K} vector of percentile levels
#' \item \code{eta}: a \code{K} by 3 matrix with point estimates, lower and upper confidence intervals
#' \item \code{chi}: a \code{K} by 3 matrix with point estimates, lower and upper confidence intervals
#' }
#' @seealso \code{\link[evd]{chiplot}} for bivariate empirical estimates of \eqn{\chi}{chi} and \eqn{\bar{\chi}}{chibar}.
#' @examples
#' \dontrun{
#' set.seed(765)
#' # Max-stable model
#' dat <- rmev(n = 1000, d = 4, param = 0.7, model = "log")
#' taildep(dat, confint = 'wald')
#' }
#' @importFrom "utils" "combn"
taildep <- function (data,
                     u = NULL,
                     nq = 40,
                     qlim = c(0.8, 0.99),
                     depmeas = c("eta","chi"),
                     method = list(eta = c("emp","betacop", "gpd", "hill"),
                                   chi = c("emp","betacop")),
                     confint = c("wald","lrt"),
                     level = 0.95,
                     trunc = TRUE,
                     empirical.transformation = TRUE,
                     ties.method = "random",
                     plot = TRUE, ...) {
  if(is.character(depmeas)){
    depmeas <- which(c("eta", "chi") %in% match.arg(depmeas, c("eta", "chi"), several.ok = TRUE))
  }
  depmeas_names <- c("eta", "chi")[depmeas]
  if((! length(depmeas) %in% c(1L, 2L)) || !all(depmeas  %in% c(1L, 2L))){
    stop("Invalid argument for \"depmeas\": must be either \"eta\",\"chi\" or both.")
  }
  if(! (is.list(method) && isTRUE(all(depmeas_names %in% names(method))))){
    stop(paste0("\"method\" should be a list with components ",
               paste(depmeas_names[depmeas], collapse = " and "),"."))
  }
  ties.method <- match.arg(ties.method, c("average", "first", "last", "random", "max", "min"))
  stopifnot(is.logical(empirical.transformation))
  data <- as.matrix(data)
  data <- na.omit(data)
  n <- nrow(data)
  D <- ncol(data)
  stopifnot(D > 1)
  if(1 %in% depmeas){
    methodeta <- match.arg(method$eta[1], choices = c("emp","betacop", "gpd", "hill"))
  } else{
   methodeta <- c()
  }
  if(2 %in% depmeas){
   methodchi <- match.arg(method$chi[1], choices = c("emp","betacop"))
  } else{
   methodchi <- c()
  }
  method <- c(methodeta, methodchi)
  if(empirical.transformation){
  # map observations to uniform scale
  if("betacop" %in% method){
    if(ties.method != "random"){
      warning("Beta copula does not allow for ties; switching to `ties.method = \"random\"`")
      ties.method <- "random"
    }
  }
   datarank <- apply(data, 2, rank, ties.method = ties.method)
   rowmin <- apply(datarank, 1, min)/(n+1)
  } else{
  if("betacop" %in% method){
    warning("Beta copula is rank-based; ignoring argument \"empirical.transformation\".")
    datarank <- apply(data, 2, rank, ties.method = ties.method)
  }
  # observations are already uniform
   rowmin <- apply(data, 1, min)
   stopifnot(min(rowmin) > 0,
             max(rowmin) < 1)
  }
  eps <- .Machine$double.eps^0.5
  qlim2 <- c(min(rowmin) + eps, max(rowmin) - eps)
  if(is.null(u)){
     if (!is.null(qlim)) {
       qlim <- sort(qlim)
        stopifnot(length(qlim) == 2L,
                  qlim[1] > 1/(n+1),
                  qlim[2] < n/(n+1))
     } else{
       stop("Invalid input: neither \"qlim\" nor \"u\" is provided.")
     }
    # Define u in terms of quantiles of the structure variable min Y
    qlims <- quantile(x = rowmin, probs = c(qlim[1], qlim[2]))
    u <- seq(from = qlims[1], to = qlims[2], length = nq)
  } else{
    u <- sort(u)

    if(min(u) < qlim2[1] || max(u) >  qlim2[2]){
      warning("upper quantile limit is too high or lower quantile limit is too low")
    }
    #u <- u[u < qlim2[2] & u > qlim2[1]]
    nq <- length(u)
  }


  confint <- match.arg(confint)
  cnst <- qnorm((1 + level)/2)
  if("emp" %in% method){
    rmin <- rowmin
    cbaru <- sapply(u, function(ui){sum(rmin > ui)})/n
    if(1 %in% depmeas && methodeta == "emp"){
    etau <- log(1 - u) / log(cbaru)
    seeta <- sqrt((((log(1 - u)^2)/(log(cbaru)^4 * cbaru)) * (1 - cbaru)) / n)
    est_eta <- cbind(eta = etau, etalow = etau - cnst * seeta, etaupp = etau + cnst * seeta)
    }
    if(2 %in% depmeas && methodchi == "emp"){
    chiu <- cbaru / (1 - u)
    sechi <- sqrt(cbaru * (1 - cbaru)/(1 - u)^2/n)
    est_chi <- cbind(chi = chiu, chilow = chiu - cnst * sechi, chiupp = chiu + cnst * sechi)
    }

  }
  if (1 %in% depmeas && methodeta == "gpd"){
    ps <- 1/(1-rowmin)
    if(confint == "lrt"){
    est_eta <- t(sapply(1/(1-u), function(th){
    # At least 15 exceedances for model fit
      if(sum(ps>th) > 15){
      fitu <- suppressWarnings(fit.gpd(ps, threshold = th))
      prof <- try(suppressWarnings(gpd.pll(param = "shape",
      dat = ps, threshold = th, mod = "prof", mle = fitu$estimate, plot = FALSE)))
      co <- try(suppressWarnings(confint(prof, prob = c((1-level)/2, (1+level)/2), print = FALSE)))
      if(inherits(x = co, what = "try-error")){
        return(c(fitu$estimate[2], rep(NA, 2)))
      } else{
        return(co)
      }
    } else{
      return(rep(NA, 3))
    }
    }))
    } else{
    est_eta <- t(sapply(1/(1-u), function(th){
      if(sum(ps>th) > 15){
        fitu <- suppressWarnings(fit.gpd(ps, threshold = th))
        if(fitu$estimate[2] < -0.5){
         c(fitu$estimate[2], NA, NA)
        } else{
          c(fitu$estimate[2], fitu$estimate[2] + c(-1,1)*fitu$std.err[2]*qnorm((1+level)/2))
        }
      } else{
        return(rep(NA, 3))
      }
    }))
    }
  } else if (1 %in% depmeas && methodeta == "hill"){
    expll <- function(x, il){
     sum(-log(il) - x/il)
    }
    es <- -log(1-rowmin)
    est_eta <- t(sapply(-log(1-u), function(th){
      samp <- es[es>th] - th
      mle <- mean(samp)
      maxll <- expll(samp, mle)
      if(confint == "lrt"){
        low <- try(uniroot(f = function(il){ 2*(expll(samp, il) - maxll) + qchisq(level, 1)}, interval = c(1e-8, mle)))
        upp <- try(uniroot(f = function(il){ 2*(expll(samp, il) - maxll) + qchisq(level, 1)}, interval = c(mle, 5)))
       return(c(mle,
                ifelse(inherits(low, what = "try-error"), NA, low),
                ifelse(inherits(upp, what = "try-error"), NA, upp)))
      } else if(confint == "wald"){
       return(c(mle, mle + c(-1,1)*mle/sqrt(length(samp))*qnorm((1+level)/2)))
      }
    }))
  }
    if("betacop" %in% method){
    # This uses the ranks, not the observations
    cbaru <- numeric(nq)
    for(i in seq_len(nq)){
      # Beta copula
      # Evaluate quantiles of beta at u
      # Model defined in terms of the survival copula Pr(V_i > u, i = 1,...,d)
      # Flip arguments u -> 1-u
      # and ranks (via n:1 instead of 1:n)
      Fu <- sapply(rev(seq_len(n)), function(r){
        suppressWarnings(pbeta(1 - u[i], r, n+1-r, log.p = TRUE))
        })
      cbaru[i] <- mean(exp(rowSums(matrix(Fu[datarank], nrow = n))))
    }
    cbaru <- pmax(cbaru, 0)
    if(2 %in% depmeas && methodchi == "betacop"){
      chiu <- cbaru / (1 - u)
      sechi <- sqrt(cbaru * (1 - cbaru)/(1 - u)^2/n)
      est_chi <- cbind(chiu, chiu - cnst * sechi, chiu + cnst * sechi)
    }
    if(1 %in% depmeas && methodeta == "betacop"){
      etau <- log(1 - u) / log(cbaru)
      seeta <- sqrt(((((log(1 - u))^2)/(log(cbaru)^4 * cbaru)) * (1 - cbaru)) / n)
      est_eta <- cbind(etau, etau - cnst * seeta, etau + cnst * seeta)
    }
  }
    if (trunc) {
      if(1 %in% depmeas){
        est_eta[est_eta > 1] <- 1
        est_eta[est_eta < 0] <- 0
      }
      if(2 %in% depmeas){
        est_chi[est_chi > 1] <- 1
        est_chi[est_chi < 0] <- 0
      }
    }
    out <- list(u = u)
    if(1 %in% depmeas){
      colnames(est_eta) <- c("coef", "lowerci", "upperci")
      out$eta <- est_eta
      out$eta_method <- methodeta
      out$eta_confint_method <- confint
    }
    if(2 %in% depmeas){
      colnames(est_chi) <- c("coef", "lowerci", "upperci")
      out$chi <- est_chi
      out$chi_method <- methodchi
    }
    out$dim <- dim(data)
    class(out) <- "mev_taildep"
    if(plot){
     plot(out, ...)
    }
    invisible(out)
}

#' @export
plot.mev_taildep <- function(x, ...){
  if(length(x$u) < 3){
    return(invisible(NULL))
  }
  ellips <- list(...)
  depmeas <- c(!is.null(x$eta), !is.null(x$chi))
  if(sum(depmeas) == 0){
    stop("Invalid argument")
  } else if(sum(depmeas) == 2){
    if(!is.null(ellips$which)){
    if(is.character(ellips$which)){
     show <- c("eta", "chi") %in% match.arg(ellips$which, choices = c("eta", "chi"), several.ok = TRUE)
    } else if(is.numeric(ellips$which)){
     show <- 1:2 %in% ellips$which
    } else if(is.logical(ellips$which)){
     show <- ellips$which
    }
    } else{
     show <- rep(TRUE, 2)
    }
  } else{
    show <- depmeas
  }

  old.par <- par(no.readonly = TRUE)
  if(!is.null(ellips$add)){
    if(!isTRUE(ellips$add)){
      if(sum(show) == 2){
        par(mfrow = c(1,2))
      }
      on.exit(par(old.par))
    }
  } else{
    if(sum(show) == 2){
      par(mfrow = c(1,2))
    }
    on.exit(par(old.par))
  }
  if (is.null(ellips$main)) {
    main <- rep("",2)[show]
  } else{
    if(length(ellips$main) != sum(show)){
      stop("Invalid input: \"main\" must be of the same length as \"which\"")
    }
    main <- ellips$main
  }
  if (is.null(ellips$xlab)) {
    xlab <- rep("u", 2)
  } else{
    xlab <- rep(ellips$xlab, length.out = 2)
  }
  if (is.null(ellips$ylab)) {
    ylab <- expression(eta, chi)
  } else{
    if(length(ellips$ylab) != sum(show)){
      stop("Invalid input: \"ylab\" must be of the same length as \"which\"")
    }
    ylab <- rep(ellips$ylab, length.out = 2)
  }
  if (is.null(ellips$col)) {
    col <- rep(1, 2)
  } else{
    col <- rep(ellips$col, length.out = 2)
  }
  if (is.null(ellips$lty)) {
    lty <- rep(1, 2)
  } else{
    lty <- rep(ellips$lty, length.out = 2)
  }
  if (is.null(ellips$pch)) {
    pch <- rep(20, 2)
  } else{
    pch <- rep(ellips$pch, length.out = 2)
  }
  if (is.null(ellips$type)) {
    type <- rep("p", 2)
  } else{
    type <- rep(ellips$type, length.out = 2)
  }
  #Arguments for confidence lines, with ci

  if (is.null(ellips$cicol)) {
    cicol <- rep("grey", 2)
  } else{
    cicol <- rep(ellips$cicol, length.out = 2)
  }
  if (is.null(ellips$cilty)) {
    cilty <- rep(2, 2)
  } else{
    cilty <- rep(ellips$cilty, length.out = 2)
  }
  if (is.null(ellips$cipch)) {
    cipch <- rep(95, 2) #corresponds to - sign
  } else{
    cipch <- rep(ellips$cipch, length.out = 2)
  }
  if (is.null(ellips$citype)) {
    citype <- rep("p", 2)
  } else{
    citype <- rep(ellips$citype, length.out = 2)
  }
  D <- x$dim
  if (show[1]) {
    matplot(x$u, x$eta, type = c(type[1], rep(citype[1], 2)), pch = c(pch[1], rep(cipch[1], 2)),
            lty = c(lty[1], rep(cilty[1], 2)), col = c(col[1], rep(cicol[1], 2)),
            xlim = c(min(x$u) - diff(range(x$u))/100, pmin(1, max(x$u) + diff(range(x$u))/100)),
            ylim = c(0, 1),
            xaxs = "i",
            yaxs ="i",
            main = main[1],
            xlab = xlab[1],
            ylab = ylab[1],
            bty = "l",
            panel.first = {
                segments(x$u[1], 0, 1, 0, lty = 5, col = "grey")
                segments(x$u[1], 1, 1, 1, lty = 5, col = "grey")
                segments(x$u[1], 1/D, 1, 1/D, lty = 5, col = "grey")
              })
  }
  if (show[2]) {
    matplot(x$u, x$chi, type = c(type[2], rep(citype[2], 2)), pch = c(pch[2], rep(cipch[2], 2)),
            lty = c(lty[2], rep(cilty[2], 2)),
            col = c(col[2], rep(cicol[2], 2)),
            xlim = c(min(x$u) - diff(range(x$u))/100, pmin(1, max(x$u) + diff(range(x$u))/100)),
            ylim = c(0, 1),
            xaxs = "i",
            yaxs ="i",
            main = main[2],
            bty = "l",
            xlab = xlab[2],
            ylab = ylab[2],
            panel.first = {
                segments(x$u[1], 0, 1, 0, lty = 5, col = "grey")
                segments(x$u[1], 1, 1, 1, lty = 5, col = "grey")
             })
  }

}
