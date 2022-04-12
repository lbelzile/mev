## ----setup, eval = TRUE, echo = FALSE-----------------------------------------
library(mev)

## ----gevpll-------------------------------------------------------------------
# Fetch data and dates (see ?maiquetia)
data(maiquetia, package = "mev")
day <- seq.Date(from = as.Date("1961-01-01"), 
                to = as.Date("1999-12-31"), 
                by = "day")
# Compute yearly maximum of daily rainfall
ymax <- tapply(maiquetia, factor(substr(day, 1, 4)), max)

## -----------------------------------------------------------------------------
# Creates plot by default
prof <- gev.pll(param = "Nmean", 
        mod = c("profile", "tem"), 
        dat = ymax[-length(ymax)],
        N = 50)
# Confidence intervals
confint(prof, print = TRUE)

## ----gp-----------------------------------------------------------------------
library(mev)
set.seed(1234)
dat <- evd::rgpd(n = 10, shape = -0.8)
fitted <- fit.gpd(dat, threshold = 0, show = TRUE)
# Empirical coefficient of variation
# Theoretical quantity defined as standard deviation/mean
sd(dat)/mean(dat)

## ----gpprofile----------------------------------------------------------------
prof_gpd_eta <- function(eta, xdat){
  ll <- -length(xdat) - sum(log(1-eta*xdat))-
    length(xdat)*log(-mean(log(1-eta*xdat))/eta)
}
# Grid value for eta, excluding a neighborhood of zero
etagrid <- c(seq(-4, -1e-8, length = 201L),
             seq(1e-8, 1/max(dat)+1e-10, length.out = 301L))
ll <- sapply(etagrid, FUN = prof_gpd_eta, xdat = dat)
# For zero, we get xi=0 and exponential model,
# whose mle for scale is the sample mean
etagrid <- c(etagrid,0)
ll <- c(ll, sum(dexp(dat, rate = 1/mean(dat), log = TRUE)))
ll <- ll[order(etagrid)]
etagrid <- sort(etagrid)
xis <- sapply(etagrid, function(eta){mean(log(1-eta*dat))})
sub <- xis > -1
xis <- xis[sub]
etagrid <- etagrid[sub]
ll <- ll[sub]
# Plot the log likelihood
mllbound <- mev::gpd.ll(dat = dat, par = c(max(dat),-1))
mll <- pmax(max(ll, na.rm = TRUE), mllbound)
par(mfrow = c(1,2), mar = c(4,4,1,1), bty = 'l')
 plot(etagrid, ll-mll, type = 'l',  ylim = c(-5,0),
      ylab = "profile log likelihood", xlab = expression(eta))
 points(1/max(dat), mllbound-mll)
# plot(xis, ll-mll, type = 'l', xlim = c(-1,1), ylim = c(-5,0),
#      ylab = "profile log likelihood", xlab = expression(xi))
# points(-1, mllbound-mll)
dat <- evd::rgpd(n = 20, scale = 1, shape = 0.1)
etagrid <- c(seq(-4, -1e-8, length = 201L),
             seq(1e-8, 1/max(dat)+1e-10, length.out = 301L))
ll <- sapply(etagrid, FUN = prof_gpd_eta, xdat = dat)
# For zero, we get xi=0 and exponential model,
# whose mle for scale is the sample mean
etagrid <- c(etagrid,0)
ll <- c(ll, sum(dexp(dat, rate = 1/mean(dat), log = TRUE)))
ll <- ll[order(etagrid)]
etagrid <- sort(etagrid)
xis <- sapply(etagrid, function(eta){mean(log(1-eta*dat))})
sub <- xis > -1
xis <- xis[sub]
etagrid <- etagrid[sub]
ll <- ll[sub]
# Plot the log likelihood
mllbound <- mev::gpd.ll(dat = dat, par = c(max(dat),-1))
mll <- pmax(max(ll, na.rm = TRUE), mllbound)
# par(mfrow = c(1,2), mar = c(4,4,1,1), bty = 'l')
 plot(etagrid, ll-mll, type = 'l',  ylim = c(-5,0),
      ylab = "profile log likelihood", xlab = expression(eta))
 points(1/max(dat), mllbound-mll)
# plot(xis, ll-mll, type = 'l', xlim = c(-1,1), ylim = c(-5,0),
#      ylab = "profile log likelihood", xlab = expression(xi))
# points(-1, mllbound-mll)s

## ----gp2----------------------------------------------------------------------
# Another example where the solution lies inside the parameter space
dat <- evd::rgpd(n = 25, shape = 0.2)
fitted <- fit.gpd(dat, threshold = 0, show = FALSE)
# Check convergence - is gradient zero?
isTRUE(all.equal(gpd.score(par = coef(fitted), dat = dat),
                 rep(0, 2)))
# Various methods are available
methods(class = "mev_gpd")

# P-P and Q-Q diagnostic plots 
par(mfrow = c(1, 2))
plot(fitted)

# Fit exponential by passing a list with a fixed parameter
reduced <- fit.gpd(dat, threshold = 0, fpar = list(shape = 0))
# The MLE is sample mean of exceedances - check this
isTRUE(coef(reduced) == mean(dat))
# Compare models using likelihood ratio test
anova(fitted, reduced)

## ----gpalternative------------------------------------------------------------
# Bayesian point estimates (based on MAP)
fit.gpd(dat, threshold = 0, 
        show = TRUE, 
        method = "zhang")
# With MCMC
fit.gpd(dat, threshold = 0, 
        show = TRUE, 
        MCMC = TRUE,
        method = "zhang")
# OBRE fit - a weight, attached to the largest
# observations is returned
fit_robust <- fit.gpd(dat, 
                      threshold = 0, 
                      show = TRUE, 
                      method = "obre")
# See fit_robust$weights

# First-order bias corrected estimates
corr_coef <- gpd.bcor(par = coef(fitted), 
                      dat = dat, 
                      corr = "firth")

# Many methods are available for these objects
# including the following `S3` classes
methods(class = "mev_gpd")

## ----ppfitmultimodal----------------------------------------------------------
set.seed(202010)
xdat <- evd::rgpd(n = 1000, shape = -0.5)
mle <- mev::fit.pp(xdat, np = 10)$param
shape_v <- seq(-0.7, 1.5, by = 0.005)
loc_v <- seq(-2,3, by = 0.005)
ll_s <- matrix(NA, nrow = length(shape_v), ncol = length(loc_v))
for(shape_i in seq_along(shape_v)){
  for(loc_i in seq_along(loc_v)){
    ll <- try(pp.ll(par = c(loc_v[loc_i], mle[2], shape_v[shape_i]),
                dat = xdat, 
                u = 0, np = 10))
    if(!inherits(ll, "try-error")){
      ll_s[shape_i, loc_i] <- ll
    }
  }
}
# Create an image with the contour curves
image(shape_v, loc_v, 
      is.finite(ll_s), 
      useRaster = TRUE, 
      col = grey.colors(n = 100,
                        start = 1, 
                        end = 0.6,
                        alpha = 1),
      xlab = expression(xi),
      ylab = expression(mu))
contour(shape_v, loc_v, z = -(max(ll_s, na.rm = TRUE)-ll_s),
      xlab = expression(xi),
      ylab = expression(mu),
      zlim = c(-1e4, 0), 
      add = TRUE)
points(mle[3], mle[1], pch = 19)
points(mle[3], mle[1], col = "white", pch = 20)
fake_mle <- mev::fit.pp(xdat, np = 10, 
                        fpar = list(scale = 0.08560669),
                        start = c(0, 1))$estimate
points(fake_mle[2], fake_mle[1], pch = 19)
points(fake_mle[2], fake_mle[1], col = "white", pch = 20)

## ----"prelimfitmaiquetia"-----------------------------------------------------
library(mev)
data("maiquetia", package = "mev")
day <- seq.Date(from = as.Date("1961-01-01"), 
                to = as.Date("1999-12-31"), by = "day")
# Keep non-zero rainfall, exclude 1999 observations
nzrain <- maiquetia[substr(day, 3, 4) < 99 & maiquetia > 0]
gpdf <- fit.gpd(nzrain, threshold = 20)
print(gpdf)

## ----checkscore---------------------------------------------------------------
isTRUE(all.equal(
  gpd.score(gpdf$estimate, dat = gpdf$exceedances),
  c(0,0), tolerance = 1e-5))

## ----biascor------------------------------------------------------------------
gpdbcor <- gpd.bcor(dat = gpdf$exceedances, par = gpdf$estimate)
#print the differences between MLE and bias-corrected estimates
gpdf$estimate - gpdbcor

