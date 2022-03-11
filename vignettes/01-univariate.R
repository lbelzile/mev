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

