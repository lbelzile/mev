## ----maxstabtest, cache = TRUE, fig.width = 10, fig.height = 6, out.width = '90%', fig.align = "center"----
library(mev)
set.seed(0)
samp <- rmev(
  n = 1000,
  d = 5,
  param = 0.1,
  model = "neglog"
)
fgev <- fit.gev(c(samp), show = FALSE)
fgev$estimate
par(mfrow = c(1, 2))
# Test of max-stability
maxstabtest(dat = samp)
# Probability-probability plot
plot(fgev, which = 1, main = "")

## ----angdensplot, cache = TRUE, fig.width = 10, fig.height = 6, out.width = '90%', fig.align = "center"----
samp <-
  rmev(
    n = 1000,
    d = 3,
    param = c(0.4, 0.6, 2.9, 0.1),
    model = "sdir"
  )
taildep(samp, method = list(eta = "betacop", chi = "betacop"))

## ----plotweights, fig.width = 8, fig.height = 6, out.width = '90%', fig.align = "center"----
# Plot the probability weights and compute the column mean
nparangmeas <- mev::angmeas(samp, th = 0.5)
plot(nparangmeas$wts, 
     ylab = "weights", 
     xlab = "observation index")
abline(h = 1 / nrow(nparangmeas$ang))
colSums(nparangmeas$wts * nparangmeas$ang)
dirangmeas <- mev::angmeasdir(samp, th = 0.5)

## ----extcoef, cache= TRUE, fig.width = 10, fig.height = 8, out.width = '90%', fig.align = "center"----
coord <- 10 * cbind(runif(50), runif(50))
di <- as.matrix(dist(coord))
dat <-
  rmev(
    n = 1000,
    d = 100,
    param = 3,
    sigma = exp(-di / 2),
    model = 'xstud'
  )
res <- extcoef(dat = dat, coord = coord)
# Extremal Student extremal coefficient function
XT.extcoeffun <- function(h, nu, corrfun, ...) {
  if (!is.function(corrfun)) {
    stop('Invalid function `corrfun`.')
  }
  h <- unique(as.vector(h))
  rhoh <- sapply(h, corrfun, ...)
  cbind(h = h, extcoef = 2 * pt(sqrt((nu + 1) * (1 - rhoh) / (1 + rhoh)), nu +
                                  1))
}
#This time, only one graph with theoretical extremal coef
plot(
  res$dist,
  res$extcoef,
  ylim = c(1, 2),
  pch = 20,
  ylab = "extremal coefficient",
  xlab = "distance"
)
extcoefxt <- XT.extcoeffun(
  seq(0, max(res$dist), by = 0.1),
  nu = 3,
  corrfun = function(x) {
    exp(-x / 2)
  }
)
lines(extcoefxt[, 'h'],
      extcoefxt[, 'extcoef'],
      type = 'l',
      col = 'blue',
      lwd = 2)


## ----simuRgpareto, cache = TRUE,fig.width = 20, fig.height= 10, out.width = '90%', fig.align = "center"----
lon <- seq(650, 720, length = 50)
lat <- seq(215, 290, length = 50)
# Create a grid
grid <- expand.grid(lon, lat)
coord <- as.matrix(grid)
dianiso <- distg(coord, 1.5, 0.5)
sgrid <- scale(grid, scale = FALSE)
# Specify marginal parameters `loc` and `scale` over grid
eta <- 26 + 0.05 * sgrid[, 1] - 0.16 * sgrid[, 2]
tau <- 9 + 0.05 * sgrid[, 1] - 0.04 * sgrid[, 2]
# Parameter matrix of Huesler--Reiss
# associated to power variogram
Lambda <- ((dianiso / 30) ^ 0.7) / 4
# Simulate generalized max-Pareto field above u=50
set.seed(345)
simu1 <- rgparp(
  n = 1,
  thresh = 50,
  shape = 0.1,
  riskf = "max",
  scale = tau,
  loc = eta,
  sigma = Lambda,
  model = "hr"
)
# The same, but conditional on an exceedance at a site
simu2 <- rgparp(
  n = 1,
  thresh = 50,
  shape = 0.1,
  riskf = "site",
  siteindex = 1225,
  scale = tau,
  loc = eta,
  sigma = Lambda,
  model = "hr"
)
#Plot the generalized max-Pareto field
par(mfrow = c(1, 2))
fields::quilt.plot(grid[, 1], grid[, 2], simu1, nx = 50, ny = 50)
SpatialExtremes::swiss(add = TRUE)
fields::quilt.plot(grid[, 1], grid[, 2], simu2, nx = 50, ny = 50)
SpatialExtremes::swiss(add = TRUE)
# Value at conditioning coordinate should be greater than 50
simu2[1225]

## ---- fitvariocloud, cache = TRUE,fig.width = 10, fig.height = 6, out.width = '90%', fig.align = "center"----
lon <- seq(650, 720, length = 10)
lat <- seq(215, 290, length = 10)
# Create a grid
grid <- expand.grid(lon, lat)
coord <- as.matrix(grid)
dianiso <- distg(coord, 1.5, 0.5)
sgrid <- scale(grid, scale = FALSE)
# Specify marginal parameters `loc` and `scale` over grid
eta <- 26 + 0.05 * sgrid[, 1] - 0.16 * sgrid[, 2]
tau <- 9 + 0.05 * sgrid[, 1] - 0.04 * sgrid[, 2]
# Parameter matrix of Huesler--Reiss
# associated to power variogram
Lambda <- ((dianiso / 30) ^ 0.7) / 4
# Simulate generalized max-Pareto field above u=50
set.seed(345)
simu1 <- rgparp(
  n = 1000,
  thresh = 50,
  shape = 0.1,
  riskf = "max",
  scale = tau,
  loc = eta,
  sigma = Lambda,
  model = "hr"
)
extdat <- extremo(
  dat = simu1,
  margp = 0.9,
  coord = coord,
  scale = 1.5,
  rho = 0.5,
  plot = TRUE
)

# Constrained optimization
# Minimize distance between extremal coefficient from fitted variogram
mindistpvario <- function(par, emp, coord) {
  alpha <-
    par[1]
  if (!isTRUE(all(alpha > 0, alpha < 2))) {
    return(1e10)
  }
  scale <- par[2]
  if (scale <= 0) {
    return(1e10)
  }
  a <- par[3]
  if (a < 1) {
    return(1e10)
  }
  rho <- par[4]
  if (abs(rho) >= pi / 2) {
    return(1e10)
  }
  semivariomat <-
    power.vario(distg(coord, a, rho), alpha = alpha, scale = scale)
  sum((2 * (1 - pnorm(
    sqrt(semivariomat[lower.tri(semivariomat)] / 2)
  )) - emp) ^ 2)
}
# constrained optimization for the parameters
hin <- function(par, ...) {
  c(1.99 - par[1],
    -1e-5 + par[1],
    -1e-5 + par[2],
    par[3] - 1,
    pi / 2 - par[4],
    par[4] + pi / 2)
}
opt <- alabama::auglag(
  par = c(0.5, 30, 1.5, 0.5),
  hin = hin,
  control.optim = list(parscale = c(0.5, 30, 1.5, 0.5)),
  fn = function(par) {
    mindistpvario(par, emp = extdat[, 'prob'], coord = coord)
  }
)
stopifnot(opt$kkt1, opt$kkt2)
# Plotting the extremogram in the deformed space
distfa <- distg(loc = coord, opt$par[3], opt$par[4])
plot(
  c(distfa[lower.tri(distfa)]),
  extdat[, 2],
  pch = 20,
  col = scales::alpha(1, 0.1),
  yaxs = "i",
  xaxs = "i",
  bty = 'l',
  xlab = "distance",
  ylab = "cond. prob. of exceedance",
  ylim = c(0, 1)
)
lines(
  x = (distvec <- seq(0, 200, length = 1000)),
  col = 2,
  lwd = 2,
  2 * (1 - pnorm(sqrt(
    power.vario(distvec, alpha = opt$par[1], scale = opt$par[2]) / 2
  )))
)

## ----snippetsimurparp, eval = FALSE-------------------------------------------
#  grid <- as.matrix(expand.grid(1:5, 1:5))
#  depmat <-
#    power.vario(h = distg(grid, scale = 1, rho = 0),
#                alpha = 1,
#                scale = 2) / 4
#  # This is where composition sampling shines!
#  samp <- rparpcs(
#    n = 1000,
#    shape = 0.1,
#    riskf = "min",
#    Lambda = depmat,
#    model = "br"
#  )
#  #rparp is computationally intensive with "min" - only for
#  samp2 <- rparp(
#    n = 1000,
#    shape = 0.1,
#    riskf = "max",
#    sigma = depmat,
#    model = "hr"
#  )

