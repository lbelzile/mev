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


