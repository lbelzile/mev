library(tinytest)

set.seed(12345)
n <- 1000L
exc <- mev::rgp(n = n, scale = 100, shape = 0.1)
expect_equal(
  stein_gp_lik(c(100, 0.1), xdat = exc, weights = rep(1, n)),
  mev::gpd.ll(par = c(100, 0.1), dat = exc)
)

opt_vanilla <- fit.wgpd(xdat = exc, threshold = 0.1, weightfun = function(n) {
  rep(1, n)
})
opt_mev <- mev::fit.gpd(xdat = exc, threshold = 0.1)
expect_equal(
  as.numeric(coef(opt_mev)),
  as.numeric(opt_vanilla$param),
  tolerance = 1e-5
)
# Check log likelihood
expect_equal(
  as.numeric(opt_mev$nllh),
  opt_vanilla$nllh
)

# Check with weighting
opt <- fit.wgpd(xdat = exc)
expect_equal(
  opt$nllh,
  -stein_gp_lik(
    pars = opt$estimate,
    xdat = opt$exceedances,
    weights = Stein_weights(opt$nat, gamma = 1)
  )
)
