set.seed(1234)
theta <- 0.2
shape <- 0.5
X <- mev::rmar1(n = 10000, theta = theta, shape = shape)
theta_hat <- mev::xdep.xindex(
  xdat = X,
  qlev = 0.95,
  estimator = "mle",
  confint = "lrt"
)
# Expect the extremal index to coincide
tinytest::expect_true((theta_hat$lower < theta) & (theta_hat$upper > theta))
# Expect marginal parameters to be as specified
fit <- mev::fit.gev(X)
dev <- -2 *
  (fit$nllh +
    mev::gev.ll(
      dat = X,
      c(1, shape, shape)
    ))
pval <- pchisq(q = dev, df = 3, lower.tail = FALSE)
tinytest::expect_true(pval > 0.05)
# Check marginal distribution (bis)
tinytest::expect_true(
  ks.test(exp(-(X^(-1 / shape))), "punif")$p.value > 0.05
)

set.seed(1234)
theta <- 0.2
X <- rmar1(n = 10000, theta = theta, shape = 0)
theta_hat <- mev::xdep.xindex(
  xdat = X,
  qlev = 0.95,
  estimator = "mle",
  confint = "lrt"
)
# Expect the extremal index to coincide
tinytest::expect_true((theta_hat$lower < theta) & (theta_hat$upper > theta))
# Expect marginal parameters to be as specified
fit <- mev::fit.gev(X)
dev <- -2 *
  (fit$nllh +
    mev::gev.ll(
      dat = X,
      c(0, 1, 0)
    ))
pval <- pchisq(q = dev, df = 3, lower.tail = FALSE)
tinytest::expect_true(pval > 0.05)
# Check marginal distribution (bis)
tinytest::expect_true(
  ks.test(mev::pgev(X, 0, 1, 0), "punif")$p.value > 0.05
)
