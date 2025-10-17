# Tests for distribution functions

set.seed(2023)
nobs <- 100L
shape <- runif(n = nobs, min = -1, max = 1)
scale <- rexp(n = nobs)
loc <- runif(n = nobs, min = -10, max = 10)
upper <- ifelse(shape < 0, loc - scale / shape, Inf)
lower <- ifelse(shape > 0, loc - scale / shape, -Inf)
samp <- rgev(n = nobs, loc = loc, scale = scale, shape = shape)

# Check that pgev and qgev are inverse function (and same length!)

u <- seq(0, 1, length.out = nobs)
expect_equal(
  mev::pgev(
    mev::qgev(u, loc = loc, scale = scale, shape = shape),
    loc = loc,
    scale = scale,
    shape = shape
  ),
  u
)
expect_equal(
  mev::pgev(
    mev::qgev(u, loc = loc, scale = scale, shape = 0),
    loc = loc,
    scale = scale,
    shape = 0
  ),
  u
)
expect_equal(
  mev::pgev(
    mev::qgev(u, loc = loc, scale = scale, shape = 1e-7),
    loc = loc,
    scale = scale,
    shape = 1e-7
  ),
  u
)
expect_equal(
  mev::pgev(
    mev::qgev(u, loc = loc, scale = scale, shape = -1e-7),
    loc = loc,
    scale = scale,
    shape = -1e-7
  ),
  u
)
# Vectorized arguments
expect_equal(
  mev::pgev(
    mev::qgev(
      u,
      loc = rep(loc, length.out = length(u)),
      scale = rep(scale, length.out = length(u)),
      shape = rep(shape, length.out = length(u))
    ),
    loc = rep(loc, length.out = length(u)),
    scale = rep(scale, length.out = length(u)),
    shape = rep(shape, length.out = length(u))
  ),
  u
)

expect_equal(
  mev::qgev(ifelse(shape < 0, 1, 0), loc = loc, scale = scale, shape = shape),
  loc - scale / shape
)


# Check that function fails when arguments are incorrect

expect_error(mev::qgev(u, loc = rep(loc, 2), scale = scale, shape = shape))
expect_equal(mev::pgev(c(-Inf, Inf, NaN, NA)), c(0, 1, NaN, NA))
expect_equal(mev::pgev(c(-Inf, Inf, NaN, NA)), c(0, 1, NaN, NA), shape = 0)
# Check that dgev evaluates to zero outside of support
expect_equal(
  mev::dgev(
    c(-Inf, loc[1] - scale[1] / shape[1] - sign(shape[1]) * 1e-4, Inf),
    loc = loc[1],
    scale = scale[1],
    shape = shape[1]
  ),
  rep(0, 3)
)

# Check that quantile function qgev evaluates to zero before lower endpoint
# then is monotone increasing
expect_true(isTRUE(all(
  diff(pgev(
    seq(-10, 10, length.out = 100),
    loc = loc[1],
    scale = scale[1],
    shape = shape[1]
  )) >=
    0
))) # only lists accepted

## REPEAT TESTS, this time for generalized Pareto

expect_equal(
  mev::pgp(
    mev::qgp(u, loc = loc, scale = scale, shape = shape),
    loc = loc,
    scale = scale,
    shape = shape
  ),
  u
)
expect_equal(
  mev::pgp(
    mev::qgp(u, loc = loc, scale = scale, shape = 0),
    loc = loc,
    scale = scale,
    shape = 0
  ),
  u
)
expect_equal(
  mev::pgp(
    mev::qgp(u, loc = loc, scale = scale, shape = 1e-8),
    loc = loc,
    scale = scale,
    shape = 1e-8
  ),
  u
)


if (isTRUE(any(shape < 0))) {
  neg <- shape < 0
  expect_equal(
    mev::qgp(
      rep(1, sum(neg)),
      loc = loc[neg],
      scale = scale[neg],
      shape = shape[neg]
    ),
    loc[neg] - scale[neg] / shape[neg]
  )
}
expect_equal(
  mev::qgp(rep(0, length(loc)), loc = loc, scale = scale, shape = shape),
  loc
)

# Check that function fails when arguments are incorrect

expect_error(mev::qgp(u, loc = rep(loc, 2), scale = scale, shape = shape))
expect_equal(mev::pgp(c(-Inf, Inf, NaN, NA)), c(0, 1, NaN, NA))
expect_equal(mev::pgp(c(-Inf, Inf, NaN, NA)), c(0, 1, NaN, NA), shape = 0)
# Check that dgev evaluates to zero outside of support

# Check that quantile function qgev evaluates to zero before lower endpoint
# then is monotone increasing
expect_true(isTRUE(all(
  diff(mev::pgp(
    seq(-10, 10, length.out = 100),
    loc = loc[1],
    scale = scale[1],
    shape = shape[1]
  )) >=
    0
))) # only lists accepted

shapeneg <- -runif(1)
expect_equal(
  mev::dgp(
    x = c(-Inf, loc[1] - scale[1] / shapeneg + 1e-4, Inf),
    loc = loc[1],
    scale = scale[1],
    shape = shapeneg
  ),
  rep(0, 3)
)
