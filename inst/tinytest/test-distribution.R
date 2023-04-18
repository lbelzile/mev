# Tests for distribution functions

set.seed(2023)
shape <- runif(n = 1, min = -1, max = 1)
scale <- rexp(n = 1)
loc <- runif(n = 1, min = -10, max = 10)
upper <- ifelse(shape < 0, loc-scale/shape, Inf)
lower <- ifelse(shape > 0, loc-scale/shape, -Inf)
samp <- mev::rgev(n = 1000, loc = loc, scale = scale, shape = shape)

# Check that pgev and qgev are inverse function (and same length!)

u <- seq(0, 1, by = 0.02)
expect_equal(mev::pgev(mev::qgev(u, loc = loc, scale = scale, shape = shape),
          loc = loc, scale = scale, shape = shape), u)
expect_equal(mev::pgev(mev::qgev(u, loc = loc, scale = scale, shape = 0),
                       loc = loc, scale = scale, shape = 0), u)
expect_equal(mev::pgev(mev::qgev(u, loc = loc, scale = scale, shape = 1e-8),
                       loc = loc, scale = scale, shape = 1e-8), u)
expect_equal(mev::pgev(mev::qgev(u, loc = loc, scale = scale, shape = -1e-8),
                       loc = loc, scale = scale, shape = -1e-8), u)
# Vectorized arguments
expect_equal(mev::pgev(mev::qgev(u, loc = rep(loc, length(u)), scale = rep(scale, length(u)), shape = shape),
                       loc = loc, scale = scale, shape = shape), u)
if(shape < 0){
  expect_equal(mev::qgev(1, loc = loc, scale = scale, shape = shape),
               loc - scale/shape)
} else if(shape > 0){
  expect_equal(mev::qgev(0, loc = loc, scale = scale, shape = shape),
               loc - scale/shape)
}

# Check that function fails when arguments are incorrect

expect_error(mev::qgev(u, loc = rep(loc,2), scale = scale, shape = shape))
expect_equal(mev::pgev(c(-Inf, Inf, NaN, NA)), c(0, 1, NaN, NA))
expect_equal(mev::pgev(c(-Inf, Inf, NaN, NA)), c(0, 1, NaN, NA), shape = 0)
# Check that dgev evaluates to zero outside of support
expect_equal(mev::dgev(c(-Inf, loc - scale/shape - sign(shape) * 1e-4, Inf),
                       loc = loc, scale = scale, shape = shape),
               rep(0, 3))

# Check that quantile function qgev evaluates to zero before lower endpoint
# then is monotone increasing
  expect_true(isTRUE(all(diff(pgev(seq(-10, 10, length.out = 100), loc = loc, scale = scale, shape = shape)) >= 0))) # only lists accepted

  ## REPEAT TESTS, this time for generalized Pareto

  expect_equal(mev::pgp(mev::qgp(u, loc = loc, scale = scale, shape = shape),
                         loc = loc, scale = scale, shape = shape), u)
  expect_equal(mev::pgp(mev::qgp(u, loc = loc, scale = scale, shape = 0),
                         loc = loc, scale = scale, shape = 0), u)
  expect_equal(mev::pgp(mev::qgp(u, loc = loc, scale = scale, shape = 1e-8),
                         loc = loc, scale = scale, shape = 1e-8), u)

  # Vectorized arguments
  expect_equal(mev::pgp(mev::qgp(u, loc = loc, scale = rep(scale, length(u)), shape = shape),
                         loc = loc, scale = scale, shape = shape), u)
  if(shape < 0){
    expect_equal(mev::qgp(1, loc = loc, scale = scale, shape = shape),
                 loc - scale/shape)
  }
    expect_equal(mev::qgp(0, loc = loc, scale = scale, shape = shape),
                 loc)

  # Check that function fails when arguments are incorrect

  expect_error(mev::qgp(u, loc = rep(loc,2), scale = scale, shape = shape))
  expect_equal(mev::pgp(c(-Inf, Inf, NaN, NA)), c(0, 1, NaN, NA))
  expect_equal(mev::pgp(c(-Inf, Inf, NaN, NA)), c(0, 1, NaN, NA), shape = 0)
  # Check that dgev evaluates to zero outside of support

  # Check that quantile function qgev evaluates to zero before lower endpoint
  # then is monotone increasing
  expect_true(isTRUE(all(diff(mev::pgp(seq(-10, 10, length.out = 100),
                                       loc = loc, scale = scale, shape = shape)) >= 0))) # only lists accepted

  shapeneg <- -runif(1)
  expect_equal(mev::dgp(c(-Inf, loc - scale/shapeneg + 1e-4, Inf),
                        loc = loc, scale = scale, shape = shapeneg),
               rep(0, 3))
