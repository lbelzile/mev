context("Maximum likelihood estimation")

set.seed(202010)
xdat <- evd::rgev(n = 1000, shape = -0.5)

test_that("`fpar` input for `fit.gev`", {
  expect_error(fit.gev(xdat, fpar = 2)) # only lists accepted
  expect_error(fit.gev(xdat, fpar = list(2))) # arguments must be named
  expect_error(fit.gev(xdat, fpar = list(scale = 1, shape=0, loc=0))) #cannot pass length 3
})

# Check user passing list/named vector / vector for start
# checked fixed param yields correct result
# check that one-D optim is not problematic
