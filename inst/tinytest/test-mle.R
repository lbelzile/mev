## Maximum likelihood estimation

set.seed(202010)
xdat <- evd::rgev(n = 1000, shape = -0.5)

# GEV: `fpar` input for `fit.gev`
  expect_error(fit.gev(xdat, fpar = 2)) # only lists accepted
  expect_error(fit.gev(xdat, fpar = list(2))) # arguments must be named
  expect_error(fit.gev(xdat, fpar = list(scale = 1, shape=0, loc=0))) #cannot pass length 3


# Check that MLE is the same for fixed arguments
#
mle <- mev::fit.gev(xdat)
# GEV: profile yields same maximum as full likelihood
expect_equal(mev::fit.gev(xdat, fpar = list(loc = mle$param[1]))$param, mle$param,tolerance = 1e-5)
expect_equal(mev::fit.gev(xdat, fpar = list(scale = mle$param[2]))$param, mle$param,tolerance = 1e-5)
expect_equal(mev::fit.gev(xdat, fpar = list(shape = mle$param[3]))$param, mle$param,tolerance = 1e-5)
expect_equal(mev::fit.gev(xdat, fpar = list(scale = mle$param[2], shape = mle$param[3]))$param, mle$param,tolerance = 1e-5)
expect_equal(mev::fit.gev(xdat, fpar = list(loc = mle$param[1], shape = mle$param[3]))$param, mle$param,tolerance = 1e-5)
expect_equal(mev::fit.gev(xdat, fpar = list(loc = mle$param[1], scale = mle$param[2]))$param, mle$param,tolerance = 1e-5)



set.seed(202010)
xdat <- evd::rgpd(n = 1000, shape = -0.5)

# GP: `fpar` input for `fit.gpd`
expect_error(fit.gpd(xdat, fpar = 2)) # only lists accepted
expect_error(fit.gpd(xdat, fpar = list(2))) # arguments must be named
expect_error(fit.gpd(xdat, fpar = list(scale = 1, shape=0, loc=0))) #cannot pass length 3


# Check that MLE is the same for fixed arguments
#
mle <- mev::fit.gpd(xdat)
# GP: Profile yields same maximum as full likelihood
expect_equal(mev::fit.gpd(xdat, fpar = list(scale = mle$param[1]))$param, mle$param,tolerance = 1e-5)
expect_equal(mev::fit.gpd(xdat, fpar = list(shape = mle$param[2]))$param, mle$param,tolerance = 1e-5)

# GP: Profile yields same maximum as full likelihood
expect_equal(mev::fit.gpd(xdat, fpar = list(scale = mle$param[1]))$param, mle$param,tolerance = 1e-5)
expect_equal(mev::fit.gpd(xdat, fpar = list(shape = mle$param[2]))$param, mle$param,tolerance = 1e-5)

# Check user passing list/named vector / vector for start
# checked fixed param yields correct result
# check that one-D optim is not problematic

mle <- mev::fit.pp(xdat, np = 10)
# IPP: profile yields same maximum as full likelihood
expect_equal(mev::fit.pp(xdat,
                           np = 10,
                           fpar = list(loc = mle$param[1]))$param,
               mle$param,
               tolerance = 1e-5)
  expect_equal(mev::fit.pp(xdat,
                           np = 10,
                           fpar = list(scale = mle$param[2]))$param,
               mle$param,
               tolerance = 1e-5)
  expect_equal(mev::fit.pp(xdat,
                           np = 10,
                           fpar = list(shape = mle$param[3]))$param,
               mle$param,
               tolerance = 1e-5)
  expect_equal(mev::fit.pp(xdat,
                           np = 10,
                           fpar = list(scale = mle$param[2],
                                       shape = mle$param[3]))$param,
               mle$param,
               tolerance = 1e-5)
  expect_equal(mev::fit.pp(xdat,
                           np = 10,
                           fpar = list(loc = mle$param[1],
                                       shape = mle$param[3]))$param,
               mle$param,
               tolerance = 1e-5)
  expect_equal(mev::fit.pp(xdat,
                           np = 10,
                           fpar = list(loc = mle$param[1],
                                       scale = mle$param[2]))$param,
               mle$param,
               tolerance = 1e-5)


# mev::pp.score(par = mle$param, dat  = xdat, u = 0, np = 10)


# Check ANOVA calls for nested / non-nested models


#FOR GP, is this ever triggered?
#
# if (temp$mle[2] < -1 && temp$conv == 0) {
#   warning("The MLE is not a solution to the score equation for \"xi < -1'")
# }
