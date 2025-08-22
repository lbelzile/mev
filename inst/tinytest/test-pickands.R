# xdat <- abs(rnorm(n = 1000))
#
# (pwm_lmom <- lmom::samlmu(x = xdat, ratio = FALSE))
# (lmom_lmom <- lmom::samlmu(x = xdat, ratios = TRUE))
# pwm(xdat)
# lmoments(xdat)
# c(1, -1) * lmom::pelgpa(bound = 0, lmom::samlmu(x = xdat))[-1]
# gpd.lmom(xdat)
# extRemes::fevd(
#   x = xdat,
#   threshold = 0,
#   type = "GP",
#   method = "Lmoments"
# )
#
# gpd.lmom(xdat, Lskew = TRUE)
#
# extRemes::fevd(
#   x = xdat,
#   threshold = 0,
#   type = "GP",
#   method = "MLE"
# )
# gpd.mle(xdat, args = c("scale", "shape"))
#
#
# xdat <- rexp(n = 1000, rate = 2)
# lmom::pelgpa(bound = 0, lmom::samlmu(x = xdat, ratios = FALSE))
# gpd.lmom(xdat)
#
# isTRUE(all.equal(
#   lmoments(xdat),
#   lmom::samlmu(x = xdat, ratios = FALSE),
#   check.attributes = FALSE
# ))
#
# test <- thselect.alrs(
#   xdat,
#   thresh = quantile(xdat, seq(0.8, 0.99, by = 0.01))
# )
# print(test)
# plot(test)
#
# test <- thselect.pickands(xdat, method = "lmom")
# plot(test$dist, x = test$thresh)
#
# test <- thselect.pickands(xdat, method = "mle")
# plot(test$dist, x = test$thresh)
#
# test <- thselect.pickands(xdat, method = "quartiles")
# plot(test$dist, x = test$thresh)
