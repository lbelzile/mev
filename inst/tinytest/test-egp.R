# library(tinytest)
# set.seed(0)
# kappa <- rexp(1)
# i1 <- integrate(mev:::degp.G1, 0, 1, kappa = kappa, shape = runif(1, -1, 1))
# i2 <- integrate(mev:::degp.G2, 0, 1, kappa = kappa)
# i3 <- integrate(mev:::degp.G3, 0, 1, kappa = kappa)
# i4 <- integrate(mev:::degp.G4, 0, 1, kappa = kappa)
# i5 <- integrate(mev:::degp.G5, 0, 1, kappa = kappa)
# i6 <- integrate(mev:::degp.G6, 0, 1, kappa = kappa)
# i7 <- integrate(mev:::degp.G7, 0, 1, kappa = kappa)
# expect_equal(1, i1$value, tolerance = 1e-4)
# expect_equal(1, i2$value, tolerance = 1e-4)
# expect_equal(1, i3$value, tolerance = 1e-4)
# expect_equal(1, i4$value, tolerance = 1e-4)
# expect_equal(1, i5$value, tolerance = 1e-4)
# expect_equal(1, i6$value, tolerance = 1e-4)
# expect_equal(1, i7$value, tolerance = 1e-4)
# x <- c(NA, seq(-0.1, 1.2, by = 0.1))
# F1 <- mev:::pegp.G1(x, kappa = kappa, shape = 0.1)
# F2 <- mev:::pegp.G2(x, kappa = kappa)
# F3 <- mev:::pegp.G3(x, kappa = kappa)
# F4 <- mev:::pegp.G4(x, kappa = kappa)
# F5 <- mev:::pegp.G5(x, kappa = kappa)
# F6 <- mev:::pegp.G6(x, kappa = kappa)
# F7 <- mev:::pegp.G7(x, kappa = kappa)
# expect_true(is.na(F1[1]))
# expect_true(isTRUE(all(diff(F1)[-1] >= 0)))
# expect_equal(tail(F1, 1), 1)
# expect_equal(F1[2], 0)
# expect_true(is.na(F2[1]))
# expect_true(isTRUE(all(diff(F2)[-1] >= 0)))
# expect_equal(tail(F2, 1), 1)
# expect_equal(F2[2], 0)
# expect_true(is.na(F3[1]))
# expect_true(isTRUE(all(diff(F3)[-1] >= 0)))
# expect_equal(tail(F3, 1), 1)
# expect_equal(F3[2], 0)
# expect_true(is.na(F4[1]))
# expect_true(isTRUE(all(diff(F4)[-1] >= 0)))
# expect_equal(tail(F4, 1), 1)
# expect_equal(F4[2], 0)
# expect_true(is.na(F5[1]))
# expect_true(isTRUE(all(diff(F5)[-1] >= 0)))
# expect_equal(tail(F5, 1), 1)
# expect_equal(F5[2], 0)
# expect_true(is.na(F6[1]))
# expect_true(isTRUE(all(diff(F6)[-1] >= 0)))
# expect_equal(tail(F6, 1), 1)
# expect_equal(F6[2], 0)
# expect_true(is.na(F7[1]))
# expect_true(isTRUE(all(diff(F7)[-1] >= 0)))
# expect_equal(tail(F7, 1), 1)
# expect_equal(F7[2], 0)
#
# x <- seq(0, 1, by = 0.1)
# qp1 <- mev:::qegp.G1(
#   mev:::pegp.G1(x = x, kappa = kappa, shape = 0.1),
#   shape = 0.1,
#   kappa = kappa
# )
# qp2 <- mev:::qegp.G2(mev:::pegp.G2(x = x, kappa = kappa), kappa = kappa)
# qp3 <- mev:::qegp.G3(mev:::pegp.G3(x = x, kappa = kappa), kappa = kappa)
# qp4 <- mev:::qegp.G4(mev:::pegp.G4(x = x, kappa = kappa), kappa = kappa)
# qp5 <- mev:::qegp.G5(mev:::pegp.G5(x = x, kappa = kappa), kappa = kappa)
# qp6 <- mev:::qegp.G6(mev:::pegp.G6(x = x, kappa = kappa), kappa = kappa)
# qp7 <- mev:::qegp.G7(mev:::pegp.G7(x = x, kappa = kappa), kappa = kappa)
# expect_equal(qp1, x)
# expect_equal(qp2, x)
# expect_equal(qp3, x)
# expect_equal(qp4, x)
# expect_equal(qp5, x)
# expect_equal(qp6, x)
# expect_equal(qp7, x)
#
# x <- runif(100)
# pq1 <- mev:::pegp.G1(
#   mev:::qegp.G1(x = x, kappa = kappa, shape = 0.1),
#   shape = 0.1,
#   kappa = kappa
# )
# pq2 <- mev:::pegp.G2(mev:::qegp.G2(x = x, kappa = kappa), kappa = kappa)
# pq3 <- mev:::pegp.G3(mev:::qegp.G3(x = x, kappa = kappa), kappa = kappa)
# pq4 <- mev:::pegp.G4(mev:::qegp.G4(x = x, kappa = kappa), kappa = kappa)
# pq5 <- mev:::pegp.G5(mev:::qegp.G5(x = x, kappa = kappa), kappa = kappa)
# pq6 <- mev:::pegp.G6(mev:::qegp.G6(x = x, kappa = kappa), kappa = kappa)
# pq7 <- mev:::pegp.G7(mev:::qegp.G7(x = x, kappa = kappa), kappa = kappa)
# expect_equal(pq1, x)
# expect_equal(pq2, x)
# expect_equal(pq3, x)
# expect_equal(pq4, x)
# expect_equal(pq5, x)
# expect_equal(pq6, x)
# expect_equal(pq7, x)
#
#
# x <- seq(0.01, 0.99, length.out = 11)
# if (requireNamespace("numDeriv", quietly = TRUE)) {
#   expect_equal(
#     numDeriv::grad(func = mev:::pegp.G1, x, shape = 0.1, kappa = kappa),
#     mev:::degp.G1(x, kappa, shape = 0.1)
#   )
#   expect_equal(
#     numDeriv::grad(func = mev:::pegp.G2, x, kappa = kappa),
#     mev:::degp.G2(x, kappa)
#   )
#   expect_equal(
#     numDeriv::grad(func = mev:::pegp.G3, x, kappa = kappa),
#     mev:::degp.G3(x, kappa)
#   )
#   expect_equal(
#     numDeriv::grad(func = mev:::pegp.G4, x, kappa = kappa),
#     mev:::degp.G4(x, kappa)
#   )
#   expect_equal(
#     numDeriv::grad(func = mev:::pegp.G5, x, kappa = kappa),
#     mev:::degp.G5(x, kappa)
#   )
#   expect_equal(
#     numDeriv::grad(func = mev:::pegp.G6, x, kappa = kappa),
#     mev:::degp.G6(x, kappa)
#   )
#   expect_equal(
#     numDeriv::grad(func = mev:::pegp.G7, x, kappa = kappa),
#     mev:::degp.G7(x, kappa)
#   )
# }
#
# expect_equal(mev:::pegp.G1(x, kappa = 1, shape = runif(1)), x)
# expect_equal(mev:::pegp.G2(x, kappa = 1), x)
# expect_equal(mev:::pegp.G3(x, kappa = 1), x)
# expect_equal(mev:::pegp.G4(x, kappa = 1), x)
# expect_equal(mev:::pegp.G5(x, kappa = 0), x)
# expect_equal(mev:::pegp.G6(x, kappa = 1), x)
# expect_equal(mev:::pegp.G7(x, kappa = 0), x)
