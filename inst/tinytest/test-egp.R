# library(tinytest)
# set.seed(0)
# kappa <- rexp(1)
# i1 <- integrate(degp.G1, 0, 1, kappa = kappa, shape = runif(1, -1, 1))
# i2 <- integrate(degp.G2, 0, 1, kappa = kappa)
# i3 <- integrate(degp.G3, 0, 1, kappa = kappa)
# i4 <- integrate(degp.G4, 0, 1, kappa = kappa)
# i5 <- integrate(degp.G5, 0, 1, kappa = kappa)
# i6 <- integrate(degp.G6, 0, 1, kappa = kappa)
# i7 <- integrate(degp.G7, 0, 1, kappa = kappa)
# expect_equal(1, i1$value, tolerance = 1e-4)
# expect_equal(1, i2$value, tolerance = 1e-4)
# expect_equal(1, i3$value, tolerance = 1e-4)
# expect_equal(1, i4$value, tolerance = 1e-4)
# expect_equal(1, i5$value, tolerance = 1e-4)
# expect_equal(1, i6$value, tolerance = 1e-4)
# expect_equal(1, i7$value, tolerance = 1e-4)
# x <- c(NA, seq(-0.1, 1.2, by = 0.1))
# F1 <- pegp.G1(x, kappa = kappa, shape = 0.1)
# F2 <- pegp.G2(x, kappa = kappa)
# F3 <- pegp.G3(x, kappa = kappa)
# F4 <- pegp.G4(x, kappa = kappa)
# F5 <- pegp.G5(x, kappa = kappa)
# F6 <- pegp.G6(x, kappa = kappa)
# F7 <- pegp.G7(x, kappa = kappa)
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
# qp1 <- qegp.G1(
#   pegp.G1(x = x, kappa = kappa, shape = 0.1),
#   shape = 0.1,
#   kappa = kappa
# )
# qp2 <- qegp.G2(pegp.G2(x = x, kappa = kappa), kappa = kappa)
# qp3 <- qegp.G3(pegp.G3(x = x, kappa = kappa), kappa = kappa)
# qp4 <- qegp.G4(pegp.G4(x = x, kappa = kappa), kappa = kappa)
# qp5 <- qegp.G5(pegp.G5(x = x, kappa = kappa), kappa = kappa)
# qp6 <- qegp.G6(pegp.G6(x = x, kappa = kappa), kappa = kappa)
# qp7 <- qegp.G7(pegp.G7(x = x, kappa = kappa), kappa = kappa)
# expect_equal(qp1, x)
# expect_equal(qp2, x)
# expect_equal(qp3, x)
# expect_equal(qp4, x)
# expect_equal(qp5, x)
# expect_equal(qp6, x)
# expect_equal(qp7, x)
#
# x <- runif(100)
# pq1 <- pegp.G1(
#   qegp.G1(x = x, kappa = kappa, shape = 0.1),
#   shape = 0.1,
#   kappa = kappa
# )
# pq2 <- pegp.G2(qegp.G2(x = x, kappa = kappa), kappa = kappa)
# pq3 <- pegp.G3(qegp.G3(x = x, kappa = kappa), kappa = kappa)
# pq4 <- pegp.G4(qegp.G4(x = x, kappa = kappa), kappa = kappa)
# pq5 <- pegp.G5(qegp.G5(x = x, kappa = kappa), kappa = kappa)
# pq6 <- pegp.G6(qegp.G6(x = x, kappa = kappa), kappa = kappa)
# pq7 <- pegp.G7(qegp.G7(x = x, kappa = kappa), kappa = kappa)
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
#     numDeriv::grad(func = pegp.G1, x, shape = 0.1, kappa = kappa),
#     degp.G1(x, kappa, shape = 0.1)
#   )
#   expect_equal(
#     numDeriv::grad(func = pegp.G2, x, kappa = kappa),
#     degp.G2(x, kappa)
#   )
#   expect_equal(
#     numDeriv::grad(func = pegp.G3, x, kappa = kappa),
#     degp.G3(x, kappa)
#   )
#   expect_equal(
#     numDeriv::grad(func = pegp.G4, x, kappa = kappa),
#     degp.G4(x, kappa)
#   )
#   expect_equal(
#     numDeriv::grad(func = pegp.G5, x, kappa = kappa),
#     degp.G5(x, kappa)
#   )
#   expect_equal(
#     numDeriv::grad(func = pegp.G6, x, kappa = kappa),
#     degp.G6(x, kappa)
#   )
#   expect_equal(
#     numDeriv::grad(func = pegp.G7, x, kappa = kappa),
#     degp.G7(x, kappa)
#   )
# }
#
# expect_equal(pegp.G1(x, kappa = 1, shape = runif(1)), x)
# expect_equal(pegp.G2(x, kappa = 1), x)
# expect_equal(pegp.G3(x, kappa = 1), x)
# expect_equal(pegp.G4(x, kappa = 1), x)
# expect_equal(pegp.G5(x, kappa = 0), x)
# expect_equal(pegp.G6(x, kappa = 1), x)
# expect_equal(pegp.G7(x, kappa = 0), x)
