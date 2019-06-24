#' Cox-Snell first order bias for the GEV distribution
#'
#' Bias vector for the GEV distribution based on an \code{n} sample.
#' Due to numerical instability, values of the information matrix and the bias
#' are linearly interpolated when the value of the shape parameter is close to zero.
#' @inheritParams gev
#' @export
#' @keywords internal
#' @seealso \code{\link{gev}}
gev.bias <- function(par, n) {
    if (length(n) > 1) {
        stop("Invalid argument for sample size")
    }
    if (length(par) != 3) {
        stop("Invalid argument for parameter vector, must be a vector of length 3.")
    }
    sigma <- par[2]
    xi <- par[3]
    if (xi < -1/3) {
        stop("Cox-Snell correction only valid if the shape is greater than -1/3")
    }
    zeta3 <- 1.20205690315959
    zeta5 <- 1.03692775514337  #gsl::zeta(5)
    k111 <- ((1 + xi)^2 * (1 + 4 * xi) * gamma(1 + 3 * xi))/sigma^3
    k11.2 <- -2 * (xi + 1)^2 * gamma(2 * xi + 1)/sigma^3
    k11.3 <- 2 * (xi + 1)^2 * psigamma(2 * xi + 1) * gamma(2 * xi + 1)/sigma^2 + 2 * (xi + 1) * gamma(2 * xi + 1)/sigma^2
    k33.2 <- 0
    euler_gamma <- -psigamma(1)

    # Numerical tolerance 1e-10 (all.equal has this)
    if (abs(xi) > 1e-10) {
        k112 <- (xi + 1) * (gamma(2 * xi + 2) - (xi + 1) * (4 * xi + 1) * gamma(3 * xi + 1))/(sigma^3 * xi)
        k12.2 <- 2 * ((xi + 1)^2 * gamma(2 * xi + 1) - gamma(xi + 2))/(sigma^3 * xi)
    } else {
        k112 <- (euler_gamma - 3)/sigma^3
        k12.2 <- -2 * (euler_gamma - 1)/sigma^3
    }


    # Numerical tolerance 1e-6 If function is not too numerically unstable
    if (abs(xi) > 1e-06) {
        k113 <- (1 + xi) * ((1 + xi) * (1 + 4 * xi) * gamma(1 + 3 * xi) - gamma(1 + 2 * xi) * (1 + 2 * xi * (2 + xi) + xi * (1 + 2 *
            xi) * psigamma(2 + 2 * xi, 0)))/(sigma^2 * xi^2)
        k122 <- ((1 - xi) * gamma(2 + xi) - gamma(3 + 2 * xi) + (1 + xi)^2 * (1 + 4 * xi) * gamma(1 + 3 * xi))/(sigma^3 * xi^2)
        k12.3 <- -(2 * (xi + 1)^2 * psigamma(2 * xi + 1) * gamma(2 * xi + 1) + 2 * (xi + 1) * gamma(2 * xi + 1) - psigamma(xi + 2) *
            gamma(xi + 2))/(sigma^2 * xi) + ((xi + 1)^2 * gamma(2 * xi + 1) - gamma(xi + 2))/(sigma^2 * xi^2)
        k13.2 <- -(((xi + 1)^2 * gamma(2 * xi + 1) - xi * ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2))/(sigma^2 * xi^2))
        k22.2 <- -2 * ((xi + 1)^2 * gamma(2 * xi + 1) - 2 * gamma(xi + 2) + 1)/(sigma^3 * xi^2)
    } else {
        # Compute the approximation for xi=0, (limit)
        k113_0 <- -1/12 * (36 * euler_gamma - 6 * euler_gamma^2 - pi^2 - 24)/sigma^2
        k122_0 <- -1/6 * (36 * euler_gamma - 6 * euler_gamma^2 - pi^2 - 24)/sigma^3
        k12.3_0 <- (6 * euler_gamma - 3 * euler_gamma^2 - 1/2 * pi^2 - 2)/(2 * sigma^2)
        k13.2_0 <- 1/12 * (12 * euler_gamma - 6 * euler_gamma^2 - pi^2)/sigma^2
        k22.2_0 <- 1/3 * (12 * euler_gamma - 6 * euler_gamma^2 - pi^2 - 6)/sigma^3
        if (isTRUE(all.equal(xi, 0))) {
            # if xi ==0 with roughly precision 1e-8, then set value to this limit
            k113 <- k113_0
            k122 <- k122_0
            k12.3 <- k12.3_0
            k13.2 <- k13.2_0
            k22.2 <- k22.2_0
        } else {
            # if xi !=0, but numerical breakdown,
            xit <- sign(xi) * 1e-06
            k113_l <- sapply(xit, function(xi) {
                (1 + xi) * ((1 + xi) * (1 + 4 * xi) * gamma(1 + 3 * xi) - gamma(1 + 2 * xi) * (1 + 2 * xi * (2 + xi) + xi * (1 + 2 *
                  xi) * psigamma(2 + 2 * xi, 0)))/(sigma^2 * xi^2)
            })
            k122_l <- sapply(xit, function(xi) {
                ((1 - xi) * gamma(2 + xi) - gamma(3 + 2 * xi) + (1 + xi)^2 * (1 + 4 * xi) * gamma(1 + 3 * xi))/(sigma^3 * xi^2)
            })
            k12.3_l <- sapply(xit, function(xi) {
                -(2 * (xi + 1)^2 * psigamma(2 * xi + 1) * gamma(2 * xi + 1) + 2 * (xi + 1) * gamma(2 * xi + 1) - psigamma(xi + 2) *
                  gamma(xi + 2))/(sigma^2 * xi) + ((xi + 1)^2 * gamma(2 * xi + 1) - gamma(xi + 2))/(sigma^2 * xi^2)
            })
            k13.2_l <- sapply(xit, function(xi) {
                -(((xi + 1)^2 * gamma(2 * xi + 1) - xi * ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2))/(sigma^2 * xi^2))
            })
            k22.2_l <- sapply(xit, function(xi) {
                -2 * ((xi + 1)^2 * gamma(2 * xi + 1) - 2 * gamma(xi + 2) + 1)/(sigma^3 * xi^2)
            })
            # use linear interpolation between two values 0 and 1e-6
            k113 <- approx(x = c(0, xit), y = c(k113_0, k113_l), xout = xi)$y
            k122 <- approx(x = c(0, xit), y = c(k122_0, k122_l), xout = xi)$y
            k12.3 <- approx(x = c(0, xit), y = c(k12.3_0, k12.3_l), xout = xi)$y
            k13.2 <- approx(x = c(0, xit), y = c(k13.2_0, k13.2_l), xout = xi)$y
            k22.2 <- approx(x = c(0, xit), y = c(k22.2_0, k22.2_l), xout = xi)$y
        }

    }


    # Numerical tolerance 1e-4 If function is not too numerically unstable
    if (abs(xi) > 1e-04) {
        k123 <- ((-gamma(2 + xi)) * (1 + 2 * xi + xi * psigamma(2 + xi, 0)) + (1 + xi) * ((-(1 + 5 * xi + 4 * xi^2)) * gamma(1 + 3 *
            xi) + gamma(1 + 2 * xi) * (2 + 7 * xi + 3 * xi^2 + xi * (1 + 2 * xi) * psigamma(2 + 2 * xi, 0))))/(sigma^2 * xi^3)
        k222 <- (1 - 3 * xi + 3 * (xi - 1) * gamma(2 + xi) + 1.5 * gamma(3 + 2 * xi) - (1 + xi)^2 * (1 + 4 * xi) * gamma(1 + 3 * xi))/(sigma^3 *
            xi^3)
        k13.3 <- -(-(2 * (xi + 1)^2 * psigamma(2 * xi + 1) * gamma(2 * xi + 1) - xi * ((xi + 1)/xi + psigamma(xi + 1)) * psigamma(xi +
            2) * gamma(xi + 2) + xi * ((xi + 1)/xi^2 - 1/xi - psigamma(xi + 1, 1)) * gamma(xi + 2) + 2 * (xi + 1) * gamma(2 * xi +
            1) - ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2))/(sigma * xi^2) + 2 * ((xi + 1)^2 * gamma(2 * xi + 1) - xi * ((xi +
            1)/xi + psigamma(xi + 1)) * gamma(xi + 2))/(sigma * xi^3))
        k22.3 <- 2 * ((xi + 1)^2 * psigamma(2 * xi + 1) * gamma(2 * xi + 1) + (xi + 1) * gamma(2 * xi + 1) - psigamma(xi + 2) * gamma(xi +
            2))/(sigma^2 * xi^2) - 2 * ((xi + 1)^2 * gamma(2 * xi + 1) - 2 * gamma(xi + 2) + 1)/(sigma^2 * xi^3)
        k23.2 <- -(-((digamma(1)) + (xi + 1)^2 * gamma(2 * xi + 1)/xi - ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2) - (gamma(xi +
            2) - 1)/xi + 1)/(sigma^2 * xi^2))
    } else {
        # Compute the approximation for xi=0, (limit)
        k123_0 <- 1/12 * (60 * euler_gamma + 6 * euler_gamma^3 - euler_gamma * pi^2 + 4 * pi^2 * (euler_gamma - 1) - 48 * euler_gamma^2 -
            4 * pi^2 + 12 * zeta3 - 12)/sigma^2
        k222_0 <- 1/4 * (48 * euler_gamma + 4 * euler_gamma^3 + 9 * euler_gamma * pi^2 - 4 * pi^2 * (2 * euler_gamma - 3) + pi^2 *
            (euler_gamma - 1) - 36 * euler_gamma^2 - 17 * pi^2 + 8 * zeta3 - 16)/sigma^3
        k13.3_0 <- (-6 * euler_gamma - 4 * euler_gamma^3 - 7/2 * euler_gamma * pi^2 + 3/2 * pi^2 * (euler_gamma - 1) + 12 * euler_gamma^2 +
            7/2 * pi^2 - 8 * zeta3)/(6 * sigma)
        k22.3_0 <- -1/6 * (12 * euler_gamma + 6 * euler_gamma^3 + 4 * euler_gamma * pi^2 - pi^2 * (euler_gamma - 1) - 18 * euler_gamma^2 -
            4 * pi^2 + 12 * zeta3)/sigma^2
        k23.2_0 <- -1/12 * (12 * euler_gamma + 6 * euler_gamma^3 - euler_gamma * pi^2 + 4 * pi^2 * (euler_gamma - 1) - 18 * euler_gamma^2 +
            pi^2 + 12 * zeta3)/sigma^2
        if (isTRUE(all.equal(xi, 0))) {
            # if xi numerically zero, use the latter
            k123 <- k123_0
            k222 <- k222_0
            k13.3 <- k13.3_0
            k22.3 <- k22.3_0
            k23.2 <- k23.2_0
        } else {
            # else if less than tol, but not zero, interpolate linearly
            xit <- sign(xi) * 1e-04
            k123_l <- sapply(xit, function(xi) {
                ((-gamma(2 + xi)) * (1 + 2 * xi + xi * psigamma(2 + xi, 0)) + (1 + xi) * ((-(1 + 5 * xi + 4 * xi^2)) * gamma(1 + 3 *
                  xi) + gamma(1 + 2 * xi) * (2 + 7 * xi + 3 * xi^2 + xi * (1 + 2 * xi) * psigamma(2 + 2 * xi, 0))))/(sigma^2 * xi^3)
            })
            k222_l <- sapply(xit, function(xi) {
                (1 - 3 * xi + 3 * (xi - 1) * gamma(2 + xi) + 1.5 * gamma(3 + 2 * xi) - (1 + xi)^2 * (1 + 4 * xi) * gamma(1 + 3 * xi))/(sigma^3 *
                  xi^3)
            })
            k13.3_l <- sapply(xit, function(xi) {
                -(-(2 * (xi + 1)^2 * psigamma(2 * xi + 1) * gamma(2 * xi + 1) - xi * ((xi + 1)/xi + psigamma(xi + 1)) * psigamma(xi +
                  2) * gamma(xi + 2) + xi * ((xi + 1)/xi^2 - 1/xi - psigamma(xi + 1, 1)) * gamma(xi + 2) + 2 * (xi + 1) * gamma(2 *
                  xi + 1) - ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2))/(sigma * xi^2) + 2 * ((xi + 1)^2 * gamma(2 * xi + 1) -
                  xi * ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2))/(sigma * xi^3))
            })
            k22.3_l <- sapply(xit, function(xi) {
                2 * ((xi + 1)^2 * psigamma(2 * xi + 1) * gamma(2 * xi + 1) + (xi + 1) * gamma(2 * xi + 1) - psigamma(xi + 2) * gamma(xi +
                  2))/(sigma^2 * xi^2) - 2 * ((xi + 1)^2 * gamma(2 * xi + 1) - 2 * gamma(xi + 2) + 1)/(sigma^2 * xi^3)
            })
            k23.2_l <- sapply(xit, function(xi) {
                -(-((digamma(1)) + (xi + 1)^2 * gamma(2 * xi + 1)/xi - ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2) - (gamma(xi +
                  2) - 1)/xi + 1)/(sigma^2 * xi^2))
            })
            # Linear interpolation
            k123 <- approx(x = c(0, xit), y = c(k123_0, k123_l), xout = xi)$y
            k222 <- approx(x = c(0, xit), y = c(k222_0, k222_l), xout = xi)$y
            k13.3 <- approx(x = c(0, xit), y = c(k13.3_0, k13.3_l), xout = xi)$y
            k22.3 <- approx(x = c(0, xit), y = c(k22.3_0, k22.3_l), xout = xi)$y
            k23.2 <- approx(x = c(0, xit), y = c(k23.2_0, k23.2_l), xout = xi)$y
        }
    }

    # Numerical tolerance 1e-3 If function is not too numerically unstable
    if (abs(xi) > 0.001) {
        k23.3 <- -((2 * (xi + 1)^2 * psigamma(2 * xi + 1) * gamma(2 * xi + 1)/xi - ((xi + 1)/xi + psigamma(xi + 1)) * psigamma(xi +
            2) * gamma(xi + 2) + ((xi + 1)/xi^2 - 1/xi - psigamma(xi + 1, 1)) * gamma(xi + 2) - (xi + 1)^2 * gamma(2 * xi + 1)/xi^2 +
            2 * (xi + 1) * gamma(2 * xi + 1)/xi - psigamma(xi + 2) * gamma(xi + 2)/xi + (gamma(xi + 2) - 1)/xi^2)/(sigma * xi^2) -
            2 * ((digamma(1)) + (xi + 1)^2 * gamma(2 * xi + 1)/xi - ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2) - (gamma(xi +
                2) - 1)/xi + 1)/(sigma * xi^3))
        k133 <- ((1 + xi)^2 * (1 + 6 * xi + 8 * xi^2) * gamma(1 + 3 * xi) - gamma(3 + 2 * xi) * (1 + 5 * xi + 3 * xi^2 + xi * (1 +
            2 * xi) * psigamma(2 + 2 * xi, 0)) + (1 + 2 * xi) * gamma(1 + xi) * (1 + 6 * xi + 5 * xi^2 + 2 * xi^3 + 2 * xi * (1 +
            3 * xi + 2 * xi^2) * psigamma(2 + xi, 0) + xi^2 * (1 + xi) * psigamma(2 + xi, 0)^2 + (xi^2) * (1 + xi) * psigamma(2 +
            xi, 1)))/(sigma * xi^4 * (1 + 2 * xi))
        k223 <- -(1 + 2 * xi + digamma(1) * xi - xi^2 - digamma(1) * xi^2 - 27 * xi^3 * gamma(3 * xi) - 12 * xi^4 * gamma(3 * xi) -
            3 * gamma(2 + xi) - 5 * xi * gamma(2 + xi) + 3 * (1 + xi) * gamma(1 + 2 * xi) + 10 * xi * (1 + xi) * gamma(1 + 2 * xi) +
            4 * xi^2 * (1 + xi) * gamma(1 + 2 * xi) - gamma(1 + 3 * xi) - 6 * xi * gamma(1 + 3 * xi) - 2 * xi * gamma(2 + xi) * psigamma(2 +
            xi, 0) + xi * (1 + xi) * (1 + 2 * xi) * gamma(1 + 2 * xi) * psigamma(2 + 2 * xi, 0))/(sigma^2 * xi^4)
        k233 <- (1 + 7 * xi + 2 * digamma(1) * xi + 4 * xi^2 + 6 * digamma(1) * xi^2 + digamma(1)^2 * xi^2 + (pi^2 * xi^2)/6 - 3 *
            xi^3 * (9 + 4 * xi) * gamma(3 * xi) + 3 * gamma(1 + 2 * xi) + 17 * xi * gamma(1 + 2 * xi) + 22 * xi^2 * gamma(1 + 2 *
            xi) + 8 * xi^3 * gamma(1 + 2 * xi) - (1 + 6 * xi) * gamma(1 + 3 * xi) + (1 + 3 * xi + 2 * xi^2) * 2 * xi * gamma(1 + 2 *
            xi) * psigamma(2 + 2 * xi, 0) - gamma(1 + xi) * (3 + 16 * xi + 13 * xi^2 + 4 * xi^3 + 2 * xi * (2 + 5 * xi + 3 * xi^2) *
            psigamma(2 + xi, 0) + xi^2 * (1 + xi) * psigamma(2 + xi, 0)^2 + xi^2 * (1 + xi) * psigamma(2 + xi, 1)))/(sigma * xi^5)
    } else {
        k23.3_0 <- -3.70965809351906/sigma
        k133_0 <- 0.10683192718888/sigma
        k223_0 <- 1/40 * (20 * euler_gamma^4 + 3 * pi^4 - 200 * euler_gamma^3 + 20 * euler_gamma^2 * (pi^2 + 18) + 60 * pi^2 - 20 *
            euler_gamma * (5 * pi^2 - 8 * zeta3 + 8) - 400 * zeta3)/sigma^2
        k233_0 <- 1/48 * (12 * euler_gamma^5 - 140 * euler_gamma^4 - 21 * pi^4 + 20 * euler_gamma^3 * (pi^2 + 16) - 4 * euler_gamma^2 *
            (35 * pi^2 - 60 * zeta3 + 48) + 8 * pi^2 * (5 * zeta3 - 4) + euler_gamma * (9 * pi^4 + 160 * pi^2 - 1120 * zeta3) + 288 *
            zeta5 + 640 * zeta3)/sigma
        if (isTRUE(all.equal(xi, 0))) {
            k23.3 <- k23.3_0
            k133 <- k133_0
            k223 <- k223_0
            k233 <- k233_0
        } else {
            xit <- sign(xi) * 0.01
            k23.3_l <- sapply(xit, function(xi) {
                -((2 * (xi + 1)^2 * psigamma(2 * xi + 1) * gamma(2 * xi + 1)/xi - ((xi + 1)/xi + psigamma(xi + 1)) * psigamma(xi +
                  2) * gamma(xi + 2) + ((xi + 1)/xi^2 - 1/xi - psigamma(xi + 1, 1)) * gamma(xi + 2) - (xi + 1)^2 * gamma(2 * xi +
                  1)/xi^2 + 2 * (xi + 1) * gamma(2 * xi + 1)/xi - psigamma(xi + 2) * gamma(xi + 2)/xi + (gamma(xi + 2) - 1)/xi^2)/(sigma *
                  xi^2) - 2 * ((digamma(1)) + (xi + 1)^2 * gamma(2 * xi + 1)/xi - ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2) -
                  (gamma(xi + 2) - 1)/xi + 1)/(sigma * xi^3))
            })
            k133_l <- sapply(xit, function(xi) {
                ((1 + xi)^2 * (1 + 6 * xi + 8 * xi^2) * gamma(1 + 3 * xi) - gamma(3 + 2 * xi) * (1 + 5 * xi + 3 * xi^2 + xi * (1 +
                  2 * xi) * psigamma(2 + 2 * xi, 0)) + (1 + 2 * xi) * gamma(1 + xi) * (1 + 6 * xi + 5 * xi^2 + 2 * xi^3 + 2 * xi *
                  (1 + 3 * xi + 2 * xi^2) * psigamma(2 + xi, 0) + xi^2 * (1 + xi) * psigamma(2 + xi, 0)^2 + (xi^2) * (1 + xi) * psigamma(2 +
                  xi, 1)))/(sigma * xi^4 * (1 + 2 * xi))
            })
            k223_l <- sapply(xit, function(xi) {
                -(1 + 2 * xi + digamma(1) * xi - xi^2 - digamma(1) * xi^2 - 27 * xi^3 * gamma(3 * xi) - 12 * xi^4 * gamma(3 * xi) -
                  3 * gamma(2 + xi) - 5 * xi * gamma(2 + xi) + 3 * (1 + xi) * gamma(1 + 2 * xi) + 10 * xi * (1 + xi) * gamma(1 + 2 *
                  xi) + 4 * xi^2 * (1 + xi) * gamma(1 + 2 * xi) - gamma(1 + 3 * xi) - 6 * xi * gamma(1 + 3 * xi) - 2 * xi * gamma(2 +
                  xi) * psigamma(2 + xi, 0) + xi * (1 + xi) * (1 + 2 * xi) * gamma(1 + 2 * xi) * psigamma(2 + 2 * xi, 0))/(sigma^2 *
                  xi^4)
            })
            k233_l <- sapply(xit, function(xi) {
                (1 + 7 * xi + 2 * digamma(1) * xi + 4 * xi^2 + 6 * digamma(1) * xi^2 + digamma(1)^2 * xi^2 + (pi^2 * xi^2)/6 - 3 *
                  xi^3 * (9 + 4 * xi) * gamma(3 * xi) + 3 * gamma(1 + 2 * xi) + 17 * xi * gamma(1 + 2 * xi) + 22 * xi^2 * gamma(1 +
                  2 * xi) + 8 * xi^3 * gamma(1 + 2 * xi) - (1 + 6 * xi) * gamma(1 + 3 * xi) + (1 + 3 * xi + 2 * xi^2) * 2 * xi * gamma(1 +
                  2 * xi) * psigamma(2 + 2 * xi, 0) - gamma(1 + xi) * (3 + 16 * xi + 13 * xi^2 + 4 * xi^3 + 2 * xi * (2 + 5 * xi +
                  3 * xi^2) * psigamma(2 + xi, 0) + xi^2 * (1 + xi) * psigamma(2 + xi, 0)^2 + xi^2 * (1 + xi) * psigamma(2 + xi, 1)))/(sigma *
                  xi^5)
            })
            k133 <- approx(x = c(0, xit), y = c(k133_0, k133_l), xout = xi)$y
            k223 <- approx(x = c(0, xit), y = c(k223_0, k223_l), xout = xi)$y
            k233 <- approx(x = c(0, xit), y = c(k233_0, k233_l), xout = xi)$y
            k23.3 <- approx(x = c(0, xit), y = c(k23.3_0, k23.3_l), xout = xi)$y
        }
    }

    # Numerical tolerance 1e-2 If function is not too numerically unstable
    if (abs(xi) > 0.01) {
        k333 <- (-3 * gamma(3 + 2 * xi) * (1 + 6 * xi + 4 * xi^2 + xi * (1 + 2 * xi) * psigamma(2 + 2 * xi, 0)) + 6 * (1 + 2 * xi) *
            gamma(1 + xi) * (1 + 8 * xi + 7 * xi^2 + 4 * xi^3 + 2 * xi * (1 + 4 * xi + 3 * xi^2) * psigamma(2 + xi, 0) + xi^2 * (1 +
            xi) * psigamma(2 + xi, 0)^2 + xi^2 * (1 + xi) * psigamma(2 + xi, 1)) + (1 + 2 * xi) * (-2 - 6 * (4 + digamma(1)) * xi -
            (30 + 42 * digamma(1) + 6 * digamma(1)^2 + pi^2) * xi^2 + 2 * (1 + xi)^2 * (1 + 4 * xi) * gamma(1 + 3 * xi) + xi^3 * (-8 -
            18 * digamma(1)^2 - 2 * digamma(1)^3 - 3 * pi^2 - digamma(1) * (24 + pi^2) + 4 * zeta3)))/(xi^6 * (2 + 4 * xi))
        k33.3 <- 2 * ((xi + 1)^2 * psigamma(2 * xi + 1) * gamma(2 * xi + 1)/xi^2 - ((xi + 1)/xi + psigamma(xi + 1)) * psigamma(xi +
            2) * gamma(xi + 2)/xi + ((xi + 1)/xi^2 - 1/xi - psigamma(xi + 1, 1)) * gamma(xi + 2)/xi - (xi + 1)^2 * gamma(2 * xi +
            1)/xi^3 + (xi + 1) * gamma(2 * xi + 1)/xi^2 + ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2)/xi^2 - (digamma(1) + 1/xi +
            1)/xi^2)/xi^2 - 1/3 * (pi^2 + 6 * (digamma(1) + 1/xi + 1)^2 + 6 * (xi + 1)^2 * gamma(2 * xi + 1)/xi^2 - 12 * ((xi + 1)/xi +
            psigamma(xi + 1)) * gamma(xi + 2)/xi)/xi^3
    } else {
        k333_0 <- -20.8076715595589
        k33.3_0 <- -5.45021409786022
        if (isTRUE(all.equal(xi, 0))) {
            k333 <- k333_0
            k33.3 <- k33.3_0
        } else {
            xit <- sign(xi) * 0.01
            k333_l <- sapply(xit, function(xi) {
                (-3 * gamma(3 + 2 * xi) * (1 + 6 * xi + 4 * xi^2 + xi * (1 + 2 * xi) * psigamma(2 + 2 * xi, 0)) + 6 * (1 + 2 * xi) *
                  gamma(1 + xi) * (1 + 8 * xi + 7 * xi^2 + 4 * xi^3 + 2 * xi * (1 + 4 * xi + 3 * xi^2) * psigamma(2 + xi, 0) + xi^2 *
                  (1 + xi) * psigamma(2 + xi, 0)^2 + xi^2 * (1 + xi) * psigamma(2 + xi, 1)) + (1 + 2 * xi) * (-2 - 6 * (4 + digamma(1)) *
                  xi - (30 + 42 * digamma(1) + 6 * digamma(1)^2 + pi^2) * xi^2 + 2 * (1 + xi)^2 * (1 + 4 * xi) * gamma(1 + 3 * xi) +
                  xi^3 * (-8 - 18 * digamma(1)^2 - 2 * digamma(1)^3 - 3 * pi^2 - digamma(1) * (24 + pi^2) + 4 * zeta3)))/(xi^6 * (2 +
                  4 * xi))
            })
            k33.3_l <- sapply(xit, function(xi) {
                2 * ((xi + 1)^2 * psigamma(2 * xi + 1) * gamma(2 * xi + 1)/xi^2 - ((xi + 1)/xi + psigamma(xi + 1)) * psigamma(xi +
                  2) * gamma(xi + 2)/xi + ((xi + 1)/xi^2 - 1/xi - psigamma(xi + 1, 1)) * gamma(xi + 2)/xi - (xi + 1)^2 * gamma(2 *
                  xi + 1)/xi^3 + (xi + 1) * gamma(2 * xi + 1)/xi^2 + ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2)/xi^2 - (digamma(1) +
                  1/xi + 1)/xi^2)/xi^2 - 1/3 * (pi^2 + 6 * (digamma(1) + 1/xi + 1)^2 + 6 * (xi + 1)^2 * gamma(2 * xi + 1)/xi^2 - 12 *
                  ((xi + 1)/xi + psigamma(xi + 1)) * gamma(xi + 2)/xi)/xi^3
            })
            k333 <- approx(x = c(0, xit), y = c(k333_0, k333_l), xout = xi)$y
            k33.3 <- approx(x = c(0, xit), y = c(k33.3_0, k33.3_l), xout = xi)$y
        }
    }
    # Derivatives of information matrix
    A1 <- 0.5 * cbind(c(k111, k112, k113), c(k112, k122, k123), c(k113, k123, k133))
    A2 <- -cbind(c(k11.2, k12.2, k13.2), c(k12.2, k22.2, k23.2), c(k13.2, k23.2, k33.2)) + 0.5 * cbind(c(k112, k122, k123), c(k122,
        k222, k223), c(k123, k223, k233))
    A3 <- -cbind(c(k11.3, k12.3, k13.3), c(k12.3, k22.3, k23.3), c(k13.3, k23.3, k33.3)) + 0.5 * cbind(c(k113, k123, k133), c(k123,
        k223, k233), c(k133, k233, k333))
    # Information matrix
    infomat <- gev.infomat(par = c(0, sigma, xi), dat = 1, method = "exp", nobs = 1)
    infoinv <- solve(infomat)

    return(infoinv %*% cbind(A1, A2, A3) %*% c(infoinv)/n)
}


#' Cox-Snell first order bias expression for the generalized Pareto distribution
#'
#' Bias vector for the GP distribution based on an \code{n} sample.
#' @inheritParams gpd
#' @references Coles, S. (2001). \emph{An Introduction to Statistical Modeling of Extreme Values}, Springer, 209 p.
#'@references Cox, D. R. and E. J. Snell (1968). A general definition of residuals, \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \strong{30}, 248--275.
#' @references Cordeiro, G. M. and R. Klein (1994). Bias correction in ARMA models, \emph{Statistics and Probability Letters}, \strong{19}(3), 169--176.
#' @references Giles, D. E., Feng, H. and R. T. Godwin (2016).  Bias-corrected maximum likelihood estimation of the  parameters of the generalized Pareto distribution, \emph{Communications in Statistics - Theory and Methods}, \strong{45}(8), 2465--2483.
#' @export
#' @keywords internal
#' @seealso \code{\link{gpd}}, \code{\link{gpd.bcor}}
gpd.bias <- function(par, n) {
    # scale, shape
    if (length(par) != 2) {
        stop("Invalid input for correction")
    }
    if (length(n) > 1) {
        stop("Invalid argument for sample size")
    }
    if (3 * par[2] < -1) {
        stop("Invalid bias correction for GPD; need shape > -1/3")
    }
    c(par[1] * (3 + 5 * par[2] + 4 * par[2]^2), -(1 + par[2]) * (3 + par[2]))/(n + 3 * n * par[2])
}

#'  Firth's modified score equation for the generalized Pareto distribution
#'
#' @inheritParams gpd
#' @references Firth, D. (1993). Bias reduction of maximum likelihood estimates, \emph{Biometrika}, \strong{80}(1), 27--38.
#' @seealso \code{\link{gpd}}, \code{\link{gpd.bcor}}
#' @export
#' @keywords internal
gpd.Fscore <- function(par, dat, method = c("obs", "exp")) {
    if (missing(method) || method != "exp") {
        method <- "obs"
    }
    gpd.score(par, dat) - gpd.infomat(par, dat, method) %*% gpd.bias(par, length(dat))
}

#'  Firth's modified score equation for the generalized extreme value distribution
#'
#' @seealso \code{\link{gev}}
#' @inheritParams gev
#' @param method string indicating whether to use the expected ('exp') or the observed ('obs' - the default) information matrix.
#' @references Firth, D. (1993). Bias reduction of maximum likelihood estimates, \emph{Biometrika}, \strong{80}(1), 27--38.
#' @export
#' @keywords internal
gev.Fscore <- function(par, dat, method = "obs") {
    if (missing(method) || method != "exp") {
        method <- "obs"
    }
    if (par[3] < 0) {
        mdat <- max(dat)
    } else if (par[3] >= 0) {
        mdat <- min(dat)
    }
    if (((mdat - par[1]) * par[3] + par[2]) < 0) {
        # warning('Error in `gev.Fscore`: data outside of range specified by parameter, yielding a zero likelihood')
        return(rep(1e+08, 3))
    } else {
        gev.score(par, dat) - gev.infomat(par, dat, method) %*% gev.bias(par, length(dat))
    }
}

#' Bias correction for GP distribution
#'
#' Bias corrected estimates for the generalized Pareto distribution using
#' Firth's modified score function or implicit bias subtraction.
#'
#'
#' Method \code{subtract} solves
#' \deqn{\tilde{\boldsymbol{\theta}} = \hat{\boldsymbol{\theta}} + b(\tilde{\boldsymbol{\theta}}}
#' for \eqn{\tilde{\boldsymbol{\theta}}}, using the first order term in the bias expansion as given by \code{\link{gpd.bias}}.
#'
#' The alternative is to use Firth's modified score and find the root of
#' \deqn{U(\tilde{\boldsymbol{\theta}})-i(\tilde{\boldsymbol{\theta}})b(\tilde{\boldsymbol{\theta}}),}
#' where \eqn{U} is the score vector, \eqn{b} is the first order bias and \eqn{i} is either the observed or Fisher information.
#'
#' The routine uses the MLE as starting value and proceeds
#' to find the solution using a root finding algorithm.
#' Since the bias-correction is not valid for \eqn{\xi < -1/3}, any solution that is unbounded
#' will return a vector of \code{NA} as the bias correction does not exist then.
#'
#' @importFrom alabama constrOptim.nl
#' @importFrom nleqslv nleqslv
#' @param par parameter vector (\code{scale}, \code{shape})
#' @param dat sample of observations
#' @param corr string indicating which correction to employ either \code{subtract} or \code{firth}
#' @param method string indicating whether to use the expected  (\code{'exp'}) or the observed (\code{'obs'} --- the default) information matrix. Used only if \code{corr='firth'}
#' @return vector of bias-corrected parameters
#' @export
#' @examples
#' set.seed(1)
#' dat <- evd::rgpd(n=40, scale=1, shape=-0.2)
#' par <- gp.fit(dat, threshold=0, show=FALSE)$estimate
#' gpd.bcor(par,dat, 'subtract')
#' gpd.bcor(par,dat, 'firth') #observed information
#' gpd.bcor(par,dat, 'firth','exp')
gpd.bcor <- function(par, dat, corr = c("subtract", "firth"), method = c("obs", "exp")) {
    corr <- match.arg(corr)
    method <- match.arg(method)
    # Basic bias correction - substract bias at MLE parbc=par-bias(par) Other bias correction - find bias corrected that solves
    # implicit eqn parbc=par-bias(parbc)
    if (length(par) != 2) {
        stop("Invalid `par` argument.")
    }
    # Basic bias correction - substract bias at MLE parbc=par-bias(par) bcor1 <- function(par, dat){ par-gpd.bias(par,length(dat))}
    # Other bias correction - find bias corrected that solves implicit eqn parbc=par-bias(parbc)
    bcor <- function(par, dat) {
        mdat <- max(dat)
        st.opt <- par
        st.opt[2] <- max(-0.25, st.opt[2])
        # Constrained optimization on L2 norm squared to find decent starting value
        bcor.st <- try(alabama::constrOptim.nl(par = st.opt, fn = function(parbc, para, dat) {
            sum((parbc - para + gpd.bias(parbc, length(dat)))^2)
        }, hin = function(x, ...) {
            c(x[1], x[2] + 1/3, -x[2] + 1, x[2] * mdat + x[1])
        }, dat = dat, para = par, control.optim = list(maxit = 1000), control.outer = list(trace = FALSE)), silent = TRUE)
        # Check if convergence, use the starting value (otherwise use par)
        if (isTRUE(try(bcor.st$value < 0.01))) {
            st <- bcor.st$par
        } else {
            st <- par
        }
        if (st[2] < -1/3 || st[2] > 1) {
            # If value is the MLE, make sure that it is a sensible starting value i.e. does not violate the constraints for the moments.
            warning("Error in gpd.bcor`: invalid starting value for the shape (less than -1/3). Aborting")
            return(rep(NA, 2))
        }
        # Root finding
        bcor.rf <- try(nleqslv::nleqslv(x = st, fn = function(parbc, par, dat) {
            parbc - par + gpd.bias(parbc, length(dat))
        }, par = par, dat = dat), silent = TRUE)
        if (!is.character(bcor.rf)) {
            if (bcor.rf$termcd == 1 || (bcor.rf$termcd == 2 && isTRUE(all.equal(bcor.rf$fvec, rep(0, 2), tolerance = 1e-06)))) {
                return(bcor.rf$x)
            }
        }
        return(rep(NA, 2))
    }
    bcorF <- function(par, dat, method = c("obs", "exp")) {
        method <- match.arg(method)
        if (par[2] < 0) {
            mdat <- max(dat)
        } else {
            mdat <- min(dat)
        }
        st.opt <- par
        st.opt[2] <- max(-0.25, st.opt[2])
        bcor.st <- try(alabama::constrOptim.nl(par = st.opt, fn = function(par, dat, met) {
            sum(gpd.Fscore(par = par, dat = dat, method = met)^2)
        }, hin = function(x, ...) {
            c(x[1], x[2] + 1/3, -x[2] + 1, x[2] * mdat + x[1])
        }, dat = dat, met = method, control.optim = list(maxit = 1000), control.outer = list(trace = FALSE)), silent = TRUE)
        if (isTRUE(try(bcor.st$value < 0.001))) {
            st <- bcor.st$par
        } else {
            st <- st.opt
        }
        # Starting values for score-based methods, depending on feasibility par.firth <-
        # try(suppressWarnings(rootSolve::multiroot(gpd.Fscore, start = st, dat = dat, method = method, positive = FALSE, maxiter =
        # 1000)), silent = TRUE)
        gpd.Fscoren <- function(par, met, dat) {
            gpd.Fscore(par = par, method = met, dat = dat)
        }
        par.firth <- try(suppressWarnings(nleqslv::nleqslv(fn = gpd.Fscoren, x = st, dat = dat, met = method, control = list(maxit = 1000))),
            silent = TRUE)
        if (!is.character(par.firth)) {
            if (par.firth$termcd == 1 || (par.firth$termcd == 2 && isTRUE(all.equal(par.firth$fvec, rep(0, 2), tolerance = 1e-06)))) {
                return(par.firth$x)
            }
        }
        return(rep(NA, 2))
    }
    # Return values
    if (corr == "subtract") {
        return(bcor(par = par, dat = dat))
    }
    if (corr == "firth") {
        return(bcorF(par = par, dat = dat, method = method))
    }
}




#' Bias correction for GEV distribution
#'
#' Bias corrected estimates for the generalized extreme value distribution using
#' Firth's modified score function or implicit bias subtraction.
#'
#' Method \code{subtract}solves
#' \deqn{\tilde{\boldsymbol{\theta}} = \hat{\boldsymbol{\theta}} + b(\tilde{\boldsymbol{\theta}}}
#' for \eqn{\tilde{\boldsymbol{\theta}}}, using the first order term in the bias expansion as given by \code{\link{gev.bias}}.
#'
#' The alternative is to use Firth's modified score and find the root of
#' \deqn{U(\tilde{\boldsymbol{\theta}})-i(\tilde{\boldsymbol{\theta}})b(\tilde{\boldsymbol{\theta}}),}
#' where \eqn{U} is the score vector, \eqn{b} is the first order bias and \eqn{i} is either the observed or Fisher information.
#'
#' The routine uses the MLE (bias-corrected) as starting values and proceeds
#' to find the solution using a root finding algorithm.
#' Since the bias-correction is not valid for \eqn{\xi < -1/3}, any solution that is unbounded
#' will return a vector of \code{NA} as the solution does not exist then.
#'
#' @param par parameter vector (\code{scale}, \code{shape})
#' @param dat sample of observations
#' @param corr string indicating which correction to employ either \code{subtract} or \code{firth}
#' @param method string indicating whether to use the expected  (\code{'exp'}) or the observed (\code{'obs'} --- the default) information matrix. Used only if \code{corr='firth'}
#' @return vector of bias-corrected parameters
#' @export
#' @examples
#' set.seed(1)
#' dat <- evd::rgev(n=40, loc = 1, scale=1, shape=-0.2)
#' par <- evd::fgev(dat)$estimate
#' gev.bcor(par,dat, 'subtract')
#' gev.bcor(par,dat, 'firth') #observed information
#' gev.bcor(par,dat, 'firth','exp')
gev.bcor <- function(par, dat, corr = c("subtract", "firth"), method = c("obs", "exp")) {
    corr <- match.arg(corr)
    # Basic bias correction - substract bias at MLE parbc=par-bias(par) bcor1 <- function(par, dat){ par-gpd.bias(par,length(dat))}
    # Other bias correction - find bias corrected that solves implicit eqn parbc=par-bias(parbc)
    bcor <- function(par, dat) {
        maxdat <- max(dat)
        mindat <- min(dat)
        st.opt <- par
        st.opt[3] <- max(-0.25, st.opt[3])
        # Constrained optimization on L2 norm squared to find decent starting value
        bcor.st <- try(alabama::constrOptim.nl(par = st.opt, fn = function(parbc, para, dat) {
            sum((parbc - para + gev.bias(parbc, length(dat)))^2)
        }, hin = function(x, ...) {
            c(x[2], x[3] + 1/3, -x[3] + 1, x[3] * (maxdat - x[1]) + x[2], x[3] * (mindat - x[1]) + x[2])
        }, dat = dat, para = par, control.optim = list(maxit = 1000), control.outer = list(trace = FALSE)), silent = TRUE)
        # Check if convergence, use the starting value (otherwise use par)
        if (isTRUE(try(bcor.st$value < 0.01))) {
            st <- bcor.st$par
        } else {
            st <- par
        }
        if (st[3] < -1/3 || st[3] > 1) {
            # If value is the MLE, make sure that it is a sensible starting value i.e. does not violate the constraints for the moments.
            warning("Error in gev.bcor`: invalid starting value for the shape (less than -1/3). Aborting")
            return(rep(NA, 3))
        }
        # Root finding

        bcor.rf <- try(nleqslv::nleqslv(x = st, fn = function(parbc, par, dat) {
            parbc - par + gev.bias(parbc, length(dat))
        }, par = par, dat = dat, control = list(maxit = 1000, xtol = 1e-10)), silent = TRUE)
        if (!is.character(bcor.rf)) {
            if (bcor.rf$termcd == 1 || (bcor.rf$termcd %in% c(2, 3) && isTRUE(all.equal(bcor.rf$fvec, rep(0, 3), tolerance = 1e-06)))) {
                return(bcor.rf$x)
            }
            bcor.rf <- try(nleqslv::nleqslv(x = st, fn = function(parbc, par, dat) {
                parbc - par + gev.bias(parbc, length(dat))
            }, global = "none", par = par, dat = dat, control = list(maxit = 1000, xtol = 1e-10)), silent = TRUE)
            if (!is.character(bcor.rf)) {
                if (bcor.rf$termcd == 1 || (bcor.rf$termcd %in% c(2, 3) && isTRUE(all.equal(bcor.rf$fvec, rep(0, 3), tolerance = 1e-06)))) {
                  return(bcor.rf$x)
                } else if (abs(bcor.rf$x[3]) < 0.025 && bcor.rf$termcd == 2 && isTRUE(all.equal(bcor.rf$fvec, rep(0, 3), tolerance = 0.01))) {
                  warning(paste0("Approximate solution for implicit bias correction - the shape is close to zero"))
                  return(bcor.rf$x)
                }
            }
        }
        return(rep(NA, 3))
    }

    bcorF <- function(par, dat, method = c("obs", "exp")) {
        method <- match.arg(method)
        maxdat <- max(dat)
        mindat <- min(dat)
        st.opt <- par
        st.opt[3] <- max(-0.25, st.opt[3])
        bcor.st <- try(alabama::constrOptim.nl(par = st.opt, fn = function(par, dat, met) {
            sum(gev.Fscore(par = par, dat = dat, method = met)^2)
        }, hin = function(x, ...) {
            c(x[2], x[3] + 0.3, -x[3] + 1, x[3] * (maxdat - x[1]) + x[2], x[3] * (mindat - x[1]) + x[2])
        }, dat = dat, met = method, control.optim = list(maxit = 1000), control.outer = list(trace = FALSE)), silent = TRUE)
        if (isTRUE(try(bcor.st$value < 0.001))) {
            st <- bcor.st$par
        } else {
            st <- st.opt
        }
        # Starting values for score-based methods, depending on feasibility par.firth <-
        # try(suppressWarnings(rootSolve::multiroot(gev.Fscore, start = st, dat = dat, method = method, positive = FALSE, maxiter =
        # 1000)), silent = TRUE)
        gev.Fscoren <- function(par, met, dat) {
            gev.Fscore(par = par, method = met, dat = dat)
        }
        par.firth <- try(suppressWarnings(nleqslv::nleqslv(fn = gev.Fscoren, x = st, dat = dat, met = method, control = list(maxit = 1000,
            xtol = 1e-10))), silent = TRUE)
        if (!is.character(par.firth)) {
            # Try finding the root with default options
            if (par.firth$termcd == 1 || (par.firth$termcd %in% c(2, 3) && isTRUE(all.equal(par.firth$fvec, rep(0, 3), tolerance = 1e-06)))) {
                return(par.firth$x)
            }
            # Sometimes, method fails for values of xi close to zero - try with a full Broyden or Newton search
            par.firth <- try(suppressWarnings(nleqslv::nleqslv(fn = gev.Fscoren, x = st, dat = dat, met = method, control = list(maxit = 1000,
                xtol = 1e-10), global = "none")), silent = TRUE)
            if (!is.character(par.firth)) {
                if (par.firth$termcd == 1 || (par.firth$termcd %in% c(2, 3) && isTRUE(all.equal(par.firth$fvec, rep(0, 3), tolerance = 1e-06)))) {
                  return(par.firth$x)
                } else if (abs(par.firth$x[3]) < 0.025 && par.firth$termcd == 2 && isTRUE(all.equal(par.firth$fvec, rep(0, 3), tolerance = 0.05))) {
                  warning(paste0("Approximate solution for Firth`s score equation with method `", method, "` - the shape is close to zero"))
                  return(par.firth$x)
                }
            }
        }
        return(rep(NA, 3))
    }
    # Return values
    if (corr == "subtract") {
        return(bcor(par = par, dat = dat))
    }
    if (corr == "firth") {
        return(bcorF(par = par, dat = dat, method = method))
    }
}

#' Posterior predictive distribution and density for the GEV distribution
#'
#' This function calculates the posterior predictive density at points x
#' based on a matrix of posterior density parameters
#'
#' @param x \code{n} vector of points
#' @param posterior \code{n} by \code{3} matrix of posterior samples
#' @param Nyr number of years to extrapolate
#' @param type string indicating whether to return the posterior \code{density} or the \code{quantile}.
#' @param vector of values for the posterior predictive density or quantile at x
#' @export
#' @keywords internal
.gev.postpred <- function(x, posterior, Nyr = 100, type = c("density", "quantile")) {
    rowMeans(cbind(apply(rbind(posterior), 1, function(par) {
        switch(type, density = evd::dgev(x = x, loc = par[1] - par[2] * (1 - Nyr^par[3])/par[3], scale = par[2] * Nyr^par[3], shape = par[3]),
            quantile = evd::qgev(p = x, loc = par[1] - par[2] * (1 - Nyr^par[3])/par[3], scale = par[2] * Nyr^par[3], shape = par[3]))
    })))
}


#' N-year return levels, median and mean estimate
#'
#' @param par vector of location, scale and shape parameters for the GEV distribution
#' @param nobs integer number of observation on which the fit is based
#' @param N integer number of observations for return level. See \strong{Details}
#' @param type string indicating the statistic to be calculated (can be abbreviated).
#' @param p probability indicating the return level, corresponding to the quantile at 1-1/p
#'
#' @details If there are \eqn{n_y} observations per year, the \code{L}-year return level is obtained by taking
#' \code{N} equal to \eqn{n_yL}.
#' @export
#' @return a list with components
#' \itemize{
#' \item{est} point estimate
#' \item{var} variance estimate based on delta-method
#' \item{type} statistic
#' }
gev.Nyr <- function(par, nobs, N, type = c("retlev", "median", "mean"), p = 1/N) {
    # Create copy of parameters
    mu <- par[1]
    sigma <- par[2]
    xi <- par[3]
    type <- match.arg(type)
    # Check whether arguments are well defined
    if (type == "retlev") {
        stopifnot(p >= 0, p < 1)
        yp <- -log(1 - p)
    } else {
        stopifnot(N >= 1)
    }

    # Euler-Masc. constant :
    emcst <- -psigamma(1)
    # Return levels, N-year median and mean for GEV
    estimate <- switch(type, retlev = ifelse(xi == 0, mu - sigma * log(yp), mu - sigma/xi * (1 - yp^(-xi))), median = ifelse(xi ==
        0, mu + sigma * (log(N) - log(log(2))), mu - sigma/xi * (1 - (N/log(2))^xi)), mean = ifelse(xi == 0, mu + sigma * (log(N) +
        emcst), mu - sigma/xi * (1 - N^xi * gamma(1 - xi))))
    if (type == "retlev") {
        if (p > 0) {
            if (xi == 0) {
                grad_retlev <- c(1, -log(yp), 0.5 * sigma * log(yp)^2)
            } else {
                grad_retlev <- c(1, -(1 - yp^(-xi))/xi, sigma * (1 - yp^(-xi))/xi^2 - sigma/xi * yp^(-xi) * log(yp))
            }
        }
        if (p == 0) {
            if (xi < 0) {
                grad_retlev <- c(1, -1/xi, sigma/xi^2)
            } else {
                stop("Invalid argument; maximum likelihood estimate of the uptyper endpoint is Inf")
            }
        }
        # Variance estimate based on delta-method
        var_retlev <- t(grad_retlev) %*% solve(gev.infomat(par = par, dat = 1, method = "exp", nobs = nobs)) %*% grad_retlev
    } else if (type == "median") {
        # Gradient of N-years maxima median
        if (xi == 0) {
            grad_Nmed <- c(1, log(N/log(2)), 0.5 * sigma * log(N/log(2))^2)
        } else {
            grad_Nmed <- c(1, ((N/log(2))^xi - 1)/xi, sigma * (N/log(2))^xi * log(N/log(2))/xi - sigma * ((N/log(2))^xi - 1)/xi^2)
        }
        # Delta-method covariance matrix
        var_Nmed <- t(grad_Nmed) %*% solve(gev.infomat(par = par, dat = 1, method = "exp", nobs = nobs)) %*% grad_Nmed
    } else if (type == "mean") {
        if (xi == 0) {
            grad_Nmean <- c(1, log(N) + emcst, 0.5 * sigma * (emcst^2 + pi^2/6 + 2 * emcst * log(N) + log(N)^2))
        } else {
            grad_Nmean <- c(1, (N^xi * gamma(-xi + 1) - 1)/xi, (N^xi * log(N) * gamma(-xi + 1) - N^xi * digamma(-xi + 1) * gamma(-xi +
                1)) * sigma/xi - (N^xi * gamma(-xi + 1) - 1) * sigma/xi^2)
        }
        var_Nmean <- t(grad_Nmean) %*% solve(gev.infomat(par = par, dat = 1, method = "exp", nobs = nobs)) %*% grad_Nmean
    }
    var_est <- switch(type, retlev = var_retlev, median = var_Nmed, mean = var_Nmean)
    return(list(est = estimate, var = var_est[1, 1], type = type))
}


#' Asymptotic bias of block maxima for fixed sample sizes
#'
#' @param shape shape parameter
#' @param rho second-order parameter, non-positive
#' @references Dombry, C. and A. Ferreira (2017). Maximum likelihood estimators based on the block maxima method. \code{https://arxiv.org/abs/1705.00465}
#' @export
#' @return a vector of length three containing the bias for location, scale and shape (in this order)
gev.abias <- function(shape, rho) {
    stopifnot(rho <= 0, shape > -0.5)
    if (shape != 0 && rho < 0) {
        bmu <- (1 + shape)/(shape * rho * (shape + rho)) * (-(shape + rho) * gamma(1 + shape) + (1 + shape) * rho * gamma(1 + 2 *
            shape) + shape * (1 - rho) * gamma(1 + shape - rho))
        bsigma <- (-shape - rho + (1 + shape) * (shape + 2 * rho) * gamma(1 + shape) - (1 + shape)^2 * rho * gamma(1 + 2 * shape) +
            shape * gamma(2 - rho) - shape * (1 + shape) * (1 - rho) * gamma(1 + shape - rho))/(shape^2 * rho * (shape + rho))
        bshape <- ((shape + rho) * (1 + shape + psigamma(1) * shape) - (shape + shape^2 * (1 + rho) + 2 * rho * (1 + shape)) * gamma(1 +
            shape) + (1 + shape)^2 * rho * gamma(1 + 2 * shape) + shape^2 * gamma(1 - rho) - shape * (1 + shape) * gamma(2 - rho) +
            shape * (1 + shape) * (1 - rho) * gamma(1 + shape - rho) - shape * rho * psigamma(2 + shape) * gamma(2 + shape) - shape^2 *
            psigamma(2 - rho) * gamma(2 - rho))/(shape^3 * rho * (shape + rho))
    } else if (rho == 0 && shape != 0) {
        bmu <- (1 + shape)/shape^2 * ((1 + shape) * gamma(1 + 2 * shape) - gamma(2 + shape) - shape * psigamma(1 + shape) * gamma(1 +
            shape))
        bsigma <- (shape^2 * gamma(1 + shape) * psigamma(1 + shape) - (shape + 1)^2 * gamma(1 + 2 * shape) + shape^2 * gamma(1 + shape) +
            shape * gamma(1 + shape) * psigamma(1 + shape) - (psigamma(1) + 1) * shape + 3 * shape * gamma(1 + shape) + 2 * gamma(1 +
            shape) - 1)/shape^3
        bshape <- ((1 + shape + shape * psigamma(1))^2 + shape^2 * pi^2/6 + (1 + shape)^2 * gamma(1 + 2 * shape) - 2 * (1 + shape) *
            ((1 + shape) * gamma(1 + shape) + shape * psigamma(1 + shape) * gamma(1 + shape)))/shape^4
    } else if (rho < 0 && shape == 0) {
        bmu <- (-1 + rho + psigamma(1) * rho + (1 - rho) * gamma(1 - rho))/rho^2
        bsigma <- (6 - 6 * rho - 6 * psigamma(1)^2 * rho - pi^2 * rho + 6 * psigamma(1) * (1 - 2 * rho) - 6 * (1 - rho) * gamma(1 -
            rho) * (1 + psigamma(1 - rho)))/(6 * rho^2)
        bshape <- -1/(12 * rho^3) * (-12 * gamma(2 - rho) * psigamma(2 - rho) + 6 * gamma(1 - rho) * (2 + 2 * (-1 + rho)^2 * psigamma(1 -
            rho) + (-1 + rho) * rho * psigamma(1 - rho)^2 + (-1 + rho) * rho * psigamma(1 - rho, deriv = 1)) + rho * (psigamma(1)^2 *
            (6 - 18 * rho) + pi^2 * (1 - 3 * rho) + 6 * (-psigamma(1)^3) * rho + 3 * (-psigamma(1)) * (-4 + (4 + pi^2) * rho) + 6 *
            rho * (-2 - 2 * psigamma(1, deriv = 2) + psigamma(2, 2))))
    } else if (rho == 0 && shape == 0) {
        # Last row of information matrix with (0, 1, shape)
        bmu <- 0.41184033042644
        bsigma <- 0.332484907160274
        bshape <- 2.42360605517703
    }
    c(solve(gev.infomat(c(0, 1, shape), dat = 1, method = "exp", nobs = 1)) %*% c(bmu, bsigma, bshape))
}

#' Asymptotic bias of threshold exceedances for k order statistics
#'
#' The formula given in de Haan and Ferreira, 2007 (Springer). Note that the latter differs from that found in Drees, Ferreira and de Haan.
#' @references Dombry, C. and A. Ferreira (2017). Maximum likelihood estimators based on the block maxima method. \code{https://arxiv.org/abs/1705.00465}
#' @param shape shape parameter
#' @param rho second-order parameter, non-positive
#' @export
#' @return a vector of length containing the bias for scale and shape (in this order)
gpd.abias <- function(shape, rho) {
    stopifnot(rho <= 0, shape > -0.5)
    c(-rho/((1 - rho) * (1 + shape - rho)), (1 + shape)/((1 - rho) * (1 + shape - rho)))
}
