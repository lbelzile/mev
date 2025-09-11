#' Robust threshold selection method via OBRE
#'
#' The function calculates the cumulative sum of the non-zero weights for the generalized Pareto model fitted via optimal B-robust estimation (OBRE). Bootstrap samples from the generalized Pareto are generated and the sum of weights for the same numbers are calculated. This exploits the fact that the weights are monotonically increasing, so downweighting (if any) occurs in the tails with increasing weights (capped at unity).
#' @param xdat vector of observations
#' @param thresh vector of candidate thresholds
#' @param nsim number of bootstrap simulations
#' @return a list with elements
#' \itemize{
#' \item \code{coef} matrix of parameter estimates for the scale and shape parameters, estimated via OBRE
#' \item \code{thresh} vector of ordered thresholds
#' \item \code{pval} p-values from the procedure
#' \item \code{ndw} number of downweighted observations
#' }
#' @references Dupuis, D.J. (1998). Exceedances over High Thresholds: A Guide to Threshold Selection,
#' \emph{Extremes}, \bold{1}(3), 251--261.
# thselect.dobre <- function(
#   xdat,
#   thresh,
#   nsim = 99L
# ) {
#   # Preprocess inputs
#   nu <- length(thresh)
#   thresh <- sort(thresh)
#   nw <- integer(nu)
#   xdat <- sort(xdat[is.finite(xdat)])
#   B <- as.integer(nsim)
#   stopifnot(B >= 19, xdat[length(xdat)] > thresh[nu])
#   pars <- matrix(nrow = nu, ncol = 3L)
#   colnames(pars) <- c("threshold", "scale", "shape")
#   pars[, 1] <- thresh
#   pvals <- list()
#   for (i in seq_along(thresh)) {
#     fit_obre <- fit.gpd(
#       xdat = xdat,
#       threshold = thresh[i],
#       method = "obre",
#       tol = 1e-5
#     )
#     dw <- which(fit_obre$weights < 1)
#     nw[i] <- length(dw)
#     pars[i, -1] <- fit_obre$estimate[1]
#     if (nw[i] >= 1) {
#       upp <- isTRUE(dw[length(dw)] == fit_obre$nat)
#       stats <- matrix(nrow = B, ncol = nw[i])
#       # OBRE weights are decreasing, but
#       for (b in seq_len(B - 1)) {
#         boot_samp <- rgp(
#           n = fit_obre$nat,
#           scale = fit_obre$estimate[1],
#           shape = fit_obre$estimate[2]
#         )
#         boot_gpfit <- fit.gpd(
#           xdat = boot_samp,
#           threshold = 0,
#           method = "obre",
#           tol = 1e-5
#         )
#         # Downweighted observations are in the upper tail
#         if (upp) {
#           stats[b, ] <- cumsum(rev(boot_gpfit$weights)[seq_len(nw[i])])
#         } else {
#           # downweight in lower tail
#           stats[b, ] <- cumsum(boot_gpfit$weights[seq_len(nw[i])])
#         }
#         if (dw[length(dw)] == fit_obre$nat) {
#           stats[B, ] <- cumsum(rev(fit_obre$weights)[seq_len(nw[i])])
#         } else {
#           # downweight in lower tail
#           stats[B, ] <- cumsum(fit_obre$weights[seq_len(nw[i])])
#         }
#       }
#       stats[B, ] <- cumsum(rev(fit_obre$weights)[seq_len(nw[i])])
#       # Small weights are extreme, so is their sum
#       # TODO check whether we get more power from looking at something different than the largest
#       pvals[[i]] <- apply(stats, 2, rank)[B, ] / B
#     } else {
#       pvals[[i]] <- NA
#     }
#   }
#   list(
#     thresh = thresh,
#     coef = pars,
#     pval = pvals,
#     ndw = nw
#   )
# }

# Quick simulation study
#
# set.seed(1234)
# shape_seq <- c(-0.3, -0.15, 0, 0.15, 0.3)
# n_seq <- c(50, 100, 200)
# nrep <- 1000
# qlev_seq <- c(0, 0.25, 0.5)
# results <- list()
# indices <- matrix(nrow = prod(c(5, 3, nrep)), ncol = 3L)
# iter <- 0L
# for (j in seq_len(nrep)) {
#   set.seed(j)
#   for (s in seq_along(shape_seq)) {
#     for (i in seq_along(n_seq)) {
#       xdat <- rgp(n = n_seq[i], scale = 1, shape = shape_seq[s])
#       thresh <- quantile(xdat, qlev_seq)
#       res <- try(D.diag(xdat = xdat, thresh = thresh, B = 199))
#       iter <- iter + 1L
#       indices[iter, ] <- c(j, s, i)
#       results[[iter]] <- res
#     }
#   }
# }
