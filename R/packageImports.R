#' @useDynLib mev, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @importFrom Rcpp evalCpp
#' @importFrom grDevices "colorRampPalette"
#' @importFrom graphics "lines" "abline" "axis" "par" "matplot" "legend" "image" "mtext" "text" "rug"
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("mev", libpath)
}

# @importFrom stats "integrate" "quantile" "cov2cor" "dnorm" "median" "nlm" "optim" "pchisq" "qnorm" "rnorm" "sd" "var" "lm"

### LIST OF R FUNCTIONS PER FILE
### bias.R
# gev.bias
# gpd.bias
# gpd.Fscore
# gev.Fscore
# gpd.bcor
# gev.bcor
# .gev.postpred [not exported]
# gev.Nyr
# gev.abias
# gpd.abias
### bivpot_methods.R
# ibvpot
### egp.R
# egp2.fit
# pegp2.G [internal]
# degp2.G [internal]
# regp2.G [internal]
# qegp2.G [internal]
# pegp2 [internal]
# degp2 [internal]
# qegp2 [internal]
# regp2 [internal]
# egp2.pwm  [not exported]
# egp2.fit.pwm [not exported]
# egp2.fit.pwm.boot [not exported]
# egp2.nll [not exported]
# egp2.fit.ml [not exported]
# egp2.fit.ml.boot [not exported]
### gp.R
# NC.diag
# .score_algebraic [not exported]
# .exp_info_algebraic [not exported]
# .gpd_1D_fit [not exported]
# .gpd_obs_info  [not exported]
# .gpd_2D_fit [not exported]
# .GP_1D_fit_nlm  [not exported]
# .gpd_grimshaw  [not exported]
# .Zhang_Stephens_posterior [not exported]
# .Zhang_posterior  [not exported]
# gp.fit
### infomat_test.R
# infomat.test
# ext.index
### multivar.R
# chibar
# angextrapo
# lambdadep
### penultimate.R
# egp.ll [internal]
# egp.ll.opt [internal]
# egp.retlev [internal]
# egp.fit
# egp.fitrange
# smith.penult
# smith.penult.fn
### profile.R
# .auglag [not exported]
# gpd.mle
# gev.mle
# confint.extprof [method]
# plot.extprof [method]
# gev.pll
# gpd.pll
# plot.fr [method]
# spline.corr
# confint.fr [method]
# gev.tem
# gpd.tem
### rmev_wrapper.R
# rmev
# rmevspec
# .mvasym.check [not exported]
# .is.CNSD [not exported]
# rparp
### specdens.R
# .returnAng [not exported]
# angmeas
# .wecdf [not exported]
# .pickands.emp [not exported]
# emplik
# angmeasdir
### univdist.R
# gpd.ll
# gpd.ll.optim
# gpd.score
# gpd.infomat
# gpd.Vfun
# gpd.phi
# gpd.dphi
# gev.ll
# gev.ll.optim
# gev.score
# gev.infomat
# gev.Vfun
# gev.phi
# gev.dphi
# gpde.ll
# gpde.ll.optim
# gpde.score
# gpde.infomat
# gpde.Vfun
# gpde.phi
# gpde.dphi
# gpdr.ll
# gpdr.ll.optim
# gpdr.score
# gpdr.infomat
# gpdr.Vfun
# gpdr.phi
# gpdr.dphi
# gev.retlev
# gevNmean
# gevN.quant
# gpdN.mean
# gpdN.quant
# gevr.ll
# gevr.ll.optim
# gevr.score
# gevr.infomat
# gevr.Vfun
# gevr.phi
# gevr.dphi
# gpdN.ll
# gpdN.ll.optim
# gpdN.score
# gpdN.infomat
# gpdN.Vfun
# gpdN.phi
# gpdN.dphi
# gevN.ll
# gevN.ll.optim
# gevN.score
# gevN.infomat
# gevN.Vfun
# gevN.phi
# gevN.dphi
### Wdiag.R
# W.diag
# .NHPP.diag
# .Expl.diag
### spunif.R
# spunif
### mle.R
# fit.gpd
# fit.gev
# fit.pp
# plot.mev_gpd
# plot.mev_gev
# plot.mev_pp
# print
