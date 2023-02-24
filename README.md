# `mev`: Modelling Extreme values

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/mev)](https://cran.r-project.org/package=mev)
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/mev?color=brightgreen)](http://www.r-pkg.org/pkg/mev)

An **R** package for the analysis of univariate, multivariate and functional extreme values. The package includes routine functions for univariate analyses multiple threshold selection diagnostics, optimization, bias-correction and tangent exponential model approximations, non-parametric spectral measure estimation using empirical likelihood methods, etc. Multivariate functionalities revolve around simulation algorithms for multivariate models, empirical likelihood, empirical dependence measures. Likelihood functions for elliptical processes and user-provided methodologies.


To install from Github, use

```R
remotes::install_github("lbelzile/mev")
```

after installing `remotes`.

# Functionalities

The functionalities of the package are sorted below by topic.

## Univariate

The package focuses on likelihood based inference for parametric models.

Log likelihood, score and information matrices for the following univariate models:

- `gpd`: generalized Pareto distribution (alternative parametrizations `gpde`, `gpdN`, `gpdr`)
- `gev`: generalized extreme value distribution (alternative parametrizations `gevN`, `gevr`)
- `pp`: inhomogeneous Poisson process for extremes
- `rlarg`: asymptotic r-largest order statistics


Fitting procedures and higher order asymptotic inference for univariate extremes

- `fit.*` for maximum likelihood estimation
- `*.bcor` for bias correction via score vectors or by subtraction
- `*.pll`: profile likelihood for objects
- `*.tem` for tangent exponential model approximation to profile likelihood

Two additional penultimate models and utilities for approximations

- `egp`: extended generalized Pareto models of Papastathopoulos and Tawn (2013)
- `extgp`: extended generalized Pareto models of Naveau et al. for rainfall
- `smith.penult`: Smith (1987) penultimate approximations to parametric models


## Threshold selection

Multiple functions can be used for threshold selection for the peaks over threshold method

- `automrl`: automatic threshold selection for mean residual life plots
- `cvselect`: threshold selection via coefficient of variation
- `tstab.egp`: threshold stability plots for `egp` models
- `infomat.test`: information matrix test for time series
- `NC.diag`: Northrop and Coleman (2014) score tests 
- `tstab.gp`: threshold stability plot for generalized Pareto distribution
- `vmetric.diag`: metric-based threshold selection of Varty et al. 
- `W.diag`: Wadsworth (2016) sequential analysis threshold diagnostics

## Multivariate 

Some functionalities (incomplete) for multivariate models. There is currently no function to optimize
multivariate threshold models, but likelihoods are provided for logistic, Brown--Resnick, Huesler--Reiss and extremal Student models

- `ibvpot`: interpretation of bivariate models (extension of `evir` for all bivariate models from `evd`)
- `likmgp`, `clikmgp`: (censored) likelihood for multivariate generalized Pareto
- `expme`: exponent measure of parametric extreme value models

Two tests, one for max-stability and the other for asymptotic independence

- `maxstabtest`: test of max-stability
- `scoreindep`: score test of asymptotic independence for bivariate logistic model

### Nonparametric

Estimation of the angular distribution using empirical estimation or empirical likelihood, with or without smoothing

- `angmeas`: rank-based estimation of the angular measure
- `angmeasdir`: Dirichlet mixture smoothing of angular measure

## Simulation

Sampling algorithms for parametric models, multivariate and spatial extreme values, angular distribution and (generalized)  risk-Pareto processes using accept-reject or composition sampling (approximate).

- `rrlarg`: simulation of $r$-largest observations from point process of extremes
- `rdir`: simulation of Dirichlet vectors
- `mvrnorm`: simulation of multivariate normal vectors
- `rmev`: exact simulation of multivariate extreme value distributions
- `rmevspec`: random samples from angular distributions of multivariate extreme value models.
- `rparp`: simulation from R-Pareto processes
- `rparpcs`: simulation from Pareto processes (max) using composition sampling
- `rparpcshr`: simulation of generalized Huesler-Reiss Pareto vectors via composition sampling
- `rgparp`: simulation from generalized R-Pareto processes


### Extremal dependence measures

Measures of tail dependence $\theta$, $\eta$, $\chi$ and $\varphi$.

- `taildep`: estimators of coefficients of tail dependence $\eta$ and tail correlation $\chi$
- `extcoef`: estimators of the extremal coefficient
- `xasym`: estimators of the extremal asymmetry coefficient
- `angextrapo`: bivariate tail dependence $\eta$ across rays
- `lambdadep`: bivariate function of Wadsworth and Tawn (2013)
- `ext.index`: extremal index estimators based on interexceedance time and gap of exceedances
- `extremo`: pairwise extremogram as a function of distance for spatial data

## Datasets

Various datasets collected here and there, (exclusively?)  for univariate peaks over threshold analysis

- `abisko`: Abisko rainfall
- `eskrain`: Eskdalemuir observatory daily rainfall
- `geomagnetic`: magnitude of geomagnetic storms
- `maiquetia`: Maiquetia daily rainfall series
- `nidd`: river Nidd daily flow
- `venice`: Venice sea level data
- `w1500m`: women 1500m track records

## Spatial

Some functionalities for fitting spatial data

- `distg`: matrix of pairwise distance with geometric anisotropy
- Variogram models (unexported functions `powerexp.cor`, `power.vario`, `schlather.vario`)
- `Lambda2cov`: conver variogram to covariance of conditional random field

## Miscellaneous

Functions used internally that could be of more general use.

- `emplik`: empirical likelihood for vector mean
- `wecdf`: weighted empirical distribution function
- `spline.corr` and `tem.corr`: corrections for Fraser--Reid objects to remove singularities nead the mode

