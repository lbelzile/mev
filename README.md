# `mev`: Modelling Extreme values

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-last-release/mev)](https://cran.r-project.org/package=mev)
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/mev?color=brightgreen)](http://www.r-pkg.org/pkg/mev)

An **R** package for the analysis of univariate, multivariate and functional extreme values. The package includes routine functions for univariate analyses multiple threshold selection diagnostics, optimization, bias-correction and tangent exponential model approximations, non-parametric spectral measure estimation using empirical likelihood methods, etc. Multivariate functionalities revolve around simulation algorithms for multivariate models, empirical likelihood, empirical dependence measures. Likelihood functions for elliptical processes and user-provided methodologies.


To install from Github, use

```R
remotes::install_github("lbelzile/mev")
```

after installing `remotes`.

# Functionalities

The functionalities of the development version of the package (GIthub) are sorted below by topic.

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

Two additional models and utilities for penultimate approximations

- `egp`: extended generalized Pareto models of Papastathopoulos and Tawn (2013), and Gamet and Jonathan (2022)
- `extgp`: extended generalized Pareto models of Naveau et al. (2017)
- `smith.penult`: Smith (1987) penultimate approximations to parametric models

## Nonparametric estimators of shape and second order regular variation

The routine `fit.shape`, or alternatively one of subroutines for real or positive (*) shape parameters.

- `shape.hill`*: Hill's estimator
- `shape.osz`: Pickands extreme U-statistic of Oorschot, Segers and Zhou
- `shape.moment`: moment estimator of Dekkers and de Haan.
- `shape.pickands`: Pickands estimator (poor performance)
- `shape.vries`*: de Vries estimator
- `shape.rbm`*: Wager's random block maxima estimator 
- `shape.genquant`*: generalized quantile
- `shape.trimhill`*: trimmed Hill estimator
- `shape.lthill`*: left-truncated Hill estimator

Note that both of the trimmed and truncated Hill estimators are not vectorized.

## Threshold selection diagnostics

Functions for automatic selection of threshold with the peaks over threshold method

- `thselect.wseq`: Wadsworth (2016) sequential analysis threshold diagnostics
- `thselect.vmetric`: metric-based threshold selection of Varty et al. (2025+)
- `thselect.ncpgp`: Northrop and Coleman (2014) piecewise generalized Pareto
- `thselect.cv`: del Castillo and Padilla (2016) coefficient of variation method
- `thselect.sdinfo`: Suveges and Davison (2010) information matrix test
- `thselect.mrl`: Langousis et al. (2016) automatization of mean residual life diagnostics
- `thselect.pickands`: Pickands (1985) goodness-of-fit threshold selection diagnostic
- `thselect.alrs`: automatic L-moments ratio selection method of Silva Lomba and Fraga Alves (2020)
- `thselect.ksmd`: Mahalanobis distance-based selection method based on L-moments of Kiran and Srivinas (2021)

Some semiparametric methods

- `thselect.bab`: Bladt, Albrecher and Beirlant (2020) minimization of AMSE for Hill estimator via lower truncated Hill
- `thselect.expgqt` Exponential generalized quantile threshold selection of Beirlant, Vynckier and Teugels (1996)
- `thselect.gbw`: Kernel-based threshold selection of Goegebeur, Beirlant and de Wet (2008)
- `thselect.rbm`: Random block maximum estimator of Wager (2014), with empirical Bayes risk minimization


## Threshold stability plots

- `tstab.gpd`: threshold stability plots for generalized Pareto
- `tstab.egp`: threshold stability plots for extended generalized Pareto
- `tstab.cv`: coefficient of variation stability plot
- `tstab.mrl`: mean residual life plot
- `tstab.hill`: Hill plot

## Multivariate 

Some functionalities (incomplete) for multivariate models. There is currently no function to optimize
multivariate threshold models, but likelihoods are provided for logistic, Brown--Resnick, Huesler--Reiss and extremal Student models

- `ibvpot`: interpretation of bivariate models (extension of `evir` for all bivariate models from `evd`)
- `likmgp`, `clikmgp`: (censored) likelihood for multivariate generalized Pareto
- `expme`: exponent measure of parametric extreme value models

Two tests, one for max-stability and the other for asymptotic independence

- `test.maxstab`: graphical test of max-stability (P-P plot)
- `test.scoreindep`: score test of asymptotic independence for bivariate logistic model

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
- `frwind`: time series of wind speeds
- `geomagnetic`: magnitude of geomagnetic storms
- `leedspollution`: multivariate air pollutant from Leeds 
- `maiquetia`: Maiquetia daily rainfall series
- `nidd`: river Nidd daily flow
- `nutrients`: interview component from NHANES on nutrients
- `pandemics`: estimated records on number of death from pandemics
- `venice`: Venice sea level data
- `w1500m`: women 1500m track records

## Spatial

Some functionalities for fitting spatial data

- `distg`: matrix of pairwise distance with geometric anisotropy
- Variogram models (unexported functions `powerexp.cor`, `power.vario`, `schlather.vario`)
- `Lambda2cov`: convert variogram to covariance of conditional random field

## Miscellaneous

Functions used internally that could be of more general use.

- `emplik`: empirical likelihood for vector mean
- `wecdf`: weighted empirical distribution function
- `spline.corr` and `tem.corr`: corrections for Fraser--Reid objects to remove singularities nead the mode

