# mev v.1.14
==============
## New:

* bivariate coefficient of extremal asymmetry
* four new families of max-stable models (pairwise beta, pairwise exponential, weighted Dirichlet and weighted exponential) for rmev, following Ballani and Schlather (2011)
* maximum likelihood estimation routines (`fit.gpd`, `fit.gev`, etc.) now accept fixed parameters
* mean residual life plots with weighted least square fit
* coefficient of variation tests for threshold selection
* Varty et al. metric based threshold selection diagnostic
* anova methods for `mev_gpd` and `mev_gev` objects
* `rmev` and `rmevspec` now accept a distance matrix in place of coordinates for spatial models.
* new datasets
* website with vignettes

## Changes

* Functions W.diag and NC.diag now have S3 plot and print methods
* Changes to arguments (backward compatible) to xdat throughout
* Many dependencies used by single functions are now listed in Suggest.

## Fixes:

* ext.coef correctly handles arrays with missing values (reported by M. Jousset)
* Optimization method in fit.gev now uses the PWM solution of Hosking (1985) as starting value
* gev.pll now returns confidence intervals for param = "quant" (reported by D. Dupuis)
* Fixes to NHPP order statistics density (returns -Inf outside of domain, also correctly evaluate for boundary case when xi=-1)
* Optimization routines fit.pp, fit.gev and fit.rlarg now return correct MLE when solution lies on boundary (xi=-1) and are more robust to failure (gradients for nlminb return large values rather than NAs which caused the algorithm to stop).
* Grimshaw's algorithm sometimes returned incorrect value because of too low tolerance for eta near zero. Set back to default settings.
* fit.gpd(..., method = "obre") now returns additional failure messages if the algorithm drifts towards infeasible values.
* rparp now correctly handles xi=0
* Extended GP model now has 'step' for discretization, and a valid distribution function that returns real arguments whenever the input is finite (#9)
* W.diag and egp.fitrange include arguments for changing 'par' (#10)
* smith.penult now computes reciprocal hazard and it's derivative on the log scale (when possible) to avoid numerical overflow.

# mev v.1.13 (Release date: 2019-12-17)
==============

## Fixes:

* rlarg.infomat incorrect sign for expected information for d=1
* rmev function did not work for `alog` and `aneglog` (reported by Michael Lalancette)
* Remove class()!= "matrix due to changes in R 4.0.0
* egp.retlev now returns invisible object (new ordering and format).

## Changes:

* New S3 methods for objects returned by "fit" routines, for use in "lax" package

# mev v.1.12 (Release date: 2019-06-24)
==============

## New:

* Function 'taildep' is multivariate equivalent to 'evd::chiplot'
* Function 'rparpcs' for simulating from elliptical Pareto processes associated with max
* 'rgparp' for simulation of generalized R Pareto processes
* Exponent measure for Brown-Resnick and extremal Student model
* 'spunif' for semi-parametric transform to uniform (tail modelled using GP)
* 'rparp' now has attributes to give acceptance rate of accept/reject procedure.
* Exponent measure estimators 'extcoef'
* New method for robust OBRE estimates of GP parameters
* Functions 'fit.gev','fit.gpd', 'fit.rlarg' and 'fit.pp' for maximum likelihood estimation
* Default printing and plot optims (p-p and q-q plots) for mev_ objects
* Changes to optimization routine for pp  in 'W.diag' function, use of expected information matrix, change to default tuning parameter
* New datasets: eskrain, maiquetia, nidd, venice, w1500m

## Fixes:

* Scaling of Brown-Resnick is now consistent with literature (half of semivariogram)
* Degrees of freedom argument in C++ code for 'rmev' with extremal Student family was incorrect
* 'rparp' now has a hard bound to ensure the vector of simulated vectors fits in memory.
* Change to Grimshaw routine to ensure shape not less than -1, constrained optimization method
* Information matrix, score for GEV now interpolated in a neighborhood of zero to preserve continuity and avoid numerical overflow.
* Information matrix of GEV: scale factor off. All scores and information matrix have been checked and limits as xi -> 0  implemented.
* ext.index number of exceedances off by one, causing Inf in weight vector for lm (thanks to @MCristinaMiranda). Warnings now silenced by default.

## Changes:

* coordinates for rmev, rparp, etc. now use `coord` instead of `loc` to avoid confusion with location parameter in `rgparp`
* Removed dependency to 'ismev', 'rootsolve' (replaced with 'nleqslv' routine) and 'numDeriv'.
* Change to plot for profile log-likelihood methods
* 'smith.penult' function has new arguments (backward compatible). The quantiles 'u' are now returned for Smith penultimate approximations with 'method = "pot"'


# mev v.1.11 (Release date: 2018-02-23)
==============

## New:

* Function 'rparp' for simulation from R-Pareto Processes via rejection sampling
* Function 'gev.pll' and 'gpd.pll' for penalized profile likelihood and tangent exponential model approximations
* New functions 'chibar', 'angextrapo' and 'lambdadep' for bivariate and multivariate model estimation, based on work of Tawn et al.
* Dirichlet mixture smoothing for empirical angular distribution of de Carvalho et al. (2013)
* Functions 'gev.mle' and 'gpd.mle' for maximum likelihood estimates of transformed parameters
* Functions 'gev.abias' and 'gpd.abias' for asymptotic bias of block maxima for fixed sample sizes or fixed thresholds

## Changes:

* Functions 'rmev', 'rmevspec', etc. now only accept variogram functions 'vario' that have distance as argument
* Simulation from 'rmev' and 'rmevspec' faster to refactoring of code
* Function 'smith.penult' now has a 'family' as argument for specifying distributions via a string
* Function 'gev.tem' and 'gpd.tem' are now a wrapper for 'gev.pll' and 'gpd.pll', respectively. Routine should be more robust
* TEM corrections now handle more options
* Clarifications in the vignette about the asymmetric negative logistic model (thanks to A. Stephenson)

## Fixes:

* Fixed incorrect scaling in 'infomat.test' (thanks to P. Northrop)
* Model "br" now simulates from stationary version only if argument 'sigma' is provided, and otherwise samples intrinsically Gaussian processes
* Display of p-value matrix for 'infomat.test'

# mev v.1.10 (Release date: 2017-02-01)
==============

## New:

* Added 'negdir' model to rmev

## Changes:

* Fixed argument matching in function 'egp2'

# mev v.1.9x
==============

## New:

* Changes to angmeas to include different weighting if the region of interest is 'max' or 'min'

## Fixes:

* Fixed bug affecting angmeas in the bivariate case that would cause the method to crash

# mev v.1.8x
==============

## New:

* Bias-correction, TEM added
* Penultimate approximations (Papastathopoulos & Tawn, 2013), Naveau et al. (2014) and Smith (1987)

## Changes:

* model "br" is now distinct from "hr"

## Fixes:

* fixed invalid random number generation from logistic model for near-independence cases


# mev v1.7 (Release date: 2016-06-07)
==============
## New:

* ext.index Extremal index estimates based on interexceedance and gap times
* infomat.test Information matrix test of Suveges and Davison (2010)

## Fixes:

* fixed an error in the acceptance rate for the `gp.fit` MCMC


# mev v1.6.1 (Release date: 2016-03-15)
==============

## Fixes: 

* fixed an error in the normal sampler (affecting version 1.5 and 1.6). All simulations of Brown-Resnick or extremal-Student were affected by the mistake

# mev v1.6 (Release date: 2016-03-08)
==============

## New:

*Empirical and Euclidean likelihood estimation of spectral measure

## Changes:

*`gp.fit` ample changes to the function, in particular a fix for the printing method, handling of errors and inclusion of the Zhang (2010) method and MCMC algorithm for the latter. This function is still preliminary and may updated in the nearby future to include further possibilities.

# mev v1.5 (Release date: 2016-02-16)
==============

## New:

* Wadsworth (2015) Technometrics's proposal for threshold selection based on NHPP superposition
* Northrop & Coleman diagnostic (2014) Extremes for shape equality and p-value path

# mev v1.4 (not on CRAN)
==============

## Fixes:

* fixed error for simulation on grids

## Changes: 

* check for marginal mean constraint for the Dirichlet mixture now has tolerance

# mev v1.3 (Release date: 2015-10-05)
==============

## Changes:

* Extremal Dirichlet model now implemented with "ef"
* Added Smith and asymmetric (negative) logistic model (differs from bivariate setting for aneglog, given that the later is not a valid DF according to Stephenson).
* rmev can now return arrays for random fields on regular grids ("hr","exstud" and "smith" models).



# mev v1.2 (Release date: 2015-08-23)
==============

## Changes:

* Added the negative bilogistic and the scaled Dirichlet models.
* Extremal Dirichlet model implemented with "sm" only.

# mev v1.1 (Release date: 2015-08-19)
==============

## Changes:

* Implementation of sampler from spectral distribution, moving rdirspec and rbilogspec to background along with other functions.
* Fixed a typo in rPextstud setting arguments of newly created arma vector to zero

mev v1.0 (Release date: 2014-08-16)
==============
