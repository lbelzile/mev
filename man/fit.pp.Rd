% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mle.R
\name{fit.pp}
\alias{fit.pp}
\title{Maximum likelihood estimation of the point process of extremes}
\usage{
fit.pp(
  xdat,
  threshold = 0,
  npp = 1,
  np = NULL,
  method = c("nlminb", "BFGS"),
  start = NULL,
  show = FALSE,
  fpar = NULL,
  warnSE = FALSE
)
}
\arguments{
\item{xdat}{a numeric vector of data to be fitted.}

\item{threshold}{the chosen threshold.}

\item{npp}{number of observation per period. See \bold{Details}}

\item{np}{number of periods of data, if \code{xdat} only contains exceedances.}

\item{method}{the method to be used. See \bold{Details}. Can be abbreviated.}

\item{start}{named list of starting values}

\item{show}{logical; if \code{TRUE} (the default), print details of the fit.}

\item{fpar}{a named list with optional fixed components \code{loc}, \code{scale} and \code{shape}}

\item{warnSE}{logical; if \code{TRUE}, a warning is printed if the standard errors cannot be returned from the observed information matrix when the shape is less than -0.5.}
}
\value{
a list containing the following components:
\itemize{
\item \code{estimate} a vector containing all parameters (optimized and fixed).
\item \code{std.err} a vector containing the standard errors.
\item \code{vcov} the variance covariance matrix, obtained as the numerical inverse of the observed information matrix.
\item \code{threshold} the threshold.
\item \code{method} the method used to fit the parameter. See details.
\item \code{nllh} the negative log-likelihood evaluated at the parameter \code{estimate}.
\item \code{nat} number of points lying above the threshold.
\item \code{pat} proportion of points lying above the threshold.
\item \code{convergence} components taken from the list returned by \code{\link[stats]{optim}}.
Values other than \code{0} indicate that the algorithm likely did not converge (in particular 1 and 50).
\item \code{counts} components taken from the list returned by \code{\link[stats]{optim}}.
}
}
\description{
Data above \code{threshold} is modelled using the limiting point process
of extremes.
}
\details{
The parameter \code{npp} controls the frequency of observations.
If data are recorded on a daily basis, using a value of \code{npp = 365.25}
yields location and scale parameters that correspond to those of the
 generalized extreme value distribution fitted to block maxima.
}
\examples{
data(eskrain)
pp_mle <- fit.pp(eskrain, threshold = 30, np = 6201)
plot(pp_mle)
}
\references{
Coles, S. (2001), An introduction to statistical modelling of extreme values. Springer : London, 208p.
}
