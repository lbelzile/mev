% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{.rmevA1}
\alias{.rmevA1}
\title{Multivariate extreme value distribution sampling algorithm via angular measure}
\usage{
.rmevA1(n, d, par, model, Sigma, loc)
}
\arguments{
\item{n}{sample size}

\item{d}{dimension of the multivariate distribution}

\item{par}{a vector of parameters}

\item{model}{integer, currently ranging from 1 to 9, corresponding respectively to
(1) \code{log}, (2) \code{neglog}, (3) \code{dirmix}, (4) \code{bilog},
(5) \code{extstud}, (6) \code{br}, (7) \code{ct} and \code{sdir}, (8) \code{smith} and (9) \code{hr}.}

\item{Sigma}{covariance matrix for Brown-Resnick, Smith and extremal student. Conditionally negative definite
matrix of parameters for the Huesler--Reiss model. Default matrix for compatibility}

\item{loc}{matrix of location for Smith model.}
}
\value{
a \code{n} by \code{d} matrix containing the sample
}
\description{
This algorithm corresponds to Algorithm 1 in Dombry, Engelke and Oesting (2016),
using the formulation of the Dirichlet mixture of Coles and Tawn (1991)
as described and derived in Boldi (2009) for the bilogistic and extremal
Dirichlet model. Models currently implemented include logistic, negative
logistic, extremal Dirichlet and bilogistic MEV.
}
\keyword{internal}
