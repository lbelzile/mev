\name{mev-package}
\alias{mev-package}
\alias{mev}
\docType{package}
\title{
Multivariate Extreme Value Distributions
}
\description{
 Exact simulation from max-stable processes and multivariate extreme value distributions for various parametric models.
}
\details{

The package allows exact generation of multivariate extreme value
vectors or max-stable processes. For the latter, the user can provide a
variogram function along with a set of locations that serve as input. Models
implemented include the 1-parameter logistic and negative logistic as described
in the article, the bilogistic and Coles and Tawn extremal Dirichlet model using
the algorithm of Boldi (2009) and the Dirichlet mixture.
The extremal Student and Husler-Reiss (Brown-Resnick) models are also implemented.

Other features of the package include threshold diagnostic tests, empirical likelihood estimation for the spectral measure,
bias-correction and tangent exponential model approximations for the GP and GEV distributions, extended generalized Pareto distributions, information matrix test and penultimate approximations. Many of these are not
found elsewhere. \code{mev} complements packages \code{evd} and \link[ismev]{ismev}.

}
\author{
Leo Belzile

Maintainer: <leo.belzile@epfl.ch>
}
\concept{multivariate extreme value}
\concept{exact simulation}
\concept{max-stable process}
\concept{logistic}
\concept{Dirichlet}
\concept{extremal Student process}
\concept{Dirichlet mixture}
\concept{Coles and Tawn model}
\concept{bilogistic}
\concept{negative logistic}
\concept{negative bilogistic}
\concept{Brown-Resnick process}
