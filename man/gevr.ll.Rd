% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univdist.R
\name{gevr.ll}
\alias{gevr.ll}
\alias{gevr.ll.optim}
\title{Negative log likelihood of the generalized extreme value distribution (return levels)}
\usage{
gevr.ll(par, dat, p)

gevr.ll.optim(par, dat, p)
}
\arguments{
\item{par}{vector of \code{retlev}, \code{scale} and \code{shape}}

\item{dat}{sample vector}

\item{p}{tail probability, corresponding to \eqn{(1-p)}th quantile for \eqn{z}}
}
\description{
Negative log likelihood of the generalized extreme value distribution (return levels)

Negative log likelihood parametrized in terms of location, log return level and shape in order to perform unconstrained optimization
}
\seealso{
\code{\link{gevr}}
}
\keyword{internal}
