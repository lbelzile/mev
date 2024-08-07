% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tailindex.R
\name{shape.hill}
\alias{shape.hill}
\title{Hill's estimator of the shape parameter}
\usage{
shape.hill(xdat, k)
}
\arguments{
\item{xdat}{vector of positive observations}

\item{kmin}{minimum number of upper order statistics (exceedances) for the estimator}

\item{kmax}{maximum number of upper order statistics (exceedances) for the estimator}
}
\value{
a data frame with the number of order statistics \code{k} and the shape parameter estimate \code{shape}, or a single numeric value if \code{kmin} equals \code{kmax}.
}
\description{
Given a sample of positive observations, calculate the tail index or
shape parameter. The shape estimate returned is positive.
}
\examples{
xdat <- mev::rgp(n = 200, loc = 1, scale = 0.5, shape = 0.5)
shape.hill(xdat)
}
\references{
Hill, B.M. (1975). \emph{A simple general approach to inference about the tail of a distribution.} Annals of Statistics, \bold{3}, 1163-1173.
}
