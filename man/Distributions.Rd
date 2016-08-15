\name{distributions}
\alias{distributions}
\alias{ddist_univ}
\alias{ddist_multi}
\alias{pdist_univ}
\alias{qdist_univ}
\alias{rdist_univ}
\alias{mdist_univ}
\alias{Score_univ}
\alias{Quantiles}
\alias{Score_multi}
\alias{rmvt_mat}
\alias{IM_univ}
\alias{rdist_multi}
\title{
	Distributions of the GAS package.
}
\description{
	Density, distribution function, quantile funcion, random generator, moments, scores and information matrix of univariate and multivariate distributions of the GAS package.
}
\usage{
ddist_univ(dY, vTheta, Dist, bLog)
pdist_univ(dQ, vTheta, Dist)
qdist_univ(dP, vTheta, Dist)
rdist_univ(vTheta, Dist)
mdist_univ(vTheta, Dist)
Score_univ(dY, vTheta, Dist)
IM_univ(vTheta, Dist)

ddist_multi(vY, vTheta, iN, Dist, bLog)
rdist_multi(vTheta, iN, Dist)
Score_multi(vY, vTheta, iN, Dist)

}
\arguments{
\item{dY, dQ}{
\code{numeric} scalar quantile.}
%
\item{dP}{
\code{numeric} scalar probability.}
%
\item{vTheta}{
\code{numeric} vector of distribution parameters. The order of the parameters is generally: \code{location, scale, skewness, shape, shape2}. When the distribution defined by \code{Dist} does not have, say, the \code{shape} parameter, this should be simply omitted. See also \link{DistInfo} for specific distributions.}
%
\item{Dist}{
\code{character}. Label of the conditional distribution, see \link{DistInfo}.}
%
\item{bLog}{
\code{logical} if \code{TRUE}, the density value \eqn{p(y)} is given as \eqn{log(p(y))}.}
%
\item{vY}{
\code{numeric} vector of quantiles.}
%
\item{iN}{
\code{numeric} integer. \code{vY} dimension.}
}
\details{
The function \code{mdist_univ} returns a vector with four elements: mean, variance, skewness and kurtosis coefficients. The functions
\code{Score_univ} and \code{IM_univ} returns the score and the Fisher information matrix for univariate distributions. The function \code{Score_multi} reuturns the score for multivariate distributions. See \link{DistInfo} for the lists of supported distributions. These functions are not vectorized.
}
\value{A \code{list} with elements: LRuc, LRcc, DQ, AD, AE.}
\author{Leopoldo Catania}
%
\examples{
# Skew Student-t distribution

# log density
vTheta = c("location" = 0, "scales" = 1, "skewness" = 1.2, "shape" = 7)

ddist_univ(dY = 0.5, vTheta, "sstd", TRUE)

# probability
pdist_univ(dQ = -1.69, vTheta, "sstd")

#quantile
qdist_univ(dP = 0.05, vTheta, "sstd")

#random generator
rdist_univ(vTheta, "sstd")

#moments
mdist_univ(vTheta, "sstd")

}

