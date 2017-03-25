\name{distributions}
\alias{distributions}
\alias{ddist_Uni}
\alias{ddist_Multi}
\alias{pdist_Uni}
\alias{qdist_Uni}
\alias{rdist_Uni}
\alias{mdist_Uni}
\alias{Score_Uni}
\alias{Quantiles}
\alias{Score_Multi}
\alias{rmvt_mat}
\alias{IM_Uni}
\alias{rdist_Multi}
\title{
	Distributions of the GAS package
}
\description{
	Density, distribution function, quantile function, random generator, moments, scores and 
	information matrix of univariate and multivariate distributions of the GAS package.
}
\usage{
ddist_Uni(y, Theta, Dist, log = FALSE)
pdist_Uni(q, Theta, Dist)
qdist_Uni(p, Theta, Dist)
rdist_Uni(Theta, Dist)
mdist_Uni(Theta, Dist)
Score_Uni(y, Theta, Dist)
IM_Uni(Theta, Dist)

ddist_Multi(y, Theta, Dist, log = FALSE)
rdist_Multi(Theta, N, Dist)
Score_Multi(y, Theta, Dist)

}
\arguments{
\item{y, q}{
\code{numeric} Scalar quantile. For \code{Score_Multi} and \code{ddist_Multi}, \code{y} is a \code{numeric} vector.}
%
\item{p}{
\code{numeric} Scalar probability.}
%
\item{Theta}{
\code{numeric} Vector of distribution parameters. The order of the parameters is generally: \code{location, scale, skewness, shape, shape2}. 
When the distribution defined by \code{Dist} does not have, say, the \code{shape} parameter, this should be simply omitted. 
See also \link{DistInfo} for specific distributions.}
%
\item{Dist}{
\code{character} Label of the conditional distribution, see \link{DistInfo}.}
%
\item{log}{
\code{logical} If \code{TRUE}, the density value \eqn{p(y)} is given as \eqn{log(p(y))}. Dy Default \code{log = FALSE}.}
%
\item{N}{
\code{numeric} Integer. cross-sectional dimension for the multivariate case.}
}
\details{
The function \code{mdist_Uni} returns a vector with four elements: \code{mean}, \code{variance}, \code{skewness} 
and \code{kurtosis} coefficients. The functions \code{Score_Uni} and \code{IM_Uni} returns the score and the Fisher information matrix 
for univariate distributions. The function \code{Score_Multi} returns the score for multivariate distributions. 
See \link{DistInfo} for the lists of supported distributions. 
These functions are not vectorized. \code{ddist_Uni} and \code{ddist_Multi} give the density, \code{pdist_Uni} gives 
the distribution function, \code{qdist_Uni} gives the quantile function, and \code{rdist_Uni} and \code{rdist_Multi} generate random deviates.
}
\value{
\enumerate{
\item{\code{numeric} scalar for: \code{ddist_Uni}, \code{pdist_Uni}, \code{qdist_Uni}, \code{rdist_Uni},}
\item{\code{numeric} vector for: \code{Score_Uni}, \code{Score_Multi} and \code{rdist_Multi}},
\item{\code{matrix} for \code{IM_Uni}.}
}
}
\author{Leopoldo Catania}
%
\examples{
# Skew Student-t distribution

# log density
Theta = c("location" = 0, "scales" = 1, "skewness" = 1.2, "shape" = 7)

ddist_Uni(y = 0.5, Theta, "sstd", TRUE)

# probability
pdist_Uni(q = -1.69, Theta, "sstd")

#quantile
qdist_Uni(p = 0.05, Theta, "sstd")

#random generator
rdist_Uni(Theta, "sstd")

#moments
mdist_Uni(Theta, "sstd")
}

