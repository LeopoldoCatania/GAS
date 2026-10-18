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
Quantiles(mTheta, Dist, vProbs)

ddist_Multi(y, Theta, Dist, log = FALSE)
rdist_Multi(Theta, N, Dist)
Score_Multi(y, Theta, Dist)

}
\arguments{
\item{mTheta}{A K by T matrix of distribution parameters, with one time point per column.}
\item{vProbs}{A nonempty numeric vector of probabilities in [0,1].}
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

\code{Quantiles} evaluates every requested probability at every time
point and returns a T by length(vProbs) matrix, even when the number
of probabilities equals T. It does not interpret equal lengths as
one probability per time point.



For \code{ast}, \code{ast1} and \code{ald}, quantiles are evaluated
analytically without a finite search interval. Their zero- and
one-probability quantiles are minus and plus infinity.

For \code{snorm} and \code{sstd}, the two-piece split has probability
\eqn{1/(1+\xi^2)}. Location and scale are the mean and standard deviation
after centering and standardizing the two-piece distribution.
The \code{sstd} parameterization requires degrees of freedom greater
than two.

For \code{skellam}, \code{Theta = c(mu, v)} denotes mean and variance,
with \eqn{v>|\mu|}. The underlying Poisson intensities are
\eqn{(v+\mu)/2} and \eqn{(v-\mu)/2}. Noninteger observations have zero
probability mass. The CDF is evaluated at the floor of its argument.
For \eqn{0<p<1}, the quantile is the smallest integer q satisfying
\eqn{F(q)\ge p}; the endpoints are minus and plus infinity.
The implementation returns \code{NA} if its bounded numerical search
cannot resolve a quantile within exactly representable integer arithmetic.

Kurtosis is reported as ordinary kurtosis, not excess kurtosis.
For \code{std}, the mean is undefined when degrees of freedom are at most
one; variance is infinite when the mean exists but degrees of freedom
are at most two. Skewness requires degrees of freedom greater than three.
For \code{std} and \code{sstd}, kurtosis is infinite when variance is
finite but degrees of freedom are at most four. Undefined moments are
reported as \code{NA}. These conventions agree with \code{ast} and
\code{ast1} at the corresponding moment-existence boundaries.


The function \code{mdist_Uni} returns a vector with four elements: \code{mean}, \code{variance}, \code{skewness} 
and \code{kurtosis} coefficients. The functions \code{Score_Uni} and \code{IM_Uni} returns the score and the Fisher information matrix 
for univariate distributions. The function \code{Score_Multi} returns the score for multivariate distributions. 
See \link{DistInfo} for the lists of supported distributions. 
The scalar distribution helpers are not vectorized; \code{Quantiles} evaluates a time-by-probability grid. \code{ddist_Uni} and \code{ddist_Multi} give the density, \code{pdist_Uni} gives
the distribution function, \code{qdist_Uni} gives the quantile function, and \code{rdist_Uni} and \code{rdist_Multi} generate random deviates.
}
\value{
\enumerate{
\item{\code{numeric} scalar for: \code{ddist_Uni}, \code{pdist_Uni}, \code{qdist_Uni}, \code{rdist_Uni},}
\item{\code{numeric} vector for: \code{Score_Uni}, \code{Score_Multi} and \code{rdist_Multi}},
\item{\code{matrix} for \code{IM_Uni} and \code{Quantiles}.}
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

