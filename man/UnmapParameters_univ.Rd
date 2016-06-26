\name{UnmapParameters_univ}
\alias{UnmapParameters_univ}
\title{
	Unmapping function for univariate distributions, i.e. inverse of \link{MapParameters_univ}
}
\description{
	Transform distribution parameters into the unrestricted parameters. The unrestricted vector of parameters is updated using the GAS recursion.
	}
\usage{
UnmapParameters_univ(vTheta, Dist, iK)
}
\arguments{
\item{vTheta}{
\code{numeric} Vector of parameters, see Details.}
%
\item{Dist}{
\code{character} Label of the conditional distribution, see \link{DistInfo}.}
%
\item{iK}{
\code{numeric} Length of \code{vTheta}.}
%
}
\details{
	The order of the parameters is generally: location, scale, skewness, shape, shape2. When the distribution defined by \code{Dist} does not have, say, the shape parameter, this should be simply omitted. See also \link{DistInfo} for specific distributions.
}
\value{A \code{numeric} vector of parameters.}
\author{Leopoldo Catania}
\note{}
\examples{
\dontrun{
# unmap parameters for the Student-t distribution
library(GAS)

Dist = "std"

vTheta = c(0.1, 1.5, 7) # vector of parameters such that,
                        # \eqn{\theta}_{t} is (0, 1.5 ,7), i.e. location = 0, scale = 1.5, degrees of freedom = 7

vTheta_tilde = UnmapParameters_univ(vTheta, Dist, iK = length(vTheta))

vTheta_tilde

}
}
