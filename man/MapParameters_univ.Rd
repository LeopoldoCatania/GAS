\name{MapParameters_univ}
\alias{MapParameters_univ}
\title{
	Mapping function for univariate distributions
}
\description{
	Map unrestricted vector of parameters into the proper space. This function transforms the parameters updated using the GAS recursion into their proper space.
	}
\usage{
MapParameters_univ(vTheta_tilde, Dist, iK)
}
\arguments{
\item{vTheta_tilde}{
\code{numeric} Vector of reparametrised parameters, see Details.}
%
\item{Dist}{
\code{character} Label of the conditional distribution, see \link{DistInfo}.}
%
\item{iK}{
\code{numeric} Length of \code{vTheta_tilde}.}
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
# map unrestricted parameters for the Student-t distribution
library(GAS)

Dist = "std"

vTheta_tilde = c(0.1, log(1.5) , log(7 - 2.01)) # vector of unconditional reparametrised parameters such that,
                                                # \eqn{\theta = (0, 1.5 ,7)'}, i.e. location = 0, scale = 1.5, degrees of freedom = 7

vTheta = MapParameters_univ(vTheta_tilde, Dist, iK = length(vTheta_tilde))

vTheta

}
}
