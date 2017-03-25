\name{UniUnmapParameters}
\alias{UniUnmapParameters}
\title{
	Unmapping function for univariate distributions, i.e. inverse of \link{UniMapParameters}
}
\description{
	Transform distribution parameters into the unrestricted parameters. 
	The unrestricted vector of parameters is updated using the GAS recursion.
	}
\usage{
UniUnmapParameters(Theta, Dist)
}
\arguments{
\item{Theta}{
\code{numeric} Vector of parameters, see Details.}
%
\item{Dist}{
\code{character} Label of the conditional distribution, see \link{DistInfo}.}
%
}
\details{
The order of the parameters is generally: \code{location}, \code{scale}, \code{skewness}, \code{shape}, \code{shape2}. 
When the distribution defined by \code{Dist} does not have, say, the \code{shape} parameter, this should be simply omitted. 
See also \link{DistInfo} for specific distributions.
}
\value{A \code{numeric} vector of parameters.}
\author{Leopoldo Catania}
\examples{
# Unmap parameters for the Student-t distribution
library("GAS")

Dist = "std"

# Vector of parameters such that,
# Theta = c(0, 1.5 ,7), i.e., location = 0, scale = 1.5,
# degrees of freedom = 7.

Theta = c(0.1, 1.5, 7)

Theta_tilde = UniUnmapParameters(Theta, Dist)

Theta_tilde

# It works.
all(abs(UniMapParameters(Theta_tilde, Dist) - Theta) < 1e-16)
}
