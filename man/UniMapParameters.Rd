\name{UniMapParameters}
\alias{UniMapParameters}
\title{
	Mapping function for univariate distributions
}
\description{
	Map unrestricted vector of parameters into the proper space. 
	This function transforms the parameters updated using the GAS recursion into their proper space.
	}
\usage{
UniMapParameters(Theta_tilde, Dist)
}
\arguments{
\item{Theta_tilde}{
\code{numeric} Vector of reparametrised parameters, see Details.}
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
# Map unrestricted parameters for the Student-t distribution.
library("GAS")

Dist = "std"

# Vector of unconditional reparametrised parameters such that,
# Theta = c(0, 1.5 ,7), i.e., location = 0, scale = 1.5,
# degrees of freedom = 7.

# LowerNu() prints the lower bound numerical parameter for the degree
# of freedom, see help(LowerNu).

Theta_tilde = c(0.1, log(1.5), log(7 - LowerNu()))

Theta = UniMapParameters(Theta_tilde, Dist)

Theta
}
