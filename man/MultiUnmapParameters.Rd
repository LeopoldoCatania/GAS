\name{MultiUnmapParameters}
\alias{MultiUnmapParameters}
\title{
	Inverse of \link{MultiMapParameters}
}
\description{
		Transform distribution parameters into the unrestricted parameters. The unrestricted vector of parameters is updated using the GAS recursion.
	}
\usage{
MultiUnmapParameters(Theta, Dist, N)
}
\arguments{
\item{Theta}{
\code{numeric} Vector parameters, see Details.}
%
\item{Dist}{
\code{character} Label of the conditional distribution, see \link{DistInfo}.}
%
\item{N}{
\code{numeric} Cross sectional dimension. Note that only \code{iN<5} is supported.}
%
}
\details{
	The order of the parameters is generally: locations, scales, correlations, shape. When the distribution defined by \code{Dist} does not have, say, the shape parameter, this should be simply omitted. See also \link{DistInfo} for specific distributions.
}
\value{A \code{numeric} vector of parameters.}
\author{Leopoldo Catania}
\examples{
\dontrun{
# unmap parameters for the Multivariate Student-t distribution with N=3
library(GAS)

N = 3

Dist = "mvt"
 # vector of location parameters (this is not transformed)
Mu  = c(0.1,0.2,0.3)

# vector of scales parameters for the firs, second and third variables.
Phi = c(1.0, 1.2, 0.3)

# This represents vec(R), where R is the correlation matrix.
# Note that is up to the user to ensure that vec(R) implies a proper correlation matrix
Rho = c(0.1,0.2,0.3)

#vector of parameters such that the degrees of freedom are 7.
Theta = c(Mu, Phi, Rho, 7)

Theta_tilde = MultiUnmapParameters(Theta, Dist, N)

Theta_tilde

}
}
