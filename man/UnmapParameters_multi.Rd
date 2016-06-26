\name{UnmapParameters_multi}
\alias{UnmapParameters_multi}
\title{
	Inverse of \link{MapParameters_multi}
}
\description{
		Transform distribution parameters into the unrestricted parameters. The unrestricted vector of parameters is updated using the GAS recursion.
	}
\usage{
UnmapParameters_multi(vTheta, Dist, iN, iK)
}
\arguments{
\item{vTheta}{
\code{numeric} Vector parameters, see Details.}
%
\item{Dist}{
\code{character} Label of the conditional distribution, see \link{DistInfo}.}
%
\item{iN}{
\code{numeric} Cross sectional dimension. Note that only \code{iN<5} is supported.}
%
\item{iK}{
\code{numeric} Length of \code{vTheta}.}
%
}
\details{
	The order of the parameters is generally: locations, scales, correlations, shape. When the distribution defined by \code{Dist} does not have, say, the shape parameter, this should be simply omitted. See also \link{DistInfo} for specific distributions.
}
\value{A \code{numeric} vector of parameters.}
\author{Leopoldo Catania}
\note{}
\examples{
\dontrun{
# unmap parameters for the Multivariate Student-t distribution with N=3
library(GAS)

iN = 3

Dist = "mvt"

vMu  = c(0.1,0.2,0.3)   # vector of location parameters (this is not transformed)
vPhi = c(1.0, 1.2, 0.3) # vector of scales parameters for the firs, second and third variables.

vRho = c(0.1,0.2,0.3)   # This represents vec(R), where R is the correlation matrix.
                        # Note that is up to the user to ensure that vec(R) implies a proper correlation matrix

vTheta = c(vMu, vPhi, vRho, 7) # vector of parameters such that the degrees of freedom are 7.

vTheta_tilde = UnmapParameters_multi(vTheta, Dist, iN, iK = length(vTheta))

vTheta_tilde

}
}
