\name{MapParameters_multi}
\alias{MapParameters_multi}
\title{
	Mapping function for univariate distributions
}
\description{
	Map unrestricted vector of parameters into the proper space. This function transforms the parameters updated using the GAS recursion into their proper space.
	}
\usage{
MapParameters_multi(vTheta_tilde, Dist, iN, iK)
}
\arguments{
\item{vTheta_tilde}{
\code{numeric} Vector of reparametrised parameters, see Details.}
%
\item{Dist}{
\code{character} Label of the conditional distribution, see \link{DistInfo}.}
%
\item{iN}{
\code{numeric} Cross sectional dimension. Note that only \code{iN<5} is supported.}
%
\item{iK}{
\code{numeric} Length of \code{vTheta_tilde}.}
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
# map unrestricted parameters for the Multivariate Student-t distribution with N=3
library(GAS)

iN = 3

Dist = "mvt"

vMu_tilde  = c(0.1,0.2,0.3) # vector of location parameters (this is not transformed)
vPhi_tilde = c(log(1.0), log(1.2), log(0.3)) # vector of unrestricted scales parameters such that
                                             # the scales will be equal to 1.0, 1.2 and 0.3, for the first, second and
                                             # third variables, respectively.

vRho_tilde = UnMapR_C(c(0.1,0.2,0.3), iN)    # The vector c(0.1,0.2,0.3) represents vec(R), where R is the correlation matrix.
                                             # Note that is up to the user to ensure that vec(R) implies a proper correlation matrix
                                             # The function UnMapR_C transforms vec(R) in a vector of unrestricted parameters. It is
                                             # the inverse of the hyperspherical coordinates transformration.

vTheta_tilde = c(vMu_tilde, vPhi_tilde , vRho_tilde, log(7 - 2.01)) # vector of unconditional reparametrised parameters such that the                                                                        # degrees of freedom are 7.


vTheta = MapParameters_multi(vTheta_tilde, Dist, iN, iK = length(vTheta_tilde))

vTheta

}
}
