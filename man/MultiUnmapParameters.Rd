\name{MultiUnmapParameters}
\alias{UnMapR_C}
\alias{MultiUnmapParameters}
\title{
	Inverse of \link{MultiMapParameters}
}
\description{
		Transform distribution parameters into the unrestricted parameters.
		The unrestricted vector of parameters is updated using the GAS recursion.
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
An integer greater than or equal to 2 giving the cross-sectional
dimension. Dimensions greater than four are supported.
}
%
}
\details{
The input vector contains, in order, N locations, N positive scales,
N(N-1)/2 distinct correlations, and, for the multivariate Student-t
distribution, one degrees-of-freedom parameter.

Correlations must be ordered as
\eqn{\rho_{12},\rho_{13},\ldots,\rho_{1N},
\rho_{23},\ldots,\rho_{2N},\ldots,\rho_{N-1,N}}.
Together with a unit diagonal, they must define a positive definite
correlation matrix.

The function returns locations, log-scales, correlation angles,
and, where applicable, transformed degrees of freedom.
The degrees of freedom must lie strictly inside the bounds used
by the package.

The inverse correlation mapping selects angles in the interval
\eqn{(0,\pi)}. Because the forward mapping is periodic, transforming
arbitrary angles to correlations and back need not recover the
original angles. Transforming admissible distribution parameters
to the unrestricted representation and back recovers the original
distribution parameters up to numerical precision.

See \code{\link{MultiMapParameters}} for the hyperspherical mapping.
}
\value{A \code{numeric} vector of parameters.}
\author{Leopoldo Catania}
\examples{
# Unmap parameters for the Multivariate Student-t distribution with N=3
library(GAS)

N = 3

Dist = "mvt"

# Vector of location parameters (this is not transformed).
Mu  = c(0.1, 0.2, 0.3)

# Vector of scales parameters for the firs, second and third variables.
Phi = c(1.0, 1.2, 0.3)

# This represents the distinct correlations in the order \(\rho_{12},\rho_{13},\rho_{23}\).
# Note that is up to the user to ensure that vec(R) implies a proper correlation matrix
Rho = c(0.1, 0.2, 0.3)

# Vector of parameters such that the degrees of freedom are 7.
Theta = c(Mu, Phi, Rho, 7)

Theta_tilde = MultiUnmapParameters(Theta, Dist, N)

Theta_tilde

# It works
all(abs(MultiMapParameters(Theta_tilde, Dist, N) - Theta) < 1e-10)
}
