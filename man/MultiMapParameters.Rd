\name{MultiMapParameters}
\alias{MultiMapParameters}
\alias{build_mR}
\alias{MapR_C}
\alias{build_vR}
\title{Mapping function for multivariate distributions}
\description{
	Map unrestricted vector of parameters into the proper space. This function transforms
	the parameters updated using the GAS recursion into their proper space.
	}
\usage{
MultiMapParameters(Theta_tilde, Dist, N)
}
\arguments{
\item{Theta_tilde}{
\code{numeric} Vector of reparametrised parameters, see Details.}
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
The input vector contains, in order, N location parameters,
N log-scale parameters, N(N-1)/2 correlation angles, and,
for the multivariate Student-t distribution, one transformed
degrees-of-freedom parameter.

Locations are unchanged and scales are obtained by exponentiation,
subject to the package's numerical safeguards. The Student-t
degrees-of-freedom parameter uses the bounded transformation
implemented by the package. Use \code{\link{MultiUnmapParameters}}
to construct the transformed vector from distribution parameters.

The correlation angles are ordered as
\eqn{\phi_{12},\phi_{13},\ldots,\phi_{1N},
\phi_{23},\ldots,\phi_{2N},\ldots,\phi_{N-1,N}}.
They define an upper-triangular matrix X through
\deqn{X_{rb}=\cos(\phi_{rb})\prod_{s<r}\sin(\phi_{sb}),\quad r<b,}
\deqn{X_{bb}=\prod_{s<b}\sin(\phi_{sb}),\qquad X_{rb}=0,\quad r>b,}
where empty products equal one. The correlation matrix is
\deqn{R=X^\top X.}

The returned correlation entries follow the same pair ordering.
They contain only the distinct off-diagonal entries, rather than
the full vectorization of R.

This construction gives a positive semidefinite matrix with unit
diagonal. It is positive definite when all diagonal entries of X
are nonzero. In particular, angles strictly between zero and pi
give a positive definite matrix in exact arithmetic. Angles near
singular configurations can cause numerical ill-conditioning.

The correlation Jacobian is computed analytically for every
supported dimension by differentiating \eqn{R=X^\top X}.
}

\value{A \code{numeric} vector of parameters.}
\author{Leopoldo Catania}
\examples{
N <- 5L
Dist <- "mvt"

# Equicorrelation matrix with correlation 0.2.
# Parameter order: locations, scales, correlations, degrees of freedom.
Theta <- c(
    rep(0, N),
    rep(1, N),
    rep(0.2, N * (N - 1L) / 2L),
    8
)

Theta_tilde <- MultiUnmapParameters(Theta, Dist, N)
Theta_back <- MultiMapParameters(Theta_tilde, Dist, N)

stopifnot(max(abs(Theta_back - Theta)) < 1e-10)
}
