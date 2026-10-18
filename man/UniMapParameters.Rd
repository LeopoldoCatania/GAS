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

For \code{Dist = "skellam"}, the public distribution parameters are
mean \eqn{\mu} and variance \eqn{v>|\mu|}. The state coordinates are
\deqn{f_1=\mu,\qquad f_2=\tfrac12\log(v^2-\mu^2).}
The forward map is
\deqn{\mu=f_1,\qquad v=\sqrt{f_1^2+\exp(2f_2)}.}
Its Jacobian, with parameter coordinates in rows and state coordinates
in columns, is
\deqn{J=\begin{pmatrix}1&0\\f_1/v&\exp(2f_2)/v\end{pmatrix}.}
States that cannot represent a finite interior parameter pair in
floating-point arithmetic return \code{NA}; the coordinates are not
clipped independently. Models estimated using the previous Skellam
mapping must be refitted.


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

# Use the inverse mapping so the distribution-specific bounds match.
Theta_tilde = UniUnmapParameters(c(0, 1.5, 7), Dist)

Theta = UniMapParameters(Theta_tilde, Dist)

Theta
}
