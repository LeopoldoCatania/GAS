\name{MultiGASSim}
\alias{MultiGASSim}
\title{
	Simulate multivariate GAS processes
}
\description{
	Simulate multivariate GAS processes.
}
\usage{
MultiGASSim(T.sim, N, kappa, A, B, Dist, ScalingType)
}
\arguments{
\item{T.sim}{
\code{numeric} Length of the simulated time series.}
%
\item{N}{
\code{numeric} Cross sectional dimension (Only N<5 supported for now).}
%
\item{kappa}{
\code{numeric} vector of unconditional level for the reparametrized vector of parameters.}
%
\item{A}{
\code{matrix} of coefficients of dimension K x K that premultiply the conditional score in the GAS updating recursion, see Details.}
%
\item{B}{
\code{matrix} of autoregressive coefficients of dimension K x K, see Details.}
%
\item{Dist}{
\code{character} Label of the conditional distribution, see \link{DistInfo}.}
%
\item{ScalingType}{
\code{character} indicating the scaling mechanism for the conditional score. Only \code{ScalingType = "Identity"} is supported for multivariate distributions, see the function \link{DistInfo}.}
%
}
\details{
	All the information regarding the supported multivariate conditional distributions can be investigated using the \link{DistInfo} function. The model is specified as: \deqn{y_{t}\sim p(y|\theta_{t})} where \eqn{\theta_{t}} is the vector of parameters for the density \eqn{p(y|.)}. Note that, \eqn{\theta_{t}} includes also those parameters that are not time-varying. The GAS recursion for \eqn{\theta_{t}} is: \deqn{\theta_{t} = \Lambda(\tilde{\theta}_{t})}\deqn{\tilde{\theta}_{t}=\kappa + A*s_{t-1} + B*\tilde{\theta}_{t-1}} where \eqn{h(.)} is the mapping function (see \link{MultiMapParameters}) and \eqn{\tilde{\theta}_{t}} is the vector of reparametrised parameters. The process is initialized at \eqn{\theta_{1}=(I - B)^{-1}\kappa}, where \eqn{\kappa} is the \code{Kappa} vector. The vector \eqn{s_{t}} is the scaled score of \eqn{p(y|.)} with respect to \eqn{\tilde{\theta}_{t}}. See Ardia et. al. (2016) for further details.
}
\value{
An object of the class \link{mGASSim}
}
\references{
Ardia D, Boudt K and Catania L (2016). "Generalized Autoregressive Score Models in R: The GAS Package." Available at SSRN: \url{http://ssrn.com/abstract=2825380}.\cr
%

Creal D, Koopman SJ, Lucas A (2013). "Generalized Autoregressive Score Models with Applications." Journal of Applied Econometrics, 28(5), 777-795. \doi{10.1002/jae.1279}.\cr
%

Harvey AC (2013). Dynamic Models for Volatility and Heavy Tails: With Applications to Financial and Economic Time Series. Cambridge University Press.\cr
}
\author{Leopoldo Catania}
\examples{
# Simulate from a GAS process with Multivariate Student-t conditional
# distribution, time-varying locations, scales, correlations
# and fixed shape parameter.
library("GAS")

set.seed(786)

T.sim = 1000 # Number of observations to simulate.
N = 3        # Trivariate series.
Dist = "mvt" # Conditional Multivariate Studen-t distribution.

# Build unconditional vector of reparametrised parameters.

Mu = c(0.1, 0.2, 0.3)  # Vector of location parameters (this is not transformed).
Phi = c(1.0, 1.2, 0.3) # Vector of scale parameters for the firs, second and third variables.

Rho = c(0.1, 0.2, 0.3) # This represents vec(R), where R is the correlation matrix.
                       # Note that is up to the user to ensure that vec(R) implies a
                       # proper correlation matrix.

Theta = c(Mu, Phi, Rho, 7) # Vector of parameters such that the degrees of freedom are 7.

kappa = MultiUnmapParameters(Theta, Dist, N)

A = matrix(0, length(kappa), length(kappa))

# Update scales and correlations, do not update locations and shape parameters.
diag(A) = c(0, 0, 0, 0.05, 0.01, 0.09, 0.01, 0.04, 0.07, 0)

B = matrix(0, length(kappa), length(kappa))

# Update scales and correlations, do not update locations and shape parameters.
diag(B) = c(0, 0, 0, 0.7, 0.7, 0.5, 0.94, 0.97, 0.92, 0)

Sim = MultiGASSim(T.sim, N, kappa, A, B, Dist, ScalingType = "Identity")

Sim

}
