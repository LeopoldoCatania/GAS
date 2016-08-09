\name{MultiGASSim}
\alias{MultiGASSim}
\title{
	Simulate Multivariate GAS processes
}
\description{
	Simulate Multivariate GAS processes.
}
\usage{
MultiGASSim(iT, N, kappa, A, B, Dist, ScalingType)
}
\arguments{
\item{iT}{
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
	All the information regarding the supported multivariate conditional distributions can be investigated using the \link{DistInfo} function. The model is specified as: \deqn{y_{t}\sim p(y|\theta_{t})} where \eqn{\theta_{t}} is the vector of parameters for the density \eqn{p(y|.)}. Note that, \eqn{\theta_{t}} includes also those parameters that are not time-varying. The GAS recursion for \eqn{\theta_{t}} is: \deqn{\theta_{t} = h(\tilde{\theta}_{t})}\deqn{\tilde{\theta}_{t}=\kappa + A*s_{t-1} + B*\tilde{\theta}_{t-1}} where \eqn{h(.)} is the mapping function (see \link{MultiMapParameters}) and \eqn{\tilde{\theta}_{t}} is the vector of reparametrised parameters. The process is initialized at \eqn{\theta_{1}=(I - B)^{-1}\kappa}, where \eqn{\kappa} is the \code{Kappa} vector. The vector \eqn{s_{t}} is the scaled score of \eqn{p(y|.)} with respect to \eqn{\tilde{\theta}_{t}}.
}
\value{
An object of the class \link{mGASSim}
}
\references{
Creal, D., Koopman, S. J., & Lucas, A. (2013). Generalized autoregressive score models with applications. Journal of Applied Econometrics, 28(5), 777-795.\cr
%

Harvey, A. C. (2013). Dynamic models for volatility and heavy tails: with applications to financial and economic time series (Vol. 52). Cambridge University Press.\cr
}
\author{Leopoldo Catania}
\examples{
# Simulate from a GAS process with Multivariate Student-t conditional
# distribution, time-varying locations, scales, correlations
# and fixed shape parameter.
library(GAS)

set.seed(786)

iT     = 1000 # number of observations to simulate
N     = 3     # trivariate series
Dist   = "mvt" # conditional Multivariate Studen-t distribution

# build unconditional vector of reparametrised parameters

Mu  = c(0.1,0.2,0.3)   # vector of location parameters (this is not transformed)
Phi = c(1.0, 1.2, 0.3) # vector of scale parameters for the firs, second and third variables.

Rho = c(0.1,0.2,0.3)   # This represents vec(R), where R is the correlation matrix.
# Note that is up to the user to ensure that vec(R) implies a proper correlation matrix

Theta = c(Mu, Phi, Rho, 7) # vector of parameters such that the degrees of freedom are 7.

kappa = MultiUnmapParameters(Theta, Dist, N)

A     = matrix(0,length(kappa),length(kappa))

# update scales and correlations, do not update locations and shape parameters
diag(A) = c(0, 0, 0, 0.05, 0.01, 0.09, 0.01, 0.04, 0.07, 0)

B     = matrix(0,length(kappa),length(kappa))

# update scales and correlations, do not update locations and shape parameters
diag(B) = c(0, 0, 0, 0.7, 0.7, 0.5, 0.94, 0.97, 0.92, 0)

Sim = MultiGASSim(iT,N, kappa, A, B, Dist, ScalingType = "Identity")

Sim

}
