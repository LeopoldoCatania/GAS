\name{MultiGASSim}
\alias{MultiGASSim}
\title{
	Simulate Multivariate GAS processes
}
\description{
	Simulate Multivariate GAS processes.
}
\usage{
MultiGASSim(iT, iN,vKappa, mA, mB, Dist, ScalingType)
}
\arguments{
\item{iT}{
\code{numeric} Length of the simulated time series.}
%
\item{iN}{
\code{numeric} Cross sectional dimension (Only iN<5 supported for now).}
%
\item{vKappa}{
\code{numeric} vector of unconditional level for the reparametrised vector of parameters.}
%

\item{mA}{
\code{matrix} of coefficients of dimension iK x iK that premultiply the conditional score in the GAS updating recursion, see Details.}
%

\item{mB}{
\code{matrix} of autoregressive coefficients of dimension iK x iK, see Details.}
%
\item{Dist}{
\code{character} Label of the conditional distribution, see \link{DistInfo}.}
%
\item{ScalingType}{
\code{character} indicating the scaling mechanism for the conditional score. Only \code{ScalingType = "Identity"} is supported for multivariate distributions, see the function \link{DistInfo}.}
%
}
\details{
	All the information regarding the supported multivariate conditional distributions can be investigated using the \link{DistInfo} function. The model is specified as \deqn{y_{t}~p(y|\theta_{t})}, where \eqn{\theta_{t}} is the vector of parameters for the density \eqn{p(y|.)}. Note that, \eqn{\theta_{t}} includes also those parameters that are not time-varying. The GAS recursion for \eqn{\theta_{t}} is \deqn{\theta_{t} = h(\tilde{\theta}_{t})},\deqn{\tilde{\theta}_{t}=\kappa(I-mB) + mA*s_{t-1} + mB*\tilde{\theta}_{t-1}}, where \eqn{h(.)} is the mapping function (this is hard codeded in the package, see \link{MapParameters_multi}) and \eqn{\tilde{\theta}_{t}} is the vector of reparametrised parameters. The process is initialised at \eqn{\theta_{1}=\kappa}, where \eqn{\kappa} is the \code{vKappa} vector. The vector \eqn{s_{t}} is the scaled score of \eqn{p(y|.)} with respect to \eqn{\tilde{\theta}_{t}}.
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
\note{}
\examples{
\dontrun{
# Simulate from a GAS process with Multivariate Student-t conditional distribution, time-varying locations, scales, correlations and fixed shape parameter.
library(GAS)

set.seed(786)

iT     = 1000 # number of observations to simulate
iN     = 3     # trivariate series
Dist   = "mvt" # conditional Multivariate Studen-t distribution

# build unconditional vector of reparametrised parameters

vMu  = c(0.1,0.2,0.3)   # vector of location parameters (this is not transformed)
vPhi = c(1.0, 1.2, 0.3) # vector of scale parameters for the firs, second and third variables.

vRho = c(0.1,0.2,0.3)   # This represents vec(R), where R is the correlation matrix.
# Note that is up to the user to ensure that vec(R) implies a proper correlation matrix

vTheta = c(vMu, vPhi, vRho, 7) # vector of parameters such that the degrees of freedom are 7.

vKappa = UnmapParameters_multi(vTheta, Dist, iN, iK = length(vTheta))

mA     = matrix(0,length(vKappa),length(vKappa))

diag(mA) = c( 0,0,0 , 0.05,0.01,0.09  , 0.01,0.04,0.07,  0) # update scales and correlations, do not update locations and shape parameters

mB     = matrix(0,length(vKappa),length(vKappa))

diag(mB) = c( 0,0,0 , 0.7,0.7,0.5  , 0.94,0.97,0.92,  0) # update scales and correlations, do not update locations and shape parameters

Sim = MultiGASSim(iT,iN, vKappa, mA, mB, Dist, ScalingType = "Identity")

Sim
}
}
