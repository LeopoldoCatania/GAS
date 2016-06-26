\name{UniGASFit}
\alias{UniGASFit}
\title{
	Estimate univariate GAS models
}
\description{
	Estimate univariate GAS models by Maximum Likelihood.
}
\usage{
UniGASFit(GASSpec,vY)
}
\arguments{
\item{GASSpec}{
An object of the class \link{uGASSpec} created using the function \link{UniGASSpec}}
%
\item{vY}{ \code{numeric} vector of length Tx1 containing the time series of observations.}
%
}
\details{
	The function prints some information during the estimation procedure.
}
\value{
An object of the class \link{uGASFit}
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
# Specify an univariate GAS model with Student-t conditional distribution and time-varying location, scale and shape parameter
library(GAS)

load(sp500ret)

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = T, scale = T, shape = T))

Fit = UniGASFit(sp500ret, GASSpec)

Fit

}
}
