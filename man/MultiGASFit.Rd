\name{MultiGASFit}
\alias{MultiGASFit}
\title{
	Estimate multivariate GAS models
}
\description{
	Estimate multivariate GAS models by Maximum Likelihood.
}
\usage{
MultiGASFit(GASSpec,mY)
}
\arguments{
\item{GASSpec}{
An object of the class \link{mGASSpec} created using the function \link{MultiGASSpec}}
%
\item{mY}{ \code{matrix} of dimension NxT containing the multivariate time series of observations.}
%
}
\details{
	The function prints some information during the estimation procedure.
}
\value{
An object of the class \link{mGASFit}
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
