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
# Specify an GAS model with multivariate Student-t conditional distribution and time-varying scales and correlations

library(GAS)

load(StockIndices)

GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity", GASPar = list(scale = TRUE, correlation = TRUE))

Fit = MultiGASFit(StockIndices, GASSpec)

Fit

}
}
