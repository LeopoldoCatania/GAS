\name{MultiGASFit}
\alias{MultiGASFit}
\title{
	Estimate multivariate GAS models
}
\description{
	Estimate multivariate GAS models by Maximum Likelihood.
}
\usage{
MultiGASFit(GASSpec,data)
}
\arguments{
\item{GASSpec}{
An object of the class \link{mGASSpec} created using the function \link{MultiGASSpec}}
%
\item{data}{ \code{matrix} (or something coercible to that using \code{as.matrix()}) of dimension TxN containing the multivariate time series of observations. It can also be an object of the class \code{ts}, \code{xts} or \code{zoo}.}
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
\examples{
# Specify an GAS model with multivariate Student-t
# conditional distribution and time-varying scales and correlations

library(GAS)

data("StockIndices")

GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity",
                       GASPar = list(scale = TRUE, correlation = TRUE))

Fit = MultiGASFit(GASSpec, StockIndices)

Fit

}
