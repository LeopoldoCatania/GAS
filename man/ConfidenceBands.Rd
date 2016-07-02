\name{ConfidenceBands}
\alias{ConfidenceBands}
\title{
	Build confidece bands for the filtered parameters
}
\description{
	Build confidece bands for the filtered parameters sampling the coefficients from the asymptotic distribution as in Blasques et al. (2016).
}
\usage{
ConfidenceBands(object, iB = 999, probs = c(0.01,0.1,0.9,0.99))
}
\arguments{
\item{object}{
An object of the class \link{uGASFit} or \link{mGASFit}}
%
\item{iB}{
\code{numeric} Number of draws from the asymptotic distributions.}
%
\item{probs}{
\code{numeric} Quantiles to returns.}
}
\details{
	This function implements the "In-Sample Simulation-Based Bands" Sec 3.3 of Blasques et al. (2016).
}
\value{An object of the class \code{array} of dimension iT+1 x iB x iK, where iT is the length of the time series, iK is the number of parameters and iB the number of draws. The first slice reports the estimated filtered parameters. Also the one step ahead prediction is reported, this why iT + 1.}
\author{Leopoldo Catania}
\references{
Blasques, F., Koopman, S. J., Łasak, K., & Lucas, A. (2016). In-sample confidence bands and out-of-sample forecast bands for time-varying parameters in observation-driven models. International Journal of Forecasting, 32(3), 875-887.
}
\examples{
\dontrun{
# show the information of all the supported distributions
library(GAS)

data("cpichg")

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = T, scale = T, shape = F))

Fit = UniGASFit(GASSpec,cpichg)

Bands = ConfidenceBands(Fit)

}
}
