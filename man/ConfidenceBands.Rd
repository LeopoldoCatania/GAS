\name{ConfidenceBands}
\alias{ConfidenceBands}
\title{
	Build confidence bands for the filtered parameters
}
\description{
	Build confidence bands for the filtered parameters sampling the coefficients 
	from the asymptotic distribution as in Blasques et al. (2016).
}
\usage{
ConfidenceBands(object, B = 10000, probs = c(0.01,0.1,0.9,0.99), ...)
}
\arguments{
\item{object}{
An object of the class \link{uGASFit} or \link{mGASFit}}
%
\item{B}{
\code{numeric} Number of draws from the asymptotic distributions.}
%
\item{probs}{
\code{numeric} Quantiles to returns.}
%
\item{...}{
Additional arguments.}
}
\details{
	This function implements the "In-Sample Simulation-Based Bands" Sec 3.3 of Blasques et al. (2016).
}
\value{An object of the class \code{array} of dimension T+1 x B x K, where T is the length of 
the time series, K is the number of parameters and B the number of draws. The first slice reports 
the estimated filtered parameters. Also the one step ahead prediction is reported, this why T + 1.}
\author{Leopoldo Catania}
\references{
Blasques F, Koopman SJ, Lasak K, and Lucas, A (2016). 
"In-sample Confidence Bands and Out-of-Sample Forecast Bands for Time-Varying Parameters in Observation-Driven Models." 
International Journal of Forecasting, 32(3), 875-887. 
\doi{10.1016/j.ijforecast.2016.04.002}.
}
\examples{
\dontrun{

# show the information of all the supported distributions
library("GAS")

data("cpichg")

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity",
                     GASPar = list(location = TRUE, scale = TRUE,
                                   shape = FALSE))

Fit = UniGASFit(GASSpec, cpichg)

Bands = ConfidenceBands(Fit)
}
}
