\name{UniGASFor}
\alias{UniGASFor}
\title{
	Forecast with univariate GAS models
}
\description{
	Forecast with univariate GAS models. One step ahead prediction of the conditional density is available in closed form. Multistep ahead prediction are perfomed by simulation as detailed in Blasques et al. (2016).
}
\usage{
UniGASFor(uGASFit, iH, Roll = F, vOut = NULL,iB = 10000,
                     vBands = c(0.1,0.15,0.85,0.9), bReturnsDraws = F)
}
\arguments{
\item{uGASFit}{ An object of the class \link{uGASFit} created using the function \link{UniGASFit}}
%
\item{iH}{ \code{numeric} Forecast horizon. Ignored if \code{Roll = TRUE}}
\item{Roll}{ \code{boolean} Forecast should be made using a rolling procedure ? Note that if \code{Roll = TRUE}, then \code{vOut} has to be specified.}
\item{vOut}{ \code{numeric} Vector of out of sample observation for rolling forecast}
\item{iB}{ \code{numeric} Number of draws from the iH-step ahead distribution if \code{Roll = TRUE}.}
\item{vBands}{ \code{numeric} Vector of probabilities representing the confidence band levels for multistep ahead parameters forecasts. Only if \code{Roll = TRUE}.}
\item{bReturnsDraws}{ \code{boolean} Return the draws from the multistep ahead predictive distribution when \code{Roll = TRUE} ?}
%
}
\details{
}
\value{
An object of the class \link{uGASFor}
}
\references{
Blasques, F., Koopman, S. J., Łasak, K., & Lucas, A. (2016). In-sample confidence bands and out-of-sample forecast bands for time-varying parameters in observation-driven models. International Journal of Forecasting, 32(3), 875-887.
}
\author{Leopoldo Catania}
\note{}
\examples{
\dontrun{
# Specify an univariate GAS model with Student-t conditional distribution and time-varying location, scale and shape parameter
# Inflation Forecast

data("cpichg")
help(cpichg)

GASSpec   = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = T, scale = T, shape = F))

# Perform H-step ahead forecast with confidence bands

Fit       = UniGASFit(GASSpec,cpichg)
Forecast  = UniGASFor(Fit, iH = 12)

Forecast

# Perform 1-Step ahead rolling forecast

InsampleData  = cpichg[1:250]
OutSampleData = cpichg[251:276]

Fit           = UniGASFit(GASSpec,InsampleData)

Forecast  = UniGASFor(Fit, Roll = TRUE, vOut = OutSampleData)

Forecast


}
}
