\name{MultiGASFor}
\alias{MultiGASFor}
\title{
	Forecast with multivariate GAS models
}
\description{
	Forecast with multivariate GAS models. One step ahead prediction of the conditional density is available in closed form. Multistep ahead prediction are perfomed by simulation as detailed in Blasques et al. (2016).
}
\usage{
MultiGASFor(mGASFit, iH, Roll = F, mOut = NULL,iB = 1000,
                     vBands = c(0.1,0.15,0.85,0.9), bReturnsDraws = F)
}
\arguments{
\item{mGASFit}{ An object of the class \link{mGASFit} created using the function \link{MultiGASFit}}
%
\item{iH}{ \code{numeric} Forecast horizon. Ignored if \code{Roll = TRUE}}
\item{Roll}{ \code{boolean} Forecast should be made using a rolling procedure ? Note that if \code{Roll = TRUE}, then \code{mOut} has to be specified.}
\item{mOut}{ \code{matrix} of out of sample observation of dimension iN x iH for rolling forecast. iN refers to the cross sectional dimension.}
\item{iB}{ \code{numeric} Number of draws from the iH-step ahead distribution if \code{Roll = TRUE}.}
\item{vBands}{ \code{numeric} Vector of probabilities representing the confidence band levels for multistep ahead parameters forecasts. Only if \code{Roll = TRUE}.}
\item{bReturnsDraws}{ \code{boolean} Return the draws from the multistep ahead predictive distribution when \code{Roll = TRUE} ?}
%
}
\details{
}
\value{
An object of the class \link{mGASFor}
}
\references{
Blasques, F., Koopman, S. J., Łasak, K., & Lucas, A. (2016). In-sample confidence bands and out-of-sample forecast bands for time-varying parameters in observation-driven models. International Journal of Forecasting, 32(3), 875-887.
}
\author{Leopoldo Catania}
\note{}
\examples{
\dontrun{
# Specify an univariate GAS model with Student-t conditional distribution and time-varying location, scale and shape parameter
# Stock returns forecast

data("StockIndex")

mY = StockIndex[1:2,]

## Specification mvt
GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity", GASPar = list(location = F, scale = T, correlation = T, shape = F))

# Perform H-step ahead forecast with confidence bands

# estimation
Fit = MultiGASFit(GASSpec, mY)

#forecast

Forecast  = MultiGASFor(Fit, iH = 50)

Forecast

# Perform 1-Step ahead rolling forecast

InSampleData  = mY[,1:1000]
OutSampleData = mY[,1001:2404]

# estimation
Fit = MultiGASFit(GASSpec, InSampleData)

Forecast  = MultiGASFor(Fit, Roll = TRUE, mOut = OutSampleData)

Forecast
}
}
