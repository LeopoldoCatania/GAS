\name{MultiGASFor}
\alias{MultiGASFor}
\title{
	Forecast with multivariate GAS models
}
\description{
	Forecast with multivariate GAS models. One step ahead prediction of the conditional density is available in closed form. Multistep ahead prediction are performed by simulation as detailed in Blasques et al. (2016).
}
\usage{
MultiGASFor(mGASFit, H, Roll = F, out = NULL, B = 10000,
                     Bands = c(0.1, 0.15, 0.85, 0.9), ReturnDraws = FALSE)
}
\arguments{
\item{mGASFit}{ An object of the class \link{mGASFit} created using the function \link{MultiGASFit}}
%
\item{H}{ \code{numeric} Forecast horizon. Ignored if \code{Roll = TRUE}}
\item{Roll}{ \code{boolean} Forecast should be made using a rolling procedure ? Note that if \code{Roll = TRUE}, then \code{out} has to be specified.}
\item{out}{ \code{matrix} of out of sample observation of dimension H x N for rolling forecast. N refers to the cross sectional dimension.}
\item{B}{ \code{numeric} Number of draws from the iH-step ahead distribution if \code{Roll = TRUE}.}
\item{Bands}{ \code{numeric} Vector of probabilities representing the confidence band levels for multistep ahead parameters forecasts. Only if \code{Roll = TRUE}.}
\item{ReturnDraws}{ \code{boolean} Return the draws from the multistep ahead predictive distribution when \code{Roll = TRUE} ?}
%
}
\value{
An object of the class \link{mGASFor}
}
\references{
Blasques, F., Koopman, S. J., Lasak, K., & Lucas, A. (2016). In-sample confidence bands and out-of-sample forecast bands for time-varying parameters in observation-driven models. International Journal of Forecasting, 32(3), 875-887.
}
\author{Leopoldo Catania}
\examples{
# Specify a GAS model with multivatiate  Student-t conditional
# distribution and time-varying scales and correlations.

# Stock returns forecast

data("StockIndices")

mY = StockIndices[, 1:2]

## Specification mvt
GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity",
                       GASPar = list(location = FALSE, scale = TRUE,
                                     correlation = TRUE, shape = FALSE))

# Perform H-step ahead forecast with confidence bands

# estimation
Fit = MultiGASFit(GASSpec, mY)

#forecast

Forecast  = MultiGASFor(Fit, H = 50)

Forecast

# Perform 1-Step ahead rolling forecast

InSampleData  = mY[1:1000, ]
OutSampleData = mY[1001:2404, ]

# estimation
Fit = MultiGASFit(GASSpec, InSampleData)

Forecast  = MultiGASFor(Fit, Roll = TRUE, out = OutSampleData)

Forecast
}
