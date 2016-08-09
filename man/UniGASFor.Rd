\name{UniGASFor}
\alias{UniGASFor}
\title{
	Forecast with univariate GAS models
}
\description{
	Forecast with univariate GAS models. One step ahead prediction of the conditional density is available in closed form. Multistep ahead prediction are performed by simulation as detailed in Blasques et al. (2016).
}
\usage{
UniGASFor(uGASFit, H, Roll = F, out = NULL,B = 1000,
                     Bands = c(0.1,0.15,0.85,0.9), ReturnDraws = FALSE)
}
\arguments{
\item{uGASFit}{ An object of the class \link{uGASFit} created using the function \link{UniGASFit}}
%
\item{H}{ \code{numeric} Forecast horizon. Ignored if \code{Roll = TRUE}}
\item{Roll}{ \code{boolean} Forecast should be made using a rolling procedure ? Note that if \code{Roll = TRUE}, then \code{vOut} has to be specified.}
\item{out}{ \code{numeric} Vector of out of sample observation for rolling forecast}
\item{B}{ \code{numeric} Number of draws from the H-step ahead distribution if \code{Roll = TRUE}.}
\item{Bands}{ \code{numeric} Vector of probabilities representing the confidence band levels for multistep ahead parameters forecasts. Only if \code{Roll = TRUE}.}
\item{ReturnDraws}{ \code{boolean} Return the draws from the multistep ahead predictive distribution when \code{Roll = TRUE} ?}
%
}
\value{
An object of the class \link{uGASFor}
}
\references{
Blasques, F., Koopman, S. J., Lasak, K., & Lucas, A. (2016). In-sample confidence bands and out-of-sample forecast bands for time-varying parameters in observation-driven models. International Journal of Forecasting, 32(3), 875-887.
}
\author{Leopoldo Catania}
\examples{
# Specify an univariate GAS model with Student-t
# conditional distribution and time-varying location, scale and shape parameter

# Inflation Forecast

data("cpichg")

GASSpec   = UniGASSpec(Dist = "std", ScalingType = "Identity",
                       GASPar = list(location = TRUE, scale = TRUE, shape = FALSE))

# Perform H-step ahead forecast with confidence bands

Fit       = UniGASFit(GASSpec,cpichg)
Forecast  = UniGASFor(Fit, H = 12)

Forecast

# Perform 1-Step ahead rolling forecast

InsampleData  = cpichg[1:250]
OutSampleData = cpichg[251:276]

Fit           = UniGASFit(GASSpec, InsampleData)

Forecast  = UniGASFor(Fit, Roll = TRUE, out = OutSampleData)

Forecast

}
