\name{MultiGASRoll}
\alias{MultiGASRoll}
\title{
	Rolling forecast with multivariate GAS models
}
\description{
 One step ahead rolling forecasts with model re-estimation. The function also reports several quantity for backtesting for point and density forecasts.
}
\usage{
MultiGASRoll(data, GASSpec, ForecastLength = 500, Nstart = NULL,
             RefitEvery = 23, RefitWindow = c("moving", "recursive"),
             cluster=NULL)
}
\arguments{
\item{data}{ \code{matrix} of dimension (T + ForecastLength) x N containing the time series of observations. }
%
\item{GASSpec}{ An object of the class \link{mGASSpec} created using the function \link{MultiGASSpec}}
%
\item{ForecastLength}{ \code{numeric} Length of the out of sample}
\item{Nstart}{ \code{numeric} Period when perform the first forecast. Ignored if \code{ForecastLength} is supplied.}
\item{RefitEvery}{ \code{numeric} Number of periods before model coefficients re-estimation.}
\item{RefitWindow}{ \code{character} Type of window. If \code{RefitWindow = "recursive"} all the observations are used when the model is re-estimated. If \code{RefitWindow = "moving"} old observations are eliminated.}
\item{cluster}{ A \code{cluster} object created calling using the \code{paralell} package. If supplied parallel processing is used to speed up the computations.}
%
}
\value{
An object of the class \link{mGASRoll}
}
\author{Leopoldo Catania}
\examples{
# Specify a GAS model with Multivariate Student-t conditional
# distribution and time-varying scale and correlation parameters

# stock returns Forecast

data("StockIndices")

mY = StockIndices[, 1:2]

## Specification mvt
GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity",
                       GASPar = list(location = FALSE, scale = TRUE,
                                     correlation = TRUE, shape = FALSE))

# Perform 1-step ahead rolling forecast with refit
library(parallel)

Roll    = MultiGASRoll(mY, GASSpec, ForecastLength = 250,
                       RefitEvery = 100, RefitWindow = c("moving"))

Roll

}
