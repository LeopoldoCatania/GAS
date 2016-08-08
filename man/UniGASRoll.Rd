\name{UniGASRoll}
\alias{UniGASRoll}
\title{
	Rolling forecast with univariate GAS models
}
\description{
 One step ahead rolling forecasts with model re-estimation. The function also reports several quantity for backtesting for point and density forecasts.
}
\usage{
UniGASRoll(data, GASSpec, ForecastLength = 500, Nstart = NULL,
           RefitEvery = 23, RefitWindow = c("moving", "recursive"),
           cluster=NULL)
}
\arguments{
\item{data}{ \code{numeric} vector containing the time series of observations.}
%
\item{GASSpec}{ An object of the class \link{uGASSpec} created using the function \link{UniGASSpec}}
%
\item{ForecastLength}{ \code{numeric} Length of the out of sample}
\item{Nstart}{ \code{numeric} Period when perform the first forecast. Ignored if \code{ForecastLength} is supplied.}
\item{RefitEvery}{ \code{numeric} Number of periods before model coefficients re-estimation.}
\item{RefitWindow}{ \code{character} Type of window. If \code{RefitWindow = "recursive"} all the observations are used when the model is re-estimated. If \code{RefitWindow = "moving"} old observations are eliminated.}
\item{cluster}{ A \code{cluster} object created calling using the \code{paralell} package. If supplied parallel processing is used to speed up the computations.}
%
}
\value{
An object of the class \link{uGASRoll}
}
\author{Leopoldo Catania}
\examples{
\dontrun{
# Specify an univariate GAS model with Student-t
# conditional distribution and time-varying location, scale and shape parameter

# Inflation Forecast

data("cpichg")
help(cpichg)

GASSpec   = UniGASSpec(Dist = "std", ScalingType = "Identity",
                       GASPar = list(location = TRUE, scale = TRUE, shape = FALSE))

# Perform 1-step ahead rolling forecast with refit
library(parallel)

cluster = makeCluster(4)

Roll = UniGASRoll(cpichg, GASSpec, ForecastLength = 50,
                  RefitEvery = 2, RefitWindow = c("moving"), cluster = cluster)

Roll


}
}
