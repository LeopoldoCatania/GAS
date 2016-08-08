\name{uGASRoll}
\Rdversion{0.1}
\docType{class}
\alias{uGASRoll}
\alias{uGASRoll-class}
\alias{LogScore,uGASRoll-method}
\alias{getForecast,uGASRoll-method}
\alias{getMoments,uGASRoll-method}
\alias{getObs,uGASRoll-method}
\alias{pit,uGASRoll-method}
\alias{plot,uGASRoll,missing-method}
\alias{quantile,uGASRoll-method}
\alias{show,uGASRoll-method}
\alias{LogScore}
\alias{getForecast}
\title{class: Class for the Univariate GAS Rolling object}
\description{
Class for the Univariate GAS Rolling object.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{Forecast}:}{Object of class \code{"list"}. Contains forecasts:
    \itemize{
\item PointForecast : \code{matrix} with parameters forecasts.
\item Moments : \code{matrix} with centered moments forecasts.
\item vLS : \code{numeric} Log Score (Predictive Log Likelihood)
\item vU  : \code{numeric} Out of sample Probability Integral Transformation
}}
%
\item{\code{Info}:}{\code{list} with forecast information.}
%
\item{\code{Data}:}{\code{list} with original data.}
%
\item{\code{Testing}:}{Statistical tests results.}
}
}
\section{Methods}{
\itemize{
\item show \code{signature(object = "uGASFor")}: Show summary.
\item plot \code{ signature(x='uGASFor',y='missing')}: Plot forecasted quantities.
\item getForecast \code{signature(object = "uGASFor")}: Extract parameters forecast.
\item getObs \code{signature(object = "uGASFor")}: Extract original observations
\item pit \code{signature(object = "uGASFor")}: Extract Probability Integral Transformation, only if \code{Roll = TRUE}
\item quantile \code{signature(object = "uGASFor")}: Extract quantile forecasts. It accepts the addional argument \code{probs} representing the vector of probabilities
\item getMoments \code{signature(object = "uGASFor")}: Extract moments forecasts.
\item LogScore \code{signature(object = "uGASFor")}: Extract Log Scores.
}
}
\author{ Leopoldo Catania}
\keyword{classes}
