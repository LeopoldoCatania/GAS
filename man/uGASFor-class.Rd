\name{uGASFor}
\Rdversion{0.1}
\docType{class}
\alias{uGASFor}
\alias{uGASFor-class}
\alias{LogScore,uGASFor-method}
\alias{getForecast,uGASFor-method}
\alias{getMoments,uGASFor-method}
\alias{getObs,uGASFor-method}
\alias{pit,uGASFor-method}
\alias{quantile,uGASFor-method}
\alias{plot,uGASFor,missing-method}
\alias{show,uGASFor-method}
\title{class: Class for the Univariate GAS Forecast object}
\description{
Class for the Univariate GAS Forecast object.
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
}
}
 \item{\code{Bands}:}{\code{array} with confidence bands parameters forecasts. Available only if \code{Roll = TRUE}.}
%
\item{\code{Draws}:}{If \code{bReturnsDraws = TRUE} it is a iH x iB \code{matrix} of draws from the predictive distribution.}
\item{\code{Info}:}{\code{list} with forecast information.}
\item{\code{Data}:}{\code{list} with original data.}
}
}
\section{Methods}{
\itemize{
\item show \code{signature(object = "uGASFor")}: Show summary.
\item plot \code{ signature(x='uGASFor',y='missing')}: Plot forecasted quantities.
\item getForecast \code{signature(object = "uGASFor")}: Extract parameters forecast.
\item getObs \code{signature(object = "uGASFor")}: Extract original observations
\item pit \code{signature(object = "uGASFor")}: Extract Probability Integral Transformation, only if \code{Roll = TRUE}
\item quantile \code{signature(object = "uGASFor")}: Extract quantile forecasts. It accepts the additional argument \code{probs} representing the vector of probabilities
\item getMoments \code{signature(object = "uGASFor")}: Extract moments forecasts.
\item LogScore \code{signature(object = "uGASFor")}: Extract Log Scores.
}
}
\author{ Leopoldo Catania}
\keyword{classes}
