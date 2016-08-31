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
\title{Class for the univariate GAS rolling object}
\description{
Class for the univariate GAS rolling object.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{Forecast}:}{Object of class \code{list}. Contains forecasts:
    \itemize{
\item \code{PointForecast}: \code{matrix} with parameters forecasts.
\item \code{Moments}: \code{matrix} with centered moments forecasts.
\item \code{vLS}: \code{numeric} Log Score (Predictive Log Likelihood).
\item \code{vU}: \code{numeric} Out-of-sample Probability Integral Transformation (PIT).
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
\item \code{show} \code{signature(object = 'uGASFor')}: Show summary.
\item \code{plot} \code{ signature(x = 'uGASFor', y = 'missing')}: Plot forecasted quantities.
\item \code{getForecast} \code{signature(object = 'uGASFor')}: Extract parameters forecast.
\item \code{getObs} \code{signature(object = 'uGASFor')}: Extract original observations.
\item \code{pit} \code{signature(object = 'uGASFor')}: Extract Probability Integral Transformation, only if \code{Roll = TRUE}
\item \code{quantile} \code{signature(object = 'uGASFor')}: Extract quantile forecasts. It accepts the addional argument \code{probs} representing the vector of probabilities.
\item \code{getMoments} \code{signature(object = 'uGASFor')}: Extract moments forecasts.
\item \code{LogScore} \code{signature(object = 'uGASFor')}: Extract Log Scores.
}
}
\author{Leopoldo Catania}
\keyword{classes}
