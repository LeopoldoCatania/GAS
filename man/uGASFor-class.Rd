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
\title{Class for the univariate GAS forecast object}
\description{
Class for the univariate GAS forecast object.
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
}
}
 \item{\code{Bands}:}{\code{array} with confidence bands parameters forecasts. Available only if \code{Roll = TRUE}.}
%
\item{\code{Draws}:}{If \code{ReturnsDraws = TRUE} it is a iH x iB \code{matrix} of draws from the predictive distribution.}
\item{\code{Info}:}{\code{list} with forecast information.}
\item{\code{Data}:}{\code{list} with original data.}
}
}
\section{Methods}{
\itemize{
\item \code{show} \code{signature(object = 'uGASFor')}: Show summary.
\item \code{plot} \code{ signature(x = 'uGASFor', y = 'missing')}: Plot forecasted quantities.
\item \code{getForecast} \code{signature(object = 'uGASFor')}: Extract parameters forecast.
\item \code{getObs} \code{signature(object = 'uGASFor')}: Extract original observations.
\item \code{pit} \code{signature(object = 'uGASFor')}: Extract Probability Integral Transformation, only if \code{Roll = TRUE}.
\item \code{quantile} \code{signature(object = 'uGASFor')}: Extract quantile forecasts. It accepts the additional argument \code{probs} representing the vector of probabilities.
\item \code{getMoments} \code{signature(object = 'uGASFor')}: Extract moments forecasts.
\item \code{LogScore} \code{signature(object = 'uGASFor')}: Extract Log Scores.
}
}
\author{Leopoldo Catania}
\keyword{classes}
