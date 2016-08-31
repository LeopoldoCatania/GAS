\name{mGASFor}
\Rdversion{0.1}
\docType{class}
\alias{mGASFor}
\alias{mGASFor-class}
\alias{LogScore,mGASFor-method}
\alias{getForecast,mGASFor-method}
\alias{getMoments,mGASFor-method}
\alias{plot,mGASFor,missing-method}
\alias{show,mGASFor-method}
\title{Class for the Multivariate GAS Forecast object}
\description{
Class for the multivariate GAS forecast object.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{Forecast}:}{Object of class \code{list}. Contains forecasts:
    \itemize{
\item \code{PointForecast}: \code{matrix} with parameters forecasts.
\item \code{Moments}: \code{list} with centered moments forecasts. The first element contains a \code{matrix} with the predicted conditional means. The second element contains an \code{array} with the predicted conditional covariances.
\item \code{vLS}: \code{numeric} Log Score (Predictive Log Likelihood).
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
\item \code{show} \code{signature(object = "uGASFor")}: Show summary.
\item \code{plot} \code{ signature(x='uGASFor',y='missing')}: Plot forecasted quantities.
\item \code{getForecast} \code{signature(object = "uGASFor")}: Extract parameters forecast.
\item \code{getObs} \code{signature(object = "uGASFor")}: Extract original observations.
\item \code{getMoments} \code{signature(object = "uGASFor")}: Extract moments forecasts.
\item \code{LogScore} \code{signature(object = "uGASFor")}: Extract Log Scores.
}
}
\author{ Leopoldo Catania}
\keyword{classes}
