\name{uGASRoll}
\Rdversion{0.1}
\docType{class}
\alias{uGASRoll}
\alias{uGASRoll-class}
\alias{residuals,uGASRoll-method}
\alias{LogScore,uGASRoll-method}
\alias{getForecast,uGASRoll-method}
\alias{getMoments,uGASRoll-method}
\alias{getObs,uGASRoll-method}
\alias{pit,uGASRoll-method}
\alias{plot,uGASRoll,missing-method}
\alias{quantile,uGASRoll-method}
\alias{ES,uGASRoll-method}
\alias{show,uGASRoll-method}
\alias{LogScore}
\alias{getForecast}
\alias{coef,uGASRoll-method}
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
\item \code{show} \code{signature(object = 'uGASRoll')}: Show summary.
\item \code{plot} \code{ signature(x = 'uGASRoll', y = 'missing')}: Plot forecasted quantities.
\item \code{getForecast} \code{signature(object = 'uGASRoll')}: Extract parameters forecast.
\item \code{getObs} \code{signature(object = 'uGASRoll')}: Extract original observations.
\item \code{pit} \code{signature(object = 'uGASRoll')}: Extract Probability Integral Transformation,
only if \code{Roll = TRUE}
\item \code{quantile} \code{signature(object = 'uGASRoll')}: Extract quantile forecasts.
It accepts the addional argument \code{probs} representing the vector of probabilities.
\item \code{ES} \code{signature(object = 'uGASRoll')}: Extract Expected Shortfall forecasts.
It accepts the addional argument \code{probs} representing the vector of probabilities.
\item \code{getMoments} \code{signature(object = 'uGASRoll')}: Extract moments forecasts.
\item \code{LogScore} \code{signature(object = 'uGASRoll')}: Extract Log Scores.
\item \code{residuals} \code{signature(object = 'uGASRoll')}: Extract the forecast errors.
Also accepts the additional logical argument \code{standardize}. If \code{standardize = TRUE},
forecast errors are standardized by the forecast standard deviation. By default \code{standardize = FALSE}.
\item \code{coef} \code{signature(object = 'uGASFit')}: Returns a matrix of estimated coefficients.
Each row of the matrix corresponds to a refit of the model during the forecast period according to the
\code{RefitEvery} argument provided in the \link{UniGASRoll} function.
Also accepts the additional logical argument \code{do.list}. If \code{do.list = TRUE}, estimated coefficients
are organized in a list of lists according according to the \code{RefitEvery} argument provided
in the \link{UniGASRoll} function. Each list is populated by three arguments:
\code{vKappa} the intercept vector, \code{mA} the A system matrix, \code{mB} the B system matrix.
By default, \code{do.list = FALSE}.
}
}
\author{Leopoldo Catania}
\keyword{classes}
