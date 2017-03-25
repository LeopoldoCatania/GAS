\name{mGASRoll}
\Rdversion{0.1}
\docType{class}
\alias{mGASRoll}
\alias{mGASRoll-class}
\alias{residuals,mGASRoll-method}
\alias{LogScore,mGASRoll-method}
\alias{getForecast,mGASRoll-method}
\alias{getMoments,mGASRoll-method}
\alias{plot,mGASRoll,missing-method}
\alias{show,mGASRoll-method}
\alias{coef,mGASRoll-method}
\title{Class for the Multivariate GAS Rolling object}
\description{
Class for the multivariate GAS rolling object.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{Forecast}:}{Object of class \code{list}. Contains forecasts:
    \itemize{
\item \code{PointForecast}: \code{matrix} with parameters forecasts.
\item \code{Moments}: \code{list} with centered moments forecasts. The first element
contains a \code{matrix} with the predicted conditional means. The second element contains
an \code{array} with the predicted conditional covariances.
\item \code{vLS}: \code{numeric} Log Score (Predictive Log Likelihood).
}
}
%
\item{\code{Info}:}{\code{list} with forecast information.}
\item{\code{Data}:}{\code{list} with original data.}
}
}
\section{Methods}{
\itemize{
\item \code{show} \code{signature(object = 'mGASRoll')}: Show summary.
\item \code{plot} \code{signature(x = 'mGASRoll',y = 'missing')}: Plot forecasted quantities.
\item \code{getForecast} \code{signature(object = 'mGASRoll')}: Extract parameters forecast.
\item \code{getObs} \code{signature(object = 'mGASRoll')}: Extract original observations.
\item \code{getMoments} \code{signature(object = 'mGASRoll')}: Extract moments forecasts.
\item \code{LogScore} \code{signature(object = 'mGASRoll')}: Extract Log Scores.
\item \code{residuals} \code{signature(object = 'mGASRoll')}: Extract the forecast errors.
Also accepts the additional logical argument \code{standardize}. If \code{standardize = TRUE},
forecast errors are standardized by cholesky of the forecast covariance matrix. By default \code{standardize = FALSE}.
\item \code{coef} \code{signature(object = 'mGASFit')}: Returns a matrix of estimated coefficients.
Each row of the matrix corresponds to a refit of the model during the forecast period according to the
\code{RefitEvery} argument provided in the \link{MultiGASRoll} function.
Also accepts the additional logical argument \code{do.list}. If \code{do.list = TRUE}, estimated coefficients
are organized in a list of lists according according to the \code{RefitEvery} argument provided
in the \link{MultiGASRoll} function. Each list is populated by three arguments:
\code{vKappa} the intercept vector, \code{mA} the A system matrix, \code{mB} the B system matrix.
By default, \code{do.list = FALSE}.
}
}
\author{Leopoldo Catania}
\keyword{classes}
