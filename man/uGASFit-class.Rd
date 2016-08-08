\name{uGASFit}
\Rdversion{0.1}
\docType{class}
\alias{uGASFit}
\alias{uGASFit-class}
\alias{coef,uGASFit-method}
\alias{getFilteredParameters,uGASFit-method}
\alias{getMoments,uGASFit-method}
\alias{getObs,uGASFit-method}
\alias{pit,uGASFit-method}
\alias{plot,uGASFit,missing-method}
\alias{quantile,uGASFit-method}
\alias{show,uGASFit-method}
\alias{getMoments}
\alias{getObs}
\alias{pit}
\alias{getFilteredParameters}
\title{class: Class for the Univariate GAS fitted object}
\description{
Class for the Univariate GAS fitted object.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{ModelInfo}:}{Object of class \code{"list"}. Contains information abount the GAS specification:
    \itemize{
\item Spec : object of the class \link{uGASSpec} containing the GAS specification.
\item iT : \code{numeric} number of observation.
\item elapsedTime: \code{numeric} elapsed Time in seconds.
}
}
 \item{\code{GASDyn}:}{Object of class \code{"list"}. Contains: The series of filtered dynamic (\code{GASDyn$mTheta}) for the time--varying parameters, The series of scaled scores (\code{GASDyn$mInnovation}), the series of unrestricted filtered parameters (\code{GASDyn$mTheta_tilde}), The series of log densities (\code{GASDyn$vLLK}), the log likelihood evaluated at its optimum value (\code{GASDyn$dLLK})}
%
\item{\code{Estimates}:}{Object of class \code{"list"}. Contains: \code{lParList} list of estimated parameters, \code{optimiser}  object delivered from the optimisation function, \code{StaticFit} ML estimates for the constant model, \code{Inference} inferential results for the estimated parameters.}
%
\item{\code{Data}:}{The user's data.}
%
\item{\code{Testing}:}{Statistical tests results.}
}
}
\section{Methods}{
\itemize{
\item show \code{signature(object = "uGASFit")}: Show summary.
\item plot \code{ signature(x='uGASFit',y='missing')}: Plot filtered dynamic and other estimated quantities.
\item getFilteredParameters \code{signature(object = "uGASFit")}: Extract filtered parameters
\item getObs \code{signature(object = "uGASFit")}: Extract original observations
\item coef \code{signature(object = "uGASFit")}: Extract estimated coefficients
\item pit \code{signature(object = "uGASFit")}: Extract Probability Integral Transformation
\item getMoments \code{signature(object = "uGASFor")}: Extract conditional moments.
}
}
\author{ Leopoldo Catania}
\keyword{classes}
