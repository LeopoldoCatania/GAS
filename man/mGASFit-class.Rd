\name{mGASFit}
\Rdversion{0.1}
\docType{class}
\alias{mGASFit}
\alias{mGASFit-class}
\alias{show,mGASFit-method}
\alias{plot,mGASFit,missing-method}
\alias{coef,mGASFit-method}
\alias{getFilteredParameters,mGASFit-method}
\alias{getMoments,mGASFit-method}
\alias{getObs,mGASFit-method}
\title{class: Class for the Multivariate GAS fitted object}
\description{
Class for the Multivariate GAS fitted object.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
  \item{\code{Data}:}{Object of class \code{"list"}. Contains the user's data.}
  %
\item{\code{Estimates}:}{Object of class \code{"list"}. Contains: \code{lParList} list of estimated parameters, \code{optimiser}  object delivered from the optimization function, \code{StaticFit} ML estimates for the constant model, \code{Inference} inferential results for the estimated parameters.}
%
%
 \item{\code{GASDyn}:}{Object of class \code{"list"}. Contains: The series of filtered dynamic (\code{GASDyn$mTheta}) for the time--varying parameters, The series of scaled scores (\code{GASDyn$mInnovation}), the series of unrestricted filtered parameters (\code{GASDyn$mTheta_tilde}), The series of log densities (\code{GASDyn$vLLK}), the log likelihood evaluated at its optimum value (\code{GASDyn$dLLK})}
 %
    \item{\code{ModelInfo}:}{Object of class \code{"list"}. Contains information about the GAS specification:
    \itemize{
\item Spec : object of the class \link{uGASSpec} containing the GAS specification.
\item iT : \code{numeric} number of observation.
\item elapsedTime: \code{numeric} elapsed Time in seconds.
}
}
%
  }
  }
\section{Methods}{
\itemize{
\item show \code{signature(object = "mGASFit")}: Show summary.
\item plot \code{ signature(x='mGASFit',y='missing')}: Plot filtered dynamic and other estimated quantities.
\item getFilteredParameters \code{signature(object = "mGASFit")}: Extract filtered parameters
\item getObs \code{signature(object = "mGASFit")}: Extract original observations
\item coef \code{signature(object = "mGASFit")}: Extract estimated coefficients
\item getMoments \code{signature(object = "mGASFit")}: Extract conditional moments
}
}
\author{ Leopoldo Catania}
\keyword{classes}
