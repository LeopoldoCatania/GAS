\name{mGASFit}
\Rdversion{0.1}
\docType{class}
\alias{mGASFit}
\alias{mGASFit-class}
\alias{residuals,mGASFit-method}
\alias{show,mGASFit-method}
\alias{summary,mGASFit-method}
\alias{plot,mGASFit,missing-method}
\alias{coef,mGASFit-method}
\alias{getFilteredParameters,mGASFit-method}
\alias{getMoments,mGASFit-method}
\alias{getObs,mGASFit-method}
\alias{convergence}
\alias{convergence,mGASFit-method}
\title{Class for the Multivariate GAS fitted object}
\description{
Class for the multivariate GAS fitted object.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
  \item{\code{Data}:}{Object of class \code{list}. Contains the user's data.}
  %
\item{\code{Estimates}:}{Object of class \code{list}. Contains: \code{lParList} list of
estimated parameters, \code{optimiser} object delivered from the optimization function, \code{StaticFit} ML
estimates for the constant model, \code{Inference} inferential results for the estimated parameters.}
%
%
 \item{\code{GASDyn}:}{Object of class \code{list}. Contains: the series of filtered dynamic (\code{GASDyn$mTheta})
 for the time--varying parameters, the series of scaled scores (\code{GASDyn$mInnovation}), the series of
 unrestricted filtered parameters (\code{GASDyn$mTheta_tilde}), the series of log densities (\code{GASDyn$vLLK}), the
 log likelihood evaluated at its optimum value (\code{GASDyn$dLLK})}
 %
    \item{\code{ModelInfo}:}{Object of class \code{list}. Contains information about the GAS specification:
    \itemize{
\item \code{Spec} Object of the class \link{uGASSpec} containing the GAS specification.
\item \code{iT} \code{numeric} Number of observation.
\item \code{elapsedTime} \code{Numeric} elapsed time in seconds.
}
}
%
  }
  }
\section{Methods}{
\itemize{
\item \code{show} \code{signature(object = 'mGASFit')}: print object information.
\item \code{summary} \code{signature(object = 'mGASFit')}: Show summary.
\item plot \code{ signature(x='mGASFit', y='missing')}: Plot filtered dynamic and other estimated quantities.
\item getFilteredParameters \code{signature(object = "mGASFit")}: Extract filtered parameters.
\item getObs \code{signature(object = "mGASFit")}: Extract original observations.
\item \code{coef} \code{signature(object = 'uGASFit')}: Returns a named vector of estimated coefficients.
Also accepts the additional logical argument \code{do.list}. If \code{do.list = TRUE}, estimated coefficients
are organized in a list with arguments: \code{vKappa} the intercept vector, \code{mA} the A system matrix,
\code{mB} the B system matrix. By default, \code{do.list = FALSE}.
\item getMoments \code{signature(object = "mGASFit")}: Extract conditional moments.
\item \code{residuals} \code{signature(object = 'mGASFit')}: Extract the residuals.
Also accepts the additional logical argument \code{standardize}. If \code{standardize = TRUE},
residuals are standardized by cholesky of the filtered covariance matrix. By default \code{standardize = FALSE}.
\item \code{convergence} \code{signature(object = 'mGASFit')}: Extract convergence information.
}
}
\author{Leopoldo Catania}
\keyword{classes}
