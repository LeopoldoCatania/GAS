\name{uGASFit}
\Rdversion{0.1}
\docType{class}
\alias{uGASFit}
\alias{uGASFit-class}
\alias{residuals}
\alias{residuals,uGASFit-method}
\alias{coef,uGASFit-method}
\alias{getFilteredParameters,uGASFit-method}
\alias{getMoments,uGASFit-method}
\alias{getObs,uGASFit-method}
\alias{pit,uGASFit-method}
\alias{plot,uGASFit,missing-method}
\alias{quantile,uGASFit-method}
\alias{ES,uGASFit-method}
\alias{show,uGASFit-method}
\alias{summary,uGASFit-method}
\alias{getMoments}
\alias{getObs}
\alias{pit}
\alias{getFilteredParameters}
\alias{convergence,uGASFit-method}
\title{Class for the univariate GAS fitted object}
\description{
Class for the univariate GAS fitted object.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{ModelInfo}:}{Object of class \code{list}. Contains information about the GAS specification:
    \itemize{
\item \code{Spec}: object of the class \link{uGASSpec} containing the GAS specification.
\item \code{iT}: \code{numeric} number of observation.
\item \code{elapsedTime}: \code{numeric} elapsed Time in seconds.
}
}
 \item{\code{GASDyn}:}{Object of class \code{list}. Contains: the series of filtered dynamic (\code{GASDyn$mTheta})
 for the time--varying parameters, the series of scaled scores (\code{GASDyn$mInnovation}), the series of
 unrestricted filtered parameters (\code{GASDyn$mTheta_tilde}), the series of log
 densities (\code{GASDyn$vLLK}), the log likelihood evaluated at its optimum value (\code{GASDyn$dLLK})}.
%
\item{\code{Estimates}:}{Object of class \code{list}. Contains: \code{lParList} list of
estimated parameters, \code{optimiser} object delivered from the optimization
function, \code{StaticFit} ML estimates for the constant model, \code{Inference}
inferential results for the estimated parameters.}
%
\item{\code{Data}:}{The user's data.}
%
\item{\code{Testing}:}{Statistical tests results.}
}
}
\section{Methods}{
\itemize{
\item \code{show} \code{signature(object = 'uGASFit')}: print object information.
\item \code{summary} \code{signature(object = 'uGASFit')}: Show summary.
\item \code{plot} \code{ signature(x = 'uGASFit',y = 'missing')}: Plot filtered dynamic and other estimated quantities.
\item \code{getFilteredParameters} \code{signature(object = 'uGASFit')}: Extract filtered parameters.
\item \code{getObs} \code{signature(object = 'uGASFit')}: Extract original observations.
\item \code{coef} \code{signature(object = 'uGASFit')}: Returns a named vector of estimated coefficients.
Also accepts the additional logical argument \code{do.list}. If \code{do.list = TRUE}, estimated coefficients
are organized in a list with arguments: \code{vKappa} the intercept vector, \code{mA} the A system matrix,
\code{mB} the B system matrix. By default, \code{do.list = FALSE}.
\item \code{pit} \code{signature(object = 'uGASFit')}: Extract Probability Integral Transformation.
\item \code{getMoments} \code{signature(object = 'uGASFit')}: Extract conditional moments.
\item \code{residuals} \code{signature(object = 'uGASFit')}: Extract the residuals.
Also accepts the additional logical argument \code{standardize}. If \code{standardize = TRUE},
residuals are standardized by the filtered standard deviation. By default \code{standardize = FALSE}.
\item \code{convergence} \code{signature(object = 'uGASFit')}: Extract convergence information.
\item \code{quantile} \code{signature(object = 'uGASSim')}: Compute quantiles of the filtered estimated density at each point in time. It accepts the additional argument \code{probs} representing the vector of probabilities.
\item \code{ES} \code{signature(object = 'uGASSim')}: Compute Expected Shortfall of the filtered estimated density at each point in time. It accepts the additional argument \code{probs} representing the vector of probabilities.
}
}
\author{Leopoldo Catania}
\keyword{classes}
