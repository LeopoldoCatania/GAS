\name{uGASSim}
\Rdversion{0.1}
\docType{class}
\alias{uGASSim}
\alias{uGASSim-class}
\alias{coef,uGASSim-method}
\alias{getFilteredParameters,uGASSim-method}
\alias{getMoments,uGASSim-method}
\alias{getObs,uGASSim-method}
\alias{plot,uGASSim,missing-method}
\alias{quantile,uGASSim-method}
\alias{show,uGASSim-method}
\title{Class for Univariate GAS Simulation}
\description{
Class for Univariate GAS model Simulation.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{ModelInfo}:}{Object of class \code{list}. Contains information about the univariate GAS specification:
    \itemize{
\item \code{iT} \code{numeric} Time length of simulated observations.
\item \code{iK} \code{numeric} Number of (possibly) time-varying parameters implied by the distributional assumption.
\item \code{vKappa} \code{numeric} Vector of unconditional level for the reparametrised vector of parameters.
\item \code{mA} \code{matrix} Of coefficients of dimension iK x iK that premultiply the conditional score in the GAS updating recursion.
\item \code{mB} \code{matrix} Of autoregressive coefficients of dimension iK x iK.
\item \code{Dist} \code{character} Label of the conditional distribution, see \link{DistInfo}
\item \code{ScalingType} \code{character} Representing the scaling mechanism for the conditional score, see \link{DistInfo}.
}
}
 \item{\code{GASDyn}:}{Object of class \code{list}. Contains: the series of simulated parameters (\code{GASDyn$mTheta}), the series of scaled scores (\code{GASDyn$mInnovation}), the series of unrestricted simulated parameters (\code{GASDyn$mTheta_tilde}), the series of log densities (\code{GASDyn$vLLK}), the log likelihood evaluated at its optimum value (\code{GASDyn$dLLK}).}
%
\item{\code{Data}:}{Object of class \code{numeric}. Vector of length \code{iT} of simulated data.}
  }
  }
\section{Methods}{
\itemize{
\item \code{show} \code{signature(object = 'uGASSim')}: Show summary.
\item \code{plot} \code{signature(x = 'uGASSim', y = 'missing')}: Plot simulated data and parameters.
\item \code{getFilteredParameters} \code{signature(object = 'uGASSim')}: Extract simulated parameters.
\item \code{getObs} \code{signature(object = 'uGASSim')}: Extract simulated observations.
\item \code{coef} \code{signature(object = 'uGASSim')}: Extract delivered coefficients.
}
}
\author{Leopoldo Catania}
\keyword{classes}
