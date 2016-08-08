\name{mGASSim}
\Rdversion{0.1}
\docType{class}
\alias{mGASSim}
\alias{mGASSim-class}
\alias{coef,mGASSim-method}
\alias{getMoments,mGASSim-method}
\alias{getObs,mGASSim-method}
\alias{plot,mGASSim,missing-method}
\alias{show,mGASSim-method}
\alias{getFilteredParameters,mGASSim-method}
\title{class: Class for Multivariate GAS Simulation}
\description{
Class for Multivariate GAS model Simulation.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{ModelInfo}:}{Object of class \code{"list"}. Contains information abount the multivariate GAS specification:
    \itemize{
\item iT: \code{numeric} Time length of simulated observations.
\item iN: \code{numeric} Cross sectional dimension.
\item iK : \code{numeric} number of (possibly) time--varying parameters implied by the distributional assumption.
\item vKappa \code{numeric} vector of unconditional level for the reparametrised vector of parameters.
\item mA \code{matrix} of coefficients of dimension iK x iK that premultiply the conditional score in the GAS updating recursion.
\item mB \code{matrix} of autoregressive coefficients of dimension iK x iK.
\item Dist \code{character} label of the conditional distribution, see \link{DistInfo}
\item ScalingType \code{character} representing the scaling mechanism for the conditional score, see \link{DistInfo}
}
}
%
\item{\code{GASDyn}:}{Object of class \code{"list"}. Contains: The series of simulated parameters (\code{GASDyn$mTheta}), The series of scaled scores (\code{GASDyn$mInnovation}), the series of unrestricted simulated parameters (\code{GASDyn$mTheta_tilde}), The series of log densities (\code{GASDyn$vLLK}), the log likelihood evaluated at its optimum value (\code{GASDyn$dLLK})}
%
\item{\code{Data}:}{Object of class \code{"matrix"}. Matrix of dimension iN x iT of simulated data}}
  }
\section{Methods}{
\itemize{
\item show \code{signature(object = "mGASSim")}: Show summary.
\item plot \code{ signature(x='mGASSim',y='missing')}: Plot simulated data and parameters.
\item getFilteredParameters \code{signature(object = "mGASSim")}: Extract simulated parameters
\item getObs \code{signature(object = "mGASSim")}: Extract simulated observations
\item coef \code{signature(object = "mGASSim")}: Extract delivered coefficients
}
}
\author{ Leopoldo Catania}
\keyword{classes}
