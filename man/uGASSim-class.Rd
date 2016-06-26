\name{uGASSim}
\Rdversion{0.1}
\docType{class}
\alias{uGASSim}
\title{class: Class for Univariate GAS Simulation}
\description{
Class for Univariate GAS model Simulation.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Extends}{

}
\section{Slots}{
  \itemize{
    \item \code{ModelInfo} Object of class \code{"list"}. Contains information abount the univariate GAS specification:
    \itemize{
\item iT: \code{numeric} Time length of simulated observations.
\item iK : \code{numeric} number of (possibly) time--varying parameters implied by the distributional assumption.
\item vKappa \code{numeric} vector of unconditional level for the reparametrised vector of parameters.
\item mA \code{matrix} of coefficients of dimension iK x iK that premultiply the conditional score in the GAS updating recursion.
\item mB \code{matrix} of autoregressive coefficients of dimension iK x iK.
\item Dist \code{character} label of the conditional distribution, see \link{DistInfo}
\item ScalingType \code{character} representing the scaling mechanism for the conditional score, see \link{DistInfo}
}
 \item \code{GASDyn} Object of class \code{"list"}. Contains: The series of simulated parameters (\code{GASDyn$mTheta}), The series of scaled scores (\code{GASDyn$mInnovation}), the series of unrestricted simulated parameters (\code{GASDyn$mTheta_tilde}), The series of log densities (\code{GASDyn$vLLK}), the log likelihood evaluated at its optimum value (\code{GASDyn$dLLK})
%
\item \code{Data} Object of class \code{"numeric"}. Vector of length iT of simulated data
  }
  }
\section{Methods}{
\itemize{
\item show \code{signature(object = "uGASSim")}: Show summary.
\item plot \code{ signature(x='uGASSim',y='missing')}: Plot simulated data and parameters.
\item getFilteredParameters \code{signature(object = "uGASSim")}: Extract simulated parameters
\item getObs \code{signature(object = "uGASSim")}: Extract simulated observations
\item coef \code{signature(object = "uGASSim")}: Extract delivered coefficients
}
}
\author{ Leopoldo Catania}
\note{ }
\keyword{classes}
