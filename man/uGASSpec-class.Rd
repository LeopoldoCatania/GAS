\name{uGASSpec}
\Rdversion{0.1}
\docType{class}
\alias{uGASSpec}
\alias{uGASSpec-class}
\alias{show,uGASSpec-method}
\title{Class for the univariate GAS model specification}
\description{
Class for the univariate GAS model specification.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{Spec}:}{Object of class \code{list}. Contains information about the univariate GAS specification:
    \itemize{
\item \code{Dist}: \code{character} containing the conditional distribution assumption.
\item \code{ScalingType}: \code{character} indicating the scaling mechanism for the conditional score.
\item \code{iK}: \code{numeric} representing the number of (possibly) time-varying parameters implied by the distributional assumption.
\item \code{GASPar} \code{list} with elements: \code{location}, \code{scale}, \code{skewness}, \code{shape}, \code{shape2}.
}
}
    }
  }
\section{Methods}{
\itemize{
\item \code{show} \code{signature(object = 'uGASSpec')}: Show summary.
}
}
\author{Leopoldo Catania}
\keyword{classes}
