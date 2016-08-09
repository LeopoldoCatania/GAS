\name{uGASSpec}
\Rdversion{0.1}
\docType{class}
\alias{uGASSpec}
\alias{uGASSpec-class}
\alias{show,uGASSpec-method}
\title{class: Class for the Univariate GAS model specification}
\description{
Class for the Univariate GAS model specification.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{Spec}:}{Object of class \code{"list"}. Contains information about the univariate GAS specification:
    \itemize{
\item Dist: \code{character} containing the conditional distribution assumption.
\item ScalingType: \code{character} indicating the scaling mechanism for the conditional score.
\item iK : \code{numeric} representing the number of (possibly) time--varying parameters implied by the distributional assumption.
\item GASPar \code{list} with elements: location, scale, skewness, shape, shape2.
}
}
    }
  }
\section{Methods}{
\itemize{
\item show \code{signature(object = "uGASSpec")}: Show summary.
}
}
\author{Leopoldo Catania}
\keyword{classes}
