\name{uGASSpec}
\Rdversion{0.1}
\docType{class}
\alias{uGASSpec}
\title{class: Class for the Univariate GAS model specification}
\description{
Class for the Univariate GAS model specification.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Extends}{

}
\section{Slots}{
  \itemize{
    \item \code{spec} Object of class \code{"list"}. Contains information abount the univariate GAS specification:
    \itemize{
\item Dist: \code{character} containing the conditional distribution assumption.
\item ScalingType: \code{character} indicating the scaling mechanism for the conditional score.
\item iK : \code{numeric} representing the number of (possibly) time--varying parameters implied by the distributional assumption.
\item GASPar \code{list} with elements: location, scale, skewness, shape, shape2.
}
    }
  }
\section{Methods}{
\itemize{
\item show \code{signature(object = "uGASSpec")}: Show summary.
}
}
\author{ Leopoldo Catania}
\note{ }
\keyword{classes}
