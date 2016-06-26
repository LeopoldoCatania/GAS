\name{mGASSpec}
\Rdversion{0.1}
\docType{class}
\alias{mGASSpec}
\title{class: Class for the Multivariate GAS model specification}
\description{
Class for the Multivariate GAS model specification.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Extends}{

}
\section{Slots}{
  \itemize{
    \item \code{spec} Object of class \code{"list"}. Contains information abount the multivariate GAS specification:
    \itemize{
\item Dist: \code{character} containing the conditional distribution assumption.
\item ScalingType: \code{character} indicating the scaling mechanism for the conditional score.
\item GASPar \code{list} with elements: location, scale, correlation, shape.
}
    }
  }
\section{Methods}{
\itemize{
\item show \code{signature(object = "mGASSpec")}: Show summary.
}
}
\author{ Leopoldo Catania}
\note{ }
\keyword{classes}
