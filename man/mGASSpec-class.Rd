\name{mGASSpec}
\Rdversion{0.1}
\docType{class}
\alias{mGASSpec}
\alias{mGASSpec-class}
\alias{show,mGASSpec-method}
\title{Class for the Multivariate GAS model specification}
\description{
Class for the Multivariate GAS model specification.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{Spec}:}{Object of class \code{list}. Contains information abount the multivariate GAS specification:
    \itemize{
\item \code{Dist}: \code{character} Containing the conditional distribution assumption.
\item \code{ScalingType}: \code{character} indicating the scaling mechanism for the conditional score.
\item \code{GASPar}: \code{list} with elements: location, scale, correlation, shape.
\item \code{ScalarParameters}: \code{logical} indicates if the parameters of the locations, scales and 
correlations dynamic have to be scalars or a diagonal matrices.
}
}
    }
  }
\section{Methods}{
\itemize{
\item \code{show} \code{signature(object = 'mGASSpec')}: Show summary.
}
}
\author{Leopoldo Catania}
\keyword{classes}
