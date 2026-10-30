\name{DistInfo}
\alias{DistInfo}
\alias{DistLabels}
\alias{DistName}
\alias{DistType}
\alias{DistParameters}
\alias{DistScalingType}
\alias{NumberParameters}
\title{
	Information for the supported distributions
}
\description{
	Print the information regarding distributions supported in the GAS package.
}
\usage{
DistInfo(DistLabel = NULL, N = 2, FULL = TRUE)
}
\arguments{
\item{DistLabel}{
\code{character} indicating the label of the conditional distribution. if \code{DistLabel = NULL} all the
supported distributions are printed. Default \code{DistLabel = NULL}}. Run \code{DistLabels()} to see
the labels of the currently implemented distributions.
%
\item{N}{
\code{numeric} Indicating the number of asset if \code{DistLabel} link to a multivariate distribution. Default \code{N = 2}.}
%
\item{FULL}{
\code{logical} If \code{TRUE} the function prints all the the information. Default \code{FULL = TRUE}}
%
}
\details{
	The information are printed in the console.
}
\value{
\code{DistInfo} prints distribution information to the console and
returns \code{NULL} invisibly.

The other functions documented on this page return:
\describe{
  \item{\code{DistLabels}}{
    A character vector of supported distribution labels.
  }
  \item{\code{DistName}}{
    A character string containing the distribution name.
  }
  \item{\code{DistType}}{
    A character string: \code{"univariate"} or \code{"multivariate"}.
  }
  \item{\code{DistParameters}}{
    A character vector naming the distribution's parameter groups.
  }
  \item{\code{DistScalingType}}{
    A character string listing the supported score scaling methods,
    separated by commas.
  }
  \item{\code{NumberParameters}}{
    A numeric scalar giving the number of distribution parameters,
    accounting for the dimension of a multivariate distribution.
  }
}
}
\author{Leopoldo Catania}
\examples{
# show the information of all the supported distributions
library("GAS")

DistInfo()
}
