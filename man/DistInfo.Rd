\name{DistInfo}
\alias{DistInfo}
\title{
	Information for the supported distributions
}
\description{
	Print the information regarding the one or more distribution supported in the GAS package.
}
\usage{
DistInfo(DistLabel = NULL, N = 2)
}
\arguments{
\item{DistLabel}{
\code{character} indicating the label of the conditional distribution. if \code{DistLabel = NULL} all the supported distributions are printed. Default \code{DistLabel = NULL}}. Run \code{DistLabels()} to see the labels of the currently implemented distributions.
%
\item{N}{
\code{numeric} indicating the number of asset if DistLabel link to a multivariate distribution. Default \code{iN = 2}.}
%
}
\details{
	The information are printed in the console.
}
\author{Leopoldo Catania}
\examples{
\dontrun{
# show the information of all the supported distributions
library(GAS)

DistInfo()
}
}
