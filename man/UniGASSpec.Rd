\name{UniGASSpec}
\alias{UniGASSpec}
\title{
	Univariate GAS specification
}
\description{
	Specify the conditional distribution, scaling mechanism and time--varying parameters for univariate GAS models.
}
\usage{
UniGASSpec(Dist = "norm", ScalingType = "Identity", GASPar = list(location = FALSE, scale = TRUE,
                                                                  skewness = FALSE, shape = FALSE, shape2 = FALSE))
}
\arguments{
\item{Dist}{
\code{character} indicating the label of the conditional distribution. Available distribution can be displayed using the function \link{DistInfo}. Default value\code{Dist = "norm"}}
%
\item{ScalingType}{
\code{character} indicating the scaling mechanism for the conditional score. Possible choices are \code{"Identity"}, \code{"Inv"},
 \code{"Inv.Sqrt"}. Note that for some distribution only \code{ScalingType = "Identity"} is supported, see the function \link{DistInfo}. When  \code{ScalingType = "Inv.Sqrt"} the inverse of the cholesky decomposition of the information matrix is used. Default value \code{ScalingType = "Identity"}}
%
\item{GASPar}{
\code{list} containing information about which parameters of the conditional distribution have to be time-varying. \code{location = TRUE} referes to the location parameter, \code{scale = TRUE} refers to the scale parameter, \code{skewness = TRUE} refers to the parameter controlling the skewness, \code{shape = TRUE} refers to the shape parameter (e.g. the degree of freedom of the Student-t distribution), \code{shape2 = TRUE} refers to the second shape parameter. If the distribution specified in the \code{Dist} argument does not have, say, a shape parameter, the condition \code{shape = TRUE} or \code{shape = FALSE} is ignored.}
%
}
\details{
	All the information regarding the supported univariate conditional distributions can be investigated using the \link{DistInfo} function.
}
\value{
An object of the class \link{uGASSpec}
}
\references{
Creal, D., Koopman, S. J., & Lucas, A. (2013). Generalized autoregressive score models with applications. Journal of Applied Econometrics, 28(5), 777-795.\cr
%

Harvey, A. C. (2013). Dynamic models for volatility and heavy tails: with applications to financial and economic time series (Vol. 52). Cambridge University Press.\cr
}
\author{Leopoldo Catania}
\note{}
\examples{
\dontrun{
# Specify an univariate GAS model with Student-t conditional distribution and time-varying location, scale and shape parameter
library(GAS)

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = TRUE, scale = TRUE, shape = TRUE))

}
}
