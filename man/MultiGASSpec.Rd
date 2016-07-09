\name{MultiGASSpec}
\alias{MultiGASSpec}
\title{
	Multivariate GAS specification
}
\description{
	Specify the conditional distribution, scaling mechanism and time-varying parameters for multivariate GAS models.
}
\usage{
MultiGASSpec(Dist = "mvnorm", ScalingType = "Identity", GASPar = list(location = FALSE, scale = TRUE,
                                                                      correlation = TRUE, shape = FALSE),
             ScalarParameters = TRUE)
}
\arguments{
\item{Dist}{
\code{character} indicating the label of the conditional distribution. Available distribution can be displayed using the function \link{DistInfo}. Default value\code{Dist = "mvnorm"}}
%
\item{ScalingType}{
\code{character} indicating the scaling mechanism for the conditional score. Only \code{ScalingType = "Identity"} is supported for multivariate distributions.}
%
\item{GASPar}{
\code{list} containing information about which parameters of the conditional distribution have to be time-varying. \code{location = TRUE} referes to the location parameters, \code{scale = TRUE} refers to the scale parameters, \code{shape = TRUE} refers to the shape parameter (e.g. the degree of freedom of the Multivariate Student-t distribution), \code{correlation = TRUE} refers to the correlations. If the distribution specified in the \code{Dist} argument does not have, say, a shape parameter, the condition \code{shape = TRUE} is ignored.}
%
\item{ScalarParameters}{
\code{logical} indicating if the parameters of the locations, scales and correlations dynamic have to be scalars or a diagonal matrices. By default \code{ScalarParameters = TRUE}}

}
\details{
	All the information regarding the supported multivariate conditional distributions can be investigated using the \link{DistInfo} function.
}
\value{
An object of the class \link{mGASSpec}
}
\references{
Creal, D., Koopman, S. J., & Lucas, A. (2012). A dynamic multivariate heavy-tailed model for time-varying volatilities and correlations. Journal of Business & Economic Statistics.
%

Creal, D., Koopman, S. J., & Lucas, A. (2013). Generalized autoregressive score models with applications. Journal of Applied Econometrics, 28(5), 777-795.\cr
%

Harvey, A. C. (2013). Dynamic models for volatility and heavy tails: with applications to financial and economic time series (Vol. 52). Cambridge University Press.\cr
}
\author{Leopoldo Catania}
\note{}
\examples{
\dontrun{
# Specify a GAS model with Multivariate Student-t conditional distribution and time-varying locations, scales and correlations  parameters but constant shape parameter.
library(GAS)

GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity", GASPar = list(location = TRUE, scale = TRUE,
                                                                             correlation = TRUE, shape = FALSE))

}
}
