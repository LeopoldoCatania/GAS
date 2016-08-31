\name{UniGASFit}
\alias{UniGASFit}
\title{
	Estimate univariate GAS models
}
\description{
	Estimate univariate GAS models by Maximum Likelihood.
}
\usage{
UniGASFit(GASSpec, data)
}
\arguments{
\item{GASSpec}{
An object of the class \link{uGASSpec} created using the function \link{UniGASSpec}}.
%
\item{data}{ \code{numeric} vector of length Tx1 containing the time series of observations. It can also be an object of the class \code{ts}, \code{xts} or \code{zoo}.}
%
}
\details{
Maximum Likelihood estimation of GAS models is an on-going research topic.
General results are reported by Blasques et al. (2014b), Blasques et al. (2014a) and Harvey
(2013), while results for specific models have been derived by Blasques et al. (2014c) and
Andres (2014).\cr
%

The optimizer used is the General Nonlinear Augmented Lagrange Multiplier method of Ye (1988) available in
the R package Rsolnp (Ghalanos and Theussl 2016).\cr
%

Starting values for the optimizer are chosen in the following way: (i) estimate the static
version of the model (i.e., with A = 0 and B = 0) and set the initial value of the intercept parameter accordingly, and (ii) perform a grid search for the coefficients contained in A and B. Further technical details are presented in Section 3.2 of Ardia et. al. (2016).\cr
%

The function prints some information during the estimation procedure.\cr
}
\value{
An object of the class \link{uGASFit}
}
\references{
Ardia D, Boudt K and Catania L (2016). "Generalized Autoregressive Score Models in R: The GAS Package." Available at SSRN: \url{http://ssrn.com/abstract=2825380}.\cr
%

Blasques F, Koopman SJ, Lucas A (2014a). "Maximum Likelihood Estimation for Correctly Specified Generalized Autoregressive Score Models: Feedback Effects, Contraction Conditions and Asymptotic Properties. techreport TI 14-074/III, Tinbergen Institute. URL \url{http://www.tinbergen.nl/discussionpaper/?paper=2332}.

%
Blasques F, Koopman SJ, Lucas A (2014b). "Maximum Likelihood Estimation for Generalized Autoregressive Score Models." techreport TI 2014-029/III, Tinbergen Institute. URL \url{http://www.tinbergen.nl/discussionpaper/?paper=2286}.\cr
%

Blasques F, Koopman SJ, Lucas A, Schaumburg J (2014c). "Spillover Dynamics for Systemic Risk Measurement using Spatial Financial Time Series Models." techreport TI 2014-103/III, Tinbergen Institute. URL \url{http://www.tinbergen.nl/discussionpaper/?paper=2369}.\cr

Creal D, Koopman SJ, Lucas A (2013). "Generalized Autoregressive Score Models with Applications." Journal of Applied Econometrics, 28(5), 777-795. \doi{10.1002/jae.1279}.\cr
%

Ghalanos A, Theussl S (2016). \code{Rsolnp}: General Non-Linear Optimization using Augmented Lagrange Multiplier Method. R package version 1.16, URL \url{https://cran.r-project.org/package=Rsolnp}.\cr
%

Harvey AC (2013). Dynamic Models for Volatility and Heavy Tails: With Applications to Financial and Economic Time Series. Cambridge University Press.\cr
%

Ye Y (1988). Interior Algorithms for Linear, Quadratic, and Linearly Constrained Convex Programming. Ph.D. thesis, Stanford University.
}
\author{Leopoldo Catania}
\examples{
\dontrun{
# Specify an univariate GAS model with Student-t
# conditional distribution and time-varying scale.
library("GAS")

data("sp500ret")

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity",
                     GASPar = list(location = FALSE, scale = TRUE,
                                   shape = FALSE))

Fit = UniGASFit(GASSpec, sp500ret)

Fit

}
}
