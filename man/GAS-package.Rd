\name{GAS-package}
\alias{GAS-package}
\alias{GAS}
\docType{package}
\title{
Generalized Autoregressive Score models in R
}
\description{
The GAS package allows us to simulate, estimate and forecast using univariate and multivariate Generalized
Autoregressive Score (GAS) models (also known as Dynamic Conditional Score (DCS) models),
see e.g., Creal et. al. (2013) and Harvey (2013) and the website \url{www.gasmodel.com}.
A detailed implementation of the package functionalities are reported in Ardia et. al. (2016a,b).
}
\details{
\tabular{ll}{
Package: \tab GAS\cr
Type: \tab Package\cr
Version: \tab 0.2.0\cr
Date: \tab 2017-05-07\cr
License: \tab GPL (>= 2)\cr
}
The authors acknowledge Google for financial support via the Google Summer of
Code 2016 project "GAS"; see \url{https://summerofcode.withgoogle.com/archive/2016/projects/4537082387103744/}.\cr
%

Current limitations:
\itemize{
\item The multivariate GAS model for N>4 does not report the exact update for the
correlation parameters since the Jacobian of the hyperspherical coordinates transformation
needs to be coded for the case N>4. The Jacobian for N>4 is replaced by the identity matrix.
}
}
\author{
Leopoldo Catania [aut,cre], Kris Boudt [ctb], David Ardia [ctb]

Maintainer: Leopoldo Catania <leopoldo.catania@uniroma2.it>
}
\references{
Ardia D, Boudt K and Catania L (2016a).
"Generalized Autoregressive Score Models in R: The GAS Package."
\url{http://ssrn.com/abstract=2825380}.\cr

Ardia D, Boudt K and Catania L (2016b).
"Value-at-Risk Prediction in R with the GAS Package."
\url{https://ssrn.com/abstract=2871444}.\cr

Creal D, Koopman SJ, Lucas A (2013).
"Generalized Autoregressive Score Models with Applications."
Journal of Applied Econometrics, 28(5), 777-795.
\doi{10.1002/jae.1279}.\cr

Harvey AC (2013).
Dynamic Models for Volatility and Heavy Tails: With Applications to Financial and Economic Time Series.
Cambridge University Press.
}
\keyword{ package }
