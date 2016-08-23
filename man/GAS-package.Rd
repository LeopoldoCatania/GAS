\name{GAS-package}
\alias{GAS-package}
\alias{GAS}
\docType{package}
\title{
Generalized Autoregressive Score models in R
}
\description{
The GAS package allows to simulate, estimate and forecast using univariate and multivariate Generalized Autoregressive Score (GAS) models (also known as Dynamic Conditional Score (DCS) models), see e.g. Creal et. al. (2013) and Harvey (2013).
}
\details{
\tabular{ll}{
Package: \tab GAS\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2016-08-23\cr
License: \tab GPL (>= 2)\cr
}
The authors acknowledge Google for financial support via the Google Summer of Code 2016 project "GAS"; see \url{https://summerofcode.withgoogle.com/projects/#4717600793690112}.\cr
%

Known bugs and limitations:
\itemize{
\item The Multivariate GAS model for N>4 does not report the exact update for the correlation parameters since the jacobian of the hyperspherical coordinate transformation needs to be coded for the case N>4. The jacobian for N>4 is replaced by the identity matrix.
}
}
\author{
Leopoldo Catania [aut,cre], Kris Boudt [ctb], David Ardia [ctb]

Maintainer: Leopoldo Catania <leopoldo.catania@uniroma2.it>
}
\references{
Creal, D., Koopman, S. J., & Lucas, A. (2013). Generalized autoregressive score models with applications. Journal of Applied Econometrics, 28(5), 777-795.\cr
%

Harvey, A. C. (2013). Dynamic models for volatility and heavy tails: with applications to financial and economic time series (Vol. 52). Cambridge University Press.\cr
}
\keyword{ package }
