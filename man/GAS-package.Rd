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
see e.g., Creal et. al. (2013) and Harvey (2013). A detailed implementation of the package functionalities are reported in Ardia et. al. (2018, 2019).
}
\details{
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
\note{
By using \code{GAS}  you agree to the following rules:
\itemize{
\item You must cite Ardia et al. (2019) in working papers and published papers that use \code{GAS}. Use \code{citation("GAS")}.
\item You must place the following URL in a footnote to help others find \code{GAS}: \url{https://CRAN.R-project.org/package=GAS}.
\item You assume all risk for the use of \code{GAS}.
}
}
\author{
Leopoldo Catania [aut,cre], Kris Boudt [ctb], David Ardia [ctb]

Maintainer: Leopoldo Catania <leopoldo.catania@econ.au.dk>
}
\references{
Ardia D, Boudt K and Catania L (2018).
"Downside Risk Evaluation with the R Package GAS."
R Journal, 10(2), 410-421.
\doi{10.32614/RJ-2018-064}.\cr

Ardia D, Boudt K and Catania L (2019).
"Generalized Autoregressive Score Models in R: The GAS Package."
Journal of Statistical Software, 88(6), 1-28.
\doi{10.18637/jss.v088.i06}.\cr

Creal D, Koopman SJ, Lucas A (2013).
"Generalized Autoregressive Score Models with Applications."
Journal of Applied Econometrics, 28(5), 777-795.
\doi{10.1002/jae.1279}.\cr

Harvey AC (2013).
Dynamic Models for Volatility and Heavy Tails: With Applications to
Financial and Economic Time Series.
Cambridge University Press.
}
\keyword{ package }
