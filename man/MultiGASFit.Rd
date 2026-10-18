\name{MultiGASFit}
\alias{MultiGASFit}
\title{
	Estimate multivariate GAS models
}
\description{
	Estimate multivariate GAS models by Maximum Likelihood.
}
\usage{
MultiGASFit(GASSpec, data, fn.optimizer = fn.optim, Compute.SE = TRUE)
}
\arguments{
\item{GASSpec}{
An object of the class \link{mGASSpec} created using the function \link{MultiGASSpec}}
%
\item{data}{ \code{matrix} (or something coercible to that using \code{as.matrix()}) of
dimension TxN containing the multivariate time series of observations.
It can also be an object of the class \code{ts}, \code{xts} or \code{zoo}.}
%
\item{fn.optimizer}{\code{function}. This is a generic optimization function that can be provided by the user.
By default \code{fn.optimizer = fn.optim} where \link{fn.optim} is
a wrapper to the \link{optim} function. See Details for user defined optimization routines.}
%
\item{Compute.SE}{\code{logical}. Should asymptotic Standard Errors be computed? By default \code{Compute.SE = TRUE}}
}
\details{
For a fully static multivariate Gaussian model, the analytical sample-mean
and sample-covariance MLE is used, with covariance divisor T. A full-rank
sample covariance is required; the supplied optimizer is not called for
this special case.

Maximum Likelihood estimation of GAS models is an on-going research topic.
General results are reported by Blasques et al. (2014b), Blasques et al. (2014a) and Harvey
(2013), while results for specific models have been derived by Blasques et al. (2014c) and
Andres (2014).\cr
%

Starting values for the optimizer are chosen in the following way: (i) estimate the static
version of the model (i.e., with A = 0 and B = 0) and set the initial value of the
intercept parameter accordingly, and (ii) perform a grid search for the
coefficients contained in A and B. Further technical details are presented in Section 3.2 of Ardia et. al. (2016a).\cr
%

The user is free to employ his/her own optimization routine via the \code{fn.optimizer} argument. \code{fn.optimizer}
accepts a \code{function} object. The user provided optimizer has to satisfy strict requirements. The arguments of the
\code{fn.optimizer} are : i) \code{par0} a vector of starting values, ii) \code{data} the data provided, iii) \code{GASSpec}
an object of the class \link{mGASSpec}, and iv) \code{FUN} the likelihood function. The output of \code{fn.optimizer} has
to be a named list containing \code{pars}, the estimated numeric parameter
vector; \code{value}, the summed negative log likelihood; and
\code{convergence}, a numeric solver status with zero indicating success.
An optional \code{hessian} element is retained as optimizer output but is not
used for inference.

When \code{Compute.SE = TRUE}, \code{numDeriv::jacobian} differentiates the
vector of per-observation log likelihood contributions at the final estimate.
The complete filter, including initialization, is rerun at each perturbed
coefficient vector. If \eqn{S} is the resulting \eqn{T \times p} score matrix,
the outer-product information matrix is \eqn{S' S} and the covariance estimate
on the unrestricted optimizer scale is \eqn{(S' S)^{-1}}. There is no additional
sample-size divisor. The delta method transforms this covariance to the
reported coefficient scale. This is the OPG covariance for correctly specified
maximum likelihood, not a sandwich covariance under misspecification.

A missing convergence flag is recorded as unknown. A zero solver flag that
fails the numerical gradient check is replaced by code 52. The original
flag is retained in \code{Estimates$optimiser$solver_convergence}.
The maximized objective is checked against the starting objective.
Standard errors, test statistics and p-values are returned as \code{NA}
when convergence or OPG checks fail. No generalized inverse or automatic
positive-definiteness repair is used. The coefficient estimates remain available.
The score matrix, summed OPG matrix and its inverse are stored in
\code{Estimates$optimiser$scores}, \code{Estimates$optimiser$opg} and
\code{Estimates$optimiser$covariance}, respectively, when available.
The \code{covariance_method} field is \code{"OPG"}; validation results are in
\code{diagnostics}, including \code{information_ok} and \code{inference_ok}.
With \code{Compute.SE = FALSE}, these matrices are not computed and
\code{inference_ok} is false.


By default, the \code{optim} optimizer with \code{method = "BFGS"} is employed.

}
\value{
An object of the class \link{mGASFit}.
}
\references{
Ardia D, Boudt K and Catania L (2016a).
"Generalized Autoregressive Score Models in R: The GAS Package."
\url{https://www.ssrn.com/abstract=2825380}.\cr

Blasques F, Koopman SJ, Lucas A (2014a).
"Maximum Likelihood Estimation for Correctly Specified Generalized Autoregressive Score Models: Feedback Effects, Contraction Conditions and Asymptotic Properties." techreport TI 14-074/III, Tinbergen Institute.
\url{https://tinbergen.nl/discussionpaper/?paper=2332}.

Blasques F, Koopman SJ, Lucas A (2014b).
"Maximum Likelihood Estimation for Generalized Autoregressive Score Models."
techreport TI 2014-029/III, Tinbergen Institute.
\url{https://tinbergen.nl/discussionpaper/?paper=2286}.\cr

Blasques F, Koopman SJ, Lucas A, Schaumburg J (2014c).
"Spillover Dynamics for Systemic Risk Measurement using Spatial Financial Time Series Models."
techreport TI 2014-103/III, Tinbergen Institute.
\url{https://tinbergen.nl/discussionpaper/?paper=2369}.\cr

Creal D, Koopman SJ, Lucas A (2013).
"Generalized Autoregressive Score Models with Applications."
Journal of Applied Econometrics, 28(5), 777-795.
\doi{10.1002/jae.1279}.\cr

Ghalanos A, Theussl S (2016).
"Rsolnp: General Non-Linear Optimization using Augmented Lagrange Multiplier Method."
\url{https://cran.r-project.org/package=Rsolnp}.\cr

Gilbert P, Varadhan R (2016).
numDeriv: Accurate Numerical Derivatives. R package 2016.8-1,
\url{https://CRAN.R-project.org/package=numDeriv}.

Harvey AC (2013).
Dynamic Models for Volatility and Heavy Tails: With Applications to Financial and Economic Time Series.
Cambridge University Press.\cr

Ye Y (1988).
Interior Algorithms for Linear, Quadratic, and Linearly Constrained Convex Programming.
Ph.D. thesis, Stanford University.
}
\author{Leopoldo Catania}
%
\examples{
\dontrun{
# Specify an GAS model with multivariate Student-t
# conditional distribution and time-varying scales and correlations

library("GAS")

data("StockIndices")

GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity",
                       GASPar = list(scale = TRUE, correlation = TRUE))

Fit = MultiGASFit(GASSpec, StockIndices)

Fit
}
}
