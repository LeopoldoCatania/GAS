\name{UniGASFit}
\alias{UniGASFit}
\title{
	Estimate univariate GAS models
}
\description{
	Estimate univariate GAS models by Maximum Likelihood.
}
\usage{
UniGASFit(GASSpec, data, fn.optimizer = fn.solnp)
}
\arguments{
\item{GASSpec}{
An object of the class \link{uGASSpec} created using the function \link{UniGASSpec}}.
%
\item{data}{\code{numeric} vector of length Tx1 containing the time series of observations.
It can also be an object of the class \code{ts}, \code{xts} or \code{zoo}.}
%
\item{fn.optimizer}{\code{function}. This is a generic optimization function that can be provided by the user.
The default optimizer is the  General Nonlinear Augmented Lagrange Multiplier method of Ye (1988) available in
the R package Rsolnp (Ghalanos and Theussl 2016). By default \code{fn.optimizer = fn.solnp} where \link{fn.solnp} is
a wrapper to the \link{solnp} function of the Rsolnp R package. See Details for user defined optimization routines.}
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
version of the model (i.e., with A = 0 and B = 0) and set the initial value of the intercept parameter
accordingly, and (ii) perform a grid search for the coefficients contained in A and B.
Further technical details are presented in Section 3.2 of Ardia et. al. (2016a).\cr
%

The user is free to employ his/her own optimization routine via the \code{fn.optimizer} argument. \code{fn.optimizer}
accepts a \code{function} object. The user provided optimizer has to satisfy strict requirements. The arguments of the
\code{fn.optimizer} are : i) \code{par0} a vector of starting values, ii) \code{data} the data provided, iii) \code{GASSpec}
an object of the class \link{uGASSpec}, and iv) \code{FUN} the likelihood function. The output of \code{fn.optimizer} have
to be an object of the class \code{list} with four named elements: i) \code{pars}: a \code{numeric} vector
where the estimated parameters are stored, ii) \code{value}: a \code{numeric} containing the value of the negative log likelihood
evaluated at its minimum, iii) \code{hessian}, a \code{numeric} matrix containing the Hessian matrix evaluated at
the minimum of the negative log likelihood, this is used for inferential purposes, and iv) \code{convergence} a \code{numeric}
or \code{character} string reporting information about the convergence of the optimization. This quantity is printed in used by the
\code{show()} and \code{summary()} methods.

The user is allowed to not include the last two elements of the output of the \code{fn.optimizer} function, that is, the values
\code{hessian = NULL} and \code{convergence = NULL} are admissible. In the case of \code{hessian = NULL}, the Hessian matrix is
evaluated numerically using the \link{hessian} function in the \code{numDeriv} package of Gilbert and Varadhan (2016). If the provided
hessian is not positive defined, a try with the hessian evaluation used by the function \link{optim} is made.

}
\value{
An object of the class \link{uGASFit}
}
\references{
Ardia D, Boudt K and Catania L (2016a).
"Generalized Autoregressive Score Models in R: The GAS Package."
\url{http://ssrn.com/abstract=2825380}.\cr

Blasques F, Koopman SJ, Lucas A (2014a).
"Maximum Likelihood Estimation for Correctly Specified Generalized Autoregressive Score Models: Feedback Effects, Contraction Conditions and Asymptotic Properties." techreport TI 14-074/III, Tinbergen Institute.
\url{http://www.tinbergen.nl/discussionpaper/?paper=2332}.

Blasques F, Koopman SJ, Lucas A (2014b).
"Maximum Likelihood Estimation for Generalized Autoregressive Score Models."
techreport TI 2014-029/III, Tinbergen Institute.
\url{http://www.tinbergen.nl/discussionpaper/?paper=2286}.\cr

Blasques F, Koopman SJ, Lucas A, Schaumburg J (2014c).
"Spillover Dynamics for Systemic Risk Measurement using Spatial Financial Time Series Models."
techreport TI 2014-103/III, Tinbergen Institute.
\url{http://www.tinbergen.nl/discussionpaper/?paper=2369}.\cr

Creal D, Koopman SJ, Lucas A (2013).
"Generalized Autoregressive Score Models with Applications."
Journal of Applied Econometrics, 28(5), 777-795.
\doi{10.1002/jae.1279}.\cr

Gilbert P, Varadhan R (2016).
numDeriv: Accurate Numerical Derivatives. R package 2016.8-1,
\url{https://CRAN.R-project.org/package=numDeriv}.

Ghalanos A, Theussl S (2016).
"Rsolnp: General Non-Linear Optimization using Augmented Lagrange Multiplier Method."
\url{https://cran.r-project.org/package=Rsolnp}.\cr

Harvey AC (2013).
Dynamic Models for Volatility and Heavy Tails: With Applications to Financial and Economic Time Series.
Cambridge University Press.\cr

Ye Y (1988).
Interior Algorithms for Linear, Quadratic, and Linearly Constrained Convex Programming.
Ph.D. thesis, Stanford University.


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

# Estimate the model with a different optimizer.
# Assume we want to use the Nelder and Mead optimization provided by
# the optim() function, we create
# the wrapper fn.optim in this way

fn.NM.optim <- function(par0, data, GASSpec, FUN) {

  optimizer = optim(par0, FUN, data = data, GASSpec = GASSpec, method = "Nelder-Mead",
                    control = list(trace = 0), hessian = TRUE)

  out = list(pars = optimizer$par,
             value = optimizer$value,
             hessian = optimizer$hessian,
             convergence = optimizer$convergence)

  return(out)

}

Fit.NM.optim = UniGASFit(GASSpec, sp500ret, fn.optimizer = fn.NM.optim )

Fit.NM.optim

# Estimate the model with a different optimizer.
# Assume we want to use the BFGS optimization provided by
# the optim() function, we create
# the wrapper fn.optim in this way

fn.BFGS.optim <- function(par0, data, GASSpec, FUN) {

  optimizer = optim(par0, FUN, data = data, GASSpec = GASSpec, method = "BFGS",
                    control = list(trace = 0), hessian = TRUE)

  out = list(pars = optimizer$par,
             value = optimizer$value,
             hessian = optimizer$hessian,
             convergence = optimizer$convergence)

  return(out)

}

Fit.BFGS.optim = UniGASFit(GASSpec, sp500ret, fn.optimizer = fn.BFGS.optim )

Fit.BFGS.optim

}
}
