\name{fn.optim}
\alias{fn.optim}
\title{
	A wrapper to the \link{optim} function.
}
\description{
	Minimize the negative log likelihood using BFGS, with a deterministic alternative solver when numerical checks fail.
}
\usage{
 fn.optim(par0, data, GASSpec, FUN, Compute.Hessian = FALSE)
}
\arguments{
\item{par0}{
\code{numeric} vector of named model coefficients.}
%
\item{data}{
\code{numeric} vector or matrix of data.}
%
\item{GASSpec}{
An object of the class \link{uGASSpec} or \link{mGASSpec}, created via the \link{UniGASSpec} and \link{MultiGASSpec} functions.}
%
\item{Compute.Hessian}{Logical, default \code{FALSE}. Whether to compute
numerical Hessians for the curvature diagnostic and final output.}
\item{FUN}{
A \link{function} to optimize.}
}
\details{
Optimization uses the average negative log likelihood to stabilize the
numerical scale across sample sizes. The primary BFGS controls are
\code{maxit = 2000L}, \code{reltol = 1e-11}, and \code{abstol = -Inf}.
The starting point is retained as a candidate. If convergence or gradient
diagnostics fail (or curvature diagnostics fail when
\code{Compute.Hessian = TRUE}), \code{stats::nlminb} and a subsequent
BFGS refinement are attempted. The refinement supplies Richardson numerical
gradients from \code{numDeriv::grad}, using \code{reltol = 1e-14}, to avoid
the coarse default finite-difference step in \code{optim}. The best finite
objective is retained; among numerically tied objectives, a candidate that
passes both solver convergence and gradient checks is preferred.
The reported value and optional Hessian refer to the summed negative log likelihood.
With \code{Compute.Hessian = FALSE}, neither numerical Hessian is computed
and the returned \code{hessian} is \code{NULL}. Convergence and gradient
checks, including gradient-based refinement when needed, remain active.
Curvature is used for optimizer diagnostics only; fitted-model standard errors
and confidence bands use the outer product of per-observation scores.
The additional \code{method} field identifies the selected candidate.
}
\value{
It returns a named list including the following four elements: i) \code{pars}: a \code{numeric} vector
where the estimated parameters are stored, ii) \code{value}: a \code{numeric} containing the value of the negative log likelihood
evaluated at its minumum, iii) \code{hessian}, a \code{numeric} matrix containing the Hessian matrix evaluated at
the minimum of the negative log likelihood (or \code{NULL} when not requested
or evaluation fails), iv) \code{convergence} a \code{numeric} element indicating the selected solver convergence result.
}
\author{Leopoldo Catania}
%
\seealso{
\code{help(optim)}
}
