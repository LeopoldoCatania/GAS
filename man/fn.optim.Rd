\name{fn.optim}
\alias{fn.optim}
\title{
	A wrapper to the \link{optim} function.
}
\description{
	Minimize the negative log likelihood using BFGS, with a deterministic alternative solver when numerical checks fail.
}
\usage{
 fn.optim(par0, data, GASSpec, FUN)
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
\item{FUN}{
A \link{function} to optimize.}
}
\details{
Optimization uses the average negative log likelihood to stabilize the
numerical scale across sample sizes. The primary BFGS controls are
\code{maxit = 2000L}, \code{reltol = 1e-11}, and \code{abstol = -Inf}.
The starting point is retained as a candidate. If convergence, gradient
or curvature diagnostics fail, \code{stats::nlminb} and a subsequent
BFGS refinement are attempted. The best finite objective is retained.
The reported value and Hessian refer to the summed negative log likelihood.
Curvature is used for optimizer diagnostics only; fitted-model standard errors
and confidence bands use the outer product of per-observation scores.
The additional \code{method} field identifies the selected candidate.
}
\value{
It returns a named list including the following four elements: i) \code{pars}: a \code{numeric} vector
where the estimated parameters are stored, ii) \code{value}: a \code{numeric} containing the value of the negative log likelihood
evaluated at its minumum, iii) \code{hessian}, a \code{numeric} matrix containing the Hessian matrix evaluated at
the minimum of the negative log likelihood, iv) \code{convergence} a \code{numeric} element indicating the selected solver convergence result.
}
\author{Leopoldo Catania}
%
\seealso{
\code{help(optim)}
}
