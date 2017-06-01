\name{fn.optim}
\alias{fn.optim}
\title{
	A wrapper to the \link{optim} function.
}
\description{
	This function is just a wrapper to the standard \link{optim}  optimizer with \code{method = "BFGS"}.
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
The following control parameters are used for \code{control}:
%
\itemize{
   \item \code{trace = 0}
   \item \code{abstol = 1e-8}
}
%
See the documentation of \link{optim}.
%
}
\value{
It returns a named list with four elements: i) \code{pars}: a \code{numeric} vector
where the estimated parameters are stored, ii) \code{value}: a \code{numeric} containing the value of the negative log likelihood
evaluated at its minumum, iii) \code{hessian}, a \code{numeric} matrix containing the Hessian matrix evaluated at
the minimum of the negative log likelihood, iv) \code{convergence} a \code{numeric} element indicating the convergence results of
\link{optim}.
}
\author{Leopoldo Catania}
%
\seealso{
\code{help(optim)}
}
