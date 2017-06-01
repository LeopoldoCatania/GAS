\name{fn.solnp}
\alias{fn.solnp}
\title{
	A wrapper to the \link{solnp} function of the \code{Rsolnp} package of Ghalanos and Theussl (2016).
}
\description{
	This function is just a wrapper to the \link{solnp} function of the \code{Rsolnp} package of Ghalanos and Theussl (2016).
}
\usage{
 fn.solnp(par0, data, GASSpec, FUN)
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
The following control parameters are used: \code{trace = 0}, \code{rho = 1}, \code{outer.iter = 400},
\code{inner.iter = 1800}, \code{delta = 1e-08}, \code{tol = 1e-08}. See the documentation of \link{solnp}.
}
\value{
It returns a named list with four elements: i) \code{pars}: a \code{numeric} vector
where the estimated parameters are stored, ii) \code{value}: a \code{numeric} containing the value of the negative log likelihood
evaluated at its minumum, iii) \code{hessian}, a \code{numeric} matrix containing the Hessian matrix evaluated at
the minimum of the negative log likelihood, and iv) \code{convergence} a \code{numeric} element indicating the convergence results of
\link{solnp}.
}
\author{Leopoldo Catania}
\references{
Alexios Ghalanos and Stefan Theussl (2015).
"Rsolnp: General Non-linear Optimization Using Augmented Lagrange
  Multiplier Method". R package version 1.16.\cr
}
\seealso{

\code{help(solnp)}

}
