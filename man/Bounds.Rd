\name{NumericalBounds}
\alias{NumericalBounds}
\alias{LowerNu}
\alias{UpperNu}
\alias{LowerA}
\alias{LowerB}
\alias{UpperA}
\alias{UpperB}
\title{
Numerical bounds imposed in parameters transformation.
}
\description{
Prints the numerical bounds.
}
\usage{
UpperNu()
LowerNu()
UpperA()
LowerA()
UpperB()
LowerB()
}
\details{
\code{UpperNu()} and \code{LowerNu()} print the numerical upper and lower bounds for the degree 
of freedom parameter of the Student--t distribution, \code{std}. (including also \code{sstd} and \code{mvt}).\cr
%
\code{UpperA()} and \code{LowerA()} print the numerical upper and lower bounds for the score parameter 
in the GAS recursion. These bounds are applied to each diagonal element of the matrix \code{A} that premultiplies the scaled score.\cr
%
\code{UpperB()} and \code{LowerB()} print the numerical upper and lower bounds for the autoregressive 
parameter in the GAS recursion. These bounds are applied to each diagonal element 
of the matrix \code{B} that premultiplies the past value of the parameters.\cr
%
}
\value{
Prints the numerical bounds.
}
\author{
Leopoldo Catania
}
\examples{
UpperNu()
LowerNu()
UpperA()
LowerA()
UpperB()
LowerB()
}
