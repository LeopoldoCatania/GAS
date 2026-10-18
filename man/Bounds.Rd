\name{NumericalBounds}
\alias{NumericalBounds}
\alias{LowerNu}
\alias{UpperNu}
\alias{LowerA}
\alias{LowerB}
\alias{UpperA}
\alias{UpperB}
\title{Numerical bounds used in parameter transformations}
\description{Returns the numerical bounds used by the package's parameter transformations.}
\usage{
UpperNu(Dist = "std")
LowerNu(Dist = "std")
UpperA()
LowerA()
UpperB()
LowerB()
}
\arguments{
  \item{Dist}{Distribution label: \code{"std"}, \code{"sstd"}, \code{"ast"},
  \code{"ast1"}, \code{"ghskt"}, or \code{"mvt"}.}
}
\details{
\code{LowerNu} and \code{UpperNu} are R wrappers around the same C++ functions
used in the forward mappings, inverse mappings and mapping Jacobians.
The default distribution is \code{"std"}; calls without arguments remain supported.
Unsupported distribution labels cause an error.

The lower bound is 0.001 for \code{std}, \code{ast}, \code{ast1}, and \code{mvt};
2 for \code{sstd}; and 4 for \code{ghskt}. Both degrees-of-freedom parameters of
\code{ast} use the same bounds. The upper bound is 50 for all six distributions.
The mapped parameter lies strictly between the bounds. The lower bound of 0.001
is a numerical floor. The higher bounds for \code{sstd} and \code{ghskt} reflect
their current mean/variance-standardized parameterizations. Bounds are not
increased merely to ensure that third or fourth moments exist.

\code{UpperA} and \code{LowerA} return bounds for the diagonal score-loading
coefficients in the GAS recursion. \code{UpperB} and \code{LowerB} return bounds
for the diagonal autoregressive coefficients.

Models fitted with earlier shape mappings should be refitted before using
this implementation: the stored unrestricted shape states have changed meaning.
}
\value{A numeric scalar giving the requested bound.}
\author{Leopoldo Catania}
\examples{
LowerNu()           # std: 0.001
LowerNu("sstd")     # 2
LowerNu("ghskt")    # 4
UpperNu("mvt")      # 50
UpperA()
LowerA()
UpperB()
LowerB()
}
