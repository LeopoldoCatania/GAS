\name{BacktestDensity}
\alias{BacktestDensity}
\title{
	Backtest a series of one-step ahead density predictions.
}
\description{
The \code{BacktestDensity()} function accepts an object of the class \link{uGASRoll}, and returns a \code{list} with
two elements: (i) the averages Negative Log Score (NLS) and  weighted Continuous Ranked Probability Score (wCRPS)
introduced by Gneiting and Ranjan (2012), and (ii) their values at each point in time. The wCRPS is
evaluated using 5 weight functions, see Details.
}
\usage{
BacktestDensity(Roll, lower, upper, K = 1000, a = NULL, b = NULL)
}
\arguments{
\item{Roll}{
an object of the class \link{uGASRoll}.}
%
\item{lower}{
\code{numeric} the lower bound used to approximate the wCRSP by Monte Carlo integration as
detailed in Gneiting and Ranjan (2012).  This coincides with \eqn{y_{l}} in Equation 16 of
Gneiting and Ranjan (2012).}
%
\item{upper}{
\code{numeric} the upper bound used to approximate the wCRSP by Monte Carlo integration as detailed
in Gneiting and Ranjan (2012).  This coincides with \eqn{y_{u}} in Equation 16 of Gneiting and Ranjan (2012).}
%
\item{K}{
\code{numeric} integer representing the number of points used to discretize the wCRPS integral.
This is \eqn{I} in Equation 16 of Gneiting and Ranjan (2012). By default \code{K = 1000.}}
%
\item{a}{
\code{numeric}. mean of the normal distribution used in the weight function. By default \code{a = NULL}, which
means that it is set equal to the empircial mean of the in sample observations.}
%
\item{b}{
\code{numeric}. standard deviation of the normal distribution used in the weight function. By default \code{b = NULL}, which
means that it is set equal to the empircial standard deviation of the in sample observations.}
%
}
\details{
The average Negative Log Score (NLS) is computed as the negative of the average of the log scores evaluated
during the out-of-sample period. The average weighted Continuous Ranked Probability Score (wCRPS) is computed as
the average of the wCRPS evaluated during the out-of-sample period, see Gneiting and Ranjan (2012).\cr
%
The wCRPS is evaluated using Equation 16 of Gneiting and Ranjan (2012). The weights functions implemented are:
\itemize{
\item \eqn{w(z) = 1}: Uniform,
\item \eqn{w(z) = \phi_{a,b}(z)}: Center,
\item \eqn{w(z) = 1 - \phi_{a,b}(z)}:  Tails,
\item \eqn{w(z) = \Phi_{a,b}(z)}: Right tail,
\item \eqn{w(z) = 1 - \Phi_{a,b}(z)}:  Left tail,
}
where \eqn{\phi_{a,b}(z)} and \eqn{\Phi_{a,b}(z)} are the pdf and cdf of a Gaussian distribution with mean \eqn{a}
and standard deviation \eqn{b}, respectively. The label "Uniform" represents the case where equal
emphasis is given to all the parts of the distribution.
}
\value{A \code{list} with elements: \code{average}, \code{series}. The element "\code{average}" is a named vector
with the averages NLS and wCRSP. The element "\code{series}" is a \code{list}: the fist element, \code{LS}, contains
the out-of-sample Log Score (not with the negative sign), the second element, \code{WCRPS}, contains a \code{matrix}
with the wCRPS series. The columns of this matrix are named: "\code{uniform}", "\code{center}", "\code{tails}",
"\code{tail_r}", "\code{tail_l}", which are associated with the wCRSP with emphasis on: Uniform, Center, Tails,
Right tail and Left tail, respectively.}
\author{Leopoldo Catania}
\references{
Gneiting T, Ranjan R (2011).
"Comparing Density Forecasts using Threshold -and Quantile-Weighted Scoring Rules."
Journal of Business & Economic Statistics, 29(3), 411-422.
\doi{10.1198/jbes.2010.08110}.
}
\examples{
\dontrun{
data("cpichg")

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity",
                       GASPar = list(location = TRUE, scale = TRUE,
                                      shape = FALSE))

Roll = UniGASRoll(cpichg, GASSpec, ForecastLength = 50,
                  RefitEvery = 10, RefitWindow = c("moving"))

BackTest = BacktestDensity(Roll, lower = -100, upper = 100)

BackTest$average
}
}
