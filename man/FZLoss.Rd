\name{FZLoss}
\alias{FZLoss}
\title{
	Fissler and Ziegel (2016) (FZ) joint loss function for Value at Risk and Expected Shortfall.
}
\description{
	This function implements Fissler and Ziegel (2016) (FZ) joint loss function for Value at Risk and Expected Shortfall.
}
\usage{
FZLoss(data, VaR, ES, alpha)
}
\arguments{
\item{data}{
\code{numeric} Vector of observations.}
%
\item{VaR}{
\code{numeric} Vector containing the VaR series.}
%
\item{ES}{
\code{numeric} Vector containing the ES series.}
%
\item{alpha}{
\code{numeric} The VaR and ES confidence level.}
%
}
\details{
This function implements Fissler and Ziegel (2016) (FZ) joint loss function for Value at Risk and Expected Shortfall. The parameterization used is that of Patton et al. (2017) and is given by:

\eqn{\frac{1}{\alpha ES_t^\alpha}I_t^\alpha(y_t - VaR_t^\alpha) + \frac{VaR_t^\alpha}{ES_t^\alpha} + \log{-ES_t^\alpha} - 1}.

See also Fissler et al. (2015).
}
\value{A \code{numeric} vector containing the joing VaR and ES loss values.}
\author{Leopoldo Catania}
\references{

Fissler, T., Ziegel, J.F., (2016).
"Higher order elicitability and Osband's principle."
The Annals of Statistics 44, 1680-1707.
\doi{10.1214/16-AOS1439}.\cr
%
Fissler, T., Ziegel, J.F., Tilmann, G. (2015).
"Expected Shortfall is jointly elicitable with Value at Risk - Implications for backtesting."
arXiv preprint arXiv:1507.00244.\cr
%
Patton, A. J., Ziegel, J.F., Chen, R. (2017).
"Dynamic semiparametric models for expected shortfall (and Value-at-Risk)."
arXiv preprint arXiv:1707.05108.

}
\examples{

data("StockIndices")

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity",
                     GASPar = list(location = FALSE, scale = TRUE,
                                   shape = FALSE))

FTSEMIB = StockIndices[, "FTSEMIB"]

InSampleData  = FTSEMIB[1:1500]
OutSampleData = FTSEMIB[1501:2404]

Fit = UniGASFit(GASSpec, InSampleData)

Forecast = UniGASFor(Fit, Roll = TRUE, out = OutSampleData)

alpha = 0.05

vVaR = quantile(Forecast, alpha)
vES  = ES(Forecast, alpha)

FZ = FZLoss(OutSampleData, vVaR, vES, alpha)


}
