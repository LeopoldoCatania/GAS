\name{BacktestVaR}
\alias{BacktestVaR}
\title{
	Backtest Value at Risk (VaR)
}
\description{
	This function implements several backtesting procedures for the Value at Risk (VaR). These are: (i) The statistical tests of Kupiec (1995), Christoffesen (1998) and Engle and Manganelli (2004), (ii) The tick loss function detailed in Gonzalez-Rivera et al. (2004), the mean and max absolute loss used by McAleer and Da Veiga (2008) and the actual over expected exceedance  ratio.
}
\usage{
BacktestVaR(data, VaR, tau, alphaTest = 0.95, Lags = 4)
}
\arguments{
\item{data}{
\code{numeric} Vector of observations.}
%
\item{VaR}{
\code{numeric} Vector containing the VaR series.}
%
\item{tau}{
\code{numeric} The VaR confidence level.}
%
\item{alphaTest}{
\code{numeric} Confidence level used in the statistical tests.}
%
\item{Lags}{
\code{numeric} Lags used in the Dynamic Quantile test of Engle and Manganelli (2004).}
}
\details{
		This function implements several backtesting procedure for the Value at Risk. The implemented statistical tests are:
  \itemize{
    \item \code{LRuc} The unconditional coverage test of Kupiec (1995).
    \item \code{LRcc} The conditional coverage test of Christoffesen (1998).
    \item \code{DQ} The Dynamic Quantile test of Engle and Manganelli (2004).
    }
    The implemented VaR backtesting quantities are:
      \itemize{
    \item \code{AD} mean and maximum absolute deviation between the observations and the quantiles as in McAleer and Da Veiga (2008).
    \item \code{Loss} Average quantile loss and quantile loss series as in Gonzalez-Rivera et al. (2004).
    \item \code{AE} Actual over Expected exceedance  ratio.
    }
}
\value{A \code{list} with elements: \code{LRuc}, \code{LRcc}, \code{DQ}, \code{AD}, \code{AE}.}
\author{Leopoldo Catania}
\references{
Christoffersen PF (1998). "Evaluating interval forecasts." International Economic Review, 39(4), 841-862, \doi{10.2307/2527341}.\cr

Engle RF and Manganelli S. (2004). "CAViaR: Conditional Autoregressive Value at Risk by Regression Quantiles." Journal of Business & Economic Statistics, 22(4), 367-381, \doi{10.1198/073500104000000370}.\cr

Gonzalez-Rivera G, Lee TH, and Mishra, S (2004). "Forecasting Volatility: A Reality Check Based on Option Pricing, Utility Function, Value-at-Risk, and Predictive Likelihood." International Journal of Forecasting, 20(4), 629-645, \doi{10.1016/j.ijforecast.2003.10.003}.\cr

Kupiec PH (1995). "Techniques for Verifying the Accuracy of Risk Measurement Models." The Journal of Derivatives, 3(2), 73-84, \doi{10.3905/jod.1995.407942}\cr

McAleer M and Da Veiga B (2008). "Forecasting Value-at-Risk with a Parsimonious Portfolio Spillover GARCH (PS-GARCH) Model." Journal of Forecasting, 27(1), 1-19, \doi{10.1002/for.1049}.\cr
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

tau = 0.05

VaR = quantile(Forecast, tau)

BackTest = BacktestVaR(OutSampleData, VaR, tau)

}
