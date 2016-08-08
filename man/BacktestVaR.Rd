\name{BacktestVaR}
\alias{BacktestVaR}
\title{
	Backtest Value at Risk (VaR)
}
\description{
	This function implements several backtesting procedure for the Value at Risk. These are: (i) The statistical tests of Kupiec (1995), Christoffesen (1998) and Engle and Manganelli (2004), (ii) The tick loss function detailed in Gonz\'{a}lez-Rivera et al. (2004), the mean and max absolute loss used by McAleer and Da Veiga (2008) and the actual over expected exceedance  ratio.
}
\usage{
BacktestVaR(data, VaR, tau, alphaTest=0.95, Lags=4)
}
\arguments{
\item{data}{
\code{numeric} vector of observations.}
%
\item{VaR}{
\code{numeric} vector containing the Value at Risk series.}
%
\item{tau}{
\code{numeric}. The VaR confidence level.}
%
\item{alphaTest}{
\code{numeric}. Confidece level used in the statistical tests.}
%
\item{Lags}{
\code{numeric}. Lags used in the Dynamic Quantile test of Engle and Manganelli (2004).}
}
\details{
		This function implements several backtesting procedure for the Value at Risk. The implemented statistical tests are:
  \itemize{
    \item \code{LRuc}: The unconditioanl coverage test of Kupiec (1995).
    \item \code{LRcc}: The conditioanl coverage test of Christoffesen (1998).
    \item \code{DQ}: The Dynamic Quantile test of Engle and Manganelli (2004).
    }
    The implemented VaR backtesting quantities are:
      \itemize{
    \item \code{AD}: mean and maximum absolute deviation between the observations and the quantiles as in McAleer and Da Veiga (2008).
    \item \code{Loss}: Average quantile loss and quantile loss series as in Gonz\'{a}lez-Rivera et al. (2004).
    \item \code{AE}: Actual over Expected exceedance  ratio.
    }
}
\value{A \code{list} with elements: LRuc, LRcc, DQ, AD, AE.}
\author{Leopoldo Catania}
\references{
Christoffersen, P. F. (1998). Evaluating interval forecasts. International economic review, 841-862.\cr

Engle, R. F., & Manganelli, S. (2004). CAViaR: Conditional autoregressive value at risk by regression quantiles. Journal of Business & Economic Statistics, 22(4), 367-381.\cr

Gonz\'{a}lez-Rivera, G., Lee, T. H., & Mishra, S. (2004). Forecasting volatility: A reality check based on option pricing, utility function, value-at-risk, and predictive likelihood. International Journal of forecasting, 20(4), 629-645.\cr

Kupiec, P. H. (1995). Techniques for verifying the accuracy of risk measurement models. The J. of Derivatives, 3(2).\cr

McAleer, M., & Da Veiga, B. (2008). Forecasting value-at-risk with a parsimonious portfolio spillover GARCH (PS-GARCH) model. Journal of Forecasting, 27(1), 1-19.\cr
}
\examples{
\dontrun{
data("StockIndex")

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity",
                     GASPar = list(location = F, scale = T, shape = F))

FTSEMIB = StockIndex["FTSEMIB",]

InSampleData  = FTSEMIB[1:1500]
OutSampleData = FTSEMIB[1501:2404]

Fit = UniGASFit(GASSpec, InSampleData)

Forecast = UniGASFor(Fit, Roll = TRUE, vOut = OutSampleData)

dTau = 0.05

vVaR = quantile(Forecast, dTau)

BackTest = BacktestVaR(OutSampleData, vVaR, dTau)

}
}
