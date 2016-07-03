\name{PIT_test}
\alias{PIT_test}
\title{
	Backtest Value at Risk (VaR)
}
\description{
	This function implements density goodness of fit procedure of Diebold et al. (1998).
}
\usage{
PIT_test(vU, iG=20, dAlpha=0.05, plot=F)
}
\arguments{
\item{vU}{
\code{numeric} vector of Probability Integral Transformation, see \link{pit}.}
%
\item{iG}{
\code{numeric} number of bins of the empirical cumulative density function of the PIT.}
%
\item{dAlpha}{
\code{numeric}. test confidence level.}
%
\item{plot}{
\code{boolean}. Plot the result ?.}
}
\details{
 This function implements density goodness of fit procedure of Diebold et al. (1998). The test relays on the result that, if the series of estimated conditional distributions is the true one, then the PIT series evaluated accordingly are iid Unif(0, 1) distributed.  The test of the iid Uniform(0, 1) assumption consists of two parts. The first part concerns the independent assumption, and it tests if all the conditional moments of the data, up to the fourth one, have been accounted for by the model, while the second part checks if the conditional distribution assumption is reliable by testing if the PITs are Uniform over the inverval (0, 1). See also Jondeau and Rockinger (2006) and Vlaar and Palm (1993).
}
\value{A \code{list} with elements: (i) Hist and (ii) IID. The fist element Hist concerns the test of the unconditional assumption of uniformity of the PIT, it is a \code{list} with elements:
  \itemize{
    \item \code{test}: statistic test.
    \item \code{crit}: The critical value of the test.
    \item \code{pvalue}: the pvalue of the test.
    \item \code{hist} : the histogram, evaluated using the \link{hist} function
    \item \code{confidence} : approximated asymptotic confidence level.
    }
The second element \code{IID} concerns the iid assumption, it is a \code{list} with elements:
  \itemize{
    \item \code{test}: a named \code{numeric} vector with elements: test1, test2, test3, test4 representing the Lagrange Multiplier test for the first four conditional moments of the PITs.
    \item \code{crit}: The critical value of the test.
    \item \code{pvalue}: a named \code{numeric} vector with elements: pvalue1, pvalue2, pvalue3, pvalue3 representing the pvalues of the Lagrange Multiplier test for the first four conditional moments of the PITs.
    }
}
\author{Leopoldo Catania}
\references{
Diebold, F. X., Gunther, T. A., & Tay, A. S. (1998). Evaluating Density Forecasts with Applications to Financial Risk Management. International Economic Review, 863-883.\cr

Jondeau, E., & Rockinger, M. (2006). The copula-garch model of conditional dependencies: An international stock market application. Journal of international money and finance, 25(5), 827-853.\cr

Vlaar, P. J., & Palm, F. C. (1993). The message in weekly exchange rates in the European monetary system: mean reversion, conditional heteroscedasticity, and jumps. Journal of Business & Economic Statistics, 11(3), 351-360.\cr
}
\examples{
\dontrun{
data("StockIndex")

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = F, scale = T, shape = F))

FTSEMIB = StockIndex["FTSEMIB",]

Fit = UniGASFit(GASSpec, FTSEMIB)

vU = pit(Fit)

Test = PIT_test(vU, iG=20, dAlpha=0.05, plot=T)

}
}
