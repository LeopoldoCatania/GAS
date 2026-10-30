\name{PIT_test}
\alias{PIT_test}
\title{
	Goodness of fit for conditional densities
}
\description{
	This function implements density goodness of fit procedure of Diebold et al. (1998).
}
\usage{
PIT_test(U, G = 20, alpha = 0.05, plot = FALSE)
}
\arguments{
\item{U}{
\code{numeric} Vector of Probability Integral Transformation (PIT).}
%
\item{G}{
Integer of at least two giving the number of equal-width histogram bins on [0,1]. The supplied value is retained for all sample sizes.}
%
\item{alpha}{
Significance level strictly between zero and one.}
%
\item{plot}{
Logical; whether to display the PIT histogram. Default is \code{FALSE}.}
}
\details{

All input values must be finite and lie in [0,1]. The histogram uses
G fixed equal-width bins over this entire interval. Its Pearson
statistic uses expected counts length(U)/G and a chi-squared reference
distribution with G-1 degrees of freedom. If an expected count is below
five, the statistic is returned but its critical value and p-value are
\code{NA}, with an explanation in \code{reason}.
The plotted count limits are approximate two-sided marginal limits
with binomial variance; they are not simultaneous limits across bins.

The dependence diagnostics regress each of the first four centered
PIT powers on an intercept and its twenty lags. They require more than
41 observations and a full-rank design. Unavailable tests are reported
as \code{NA}, with a \code{reason}, and do not abort model fitting.

For discrete distributions, ordinary CDF values F(Y) are not uniform.
The package computes randomized PITs for \code{ber}, \code{poi},
\code{negbin} and \code{skellam}:
\deqn{U_t=F_t(Y_t-1)+V_t\{F_t(Y_t)-F_t(Y_t-1)\},}
where the V variables are independent Uniform(0,1) random variables.
Set \code{set.seed()} before fitting or forecasting for reproducibility.
This consumes random numbers. Calls to \code{pit()} extract the stored
realization without drawing new random variables.

The reference distributions are asymptotic diagnostics. Estimated
in-sample parameters and dependence can affect their finite-sample
calibration; randomized PITs correct the discrete-CDF issue but do not
remove parameter-estimation effects.


This function implements density goodness of fit procedure of Diebold et al. (1998). The test relays on the
result that, if the series of estimated conditional distributions is the true one, then the PIT series
evaluated accordingly are iid Unif(0, 1) distributed.  The test of the iid Uniform(0, 1) assumption
consists of two parts. The first part concerns the independent assumption, and it tests if all the
conditional moments of the data, up to the fourth one, have been accounted for by the model, while
the second part checks if the conditional distribution assumption is reliable by testing if the
PITs are Uniform over the interval (0, 1). See also Jondeau and Rockinger (2006) and Vlaar and Palm (1993).
}
\value{A \code{list} with elements: (i) \code{Hist} and (ii) \code{IID}. The fist element \code{Hist} concerns
the test of the unconditional assumption of uniformity of the PIT, it is a \code{list} with elements:
  \itemize{
    \item \code{test} Statistic test.
    \item \code{crit} The critical value of the test.
    \item \code{pvalue} The pvalue of the test.
    \item \code{hist} The histogram, evaluated using the \link{hist} function.
    \item \code{confidence} Approximated asymptotic confidence level.
    }
The second element \code{IID} concerns the iid assumption, it is a \code{list} with elements:
  \itemize{
    \item \code{test} A named \code{numeric} vector with elements: \code{test1}, \code{test2}, \code{test3}, \code{test4}
    representing the Lagrange Multiplier test for the first four conditional moments of the PITs.
    \item \code{crit} The critical value of the test.
    \item \code{pvalue} A named \code{numeric} vector with elements: \code{pvalue1}, \code{pvalue2}, \code{pvalue3}, \code{pvalue4}
    representing the pvalues of the Lagrange Multiplier test for the first four conditional moments of the PITs.
    }
}
\author{Leopoldo Catania}
\references{
Diebold FX, Gunther TA and Tay AS (1998).
"Evaluating Density Forecasts with Applications to Financial Risk Management."
International Economic Review, 39(4), 863-883.\cr

Jondeau E and Rockinger M (2006).
"The Copula-Garch Model of Conditional Dependencies: An International Stock Market Application."
Journal of International Money and Finance, 25(5), 827-853.
\doi{10.1016/j.jimonfin.2006.04.007}.\cr

Vlaar PJ and Palm FC (1993).
"The Message in Weekly Exchange Rates in the European Monetary System: Mean Reversion, Conditional Heteroscedasticity, and Jumps."
Journal of Business & Economic Statistics, 11(3), 351-360.
\doi{10.1080/07350015.1993.10509963}.
}
\examples{
\donttest{
data("StockIndices")

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity",
                     GASPar = list(location = FALSE, scale = TRUE,
                                   shape = FALSE))

FTSEMIB = StockIndices[, "FTSEMIB"]

Fit = UniGASFit(GASSpec, FTSEMIB)

U = pit(Fit)

Test = PIT_test(U, G = 20, alpha = 0.05, plot = TRUE)
}
}
