\name{ConfidenceBands}
\alias{ConfidenceBands}
\title{Pointwise confidence bands for time-varying parameters}

\description{
Construct pointwise confidence bands for time-varying distribution
parameters by sampling model coefficients from their estimated
asymptotic normal distribution and filtering the observed data
for each draw.
}

\usage{
ConfidenceBands(object, B = 10000L,
                probs = c(0.01, 0.1, 0.9, 0.99), ...)
}

\arguments{
\item{object}{
A fitted object of class \code{\link{uGASFit}} or
\code{\link{mGASFit}}.
}
\item{B}{
A positive integer giving the number of coefficient draws.
}
\item{probs}{
A nonempty numeric vector of probabilities between zero and one.
The returned quantiles follow the supplied order.
}
\item{...}{
Currently unused.
}
}

\details{
The function implements the in-sample simulation-based approach
described in Section 3.3 of Blasques et al. (2016).

Model coefficients are sampled on the unrestricted estimation
scale. The covariance estimate is constructed from the Hessian
stored in the fitted object's optimiser output. Reliable use of
the normal approximation requires a finite, positive definite
Hessian at an interior optimum. The function rejects fits failing the
convergence or gradient checks and rejects invalid or numerically singular
Hessians. It uses the inverse of the validated Hessian without generalized
inversion or automatic positive-definiteness repair.

Each sampled coefficient vector is transformed using the mapping
appropriate to the fitted model. For multivariate models, the
fitted \code{ScalarParameters} specification is respected.
The observed data are then filtered under each sampled vector.

Quantiles are computed separately for each time point and each
distribution parameter, using the B simulated paths. The fitted
path is not included as an additional simulation draw.

These are pointwise parameter confidence bands. They are not
simultaneous bands over the full path or prediction intervals
for future observations.

For a central 95 percent pointwise interval, use
\code{probs = c(0.025, 0.975)}.
}

\value{
A numeric array with dimensions
\code{c(T + 1L, length(probs), K)}, where T is the number of
observations and K is the number of distribution parameters.

The first dimension indexes time. Its first T entries correspond
to the in-sample parameter path, and its final entry corresponds
to the one-step-ahead parameter prediction.

The second dimension indexes the requested quantile probabilities.
The third dimension indexes the distribution parameters, in the
order returned by \code{\link{getFilteredParameters}}.

The array contains quantiles only. The fitted parameter path is
available separately through
\code{getFilteredParameters(object)}.
}
\author{Leopoldo Catania}
\references{
Blasques F, Koopman SJ, Lasak K, and Lucas, A (2016).
"In-sample Confidence Bands and Out-of-Sample Forecast Bands for Time-Varying Parameters in Observation-Driven Models."
International Journal of Forecasting, 32(3), 875-887.
\doi{10.1016/j.ijforecast.2016.04.002}.
}
\examples{
\dontrun{

library("GAS")

data("cpichg")

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity",
                     GASPar = list(location = TRUE, scale = TRUE,
                                   shape = FALSE))

Fit = UniGASFit(GASSpec, cpichg)

Bands <- ConfidenceBands(
    Fit,
    B = 1000L,
    probs = c(0.025, 0.975)
)

lower <- Bands[, 1L, ]
upper <- Bands[, 2L, ]

# Fitted path, returned separately:
fitted_path <- getFilteredParameters(Fit)
}
}
