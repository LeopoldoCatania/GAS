# OPG inference update

This complete GAS 0.4.1 source package includes the earlier distribution-specific
nu bounds and moment-existence changes, plus the OPG inference update.

## Install

Restart R, then select this complete source tarball:

```r
install.packages(file.choose(), repos = NULL, type = "source")
library(GAS)
```

Refit existing models to obtain OPG standard errors and confidence bands.
Do not replace the package's src directory with a partial set of changed files.

## Calculation

At the final unrestricted coefficient vector p, numDeriv::jacobian evaluates
all T log likelihood contributions. Each evaluation reruns the complete filter,
including initialization. Rows of the T-by-p matrix S are observation scores.
For multivariate models, one contribution is the joint log density of the
observation vector, rather than a separate contribution for each series.

The information matrix is crossprod(S). Its inverse estimates Var(p_hat), with
no extra division by T. The corresponding asymptotic covariance of
sqrt(T)*(p_hat-p_true) is T*solve(crossprod(S)). The existing delta method maps
the coefficient covariance to the reported parameter scale. This is the OPG
estimator for correctly specified ML, not a misspecification-robust sandwich.

UniGASFit, MultiGASFit and ConfidenceBands now use this information matrix.
Optimizer Hessians may still be calculated for optimization diagnostics and
retained as optimizer output. They are not used to compute standard errors or
confidence bands. No additional optimization is performed for OPG evaluation.

Inspect the result with:

```r
Fit@Estimates$optimiser$scores
Fit@Estimates$optimiser$opg
Fit@Estimates$optimiser$covariance
Fit@Estimates$optimiser$covariance_method
Fit@Estimates$optimiser$diagnostics
```

The covariance is on the unrestricted optimizer scale; reported standard errors
also include the delta-method transformation. A finite covariance is only used
for inference when convergence and gradient checks pass. Non-finite or
numerically singular information matrices produce NA standard errors and a
warning; no generalized inverse or automatic eigenvalue repair is used.
Compute.SE = FALSE skips the score Jacobian and OPG calculation.

Moment-existence warnings are unchanged. In particular, a Student-t estimate
with nu between 2 and 3 has a finite variance but undefined skewness and kurtosis.

## Source changes

- R/Likelihood.R: strict per-observation contribution evaluator.
- R/EstimationChecks.R: numerical score Jacobian, OPG and validation.
- R/Estimate.R and R/Inference.R: use OPG for fitted-model inference.
- R/ConfidenceBands.R: draw coefficients using inverse OPG covariance.
- R/OptimiserWrapper.R: clarify the retained optimizer Hessian's role.
- man/: update fit, optimizer and confidence-band documentation.
- tests/testthat/test_OPG.R: regression tests for analytic scores, full-recursion
  derivatives, scaling, multivariate contributions, delta method, confidence
  bands, invalid optimizer Hessians and Compute.SE = FALSE.

No C++ code is changed by this OPG update. Package version remains 0.4.1.

## Display update: average filtered parameters

The show and summary methods for univariate and multivariate fits now report
Average Filtered Parameters. These are the row means of GASDyn$mTheta over
columns 1 through T, excluding the next-period forecast in column T+1.
This replaces the mapped fixed point previously labeled Unconditional Parameters.
It is computed from the stored filtered path, so existing fitted objects
also display the new summary after loading this package version.

## Conditional inference for nearly inactive dynamics

UniGASFit and MultiGASFit accept SE.A.tol = 1e-5. Only if the full OPG
fails, estimated a coefficients with absolute value <= SE.A.tol on the reported
scale are excluded together with their matching b coefficients. Set SE.A.tol=0
to disable this fallback. The principal OPG submatrix of the other coefficients
must pass the same information checks before it is inverted. Point estimates
and log likelihood are unchanged; excluded coefficients have NA standard
errors, test statistics and p-values. Other standard errors are conditional
on the excluded fitted coefficients, not a full joint covariance estimator.
No reparameterization, refit, generalized inverse or eigenvalue repair is used.

The optimiser retains the full scores and opg. Its inference_parameters and
excluded_parameters record the selection. covariance is a full-sized matrix
with the reduced inverse in the retained block and NA in excluded rows/columns.
Diagnostics distinguish full_information_ok from information_ok and record
conditional_inference and se_a_tol. Convergence and gradient checks still apply.
Fit displays explicitly mark conditional standard errors. ConfidenceBands uses
the same retained block and holds excluded coefficients fixed, with a warning.

## Multivariate inverse scaling

Both `mvnorm` and `mvt` now implement `Inv` and `InvSqrt` using the full
expected Fisher information, including Student-t shape cross terms. These
options were previously silently treated as Identity. Refit old multivariate
objects bearing these labels. See `Multi-scaling.md` for formulas and the
Cholesky convention. Coefficient inference continues to use OPG scores.

## Accurate optimizer refinement

The final BFGS retry now supplies Richardson numerical derivatives through
`numDeriv::grad` and uses `reltol = 1e-14`. The default `optim` difference
step of 0.001 was too coarse for small intercepts in persistent multivariate
models. Among numerically tied likelihoods, candidate selection now prefers
a solution that passes both the solver and gradient checks. The final
gradient tolerance remains 1e-5; standard errors still require a valid OPG.

## Readable time axes

Shared time-series plot helpers use sparse, calendar-aligned date ticks
selected for the visible panel width. This removes the daily and fractional-day
tick marks that formed a black band when plotting long xts series. Numeric
observation indices use sparse numeric ticks. The change applies to fitted
series, filtered parameters and comparison plots, and to forecast plots using
the same helpers. Existing fit objects can be plotted after installing the
update; no re-estimation is needed.
