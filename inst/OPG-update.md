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
