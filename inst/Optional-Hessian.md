# Optional Hessian computation

`UniGASFit()` and `MultiGASFit()` now accept `Compute.Hessian = FALSE`.
The argument is appended after `SE.A.tol` to preserve positional calls.

```r
# Default: omit numerical Hessians, retain OPG standard errors.
Fit <- UniGASFit(GASSpec, data, Compute.Hessian = FALSE)

# Restore numerical curvature diagnostics and the final Hessian.
Fit <- UniGASFit(GASSpec, data, Compute.Hessian = TRUE)

# Omit both Hessians and OPG standard errors.
Fit <- UniGASFit(GASSpec, data, Compute.Hessian = FALSE, Compute.SE = FALSE)
```

The same options apply to `MultiGASFit`. Rolling fits forward
`Compute.Hessian` through their existing `...` argument.

With the default `fn.optim`, FALSE skips both `stats::optimHess()` calls:
the intermediate curvature check and the final numerical Hessian. The
returned `hessian` is NULL. Both the static initialization and final dynamic
optimization receive the setting. Convergence and gradient checks, together
with retries/refinement triggered by those checks, remain active. BFGS still
maintains its internal inverse-Hessian approximation as part of optimization.

TRUE restores the previous numerical curvature checks and returned Hessian.
Removing the curvature-based retry trigger can change the optimization path
when a candidate passes gradient checks but has non-positive curvature.
Enable this option when that extra diagnostic is wanted.

`Compute.SE` remains independent: OPG standard errors and confidence bands
do not require a Hessian. The Hessian, when requested, refers to the summed
negative log likelihood on the unrestricted optimizer parameter scale and
is available at `Fit@Estimates$optimiser$hessian`. A failed numerical Hessian
evaluation yields NULL.

Direct calls to `fn.optim` also accept `Compute.Hessian = FALSE`.
For `fn.solnp`, the argument only controls whether SOLNP's internally
maintained approximation is returned; its internal quasi-Newton updates
cannot be disabled by this option. No extra numerical Hessian is computed.

Custom optimizers that explicitly declare `Compute.Hessian` receive the
setting. Original four-argument custom optimizers continue to work and
retain control over their own internal calculations.
