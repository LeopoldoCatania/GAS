# Multivariate score scaling

`MultiGASSpec` supports `Identity`, `Inv`, and `InvSqrt` for both `mvnorm`
and `mvt`. Filtering, simulation and forecasting use the same C++ update.
Invalid scaling names are rejected instead of silently using Identity.

Let `f` be the unconstrained state, `theta = h(f)`, `J = dh/df`, and
`s = J' d log p(y;theta)/dtheta`. The expected information in state coordinates
is `I_f = J' I_theta J`. The innovations are:

* Identity: `s`.
* Inv: `solve(I_f, s)`.
* InvSqrt: `solve(t(U), s)`, where `U = chol(I_f)` and `I_f = t(U) U`.

The last convention is Cholesky whitening, matching the univariate code; it
is not the symmetric matrix square root. 

## Expected information

Write `Sigma = diag(phi) R diag(phi)` for the dispersion matrix, `n` for the
number of series, `D_j = d Sigma / d theta_j` for each scale/correlation
parameter, and `W_j = solve(Sigma, D_j)`.

For the Gaussian model:

* `I_mu,mu = solve(Sigma)`.
* `I_j,k = tr(W_j W_k)/2`.

For Student-t:

* `I_mu,mu = (nu+n)/(nu+n+2) solve(Sigma)`.
* `I_j,k = ((nu+n) tr(W_j W_k) - tr(W_j) tr(W_k))/(2 (nu+n+2))`.
* `I_j,nu = -tr(W_j)/((nu+n) (nu+n+2))`.
* `I_nu,nu = (trigamma(nu/2)-trigamma((nu+n)/2))/4
  - n (nu+n+4)/(2 nu (nu+n) (nu+n+2))`.

Location is information-orthogonal to scales, correlations and shape.
The formulas apply for `nu > 0`; existence of the data covariance is not
required. The information must be finite and positive definite;
failed factorizations are rejected.

Reference: Lange, Little and Taylor (1989), *Robust Statistical Modeling
Using the t Distribution*, JASA 84, 881–896, Appendix B,
https://doi.org/10.1080/01621459.1989.10478852.
