# Diagnostics for the summed negative log likelihood on the optimizer scale.
GAS_sample_size <- function(data, GASSpec) {
  if (is.list(GASSpec) && !is.null(GASSpec$iT)) return(GASSpec$iT)
  if (is(GASSpec, "mGASSpec")) return(ncol(data))
  length(data)
}

GAS_hessian_info <- function(H, p = if (is.matrix(H)) nrow(H) else 0L,
                             tol = 1e-8) {
  bad <- function(reason) list(ok = FALSE, reason = reason, matrix = H)
  if (!is.matrix(H) || !is.numeric(H) || !identical(dim(H), c(as.integer(p), as.integer(p))) ||
      p < 1L || any(!is.finite(H))) return(bad("missing, malformed or non-finite Hessian"))
  if (max(abs(H-t(H))) > 1e-7*max(1, max(abs(H))))
    return(bad("Hessian is not numerically symmetric"))
  H <- (H+t(H))/2
  if (any(diag(H) <= 0)) return(bad("Hessian has nonpositive diagonal entries"))
  # Assess conditioning after removing arbitrary diagonal parameter scales.
  d <- sqrt(diag(H)); C <- H / outer(d,d)
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= tol*max(ev)) return(bad("Hessian is indefinite or numerically singular"))
  U <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(U)) return(bad("Hessian Cholesky factorization failed"))
  list(ok = TRUE, reason = "", matrix = H, chol = U,
       scaled_condition = max(ev)/min(ev))
}

GAS_checked_covariance <- function(H, p) {
  info <- GAS_hessian_info(H, p)
  if (!info$ok) return(NULL)
  V <- chol2inv(info$chol)
  if (any(!is.finite(V))) return(NULL)
  (V+t(V))/2
}

GAS_objective_gradient <- function(FUN, par, data, GASSpec) {
  nm <- names(par); n <- GAS_sample_size(data, GASSpec)
  fn <- function(x) {
    names(x) <- nm
    FUN(x, data = data, GASSpec = GASSpec)/n
  }
  tryCatch(as.numeric(numDeriv::grad(fn, par)),
           error = function(e) rep(NA_real_, length(par)))
}

GAS_finalize_optimiser <- function(optimiser, par0, data, GASSpec, FUN,
                                  Compute.SE, grad_tol = 1e-5) {
  flag <- optimiser$convergence
  if (is.null(flag)) {
    flag <- NA_real_
    warning("The optimizer did not report convergence; status is unknown.", call. = FALSE)
  } else if (!is.numeric(flag) || length(flag) != 1L || !is.finite(flag)) {
    stop("The optimizer's convergence flag must be NULL or one finite number.")
  }
  p <- optimiser$pars
  if (!is.numeric(p) || length(p) != length(par0) || any(!is.finite(p)))
    stop("The optimizer returned invalid parameters.")
  # The optimizer contract preserves order; restore names if it discarded them.
  if (!is.null(names(p)) && !identical(names(p), names(par0)))
    stop("The optimizer changed the parameter names or order.")
  names(p) <- names(par0)
  f <- FUN(p, data = data, GASSpec = GASSpec)
  f0 <- FUN(par0, data = data, GASSpec = GASSpec)
  if (length(f) != 1L || !is.finite(f) || f == 1e10)
    stop("The returned solution has an invalid or penalized likelihood.")
  if (is.finite(f0) && f > f0 + 1e-8*(1+abs(f0)))
    stop("The optimizer returned a worse likelihood than the starting values.")
  g <- GAS_objective_gradient(FUN, p, data, GASSpec)
  stationary <- all(is.finite(g)) && max(abs(g)) <= grad_tol
  optimiser$solver_convergence <- flag
  if (isTRUE(flag == 0) && !stationary) flag <- 52L
  if (!isTRUE(flag == 0))
    warning("Estimation did not satisfy the convergence/gradient checks; inference is unavailable.", call. = FALSE)
  optimiser$convergence <- flag
  optimiser$pars <- p
  optimiser$value <- f
  H <- optimiser$hessian
  if (Compute.SE) {
    # Recompute curvature at the FINAL estimate, never at a one-step refit.
    obj <- function(x) {names(x) <- names(p); FUN(x, data = data, GASSpec = GASSpec)}
    H <- tryCatch(numDeriv::hessian(obj, p), error = function(e) NULL)
    if (is.null(H) || any(!is.finite(H)))
      H <- tryCatch(stats::optimHess(p, obj), error = function(e) NULL)
  }
  info <- GAS_hessian_info(H, length(p))
  optimiser$hessian <- info$matrix
  optimiser$diagnostics <- list(gradient_per_observation = g,
    max_gradient_per_observation = if (all(is.finite(g))) max(abs(g)) else NA_real_,
    gradient_ok = stationary, hessian_ok = info$ok, hessian_reason = info$reason,
    inference_ok = isTRUE(flag == 0) && stationary && info$ok)
  if (Compute.SE && !info$ok)
    warning(paste0("Standard errors are unavailable: ", info$reason, "."), call. = FALSE)
  optimiser
}
