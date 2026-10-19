GAS_validate_compute_hessian <- function(Compute.Hessian) {
  if (!is.logical(Compute.Hessian) || length(Compute.Hessian) != 1L ||
      is.na(Compute.Hessian))
    stop("Compute.Hessian must be TRUE or FALSE.", call. = FALSE)
  invisible(Compute.Hessian)
}

# Bind once for both initialization and the dynamic fit. Legacy optimizers
# retain their four-argument interface; opt-in optimizers declare the option.
GAS_optimizer_hessian <- function(fn.optimizer, Compute.Hessian) {
  GAS_validate_compute_hessian(Compute.Hessian)
  if (!"Compute.Hessian" %in% names(formals(fn.optimizer))) return(fn.optimizer)
  force(fn.optimizer)
  force(Compute.Hessian)
  function(par0, data, GASSpec, FUN)
    fn.optimizer(par0 = par0, data = data, GASSpec = GASSpec, FUN = FUN,
                 Compute.Hessian = Compute.Hessian)
}

fn.solnp <- function(par0, data, GASSpec, FUN, Compute.Hessian = FALSE) {

  GAS_validate_compute_hessian(Compute.Hessian)

  solver.ctr <- list(trace = 0, rho = 1, outer.iter = 400, inner.iter = 800, delta = 1e-07, tol = 1e-08)

  optimiser = suppressWarnings(solnp(par0, FUN, data = data, GASSpec = GASSpec, control = solver.ctr))

  out = list(pars = optimiser$pars,
             value = tail(optimiser$values, 1),
             # SOLNP maintains this approximation internally in either mode.
             hessian = if (Compute.Hessian) optimiser$hessian else NULL,
             convergence = optimiser$convergence)

  return(out)

}

fn.optim <- function(par0, data, GASSpec, FUN, Compute.Hessian = FALSE) {
  GAS_validate_compute_hessian(Compute.Hessian)
  nm <- names(par0)
  n <- GAS_sample_size(data, GASSpec)
  objective <- function(x) {
    names(x) <- nm
    FUN(x, data = data, GASSpec = GASSpec)/n
  }
  value0 <- objective(par0)
  if (!is.finite(value0)) stop("The starting objective is non-finite.")
  # Retain the start as a safeguard, without labeling it converged.
  candidates <- list(list(par = par0, value = value0, convergence = 1L,
                          method = "starting values"))
  add <- function(z, method) {
    if (!is.null(z) && length(z$par) == length(par0) && all(is.finite(z$par))) {
      val <- tryCatch(objective(z$par), error = function(e) Inf)
      if (is.finite(val)) candidates[[length(candidates)+1L]] <<-
        list(par = z$par, value = val, convergence = z$convergence, method = method)
    }
  }
  first <- tryCatch(stats::optim(par0, objective, method = "BFGS",
    control = list(maxit = 2000L, reltol = 1e-11, abstol = -Inf)), error = function(e) NULL)
  add(first, "BFGS")
  best <- candidates[[which.min(vapply(candidates, function(x) x$value, numeric(1)))]]
  g <- GAS_objective_gradient(FUN, stats::setNames(best$par, nm), data, GASSpec)
  curvature_failed <- FALSE
  if (Compute.Hessian) {
    Hmean <- tryCatch(stats::optimHess(best$par, objective), error = function(e) NULL)
    curvature_failed <- !GAS_hessian_info(Hmean, length(par0))$ok
  }
  needs_retry <- !isTRUE(best$convergence == 0) || any(!is.finite(g)) ||
    max(abs(g)) > 1e-5 || curvature_failed
  if (needs_retry) {
    # A deterministic alternative solver; always retain the best objective.
    second <- tryCatch(stats::nlminb(best$par, objective,
      control = list(iter.max = 2000L, eval.max = 5000L, rel.tol = 1e-11)),
      error = function(e) NULL)
    add(second, "nlminb")
    if (!is.null(second)) {
      # optim's default 1e-3 differences can be too coarse for small GAS
      # intercepts, especially with high persistence. Refine with the same
      # Richardson derivatives used in the final stationarity check.
      gradient <- function(x) GAS_objective_gradient(FUN,
        stats::setNames(x, nm), data, GASSpec)
      third <- tryCatch(stats::optim(second$par, objective, gr = gradient,
        method = "BFGS", control = list(maxit = 2000L, reltol = 1e-14,
                                       abstol = -Inf)), error = function(e) NULL)
      add(third, "BFGS gradient refinement")
    }
  }
  vals <- vapply(candidates, function(x) x$value, numeric(1))
  # Tiny likelihood rounding differences must not discard a stationary
  # solution in favor of one that only has a successful solver return code.
  ix <- which(vals <= min(vals) + 1e-12*(1+abs(min(vals))))
  converged <- ix[vapply(candidates[ix], function(x) isTRUE(x$convergence == 0), logical(1))]
  if (length(converged)) {
    ix <- converged
    stationary <- vapply(candidates[ix], function(x) {
      g <- GAS_objective_gradient(FUN, stats::setNames(x$par, nm), data, GASSpec)
      all(is.finite(g)) && max(abs(g)) <= 1e-5
    }, logical(1))
    if (any(stationary)) ix <- ix[stationary]
  }
  best <- candidates[[ix[which.min(vals[ix])]]]
  names(best$par) <- nm
  # Return TOTAL objective and optimizer curvature for compatibility.
  # Fitted-model standard errors are constructed separately from the OPG.
  H <- if (Compute.Hessian)
    tryCatch(n*stats::optimHess(best$par, objective), error = function(e) NULL) else NULL
  list(pars = best$par, value = n*best$value, hessian = H,
       convergence = best$convergence, method = best$method)
}
