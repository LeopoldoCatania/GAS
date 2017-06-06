fn.solnp <- function(par0, data, GASSpec, FUN) {

  solver.ctr <- list(trace = 0, rho = 1, outer.iter = 400, inner.iter = 800, delta = 1e-07, tol = 1e-08)

  # optimiser = suppressWarnings(solnp(par0, FUN, data = data, GASSpec = GASSpec, control = solver.ctr))
  optimiser = suppressWarnings(solnp(par0, FUN, data = data, GASSpec = GASSpec))

  out = list(pars = optimiser$pars,
             value = tail(optimiser$values, 1),
             hessian = optimiser$hessian,
             convergence = optimiser$convergence)

  return(out)

}

fn.optim <- function(par0, data, GASSpec, FUN) {

  solver.ctr <- list(trace = 0, abstol = 1e-8, reltol = 1e-8)

  optimiser = suppressWarnings(optim(par0, FUN, data = data, GASSpec = GASSpec,
                                     method = "BFGS",
                                     control = solver.ctr,  hessian = TRUE))

  out = list(pars = optimiser$par,
             value = optimiser$value,
             hessian = optimiser$hessian,
             convergence = optimiser$convergence)

  return(out)

}
