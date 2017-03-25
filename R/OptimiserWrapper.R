fn.solnp <- function(par0, data, GASSpec, FUN) {

  optimiser = suppressWarnings(solnp(par0, FUN, data = data, GASSpec = GASSpec, control = list(trace = 0)))

  out = list(pars = optimiser$pars,
             value = tail(optimiser$values, 1),
             hessian = optimiser$hessian,
             convergence = optimiser$convergence)

  return(out)

}
