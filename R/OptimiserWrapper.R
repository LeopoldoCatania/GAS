fn.solnp <- function(par0, data, GASSpec, FUN) {

  optimiser = suppressWarnings(solnp(par0, FUN, data = data, GASSpec = GASSpec, control = list(trace = 0)))

  out = list(pars = optimiser$pars,
             value = tail(optimiser$values, 1),
             hessian = optimiser$hessian,
             convergence = optimiser$convergence)

  return(out)

}

fn.nloptr = function(par0, data, GASSpec, FUN) {
  ctr = list(algorithm = "NLOPT_LN_SBPLX",
             local_opts = list(algorithm = NULL,
                               ftol_abs = 1e-8,
                               xtol_rel = 1e-6,
                               maxeval = 25000,
                               print_level = 1))
  fit = nloptr::nloptr(x0 = par0, eval_f = FUN, data = data, GASSpec = GASSpec, opts = ctr)
  out = list(pars = fit$solution, value = fit$objective, hessian = NULL, convergence = fit$status)
  return(out)
}
