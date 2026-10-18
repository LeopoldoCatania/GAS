# Diagnostics for the summed negative log likelihood on the optimizer scale.
GAS_sample_size <- function(data, GASSpec) {
  if (is.list(GASSpec) && !is.null(GASSpec$iT)) return(GASSpec$iT)
  if (is(GASSpec, "mGASSpec")) return(ncol(data))
  length(data)
}

GAS_information_info <- function(H, p = if (is.matrix(H)) nrow(H) else 0L,
                                 tol = 1e-8, label = "OPG information matrix") {
  bad <- function(reason) list(ok = FALSE, reason = paste(label, reason), matrix = H)
  if (!is.matrix(H) || !is.numeric(H) || !identical(dim(H), c(as.integer(p), as.integer(p))) ||
      p < 1L || any(!is.finite(H))) return(bad("is missing, malformed or non-finite"))
  if (max(abs(H-t(H))) > 1e-7*max(1, max(abs(H))))
    return(bad("is not numerically symmetric"))
  H <- (H+t(H))/2
  if (any(diag(H) <= 0)) return(bad("has nonpositive diagonal entries"))
  # Assess conditioning after removing arbitrary diagonal parameter scales.
  d <- sqrt(diag(H)); C <- H / outer(d,d)
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= tol*max(ev)) return(bad("is indefinite or numerically singular"))
  U <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(U)) return(bad("Cholesky factorization failed"))
  list(ok = TRUE, reason = "", matrix = H, chol = U,
       scaled_condition = max(ev)/min(ev))
}

# Curvature is still used by fn.optim to decide whether to retry optimization;
# it is never used to construct the fitted model's covariance estimate.
GAS_hessian_info <- function(H, p = if (is.matrix(H)) nrow(H) else 0L,
                             tol = 1e-8) {
  GAS_information_info(H, p, tol, label = "Hessian")
}

GAS_checked_covariance <- function(H, p) {
  info <- GAS_information_info(H, p)
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

# Rows of S are observation scores and columns are optimizer parameters.
# crossprod(S) is SUM_t s_t s_t'; its inverse estimates Var(p_hat), with
# no additional factor of T. The asymptotic covariance of sqrt(T)(p_hat-p)
# is instead the inverse of crossprod(S)/T.
GAS_opg <- function(par, data, GASSpec) {
  nm <- names(par)
  contributions <- function(x) {
    names(x) <- nm
    GAS_loglik_contributions(x, data, GASSpec)
  }
  S <- numDeriv::jacobian(contributions, par)
  n <- GAS_sample_size(data, GASSpec)
  if (!is.matrix(S) || !identical(dim(S), c(as.integer(n), as.integer(length(par)))) ||
      any(!is.finite(S)))
    stop("The likelihood score Jacobian is malformed or non-finite.")
  colnames(S) <- nm
  list(scores = S, opg = crossprod(S))
}

GAS_validate_se_a_tol <- function(tol) {
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("SE.A.tol must be one finite, nonnegative number.")
  invisible(tol)
}

# Fall back only when the full OPG fails. Fix near-zero score loadings and
# their matching persistence coefficients at their fitted values. Invert the
# principal submatrix for the remaining parameters, not a generalized inverse.
GAS_opg_inference <- function(opg, par, GASSpec, SE.A.tol = 1e-5) {
  GAS_validate_se_a_tol(SE.A.tol)
  nm <- names(par)
  full_info <- GAS_information_info(opg, length(par))
  info <- full_info
  included <- nm
  excluded <- character()
  if (!full_info$ok && SE.A.tol > 0 && is.matrix(opg) &&
      identical(dim(opg), c(length(par), length(par))) && all(is.finite(opg))) {
    Dist <- getDist(GASSpec)
    if (is(GASSpec, "uGASSpec")) {
      an <- paste0("a", seq_len(NumberParameters(Dist)))
    } else if (is(GASSpec, "mGASSpec")) {
      # Names encode scalar or diagonal grouping, without expanding state indices.
      an <- grep("^a\\.", nm, value = TRUE)
    } else {
      stop("OPG inference requires a univariate or multivariate GAS specification.")
    }
    an <- intersect(an, nm)
    if (length(an)) {
      a <- as.numeric(Map_Vec(par[an], LowerA(), UpperA()))
      near_zero <- an[is.finite(a) & abs(a) <= SE.A.tol]
      excluded <- nm[nm %in% c(near_zero, sub("^a", "b", near_zero))]
      included <- nm[!nm %in% excluded]
      if (length(excluded) && length(included)) {
        ix <- match(included, nm)
        info <- GAS_information_info(opg[ix, ix, drop = FALSE], length(ix),
                                      label = "Reduced OPG information matrix")
      }
    }
  }
  covariance <- NULL
  if (info$ok) {
    ix <- match(included, nm)
    reduced_covariance <- GAS_checked_covariance(opg[ix, ix, drop = FALSE], length(ix))
    if (is.null(reduced_covariance)) {
      info <- list(ok = FALSE, reason = "The inverse OPG is non-finite.")
    } else {
      covariance <- matrix(NA_real_, length(par), length(par), dimnames = list(nm, nm))
      covariance[ix, ix] <- reduced_covariance
    }
  }
  list(info = info, full_info = full_info, covariance = covariance,
       included = included, excluded = excluded,
       conditional = info$ok && length(excluded) > 0L)
}

GAS_finalize_optimiser <- function(optimiser, par0, data, GASSpec, FUN,
                                  Compute.SE, grad_tol = 1e-5, SE.A.tol = 1e-5) {
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
  # A custom optimizer may provide a Hessian; retain it as optimizer output
  # only. Inference and confidence bands use the OPG computed here.
  optimiser$scores <- NULL
  optimiser$opg <- NULL
  optimiser$covariance <- NULL
  optimiser$covariance_method <- "OPG"
  optimiser$inference_parameters <- names(p)
  optimiser$excluded_parameters <- character()
  info <- list(ok = FALSE, reason = "Standard errors were not requested.")
  full_info <- info
  conditional <- FALSE
  if (Compute.SE) {
    result <- tryCatch(GAS_opg(p, data, GASSpec), error = function(e) e)
    if (inherits(result, "error")) {
      info <- list(ok = FALSE, reason = paste("OPG calculation failed:", conditionMessage(result)))
      full_info <- info
    } else {
      optimiser$scores <- result$scores
      optimiser$opg <- result$opg
      selection <- GAS_opg_inference(result$opg, p, GASSpec, SE.A.tol)
      info <- selection$info
      full_info <- selection$full_info
      conditional <- selection$conditional
      optimiser$covariance <- selection$covariance
      optimiser$inference_parameters <- selection$included
      optimiser$excluded_parameters <- selection$excluded
    }
  }
  optimiser$diagnostics <- list(gradient_per_observation = g,
    max_gradient_per_observation = if (all(is.finite(g))) max(abs(g)) else NA_real_,
    gradient_ok = stationary, information_ok = info$ok, information_reason = info$reason,
    full_information_ok = full_info$ok, full_information_reason = full_info$reason,
    conditional_inference = conditional, se_a_tol = SE.A.tol,
    inference_ok = Compute.SE && isTRUE(flag == 0) && stationary && info$ok)
  if (isTRUE(optimiser$diagnostics$inference_ok) && conditional)
    warning(paste0("Conditional OPG standard errors: holding ",
      paste(optimiser$excluded_parameters, collapse = ", "),
      " at their fitted values; their standard errors are NA (SE.A.tol = ",
      format(SE.A.tol), ")."), call. = FALSE)
  if (Compute.SE && !info$ok)
    warning(paste0("Standard errors are unavailable: ", info$reason, "."), call. = FALSE)
  optimiser
}
