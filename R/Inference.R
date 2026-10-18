# The retained coefficients are differentiated with excluded coefficients held
# fixed. Excluded rows retain their estimates but have NA inferential entries.
GAS_inference_table <- function(mInformation, vPw, map, Compute.SE,
                                 inference_parameters) {
  vPn <- map(vPw)
  out <- matrix(NA_real_, length(vPn), 4L,
    dimnames = list(names(vPn), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
  out[, "Estimate"] <- vPn
  if (!Compute.SE) return(out)
  ix <- match(inference_parameters, names(vPw))
  iy <- which(names(vPn) %in% inference_parameters)
  if (!length(ix) || anyNA(ix) || !length(iy)) return(out)
  V <- GAS_checked_covariance(mInformation[ix, ix, drop = FALSE], length(ix))
  if (is.null(V)) {
    warning("OPG information matrix is unsuitable for inference; returning estimates with NA standard errors.", call. = FALSE)
    return(out)
  }
  transform <- function(x) {
    p <- vPw
    p[ix] <- x
    map(p)[iy]
  }
  J <- numDeriv::jacobian(transform, vPw[ix])
  covariance <- J %*% V %*% t(J)
  if (any(!is.finite(covariance)) || any(diag(covariance) <= 0)) {
    warning("Delta-method covariance is unsuitable for inference; returning NA standard errors.", call. = FALSE)
    return(out)
  }
  se <- sqrt(diag(covariance))
  statistic <- vPn[iy]/se
  out[iy, "Std. Error"] <- se
  out[iy, "t value"] <- statistic
  out[iy, "Pr(>|t|)"] <- 2*pnorm(abs(statistic), lower.tail = FALSE)
  out
}

InferenceFun_Uni <- function(mInformation, vPw, iK, Compute.SE = TRUE,
                              inference_parameters = names(vPw)) {
  GAS_inference_table(mInformation, vPw,
    function(p) vPw2vPn_Uni(p, iK), Compute.SE, inference_parameters)
}

InferenceFun_Multi <- function(mInformation, Dist, vPw, iK, iN, ScalarParameters,
                                Compute.SE = TRUE,
                                inference_parameters = names(vPw)) {
  GAS_inference_table(mInformation, vPw,
    function(p) vPw2vPn_Multi(p, Dist, iK, iN, ScalarParameters),
    Compute.SE, inference_parameters)
}

ICfun <- function(llk, np, iT) {
    AIC <- -2 * (llk - np)
    BIC <- -2 * llk + np * log(iT)

    IC = c(AIC = AIC, BIC = BIC, np = np, llk = llk)
    return(IC)
}
