
InferenceFun_Uni <- function(mHessian, vPw, iK, Compute.SE = TRUE) {

    vPn = vPw2vPn_Uni(vPw, iK)
    iK_s = length(vPn)

    out = matrix(NA, iK_s, 4L, dimnames = list(names(vPn), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))

    out[, "Estimate"] = vPn

    if (Compute.SE) {
      mJacob = jacobian(vPw2vPn_Uni, vPw, iK = iK)
      mInvHessian = GAS_checked_covariance(mHessian, length(vPw))
      if (is.null(mInvHessian)) {
        warning("Hessian is unsuitable for inference; returning estimates with NA standard errors.", call. = FALSE)
        return(out)
      }
      mSandwitch = mJacob %*% mInvHessian %*% t(mJacob)

      if (any(!is.finite(mSandwitch)) || any(diag(mSandwitch) <= 0)) {
        warning("Delta-method covariance is unsuitable for inference; returning NA standard errors.", call. = FALSE)
        return(out)
      }
      vSE = sqrt(diag(mSandwitch))
      vTest = vPn/vSE
      vPvalues = 2 * pnorm(abs(vTest), lower.tail = FALSE)


      out[, "Std. Error"] = vSE
      out[, "t value"] = vTest
      out[, "Pr(>|t|)"] = vPvalues
    }

    return(out)
}

InferenceFun_Multi <- function(mHessian, Dist, vPw, iK, iN, ScalarParameters, Compute.SE = TRUE) {

    vPn = vPw2vPn_Multi(vPw, Dist, iK, iN, ScalarParameters)
    iK_s = length(vPn)

    out = matrix(NA, iK_s, 4L, dimnames = list(names(vPn), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))

    out[, "Estimate"] = vPn

    if (Compute.SE) {

      mJacob = jacobian(vPw2vPn_Multi, vPw, iK = iK, iN = iN, Dist = Dist, ScalarParameters = ScalarParameters)
      mInvHessian = GAS_checked_covariance(mHessian, length(vPw))
      if (is.null(mInvHessian)) {
        warning("Hessian is unsuitable for inference; returning estimates with NA standard errors.", call. = FALSE)
        return(out)
      }
      mSandwitch = mJacob %*% mInvHessian %*% t(mJacob)

      if (any(!is.finite(mSandwitch)) || any(diag(mSandwitch) <= 0)) {
        warning("Delta-method covariance is unsuitable for inference; returning NA standard errors.", call. = FALSE)
        return(out)
      }
      vSE = sqrt(diag(mSandwitch))
      vTest = vPn/vSE
      vPvalues = 2 * pnorm(abs(vTest), lower.tail = FALSE)

      out[, "Pr(>|t|)"] = vPvalues
      out[, "Std. Error"] = vSE
      out[, "t value"] = vTest

    }

    return(out)
}

ICfun <- function(llk, np, iT) {
    AIC <- -2 * (llk - np)
    BIC <- -2 * llk + np * log(iT)

    IC = c(AIC = AIC, BIC = BIC, np = np, llk = llk)
    return(IC)
}
