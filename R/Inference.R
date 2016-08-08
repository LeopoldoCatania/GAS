
InferenceFun_Uni <- function(mHessian, vPw, iK) {

    vPn = vPw2vPn_Uni(vPw, iK)
    iK_s = length(vPn)

    out = matrix(NA, iK_s, 4, dimnames = list(names(vPn), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))

    mJacob = jacobian(vPw2vPn_Uni, vPw, iK = iK)
    mInvHessian = ginv(mHessian)
    mSandwitch = t(mJacob) %*% mInvHessian %*% mJacob

    vSE = sqrt(diag(mSandwitch))
    vTest = vPn/vSE
    vPvalues = 1 - pnorm(abs(vTest))

    out[, "Estimate"] = vPn
    out[, "Std. Error"] = vSE
    out[, "t value"] = vTest
    out[, "Pr(>|t|)"] = vPvalues

    return(out)
}

InferenceFun_Multi <- function(mHessian, Dist, vPw, iK, iN, ScalarParameters) {

    vPn = vPw2vPn_Multi(vPw, Dist, iK, iN, ScalarParameters)
    iK_s = length(vPn)

    out = matrix(NA, iK_s, 4, dimnames = list(names(vPn), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))

    mJacob = jacobian(vPw2vPn_Multi, vPw, iK = iK, iN = iN, Dist = Dist, ScalarParameters = ScalarParameters)
    mInvHessian = ginv(mHessian)
    mSandwitch = t(mJacob) %*% mInvHessian %*% mJacob

    vSE = sqrt(diag(mSandwitch))
    vTest = vPn/vSE
    vPvalues = 1 - pnorm(abs(vTest))

    out[, "Estimate"] = vPn
    out[, "Std. Error"] = vSE
    out[, "t value"] = vTest
    out[, "Pr(>|t|)"] = vPvalues

    return(out)
}

ICfun <- function(llk, np, iT) {
    AIC <- -2 * (llk - np)
    BIC <- -2 * llk + np * log(iT)

    IC = c(AIC = AIC, BIC = BIC, np = np, llk = llk)
    return(IC)
}
