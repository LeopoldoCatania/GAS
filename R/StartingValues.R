StartingNu <- function(vY) {
  dStart = c(unmapVec_C(12, LowerNu(), UpperNu()))
    dNu = try(exp(solnp(dStart, function(x, vY) {
        dNu =  c(Map_Vec(x, LowerNu(), UpperNu()))
        - sum(dt(vY, dNu, log = TRUE))
    }, vY = vY, control = list(trace = 0))$pars) + LowerNu(), silent = TRUE)

    if (is(dNu, "try-error"))
        dNu = 7

    if (dNu > UpperNu())
        dNu = UpperNu() - 1
    if (dNu <= LowerNu())
        dNu = LowerNu() + 0.02

    return(dNu)
}

StaticStarting_Uni <- function(vY, Dist, iK) {

    if (Dist == "std") {

        dMu = mean(vY)
        dNu = 8
        dPhi2 = var(vY) * (dNu - 2)/dNu

        vTheta = c(dMu, dPhi2, dNu)

    }

    if (Dist == "sstd") {

        dMu = mean(vY)
        dNu = 8
        dSigma = sd(vY)
        dXi = 1

        vTheta = c(dMu, dSigma, dXi, dNu)

    }

    if (Dist == "norm") {

        dMu = mean(vY)
        dSigma2 = var(vY)

        vTheta = c(dMu, dSigma2)

    }
    if (Dist == "snorm") {

        dMu = mean(vY)
        dSigma = sd(vY)
        dXi = 1

        vTheta = c(dMu, dSigma, dXi)

    }
    if (Dist == "ast" | Dist == "ast1") {

        dNu1 = StartingNu(vY)

        if (Dist == "ast") {

            dNu2 = dNu1 * 1.5

            if (dNu2 > UpperNu()) {
              dNu2 = UpperNu() - 1
            }

        } else {

            dNu2 = dNu1

        }

        dAlpha = 0.5
        dMu = mean(vY)
        dSigma = sd(vY)

        if (Dist == "ast") {

            vTheta = c(dMu, dSigma, dAlpha, dNu1, dNu2)

        } else {

            vTheta = c(dMu, dSigma, dAlpha, dNu1)

        }
    }
    if (Dist == "poi") {

        vTheta = mean(vY)

    }
    if (Dist == "ber") {

        dPi = sum(vY)/length(vY)

        vTheta = c(dPi)

    }
    if (Dist == "gamma") {

        dMean = mean(vY)
        dSigma2 = var(vY)

        dBeta = dMean/dSigma2
        dAlpha = dMean^2/dSigma2

        vTheta = c(dAlpha, dBeta)

    }

    if (Dist == "exp") {

        vTheta = c(1/mean(vY))

    }
    if (Dist == "beta") {

        dMean = mean(vY)
        dSigma2 = var(vY)

        dAlpha = dMean * (dMean * (1 - dMean)/dSigma2 - 1)
        dBeta = (1 - dMean) * (dMean * (1 - dMean)/dSigma2 - 1)

        vTheta = c(dAlpha, dBeta)
    }
    if (Dist == "ald") {

        dTheta = mean(vY)
        dSigma = sd(vY)
        dKappa = 1

        vTheta = c(dTheta, dSigma, dKappa)

    }

    vTheta_tilde = as.numeric(UnmapParameters_univ(vTheta, Dist, iK = iK))

    names(vTheta_tilde) = FullNamesUni(Dist)

    return(vTheta_tilde)
}

StaticStarting_Multi <- function(mY, Dist, iN) {

    vEmpRho = build_vR(cor(t(mY)), iN)
    vEmpPhi = UnMapR_C(vEmpRho, iN)

    vEmpMu = apply(mY, 1, mean)
    vEmpSigma = apply(mY, 1, sd)

    if (Dist == "mvnorm") {
        vMu_tilde = vEmpMu
        vSigma_tilde = log(vEmpSigma)
        vPw = c(vMu_tilde, vSigma_tilde, vEmpPhi)
    }
    if (Dist == "mvt") {

        dNu = StartingNu(c(mY))

        vMu_tilde = vEmpMu
        vSigma_tilde = log(vEmpSigma * (dNu - 2)/dNu)
        vPw = c(vMu_tilde, vSigma_tilde, vEmpPhi, log(dNu - LowerNu()))
    }

    return(as.numeric(vPw))

}

UniGAS_Starting <- function(vY, iT, iK, Dist, ScalingType, GASPar) {

    StaticFit = StaticMLFIT(vY, Dist)
    vUncValues = StaticFit$optimiser$pars
    names(vUncValues) = paste("kappa", 1:iK, sep = "")

    if (iK > 1) {
        vA = starting_vA_Uni(vY, vUncValues, mB = diag(rep(0.9, iK)), dA_foo = 1e-06, iT, iK, Dist,
            ScalingType = ScalingType, GASPar)
        vB = starting_vB_Uni(vY, vUncValues, dB_foo = 0.9, mA = diag(vA), iT, iK, Dist, ScalingType = ScalingType,
            GASPar)
        vKappa = (diag(iK) - diag(vB)) %*% vUncValues
        names(vKappa) = paste("kappa", 1:iK, sep = "")
    } else {
        vA = starting_vA_Uni(vY, vUncValues, mB = matrix(0.9, iK, iK), dA_foo = 1e-06, iT, iK, Dist,
            ScalingType = ScalingType, GASPar)
        vB = starting_vB_Uni(vY, vUncValues, dB_foo = 0.9, mA = matrix(vA, iK, iK), iT, iK, Dist, ScalingType = ScalingType,
            GASPar)
        vKappa = (1 - vB) * vUncValues
        names(vKappa) = paste("kappa", 1:iK, sep = "")
    }
    vA = unmapVec_C(vA, LowerA(), UpperA())
    names(vA) = paste("a", 1:iK, sep = "")
    vB = unmapVec_C(vB, LowerB(), UpperB())
    names(vB) = paste("b", 1:iK, sep = "")

    return(list(vPw = c(vKappa, vA, vB), StaticFit = StaticFit))

}

starting_vA_Uni <- function(vY, vUncValues, mB, dA_foo, iT, iK, Dist, ScalingType, GASPar) {

    seq_alpha = c(seq(1e-04, 1.5, length.out = 30))

    vKappa = (diag(iK) - mB) %*% vUncValues

    mA = matrix(0, iK, iK)

    diag(mA)[unlist(GASPar)] = dA_foo

    dAlpha_best = dA_foo

    for (i in 1:iK) {
        if (GASPar[[i]]) {
            for (l in 1:length(seq_alpha)) {
                dLLK_foo = try(GASFilter_univ(vY, vKappa, mA, mB, iT, iK, Dist, ScalingType)$dLLK, silent = TRUE)
                if (is.numeric(dLLK_foo) & !is.nan(dLLK_foo)) {
                  mA[i, i] = seq_alpha[l]
                  dLLK_post = try(GASFilter_univ(vY, vKappa, mA, mB, iT, iK, Dist, ScalingType)$dLLK,
                    silent = TRUE)
                  if (is.numeric(dLLK_post) & !is.nan(dLLK_post)) {
                    if (dLLK_post > dLLK_foo)
                      dAlpha_best = seq_alpha[l]
                  }
                  mA[i, i] = dAlpha_best
                }
            }
        }
    }

    return(diag(mA))
}

starting_vB_Uni <- function(vY, vUncValues, dB_foo, mA, iT, iK, Dist, ScalingType, GASPar) {

    seq_beta = c(seq(0.5, 0.98, length.out = 30))

    mB = matrix(0, iK, iK)

    diag(mB) = dB_foo

    vKappa = (diag(iK) - mB) %*% vUncValues

    dB_best = dB_foo

    for (i in 1:iK) {
        if (GASPar[[i]]) {
            for (l in 1:length(seq_beta)) {
                dLLK_foo = try(GASFilter_univ(vY, vKappa, mA, mB, iT, iK, Dist, ScalingType)$dLLK, silent = TRUE)
                if (is.numeric(dLLK_foo) & !is.nan(dLLK_foo)) {
                  mB[i, i] = seq_beta[l]
                  vKappa = (diag(iK) - mB) %*% vUncValues
                  dLLK_post = try(GASFilter_univ(vY, vKappa, mA, mB, iT, iK, Dist, ScalingType)$dLLK,
                    silent = TRUE)
                  if (is.numeric(dLLK_post) & !is.nan(dLLK_post)) {
                    if (dLLK_post > dLLK_foo)
                      dB_best = seq_beta[l]
                  }
                  mB[i, i] = dB_best
                  vKappa = (diag(iK) - mB) %*% vUncValues
                }
            }
        }
    }

    return(diag(mB))
}

StartingValues_mvnorm <- function(mY, iT, iN, iK, GASPar, ScalingType, ScalarParameters) {

    vEmpRho = build_vR(cor(t(mY)), iN)
    vEmpPhi = UnMapR_C(vEmpRho, iN)

    vEmpMu = apply(mY, 1, mean)
    vEmpSigma = apply(mY, 1, sd)

    vUncValues = c(vEmpMu, log(vEmpSigma), vEmpPhi)
    names(vUncValues) = paste("kappa.", mvnormParNames(iN), sep = "")
    #

    if (ScalarParameters) {

        mA = starting_mA_Multi_Scalars(mY, vUncValues, mB = diag(rep(0.9, iK)), dA_foo = 0.001, iT,
            iK, iN, "mvnorm", ScalingType, GASPar)
        mB = starting_mB_Multi_Scalars(mY, vUncValues, dB_foo = 0.9, mA, iT, iK, iN, "mvnorm", ScalingType,
            GASPar)

        vKappa = (diag(iK) - mB) %*% vUncValues
        names(vKappa) = paste("kappa.", mvnormParNames(iN), sep = "")

        vA = diag(mA)[c(1, iN + 1, 2 * iN + 1, 2 * iN + iN * (iN - 1)/2 + 1)]
        vB = diag(mB)[c(1, iN + 1, 2 * iN + 1, 2 * iN + iN * (iN - 1)/2 + 1)]

        vA = vA[!is.na(vA)]  # nas arise when the distribution do not have shape pars
        vB = vB[!is.na(vB)]  # nas arise when the distribution do not have shape pars

    } else {

        vA = starting_vA_Multi(mY, vUncValues, mB = diag(rep(0.9, iK)), dA_foo = 0.01, iT, iK, iN, "mvnorm",
            ScalingType, GASPar)
        vB = starting_vB_Multi(mY, vUncValues, dB_foo = 0.9, mA = diag(vA), iT, iK, iN, "mvnorm", ScalingType,
            GASPar)

        vKappa = (diag(iK) - diag(vB)) %*% vUncValues
        names(vKappa) = paste("kappa.", mvnormParNames(iN), sep = "")
    }

    pw = c(vKappa, vA, vB)
    return(pw)
}

StartingValues_mvt <- function(mY, iT, iN, iK, GASPar, ScalingType, ScalarParameters) {

    StaticFit = StaticMLFIT_Multiv(mY, "mvt")

    vUncValues = StaticFit$optimiser$pars
    names(vUncValues) = paste("kappa.", mvtParNames(iN), sep = "")

    if (ScalarParameters) {

        mA = starting_mA_Multi_Scalars(mY, vUncValues, mB = diag(rep(0.9, iK)), dA_foo = 0.001, iT,
            iK, iN, "mvt", ScalingType, GASPar)
        mB = starting_mB_Multi_Scalars(mY, vUncValues, dB_foo = 0.9, mA, iT, iK, iN, "mvt", ScalingType,
            GASPar)

        vKappa = (diag(iK) - mB) %*% vUncValues
        names(vKappa) = paste("kappa.", mvtParNames(iN), sep = "")

        vA = diag(mA)[c(1, iN + 1, 2 * iN + 1, 2 * iN + iN * (iN - 1)/2 + 1)]
        vB = diag(mB)[c(1, iN + 1, 2 * iN + 1, 2 * iN + iN * (iN - 1)/2 + 1)]

        vA = vA[!is.na(vA)]  # nas arise when the distribution do not have shape pars
        vB = vB[!is.na(vB)]  # nas arise when the distribution do not have shape pars

    } else {

        vA = starting_vA_Multi(mY, vUncValues, mB = diag(rep(0.9, iK)), dA_foo = 0.01, iT, iK, iN, "mvt",
            ScalingType, GASPar)
        vB = starting_vB_Multi(mY, vUncValues, dB_foo = 0.9, mA = diag(vA), iT, iK, iN, "mvt", ScalingType,
            GASPar)

        vKappa = (diag(iK) - diag(vB)) %*% vUncValues
        names(vKappa) = paste("kappa.", mvtParNames(iN), sep = "")
    }

    vA = unmapVec_C(vA, LowerA(), UpperA())
    names(vA) = paste("a.", mvtParNames(iN, ScalarParameters), sep = "")
    vB = unmapVec_C(vB, LowerB(), UpperB())
    names(vB) = paste("b.", mvtParNames(iN, ScalarParameters), sep = "")

    pw = c(vKappa, vA, vB)

    return(pw)
}

MultiGAS_Starting <- function(mY, iT, iN, iK, Dist, GASPar, ScalingType, ScalarParameters) {

    if (Dist == "mvnorm")
        vPw = StartingValues_mvnorm(mY, iT, iN, iK, GASPar, ScalingType, ScalarParameters)
    if (Dist == "mvt")
        vPw = StartingValues_mvt(mY, iT, iN, iK, GASPar, ScalingType, ScalarParameters)

    return(vPw)
}

starting_vA_Multi <- function(mY, vUncValues, mB, dA_foo, iT, iK, iN, Dist, ScalingType, GASPar) {

    seq_alpha = c(seq(1e-04, 0.5, length.out = 30))

    vKappa = (diag(iK) - mB) %*% vUncValues

    vBool = FixedDynamicPar_Multi(Dist, iN, GASPar)

    mA = matrix(0, iK, iK)
    diag(mA)[vBool] = dA_foo

    dAlpha_best = dA_foo

    for (i in 1:iK) {

        if (vBool[i]) {

            for (l in 1:length(seq_alpha)) {

                dLLK_foo = try(GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK,
                  silent = TRUE)

                if (is.numeric(dLLK_foo) & !is.nan(dLLK_foo)) {

                  mA[i, i] = seq_alpha[l]

                  dLLK_post = try(GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK,
                    silent = TRUE)

                  if (is.numeric(dLLK_post) & !is.nan(dLLK_post)) {

                    if (dLLK_post > dLLK_foo)
                      dAlpha_best = seq_alpha[l]

                  }

                  mA[i, i] = dAlpha_best

                }
            }
        }
    }

    return(diag(mA))
}

starting_mA_Multi_Scalars <- function(mY, vUncValues, mB, dA_foo, iT, iK, iN, Dist, ScalingType, GASPar) {

    seq_alpha = c(seq(1e-04, 0.5, length.out = 30))

    vKappa = (diag(iK) - mB) %*% vUncValues

    vBool = unlist(GASPar)
    iK2 = length(GASPar)

    mA = matrix(0, iK, iK)
    diag(mA)[FixedDynamicPar_Multi(Dist, iN, GASPar)] = dA_foo

    dAlpha_best = dA_foo

    lStructure = MatrixCoefficientStructure_Multi(Dist, iN, iK, GASPar)

    for (i in 1:iK2) {
        if (vBool[i]) {
            for (l in 1:length(seq_alpha)) {
                dLLK_foo = try(GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK,
                  silent = TRUE)
                if (is.numeric(dLLK_foo) & !is.nan(dLLK_foo)) {
                  mA[lStructure[[i]]] = seq_alpha[l]
                  dLLK_post = try(GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK,
                    silent = TRUE)
                  if (is.numeric(dLLK_post) & !is.nan(dLLK_post)) {
                    if (dLLK_post > dLLK_foo)
                      dAlpha_best = seq_alpha[l]
                  }
                  mA[lStructure[[i]]] = dAlpha_best
                }
            }
        }
    }



    return(mA)
}

starting_vB_Multi <- function(mY, vUncValues, dB_foo, mA, iT, iK, iN, Dist, ScalingType, GASPar) {

    seq_beta = c(seq(0.5, 0.98, length.out = 30))
    vBool = FixedDynamicPar_Multi(Dist, iN, GASPar)

    mB = matrix(0, iK, iK)
    diag(mB) = dB_foo

    vKappa = (diag(iK) - mB) %*% vUncValues

    dB_best = dB_foo

    for (i in 1:iK) {
        if (vBool[i]) {
            for (l in 1:length(seq_beta)) {
                dLLK_foo = try(GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK,
                  silent = TRUE)
                if (is.numeric(dLLK_foo) & !is.nan(dLLK_foo)) {
                  mB[i, i] = seq_beta[l]
                  vKappa = (diag(iK) - mB) %*% vUncValues
                  dLLK_post = try(GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK,
                    silent = TRUE)
                  if (is.numeric(dLLK_post) & !is.nan(dLLK_post)) {
                    if (dLLK_post > dLLK_foo)
                      dB_best = seq_beta[l]
                  }
                  mB[i, i] = dB_best
                  vKappa = (diag(iK) - mB) %*% vUncValues
                }
            }
        }
    }

    return(diag(mB))
}

starting_mB_Multi_Scalars <- function(mY, vUncValues, dB_foo, mA, iT, iK, iN, Dist, ScalingType, GASPar) {

    seq_beta = c(seq(0.5, 0.98, length.out = 30))

    vBool = unlist(GASPar)
    iK2 = length(GASPar)

    mB = matrix(0, iK, iK)
    diag(mB)[FixedDynamicPar_Multi(Dist, iN, GASPar)] = dB_foo

    vKappa = (diag(iK) - mB) %*% vUncValues

    dBeta_best = dB_foo

    lStructure = MatrixCoefficientStructure_Multi(Dist, iN, iK, GASPar)

    for (i in 1:iK2) {
        if (vBool[i]) {
            for (l in 1:length(seq_beta)) {
                dLLK_foo = try(GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK,
                  silent = TRUE)
                if (is.numeric(dLLK_foo) & !is.nan(dLLK_foo)) {
                  mB[lStructure[[i]]] = seq_beta[l]
                  vKappa = (diag(iK) - mB) %*% vUncValues
                  dLLK_post = try(GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK,
                    silent = TRUE)
                  if (is.numeric(dLLK_post) & !is.nan(dLLK_post)) {
                    if (dLLK_post > dLLK_foo)
                      dBeta_best = seq_beta[l]
                  }
                  mB[lStructure[[i]]] = dBeta_best
                  vKappa = (diag(iK) - mB) %*% vUncValues
                }
            }
        }
    }

    return(mB)
}

