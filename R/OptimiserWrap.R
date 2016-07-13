StaticLLKoptimizer_Uni <- function(vTheta_tilde, vY, Dist, iT, iK) {
    vTheta = MapParameters_univ(vTheta_tilde, Dist, iK)
    dLLK = StaticLLK_Univ(vY, vTheta, iT, Dist)
    
    if (is.na(dLLK)) {
        dLLK = -1e+50
    }
    return(-dLLK)
}
StaticLLKoptimizer_Multi <- function(vTheta_tilde, mY, Dist, iT, iN, iK) {
    vTheta = MapParameters_multi(vTheta_tilde, Dist, iN, iK)
    dLLK = StaticLLK_Multi(mY, vTheta, iT, iN, Dist)
    
    if (is.na(dLLK)) {
        dLLK = -1e+50
    }
    return(-dLLK)
}
UniGASOptimiser <- function(vPw, vY, Dist, ScalingType, iT, iK) {
    lParList = vPw2lPn_Uni(vPw, iK)
    
    lParList = AddFixedPar(lParList)
    
    dLLK = try(GASFilter_univ(vY, lParList$vKappa, lParList$mA, lParList$mB, iT, iK, Dist, ScalingType)$dLLK, silent = T)
    
    if (!is(dLLK, "try-error")) {
        dMLLK = -dLLK
    } else {
        dMLLK = 1e+50
    }
    return(dMLLK)
}
MultiGASOptimiser <- function(vPw, mY, Dist, ScalingType, iT, iN, iK, ScalarParameters) {
    
    lParList = vPw2lPn_Multi(vPw, Dist, iK, iN, ScalarParameters)
    lParList = AddFixedPar(lParList)
    
    dLLK = try(GASFilter_multi(mY, lParList$vKappa, lParList$mA, lParList$mB, iT, iN, iK, Dist, ScalingType)$dLLK, silent = T)
    
    if (!is(dLLK, "try-error")) {
        dMLLK = -dLLK
    } else {
        dMLLK = 1e+50
    }
    return(dMLLK)
}
