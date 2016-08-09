
UniGASSim <- function(iT, kappa, A, B, Dist, ScalingType) {
    
    vKappa = kappa
    mA = A
    mB = B
    
    lSim = SimulateGAS_univ(iT, vKappa, mA, mB, Dist, ScalingType)
    
    iK = NumberParameters(Dist)
    
    mMoments = EvalMoments_univ(lSim$mTheta, Dist)
    
    Sim <- new("uGASSim", ModelInfo = list(iT = iT, iK = iK, vKappa = vKappa, mA = mA, mB = mB, Dist = Dist, 
        ScalingType = ScalingType), GASDyn = lSim, Data = list(vY = lSim$vY, Moments = mMoments))
    return(Sim)
}

MultiGASSim <- function(iT, N, kappa, A, B, Dist, ScalingType) {
    
    vKappa = kappa
    mA = A
    mB = B
    iN = N
    
    lSim = SimulateGAS_multi(iT, iN, vKappa, mA, mB, Dist, ScalingType)
    
    iK = NumberParameters(Dist, iN)
    mY = lSim$mY
    rownames(mY) = paste("Series", 1:iN)
    
    lMoments = EvalMoments_multi(lSim$mTheta, Dist, iN)
    
    Sim <- new("mGASSim", ModelInfo = list(iT = iT, iN = iN, iK = iK, vKappa = vKappa, mA = mA, mB = mB, 
        Dist = Dist, ScalingType = ScalingType), GASDyn = lSim, Data = list(mY = mY, Moments = lMoments))
    return(Sim)
}




