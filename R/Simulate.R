
UniGASSim<-function(iT, vKappa, mA, mB, Dist, ScalingType){
  lSim = SimulateGAS_univ(iT, vKappa, mA, mB, Dist, ScalingType)

  iK = NumberParameters(Dist)

  mMoments = EvalMoments(lSim$mTheta,Dist)

  Sim <- new("uGASSim", ModelInfo=list(iT = iT, iK = iK,vKappa=vKappa, mA = mA, mB = mB, Dist = Dist,ScalingType=ScalingType),
             GASDyn = lSim, Data = list(vY = lSim$vY, Moments=mMoments))
  return(Sim)
}

MultiGASSim<-function(iT, iN,vKappa, mA, mB, Dist, ScalingType){
  lSim = SimulateGAS_multi(iT, iN, vKappa, mA, mB, Dist, ScalingType)

  iK = NumberParameters(Dist,iN)
  mY = lSim$mY ; rownames(mY) = paste("Series",1:iN)

  Sim <- new("mGASSim", ModelInfo=list(iT = iT, iN=iN, iK = iK, vKappa=vKappa, mA = mA, mB = mB, Dist = Dist,ScalingType=ScalingType),
             GASDyn = lSim, Data = list(mY = mY))
  return(Sim)
}




