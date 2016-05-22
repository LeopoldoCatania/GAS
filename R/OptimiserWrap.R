StaticLLKoptimizer<-function(vTheta_tilde, vY,Dist, iT, iK){
  vTheta = MapParameters(vTheta_tilde, Dist, iK)
  dLLK = StaticLLK_Univ(vY, vTheta, iT, Dist)
  if(is.na(dLLK)){
    dLLK = -1e50
  }
  return(-dLLK)
}
UnivGASOptimiser<-function(vPw, vY, Dist, ScalingType, iT, iK){
  lParList = vPw2lPn_Univ(vPw,iK)

  lParList = AddFixedPar(lParList)

  # print(vPw)
  dLLK = try(GASFilter_univ(vY, lParList$vKappa, lParList$mA, lParList$mB, iT, iK, Dist, ScalingType)$dLLK,silent = T)

  if(!is(dLLK,"try-error")){
    dMLLK = -dLLK
  }else{
    dMLLK = 1e50
  }
  return(dMLLK)
}
