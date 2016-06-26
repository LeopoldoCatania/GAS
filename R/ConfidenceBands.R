ConfidenceBands<-function(object, iB = 999, probs = c(0.01,0.05,0.1,0.9,0.95,0.99)){

  vPw  = object@Estimates$optimiser$pars
  iT   = object@ModelInfo$iT
  iK   = object@ModelInfo$iK
  Data = getObs(object)
  Dist = getDist(object)
  DType = DistType(Dist)
  ScalingType = getScalingType(object)

  vCov = ginv(object@Estimates$optimiser$hessian)/iT

  mPw = rmvnorm_mat(iB, vPw,vCov)

  mTheta = getFilteredParameters(object)

  cTheta      = array(,dim = c(dim(mTheta),iB+1),dimnames = list(1:(iT+1), colnames(mTheta),1:(iB+1)))
  cTheta[,,1] = mTheta

  for(b in 2:(iB+1)){
    vPw_foo  = mPw[b-1,] ; names(vPw_foo) = names(vPw)

    lParList = vPw2lPn_Uni(vPw_foo,iK)
    lParList = AddFixedPar(lParList)

    if(DType=="univariate") cTheta[,,b] = t(GASFilter_univ(Data, lParList$vKappa, lParList$mA, lParList$mB, iT, iK, Dist, ScalingType)$mTheta)
    if(DType=="multivariate") cTheta[,,b] = t(GASFilter_multi(Data, lParList$vKappa, lParList$mA, lParList$mB, iT, nrow(Data),iK, Dist, ScalingType)$mTheta)

  }

  cQuantile = array(,dim = c(dim(mTheta),length(probs)),dimnames = list(1:(iT+1), colnames(mTheta),paste("q",probs,sep = "")))

  for(k in 1:iK){
    cQuantile[,k,] = t(apply(cTheta[,k,],1,quantile,probs = probs))
  }
  return(cQuantile)
}




