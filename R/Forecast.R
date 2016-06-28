UniGASForc<-function(uGASFit, iH, Roll = F, vOut = NULL){

  if(Roll) {
    iH=length(vOut)
  }else{iH=1}

  vY          = getObs(uGASFit)
  iT          = length(vY)
  iK          = uGASFit@ModelInfo$iK
  ScalingType = getScalingType(uGASFit)
  lParList    = coef(uGASFit)$lCoef
  Dist        = getDist(uGASFit)

  if(Roll) {
    vY = c(vY, vOut)
    iT = iT + iH
    mTheta = GASFilter_univ(vY, lParList$vKappa, lParList$mA, lParList$mB, iT, iK, Dist, ScalingType)$mTheta
    mForc  = mTheta[,(iT-iH+1):(iT)]

    return(mForc)
  }

}

UniGASRoll<-function(vY,GASSpec,ForecastLength = 500, Nstart = NULL, RefitEvery = 23, RefitWindow = c("recursive", "moving"),
                     cluster=NULL){

  StartTime = Sys.time()

  iT = length(vY)

  Dist        = getDist(GASSpec)
  ScalingType = getScalingType(GASSpec)
  GASPar      = getGASPar(GASSpec)
  iK          = NumberParameters(Dist)

  if(!is.null(cluster)) clusterEvalQ(cluster,{library(GAS)})
  if(!is.null(Nstart)) {
    iStart=Nstart
  }else{
    iStart=iT-ForecastLength
  }

  FitIndex = seq(iStart,iT,RefitEvery)
  if(tail(FitIndex,1)==iT) FitIndex=FitIndex[-length(FitIndex)]

  lFits=list()
  lForecasts=list()
  lData=list()
  lOut=list()

  if(RefitWindow[1]=="recursive") {
    for(i in 1:length(FitIndex)) {
      lData[[i]]=vY[1:FitIndex[i]]
    }
  }
  if(RefitWindow[1]=="moving") {
    for(i in 1:length(FitIndex)){
      lData[[i]]=vY[(FitIndex[i]-iStart+1):FitIndex[i]]
    }
  }
  #
  #fits
  if(is.null(cluster))  lFits=lapply(lData,UniGASFit, GASSpec=GASSpec)
  if(!is.null(cluster)) lFits=parLapply(cluster,lData,UniGASFit,GASSpec=GASSpec)

  #coef
  lCoef=lapply(lFits, coef)

  if(RefitEvery==1){
    mForc = do.call(rbind, lapply(lFits, function(x) tail(getFilteredParameters(x),1)))
  }else{
    for(i in 1:length(FitIndex)){
      if(i!=length(FitIndex)){
        lOut[[i]]=vY[(FitIndex[i]+1):(FitIndex[i+1])]
      }else{
        lOut[[i]]=vY[(FitIndex[i]+1):iT]
      }
    }
    lForCluster=list()
    for(i in 1:length(lOut)){
      lForCluster[[i]]=list(uGASFit=lFits[[i]],vOut=lOut[[i]])
    }
    if(!is.null(cluster)) {
      lForecasts=parLapply(cluster,lForCluster,function(x){
        forc=UniGASForc(uGASFit=x$uGASFit,vOut=x$vOut,Roll=T)
        return(forc)
      })
    }
    if(is.null(cluster)){
      lForecasts=lapply(lForCluster,function(x){
        forc=UniGASForc(uGASFit=x$uGASFit,vOut=x$vOut,Roll=T)
        return(forc)
      })
    }
    mForc = t(do.call(cbind,lForecasts))
  }

  vU = EvaluatePit_Univ(t(mForc), tail(vY,ForecastLength), Dist, ForecastLength)

  return(list(vU = vU, mForc=mForc, lFits=lFits, lCoef=lCoef))
}
