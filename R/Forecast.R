UniGASFor<-function(uGASFit, iH, Roll = F, vOut = NULL,iB = 10000,
                     vBands = c(0.1,0.15,0.85,0.9), bReturnsDraws = F){

  iK          = uGASFit@ModelInfo$iK
  ScalingType = getScalingType(uGASFit)
  lParList    = coef(uGASFit)$lCoef
  Dist        = getDist(uGASFit)
  GASPar      = getGASPar(uGASFit)
  vY          = getObs(uGASFit)

  FilteredParameters = getFilteredParameters(uGASFit)

  if(Roll) {

    iH  = length(vOut)
    iT  = length(vY)
    vYf = c(vY, vOut)
    iT  = iT + iH

    mTheta  = GASFilter_univ(vYf, lParList$vKappa, lParList$mA, lParList$mB, iT, iK, Dist, ScalingType)$mTheta

    PointForecast = t(mTheta[,(iT-iH+1):(iT)])
    cBands        = array(0,dim = c(1,1,1))
    mY            = matrix(0,1,1)
    vU            = EvaluatePit_Univ(t(PointForecast), vOut, Dist, iH)
    vLS           = EvaluateLogScore_Univ(t(PointForecast), vOut, Dist, iH)
  }else{

    vTheta_tp1 = tail(getFilteredParameters(uGASFit),1)
    Forc       = uGASMultiForcast(vTheta_tp1, lParList$vKappa, lParList$mA, lParList$mB, iH, iB, iK, Dist, ScalingType, bReturnsDraws)

    PointForecast = matrix(0,iH,iK)
    cBands        = array(0,dim = c(iH,length(vBands),iK), dimnames = list(1:iH, paste("q.", vBands, sep=""), colnames(vTheta_tp1)))

    for(k in 1:iK){
      PointForecast[,k] = apply(Forc$cTheta[k,,],1,function(x) median(na.omit(x)))
      for(q in vBands){
        cBands[,paste("q.", q, sep=""),k] = apply(Forc$cTheta[k,,],1,function(x,q) quantile(na.omit(x), probs = q), q=q)
      }
    }
    if(bReturnsDraws){mY = Forc$mY}else{mY = matrix(0,1,1)}

    vU = NULL
  }

  colnames(PointForecast) = names(GASPar)
  rownames(PointForecast) = paste("T+",1:iH,sep="")

  mMoments = EvalMoments(t(PointForecast), Dist)

  Out <- new("uGASFor",
             Forecast = list(PointForecast = PointForecast, Moments = mMoments, vU = vU, vLS=vLS),
             Bands    = cBands,
             Draws    = mY,
             Info     = list(iH = iH, Roll = Roll, iB = iB, vBands = vBands,
                              bReturnsDraws = bReturnsDraws, GASPar=GASPar,
                              Dist = Dist, ScalingType = ScalingType, iK = iK),
             Data     = list(vY=vY, FilteredParameters = FilteredParameters, vOut=vOut)
             )

  return(Out)

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
    vU    = EvaluatePit_Univ(t(mForc), tail(vY,ForecastLength), Dist, ForecastLength)
    vLS   = EvaluateLogScore_Univ(t(mForc), tail(vY,ForecastLength), Dist, ForecastLength)
    Moments = EvalMoments(t(mForc), Dist)
  }else{

    for(i in 1:length(FitIndex)){
      if(i != length(FitIndex)){
        lOut[[i]] = vY[(FitIndex[i]+1):(FitIndex[i+1])]
      }else{
        lOut[[i]] = vY[(FitIndex[i]+1):iT]
      }
    }

    lForecasts = lapply(1:length(lOut),function(i, lFits, lOut){

      UniGASFor(uGASFit=lFits[[i]],vOut=lOut[[i]],Roll=T)

    }, lFits = lFits, lOut = lOut)

    mForc = do.call(rbind,lapply(lForecasts, getForecast))
    vU    = do.call(rbind,lapply(lForecasts, pit))
    vLS   = do.call(rbind,lapply(lForecasts, LogScore))
    Moments = do.call(rbind,lapply(lForecasts, getMoments))
  }

  elapsedTime =  Sys.time() - StartTime

  Out <- new("uGASRoll",
             Forecast = list(PointForecast = mForc, vU = vU, vLS = vLS, Moments=Moments),
             Info = list(GASSpec = GASSpec, ForecastLength = ForecastLength,
                         RefitEvery = RefitEvery, RefitWindow = RefitWindow[1],
                         iT = iT, iK = iK, elapsedTime = elapsedTime),
             Data = list(vY = vY))
  return(Out)

}
