getSpec<-function(object) object@ModelInfo$Spec
getDist<-function(object){
  if(is(object,"uGASSpec") | is(object,"mGASSpec")) Dist = object@Spec$Dist
  if(is(object,"uGASFit" ) | is(object,"mGASFit" )) Dist = object@ModelInfo$Spec@Spec$Dist
  if(is(object,"uGASSim" ) | is(object,"mGASSim" )) Dist = object@ModelInfo$Dist
  return(Dist)
}
getScalingType<-function(object){
  if(is(object,"uGASSpec") | is(object,"mGASSpec")) ScalingType = object@Spec$ScalingType
  if(is(object,"uGASFit" ) | is(object,"mGASFit" )) ScalingType = object@ModelInfo$Spec@Spec$ScalingType
  if(is(object,"uGASSim" ) | is(object,"mGASSim" )) ScalingType = object@ModelInfo$ScalingType
  return(ScalingType)
}
getGASPar<-function(object){
  if(is(object,"uGASSpec") | is(object,"mGASSpec")) GASPar = object@Spec$GASPar
  if(is(object,"uGASFit" ) | is(object,"mGASFit" )) GASPar = object@ModelInfo$Spec@Spec$GASPar
  return(GASPar)
}
getFilteredParameters<-function(object){
  if(is(object,"uGASFit" ) | is(object,"mGASFit" )) mTheta = object@GASDyn$mTheta
  if(is(object,"uGASSim" ) | is(object,"mGASSim" )) mTheta = object@GASDyn$mTheta

  mTheta = t(mTheta)
  parNames = getParNames(object)
  colnames(mTheta) = parNames

  return(mTheta)
}
getObs<-function(object){
  if(is(object,"uGASFit" )) Data = object@Data$vY
  if(is(object,"mGASFit" )) Data = object@Data$mY

  if(is(object,"uGASSim" )) Data = object@Data$vY
  if(is(object,"mGASSim" )) Data = object@Data$mY

  return(Data)
}
getIC<-function(object){
  return(object@Estimates$IC)
}
