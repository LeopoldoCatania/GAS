getSpec<-function(object) {
  if(is(object,"uGASFit" ) | is(object,"mGASFit" ))   Spec = object@ModelInfo$Spec
  if(is(object,"uGASRoll" ) | is(object,"mGASRoll" )) Spec = object@Info$GASSpec

  return(Spec)
}
getDist<-function(object){
  if(is(object,"uGASSpec") | is(object,"mGASSpec")) Dist = object@Spec$Dist
  if(is(object,"uGASFit" ) | is(object,"mGASFit" )) Dist = object@ModelInfo$Spec@Spec$Dist
  if(is(object,"uGASSim" ) | is(object,"mGASSim" )) Dist = object@ModelInfo$Dist
  if(is(object,"uGASFor") | is(object,"mGASFor"))   Dist = object@Info$Dist
  if(is(object,"uGASRoll") | is(object,"mGASRoll")) Dist = getDist(getSpec(object))

  return(Dist)
}
getScalingType<-function(object){
  if(is(object,"uGASSpec") | is(object,"mGASSpec")) ScalingType = object@Spec$ScalingType
  if(is(object,"uGASFit" ) | is(object,"mGASFit" )) ScalingType = object@ModelInfo$Spec@Spec$ScalingType
  if(is(object,"uGASSim" ) | is(object,"mGASSim" )) ScalingType = object@ModelInfo$ScalingType
  if(is(object,"uGASFor") | is(object,"mGASFor"))   ScalingType = object@Info$ScalingType
  if(is(object,"uGASRoll") | is(object,"mGASRoll")) ScalingType = getScalingType(getSpec(object))

  return(ScalingType)
}
getGASPar<-function(object){
  if(is(object,"uGASSpec") | is(object,"mGASSpec")) GASPar = object@Spec$GASPar
  if(is(object,"uGASFit" ) | is(object,"mGASFit" )) GASPar = object@ModelInfo$Spec@Spec$GASPar
  if(is(object,"uGASFor") | is(object,"mGASFor" ))  GASPar = object@Info$GASPar

  return(GASPar)
}
##
getIC<-function(object){
  return(object@Estimates$IC)
}
