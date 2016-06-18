#' specugas
#'
#' Univariate Generalised Autoregressive Score models specifications
#'
#' @param  Dist The conditional distribution assumed for the data
#' @return An SpecUGAS object
#' @export
UnivGASSpec<-function(Dist = "norm", ScalingType = "Identity", GASPar = list(location = F, scale = T, skewness = F, shape = F, shape2 = F)){

  iK = NumberParameters(Dist)

  return (list(Dist = Dist, ScalingType=ScalingType, iK = iK, GASPar=GASPar))
}
#' specmgas
#'
#' Multivariate Generalised Autoregressive Score models specifications
#'
#' @param  Dist The conditional distribution assumed for the data. Dist = "mvnorm" and Dist = "mvt" are currently supported
#' @return An SpecMGAS object
#' @export
MultiGASSpec<-function(Dist = "norm", ScalingType = "Identity", GASPar = list(location = F, scale = T, correlation = F, shape = F)){

  # iK = NumberParameters(Dist)

  return (list(Dist = Dist, ScalingType=ScalingType,  GASPar=GASPar))
}


getDist<-function(GASSpec){
  GASSpec$Dist
}
getScalingType<-function(GASSpec){
  GASSpec$ScalingType
}
getGASPar<-function(GASSpec){
  GASSpec$GASPar
}
