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
