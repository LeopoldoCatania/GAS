
UniGASSpec<-function(Dist = "norm", ScalingType = "Identity", GASPar = list(location = F, scale = T, skewness = F, shape = F, shape2 = F)){

  iK = NumberParameters(Dist)

  if(is.null(GASPar$location)) GASPar$location = FALSE
  if(is.null(GASPar$scale))    GASPar$scale    = FALSE
  if(is.null(GASPar$skewness)) GASPar$skewness = FALSE
  if(is.null(GASPar$shape))    GASPar$shape    = FALSE
  if(is.null(GASPar$shape2))   GASPar$shape2   = FALSE

  DistPar = DistParameters(Dist)

  GASPar = GASPar[DistPar]

  Spec <- new( "uGASSpec" ,Spec = list(Dist = Dist, ScalingType=ScalingType, iK = iK, GASPar=GASPar))
  return(Spec)
}

MultiGASSpec<-function(Dist = "mvnorm", ScalingType = "Identity", GASPar = list(location = F, scale = T, correlation = F, shape = F)){

  if(is.null(GASPar$location))    GASPar$location    = FALSE
  if(is.null(GASPar$scale))       GASPar$scale       = FALSE
  if(is.null(GASPar$correlation)) GASPar$correlation = FALSE
  if(is.null(GASPar$shape))       GASPar$shape       = FALSE

  Spec <- new( "mGASSpec" ,Spec = list(Dist = Dist, ScalingType=ScalingType, GASPar=GASPar))
  return(Spec)
}


