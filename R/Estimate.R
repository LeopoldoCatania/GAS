
StaticMLFIT<-function(vY,Dist){

  iT = length(vY)
  iK = NumberParameters(Dist)

  vTheta_tilde  = StaticStarting_Univ(vY,Dist,iK)

  optimiser     = solnp(vTheta_tilde, StaticLLKoptimizer,vY=vY,Dist=Dist, iT=iT,iK = iK, control = list(trace=0))

  vTheta_tilde  = optimiser$pars

  vTheta        = as.numeric(MapParameters(vTheta_tilde, Dist, iK))
  names(vTheta) = getParNames(Dist)

  out=list(vTheta=vTheta,dLLK=-tail(optimiser$values,1),optimiser=optimiser)

  return(out)
}

UniGASFit<-function(GASSpec,vY){

  iT = length(vY)
  iK = GASSpec$iK

  Dist        = GASSpec$Dist
  ScalingType = GASSpec$ScalingType
  GASPar      = GASSpec$GASPar

  # starting par
  lStarting = UnivGAS_Starting(vY,iT,iK,Dist,ScalingType)
  vPw       = lStarting$vPw
  StaticFit = lStarting$StaticFit

  # fixed par
  FixedPar = GetFixedPar_Univ(Dist,GASPar)
  vPw = RemoveFixedPar(vPw, FixedPar)

  #optimise
  optimiser = solnp(vPw, UnivGASOptimiser, vY=vY, Dist=Dist, ScalingType=ScalingType, iT=iT, iK=iK)

  lParList = vPw2lPn_Univ(vPw,iK)
  lParList = AddFixedPar(lParList)

  GASDyn = GASFilter_univ(vY, lParList$vKappa, lParList$mA, lParList$mB, iT, iK, Dist, ScalingType)

  return(list(Spec = GASSpec, GASDyn = GASDyn, Estimates = list(lParList=lParList, optimiser=optimiser, StaticFit=StaticFit )))
}
