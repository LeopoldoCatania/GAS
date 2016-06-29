
StaticMLFIT<-function(vY,Dist){

  iT = length(vY)
  iK = NumberParameters(Dist)

  vTheta_tilde  = StaticStarting_Uni(vY,Dist,iK)

  optimiser     = solnp(vTheta_tilde, StaticLLKoptimizer,vY=vY,Dist=Dist, iT=iT,iK = iK, control = list(trace=0))

  vTheta_tilde  = optimiser$pars

  vTheta        = as.numeric(MapParameters_univ(vTheta_tilde, Dist, iK))
  names(vTheta) = FullNamesUni(Dist)

  out=list(vTheta=vTheta,dLLK=-tail(optimiser$values,1),optimiser=optimiser)

  return(out)
}

UniGASFit<-function(GASSpec,vY){

  Start = Sys.time()

  iT = length(vY)

  Dist        = getDist(GASSpec)
  ScalingType = getScalingType(GASSpec)
  GASPar      = getGASPar(GASSpec)
  iK          = NumberParameters(Dist)

  # starting par
  lStarting = UniGAS_Starting(vY,iT,iK,Dist,ScalingType)
  vPw       = lStarting$vPw
  StaticFit = lStarting$StaticFit

  # fixed par
  FixedPar = GetFixedPar_Uni(Dist,GASPar)
  vPw = RemoveFixedPar(vPw, FixedPar)

  #optimise
  optimiser = solnp(vPw, UniGASOptimiser, vY=vY, Dist=Dist, ScalingType=ScalingType, iT=iT, iK=iK)

  vPw = optimiser$pars

  lParList = vPw2lPn_Uni(vPw,iK)
  lParList = AddFixedPar(lParList)

  Inference = InferenceFun_Uni(optimiser$hessian,vPw, iK)

  GASDyn = GASFilter_univ(vY, lParList$vKappa, lParList$mA, lParList$mB, iT, iK, Dist, ScalingType)

  IC = ICfun(-tail(optimiser$values,1),length(optimiser$pars),iT)

  vU = EvaluatePit_Univ(GASDyn$mTheta, vY, Dist, iT)

  mMoments = EvalMoments(GASDyn$mTheta,Dist)

  elapsedTime =  Sys.time() - Start

  Out <- new("uGASFit", ModelInfo = list(Spec = GASSpec, iT = iT, iK = iK, elapsedTime = elapsedTime),
             GASDyn = GASDyn,
             Estimates = list(lParList=lParList, optimiser=optimiser,
                              StaticFit=StaticFit,
                              Inference = Inference,IC=IC,vU=vU,
                              Moments   = mMoments),
             Data = list(vY = vY))

  return(Out)
}

# mY is NxT
MultiGASFit<-function(GASSpec,mY){
  Start = Sys.time()
  # getInfo
  Dist        = getDist(GASSpec)
  ScalingType = getScalingType(GASSpec)
  GASPar      = getGASPar(GASSpec)

  #dimension par
  iT = ncol(mY)
  iN = nrow(mY)
  iK = NumberParameters(Dist,iN)

  if(is.null(rownames(mY))) rownames(mY) = paste("Series",1:iN)

  # starting par
  vPw = MultiGAS_Starting(mY,iN,Dist)

  # fixed par
  FixedPar = GetFixedPar_Multi(Dist,GASPar,iN)
  vPw      = RemoveFixedPar(vPw, FixedPar)

  #optimise
  optimiser = solnp(vPw, MultiGASOptimiser, mY=mY, Dist=Dist, ScalingType=ScalingType, iT=iT, iN=iN, iK=iK)

  vPw = optimiser$pars

  lParList = vPw2lPn_Multi(vPw,Dist,iK,iN)
  lParList = AddFixedPar(lParList)

  Inference = InferenceFun_Multi(optimiser$hessian, Dist, vPw, iK, iN)

  GASDyn = GASFilter_multi(mY, lParList$vKappa, lParList$mA, lParList$mB, iT, iN, iK, Dist, ScalingType)

  IC = ICfun(-tail(optimiser$values,1),length(optimiser$pars),iT)

  elapsedTime =  Sys.time() - Start

  Out <- new("mGASFit", ModelInfo = list(Spec = GASSpec, iT = iT, iN = iN, iK = iK, elapsedTime = elapsedTime),
             GASDyn = GASDyn,
             Estimates = list(lParList=lParList, optimiser=optimiser,
                              Inference = Inference,IC = IC ),
             Data = list(mY = mY))

}
