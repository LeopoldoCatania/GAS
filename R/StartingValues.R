StaticStarting_Uni<-function(vY,Dist,iK){

  if(Dist=="std"){
    mu=mean(vY)
    df=8
    phi2=var(vY)*(df-2)/df
    vTheta=c(mu,phi2,df)
  }

  if(Dist=="norm"){
    mu=mean(vY)
    sigma2=var(vY)
    vTheta=c(mu,sigma2)
  }
  if(Dist == "ast" | Dist == "ast1" ){

    empmean=mean(vY)
    empsd=sd(vY)

    st_vY=(vY-empmean)/empsd

    df1=6
    if(Dist=="ast"){df2=8}else{df2=df1}
    alpha=0.5
    mu = mean(vY)
    sigma = sd(vY)

    if(Dist=="ast"){
      vTheta=c(mu,sigma,alpha,df1,df2)
    }else{
      vTheta=c(mu,sigma,alpha,df1)
    }
  }
  if(Dist=="poi"){
    mu=mean(vY)
    vTheta=c(mu)
  }
  if(Dist=="gamma"){
    dMean   = mean(vY)
    dSigma2 = var(vY)

    dBeta  = dMean/dSigma2
    dAlpha = dMean^2 / dSigma2

    vTheta=c(dAlpha, dBeta)
  }

  if(Dist=="exp"){
    vTheta=c(1.0/mean(vY))
  }
  if(Dist=="beta"){
    dMean   = mean(vY)
    dSigma2 = var(vY)

    dAlpha = dMean*(dMean*(1.0 - dMean)/dSigma2 - 1)
    dBeta = (1.0 - dMean)*(dMean*(1.0 - dMean)/dSigma2 - 1)

    vTheta=c(dAlpha, dBeta)
  }
  if(Dist == "ald"){

    dTheta = mean(vY)
    dSigma = sd(vY)
    dKappa = 1.0

    vTheta=c(dTheta, dSigma, dKappa)
  }

  vTheta_tilde = as.numeric(UnmapParameters_univ(vTheta, Dist, iK = iK))
  names(vTheta_tilde) = FullNamesUni(Dist)
  return(vTheta_tilde)
}

StaticStarting_Multi<-function(mY, Dist, iN){

  vEmpRho    = build_vR(cor(t(mY)),iN)
  vEmpPhi    = UnMapR_C(vEmpRho, iN)

  vEmpMu     = apply(mY,1,mean)
  vEmpSigma  = apply(mY,1,sd)

  if(Dist=="mvnorm"){
    vMu_tilde    = vEmpMu
    vSigma_tilde = log(vEmpSigma)
    vPw          = c(vMu_tilde, vSigma_tilde, vEmpPhi)
  }
  if(Dist=="mvt"){

    dNu = try(exp(solnp(log(5),
                        function(x,vY){
                          dNu = exp(x) + LowerNu()
                          -sum(dt(vY,x,log=TRUE))
                          },
                        vY = c(mY), control = list(trace = 0))$pars)
              + LowerNu(),
              silent=F)

    if(is(dNu, "try-error")) dNu = 7.0

    vMu_tilde    = vEmpMu
    vSigma_tilde = log(vEmpSigma*(dNu - 2.0)/dNu)
    vPw          = c(vMu_tilde, vSigma_tilde, vEmpPhi, log(dNu - LowerNu()))
  }

  return(as.numeric(vPw))

}

UniGAS_Starting<-function(vY,iT,iK,Dist,ScalingType, GASPar){
  StaticFit = StaticMLFIT(vY,Dist)
  vKappa = StaticFit$optimiser$pars; names(vKappa) = paste("kappa",1:iK,sep="")

  if(iK>1){
    vA     = starting_vA_Uni(vY, vKappa, mB=diag(rep(9e-1,iK)), dA_foo=1e-15, iT, iK, Dist, ScalingType=ScalingType, GASPar)
    vB     = starting_vB_Uni(vY, vKappa, dB_foo=0.9, mA=diag(vA), iT, iK, Dist,ScalingType=ScalingType, GASPar)
  }else{
    vA     = starting_vA_Uni(vY, vKappa, mB=matrix(9e-1,iK,iK), dA_foo=1e-15, iT, iK, Dist, ScalingType=ScalingType, GASPar)
    vB     = starting_vB_Uni(vY, vKappa, dB_foo=0.9, mA=matrix(vA,iK,iK), iT, iK, Dist,ScalingType=ScalingType, GASPar)
  }
  vA = unmapVec_C(vA, LowerA(), UpperA()); names(vA) = paste("a",1:iK,sep="")
  vB = unmapVec_C(vB, LowerB(), UpperB()); names(vB) = paste("b",1:iK,sep="")

  return(list(vPw = c(vKappa,vA,vB),StaticFit=StaticFit))

}
starting_vA_Uni<-function(vY, vKappa, mB, dA_foo, iT, iK, Dist,ScalingType, GASPar){

  seq_alpha = c(seq(1e-5,8.5,length.out = 30))

  mA = matrix(0,iK,iK)
  diag(mA)[unlist(GASPar)] = dA_foo

  dAlpha_best = dA_foo

  for(i in 1:iK){
    if(GASPar[[i]]){
      for(l in 1:length(seq_alpha)){
        dLLK_foo = try(GASFilter_univ(vY, vKappa, mA, mB, iT, iK, Dist, ScalingType)$dLLK,silent = T)
        if(is.numeric(dLLK_foo) & !is.nan(dLLK_foo)){
          mA[i,i]  = seq_alpha[l]
          dLLK_post = try(GASFilter_univ(vY, vKappa, mA, mB, iT, iK, Dist, ScalingType)$dLLK,silent = T)
          if(is.numeric(dLLK_post) & !is.nan(dLLK_post)){
            if(dLLK_post>dLLK_foo){
              dAlpha_best = seq_alpha[l]
            }
          }
          mA[i,i]  = dAlpha_best
        }
      }
    }
  }

  return(diag(mA))
}
starting_vB_Uni<-function(vY, vKappa, dB_foo, mA, iT, iK, Dist,ScalingType, GASPar){

  seq_beta = c(seq(0.5,0.98,length.out = 30))

  mB = matrix(0,iK,iK)
  diag(mB) = dB_foo

  dB_best = dB_foo

  for(i in 1:iK){
    if(GASPar[[i]]){
      for(l in 1:length(seq_beta)){
        dLLK_foo = try(GASFilter_univ(vY, vKappa, mA, mB, iT, iK, Dist, ScalingType)$dLLK,silent = T)
        if(is.numeric(dLLK_foo) & !is.nan(dLLK_foo)){
          mB[i,i]  = seq_beta[l]
          dLLK_post = try(GASFilter_univ(vY, vKappa, mA, mB, iT, iK, Dist, ScalingType)$dLLK,silent = T)
          if(is.numeric(dLLK_post) & !is.nan(dLLK_post)){
            if(dLLK_post>dLLK_foo){
              dB_best = seq_beta[l]
            }
          }
          mB[i,i]  = dB_best
        }
      }
    }
  }

  return(diag(mB))
}


StartingValues_mvnorm<-function(mY,iT, iN, iK,GASPar, ScalingType){

  vEmpRho    = build_vR(cor(t(mY)),iN)
  vEmpPhi    = UnMapR_C(vEmpRho, iN)

  vEmpMu     = apply(mY,1,mean)
  vEmpSigma  = apply(mY,1,sd)
  #
  # vA     = c(unmapVec_C(c(rep(0.000001,iN),rep(0.001,iN),
  #                       rep(0.0001,iN*(iN-1)/2)),
  #                     LowerA(), UpperA() ) )              ; names(vA) = paste("a.",mvnormParNames(iN),sep="")
  # vB     = c(unmapVec_C(rep(0.90,iK) , LowerB(), UpperB())) ; names(vB) = paste("b.",mvnormParNames(iN),sep="")
  vKappa = c(vEmpMu,log(vEmpSigma),vEmpPhi)               ; names(vKappa) = paste("kappa.",mvnormParNames(iN),sep="")
  #

  vA     = starting_vA_Multi(mY, vKappa, mB=diag(rep(0.9,iK)), dA_foo=1e-2, iT, iK, iN, "mvnorm", ScalingType=ScalingType, GASPar)
  vB     = starting_vB_Multi(mY, vKappa, dB_foo=0.90, mA=diag(vA), iT, iK, iN, "mvnorm", ScalingType=ScalingType, GASPar)

  vA = unmapVec_C(vA, LowerA(), UpperA()); names(vA) = paste("a.",mvnormParNames(iN),sep="")
  vB = unmapVec_C(vB, LowerB(), UpperB()); names(vB) = paste("b.",mvnormParNames(iN),sep="")

  pw = c(vKappa,vA,vB)
  return(pw)
}
StartingValues_mvt<-function(mY,iT, iN, iK,GASPar, ScalingType){

  # iT = ncol(mY)
  # iK = NumberParameters("mvt",iN)
  #
  # vEmpRho    = build_vR(cor(t(mY)),iN)
  # vEmpPhi    = UnMapR_C(vEmpRho, iN)
  #
  # dNu        = 5.0
  # vEmpMu     = apply(mY,1,mean)
  # vEmpSigma  = apply(mY,1,sd)*sqrt((dNu-2.0)/dNu)

  # vKappa = c(vEmpMu,log(vEmpSigma),vEmpPhi, log(dNu - LowerNu()))                         ; names(vKappa) = paste("kappa.",mvtParNames(iN),sep="")
  # vB     = c(unmapVec_C(rep(0.90,iK), LowerB(), UpperB()))                             ; names(vB) = paste("b.",mvtParNames(iN),sep="")
  # vA     = c(unmapVec_C(c(rep(0.000001,iN),rep(0.001,iN),
  #                       rep(0.0001,iN*(iN-1)/2), 0.00001),
  #                     LowerA(), UpperA()) )                                          ; names(vA) = paste("a.",mvtParNames(iN),sep="")

  StaticFit = StaticMLFIT_Multiv(mY, "mvt")

  vKappa = StaticFit$optimiser$pars; names(vKappa) =  paste("kappa.",mvtParNames(iN),sep="")

  vA     = starting_vA_Multi(mY, vKappa, mB=diag(rep(0.9,iK)), dA_foo=1e-2, iT, iK, iN, "mvt", ScalingType=ScalingType, GASPar)
  vB     = starting_vB_Multi(mY, vKappa, dB_foo=0.90, mA=diag(vA), iT, iK, iN, "mvt", ScalingType=ScalingType, GASPar)

  # uGASSpec_foo = UniGASSpec("std", GASPar = list( location = GASPar$location, scale = GASPar$scale))
  # uGASFit = lapply(1:iN, function(i,mY,uGASSpec_foo) coef(UniGASFit(uGASSpec_foo, mY[i,])), mY = mY, uGASSpec_foo=uGASSpec_foo  )
  #
  # vKappa_Mu  = sapply(uGASFit, function(x) x$lCoef$vKappa[1]) ; names(vKappa_Mu)  = paste("kappa.mu",1:iN,sep="")
  # vKappa_Phi = sapply(uGASFit, function(x) x$lCoef$vKappa[2]) ; names(vKappa_Phi) = paste("kappa.phi",1:iN,sep="")
  # dKappa_Nu  = mean(sapply(uGASFit, function(x) x$lCoef$vKappa[3])) ; names(dKappa_Nu)  = "kappa.nu"
  #
  # vA_Mu      = sapply(uGASFit, function(x) x$lCoef$mA[1,1]) ; names(vA_Mu)  = paste("a.mu",1:iN,sep="")
  # vA_Phi     = sapply(uGASFit, function(x) x$lCoef$mA[2,2]) ; names(vA_Phi) = paste("a.phi",1:iN,sep="")
  #
  # vB_Mu      = sapply(uGASFit, function(x) x$lCoef$mB[1,1]) ; names(vB_Mu)  = paste("b.mu",1:iN,sep="")
  # vB_Phi     = sapply(uGASFit, function(x) x$lCoef$mB[2,2]) ; names(vB_Phi) = paste("b.phi",1:iN,sep="")
  #
  # vA     = numeric(iK) ; names(vA)     = paste("a.", mvtParNames(iN),sep="")
  # vB     = numeric(iK) ; names(vB)     = paste("b.", mvtParNames(iN),sep="")
  # vKappa = numeric(iK) ; names(vKappa) = paste("kappa.", mvtParNames(iN),sep="")
  #
  # vA[c(paste("a.mu",1:iN,sep=""), paste("a.phi",1:iN,sep=""))] = c(vA_Mu, vA_Phi)
  # vB[c(paste("b.mu",1:iN,sep=""), paste("b.phi",1:iN,sep=""))] = c(vB_Mu, vB_Phi)
  # vKappa[c(paste("kappa.mu",1:iN,sep=""), paste("kappa.phi",1:iN,sep=""))] = c(vKappa_Mu, vKappa_Phi)
  #
  # vKappa[paste("kappa.",paste(RhoNames(iN)),sep="")] = vEmpPhi
  #
  # ### find vA and vB for correlations
  #
  # vA[paste("a.",paste(RhoNames(iN)),sep="")] = 1e-3
  # vB[paste("b.",paste(RhoNames(iN)),sep="")] = 0.9
  #
  # vA["a.nu"] = 1e-3
  # vB["b.nu"] = 0.9

  # if(GASPar$correlation | GASPar$shape){
  #
  #   GASPar_new = GASPar
  #   GASPar_new$location = FALSE
  #   GASPar_new$shape    = FALSE
  #
  #   vA_other = starting_vA_Multi(mY, vKappa, mB=diag(rep(9e-1,iK)), dA_foo=1e-15, iT, iK, iN, "mvt", ScalingType=ScalingType, GASPar)
  #   vB_other = starting_vB_Multi(mY, vKappa, dB_foo=0.9, mA=diag(vA), iT, iK, iN, "mvt", ScalingType=ScalingType, GASPar)
  #   vA       = unmapVec_C(vA, LowerA(), UpperA()); names(vA) = paste("a.",mvtParNames(iN),sep="")
  #   vB       = unmapVec_C(vB, LowerB(), UpperB()); names(vB) = paste("b.",mvtParNames(iN),sep="")
  #
  #
  # }

  vA = unmapVec_C(vA, LowerA(), UpperA()); names(vA) = paste("a.",mvtParNames(iN),sep="")
  vB = unmapVec_C(vB, LowerB(), UpperB()); names(vB) = paste("b.",mvtParNames(iN),sep="")

  pw = c(vKappa,vA,vB)
  return(pw)
}

MultiGAS_Starting<-function(mY,iT, iN, iK,Dist, GASPar, ScalingType){
  if(Dist=="mvnorm") vPw = StartingValues_mvnorm(mY,iT, iN, iK,  GASPar, ScalingType)
  if(Dist=="mvt")    vPw = StartingValues_mvt(mY,iT, iN, iK,  GASPar, ScalingType)

  return(vPw)
}


starting_vA_Multi<-function(mY, vKappa, mB, dA_foo, iT, iK, iN, Dist,ScalingType, GASPar){

  seq_alpha = c(seq(1e-4,0.5,length.out = 30))

  vBool = FixedDynamicPar_Multi(Dist, iN, GASPar)

  mA = matrix(0,iK,iK)
  # diag(mA)[vBool] = dA_foo

  dAlpha_best = dA_foo

  for(i in 1:iK){
    if(vBool[i]){
      for(l in 1:length(seq_alpha)){
        dLLK_foo =
          # try(
            GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK
            # ,silent = F)

        if(is.numeric(dLLK_foo) & !is.nan(dLLK_foo)){
          mA[i,i]  = seq_alpha[l]
          dLLK_post =
            # try(
              GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK
              # ,silent = F)
          if(is.numeric(dLLK_post) & !is.nan(dLLK_post)){
            if(dLLK_post>dLLK_foo){
              dAlpha_best = seq_alpha[l]
            }
          }
          mA[i,i]  = dAlpha_best
        }
      }
    }
  }

  return(diag(mA))
}

starting_vB_Multi<-function(mY, vKappa, dB_foo, mA, iT, iK, iN, Dist,ScalingType, GASPar){

  seq_beta = c(seq(0.5,0.98,length.out = 30))
  vBool    = FixedDynamicPar_Multi(Dist, iN, GASPar)

  mB = matrix(0,iK,iK)
  diag(mB) = dB_foo

  dB_best = dB_foo

  for(i in 1:iK){
    if(vBool[i]){
      for(l in 1:length(seq_beta)){
        dLLK_foo = try(GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK,silent = T)
        if(is.numeric(dLLK_foo) & !is.nan(dLLK_foo)){
          mB[i,i]  = seq_beta[l]
          dLLK_post = try(GASFilter_multi(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType)$dLLK,silent = T)
          if(is.numeric(dLLK_post) & !is.nan(dLLK_post)){
            if(dLLK_post>dLLK_foo){
              dB_best = seq_beta[l]
            }
          }
          mB[i,i]  = dB_best
        }
      }
    }
  }

  return(diag(mB))
}
