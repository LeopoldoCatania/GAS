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

UniGAS_Starting<-function(vY,iT,iK,Dist,ScalingType, GASPar){
  StaticFit = StaticMLFIT(vY,Dist)
  vKappa = StaticFit$optimiser$pars; names(vKappa) = paste("kappa",1:iK,sep="")

  if(iK>1){
    vA     = starting_vA(vY, vKappa, mB=diag(rep(9e-1,iK)), dA_foo=1e-15, iT, iK, Dist, ScalingType=ScalingType, GASPar)
    vB     = starting_vB(vY, vKappa, dB_foo=0.9, mA=diag(vA), iT, iK, Dist,ScalingType=ScalingType)
  }else{
    vA     = starting_vA(vY, vKappa, mB=matrix(9e-1,iK,iK), dA_foo=1e-15, iT, iK, Dist, ScalingType=ScalingType, GASPar)
    vB     = starting_vB(vY, vKappa, dB_foo=0.9, mA=matrix(vA,iK,iK), iT, iK, Dist,ScalingType=ScalingType)
  }
  vA = unmapVec_C(vA, LowerA(), UpperA()); names(vA) = paste("a",1:iK,sep="")
  vB = unmapVec_C(vB, LowerB(), UpperB()); names(vB) = paste("b",1:iK,sep="")

  return(list(vPw = c(vKappa,vA,vB),StaticFit=StaticFit))

}
starting_vA<-function(vY, vKappa, mB, dA_foo, iT, iK, Dist,ScalingType, GASPar){

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
starting_vB<-function(vY, vKappa, dB_foo, mA, iT, iK, Dist,ScalingType){

  seq_beta = c(seq(0.5,0.98,length.out = 30))

  mB = matrix(0,iK,iK)
  diag(mB) = dB_foo

  dB_best = dB_foo

  for(i in 1:iK){
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

  return(diag(mB))
}


StartingValues_mvnorm<-function(mY,iN){

  iT = ncol(mY)
  iK = NumberParameters("mvnorm",iN)

  vEmpRho    = build_vR(cor(t(mY)),iN)
  vEmpPhi    = UnMapR_C(vEmpRho, iN)

  vEmpMu     = apply(mY,1,mean)
  vEmpSigma  = apply(mY,1,sd)

  vA     = c(unmapVec_C(c(rep(0.000001,iN),rep(0.001,iN),
                        rep(0.0001,iN*(iN-1)/2)),
                      LowerA(), UpperA() ) )              ; names(vA) = paste("a.",mvnormParNames(iN),sep="")
  vB     = c(unmapVec_C(rep(0.90,iK) , LowerB(), UpperB())) ; names(vB) = paste("b.",mvnormParNames(iN),sep="")
  vKappa = c(vEmpMu,log(vEmpSigma),vEmpPhi)               ; names(vKappa) = paste("kappa.",mvnormParNames(iN),sep="")

  pw = c(vKappa,vA,vB)
  return(pw)
}
StartingValues_mvt<-function(mY,iN){

  iT = ncol(mY)
  iK = NumberParameters("mvt",iN)

  vEmpRho    = build_vR(cor(t(mY)),iN)
  vEmpPhi    = UnMapR_C(vEmpRho, iN)

  dNu        = 5.0
  vEmpMu     = apply(mY,1,mean)
  vEmpSigma  = apply(mY,1,sd)*sqrt((dNu-2.0)/dNu)

  vKappa = c(vEmpMu,log(vEmpSigma),vEmpPhi, log(dNu - LowerNu()))                         ; names(vKappa) = paste("kappa.",mvtParNames(iN),sep="")
  vB     = c(unmapVec_C(rep(0.90,iK), LowerB(), UpperB()))                             ; names(vB) = paste("b.",mvtParNames(iN),sep="")
  vA     = c(unmapVec_C(c(rep(0.000001,iN),rep(0.001,iN),
                        rep(0.0001,iN*(iN-1)/2), 0.00001),
                      LowerA(), UpperA()) )                                          ; names(vA) = paste("a.",mvtParNames(iN),sep="")

  pw = c(vKappa,vA,vB)
  return(pw)
}

MultiGAS_Starting<-function(mY,iN,Dist){
  if(Dist=="mvnorm") vPw = StartingValues_mvnorm(mY,iN)
  if(Dist=="mvt")    vPw = StartingValues_mvt(mY,iN)

  return(vPw)
}

