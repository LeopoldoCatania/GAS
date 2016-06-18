StaticStarting_Univ<-function(vY,Dist,iK){

  if(Dist=="std"){
    mu=mean(vY)
    df=8
    phi=sd(vY)*sqrt((df-2)/df)
    vTheta=c(mu,phi,df)
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

    K1=Kast(df1)
    K2=Kast(df2)

    alpha_star=alpha*K1/(alpha*K1+(1-alpha)*K2)
    B=alpha*K1/alpha_star

    EY=4*B*(-alpha_star^2*(df1/(df1-1) + (1-alpha_star)^2*df2/(df2-1)))
    VARY=4*(alpha*alpha_star^2*df1/(df1-2) + (1-alpha)*(1-alpha_star)^2*df2/(df2-2)) - 16*B^2*(-alpha_star^2*df1/(df1-1) +
                                                                                                 (1-alpha_star)^2*df2/(df2-1))^2

    sigma=empsd/sqrt(VARY)
    mu=empmean-sigma*EY

    if(Dist=="ast"){
      vTheta=c(mu,sigma,alpha,df1,df2)
    }else{
      vTheta=c(mu,sigma,alpha,df1)
    }
  }
  vTheta_tilde = as.numeric(UnmapParameters(vTheta, Dist, iK = iK))
  names(vTheta_tilde) = getParNames(Dist)
  return(vTheta_tilde)
}

UnivGAS_Starting<-function(vY,iT,iK,Dist,ScalingType){
  StaticFit = StaticMLFIT(vY,Dist)
  vKappa = StaticFit$optimiser$pars; names(vKappa) = paste("kappa",1:iK,sep="")

  vA     = starting_vA(vY, vKappa, mB=diag(rep(9e-1,iK)), dA_foo=1e-15, iT, iK, Dist, ScalingType=ScalingType)
  vB     = starting_vB(vY, vKappa, dB_foo=0.9, mA=diag(vA), iT, iK, Dist,ScalingType=ScalingType)

  vA = unmapVec_C(vA, LowerA(), UpperA()); names(vA) = paste("a",1:iK,sep="")
  vB = unmapVec_C(vB, LowerB(), UpperB()); names(vB) = paste("b",1:iK,sep="")

  return(list(vPw = c(vKappa,vA,vB),StaticFit=StaticFit))

}
starting_vA<-function(vY, vKappa, mB, dA_foo, iT, iK, Dist,ScalingType){

  seq_alpha = c(seq(1e-15,4.5,length.out = 30))

  mA = matrix(0,iK,iK)
  diag(mA) = dA_foo

  dAlpha_best = dA_foo

  for(i in 1:iK){
    for(l in 1:length(seq_alpha)){
      dLLK_foo = try(GASFilter_univ(vY, vKappa, mA, mB, iT, iK, Dist, ScalingType)$dLLK,silent = T)
      if(is.numeric(dLLK_foo)){
        mA[i,i]  = seq_alpha[l]
        dLLK_post = try(GASFilter_univ(vY, vKappa, mA, mB, iT, iK, Dist, ScalingType)$dLLK,silent = T)
        if(is.numeric(dLLK_post)){
          if(dLLK_post>dLLK_foo){
            dAlpha_best = seq_alpha[l]
          }
        }
        mA[i,i]  = dAlpha_best
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
      if(is.numeric(dLLK_foo)){
        mB[i,i]  = seq_beta[l]
        dLLK_post = try(GASFilter_univ(vY, vKappa, mA, mB, iT, iK, Dist, ScalingType)$dLLK,silent = T)
        if(is.numeric(dLLK_post)){
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

  vKappa = c(vEmpMu,log(vEmpSigma),vEmpPhi, log(dNu - 2.01))                         ; names(vKappa) = paste("kappa.",mvtParNames(iN),sep="")
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

