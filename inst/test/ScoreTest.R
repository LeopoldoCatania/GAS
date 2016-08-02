library(GAS)
library(numDeriv)

### norm score
dY      = 0
dMu     = 1
dSigma2 = 4

vTheta = c(dMu,dSigma2)

grad(function(vTheta,dY){
  ddist_univ(dY, vTheta, "norm", bLog=TRUE)
},vTheta, dY = dY)

Score_univ(dY, vTheta, "norm")

## std score

dY   = 0
dMu  = 1
dPhi2 =5
dNu  = 5

vTheta = c(dMu,dPhi2,dNu)

grad(function(vTheta,dY){
  ddist_univ(dY, vTheta, "std", bLog=TRUE)
},vTheta, dY = dY)

Score_univ(dY, vTheta, "std")

ddist("std",dY,(dPhi2)*(dNu-2)/dNu,shape = dNu)

ddist_univ(dY, vTheta, "std", bLog=FALSE)
dt( (dY-vTheta[1])/sqrt(vTheta[2]),vTheta[3])/sqrt(vTheta[2])

## ast score

dY   = 0
dMu  = 1
dSigma = 4
dNu1  = 5
dNu2  = 6
dAlpha = 0.4

vTheta = c(dMu,dSigma,dAlpha,dNu1,dNu2)

grad(function(vTheta,dY){
  ddist_univ(dY, vTheta, "ast", bLog=TRUE)
},vTheta, dY = dY)

Score_univ(dY, vTheta, "ast")

## ast1 score

dY   = 0
dMu  = 1
dSigma = 4
dNu  = 5
dAlpha = 0.4

vTheta = c(dMu,dSigma,dAlpha,dNu)

grad(function(vTheta,dY){
  ddist_univ(dY, vTheta, "ast1", bLog=TRUE)
},vTheta, dY = dY)

Score_univ(dY, vTheta, "ast1")


# snorm

dY      =  0.9
dMu     = 1
dSigma2 = 4
dDelta  = 1.2

vTheta = c(dMu,dSigma2,dDelta)

grad(function(vTheta,dY){
  ddist_univ(dY, vTheta, "snorm", bLog=TRUE)
},vTheta, dY = dY)

Score_univ(dY, vTheta, "snorm")


## multivariate

## mvt

###### mvt

library(GAS)
library(numDeriv)

iN = 3

vPhi  =  rnorm(iN*(iN-1)/2)
mR = MapR_C(vPhi, iN)

vD = exp(rnorm(iN))/10
mD  = diag(vD)
vMu = rnorm(iN)
dNu = 4

mSigma = mD%*%mR%*%mD

vY = rnorm(iN)

dmvt(vY, vMu, mSigma,dNu,TRUE)

adaptIntegrate(dmvt,lowerLimit=rep(-5,iN),upperLimit=rep(5,iN),vMu=vMu, mSigma=mSigma,dNu=dNu,bLog=FALSE)

vR = build_vR(mR,iN)

grad(function(vR,mD,vY, iN,vMu, dNu){
  mR     =  build_mR(vR, iN)
  mSigma = mD%*%mR%*%mD
  dmvt(vY,vMu,mSigma,dNu,TRUE)
},vR,mD = mD, vY=vY, iN=iN,vMu=vMu,dNu=dNu)

RhoScore_mvt(vR, mD, vY, vMu, dNu, iN)

grad(function(vMu,mR, mD,vY,dNu, iN){
  mSigma = mD%*%mR%*%mD
  dmvt(vY,vMu,mSigma,dNu,TRUE)
},vMu,mR = mR,mD=mD, vY=vY, iN=iN, dNu=dNu)

MuScore_mvt(vMu, mD,mR, vY, dNu,iN)


grad(function(vD,mR,vY, iN,dNu, vMu){
  mD = diag(vD)
  mSigma = mD%*%mR%*%mD
  dmvt(vY,vMu,mSigma,dNu,TRUE)
},vD,mR = mR, vY=vY, iN=iN,vMu=vMu,dNu=dNu)

DScore_mvt(mD, mR, vY, vMu,dNu, iN)

grad(function(dNu,mR,vY, iN, vD, vMu){
  mD = diag(vD)
  mSigma = mD%*%mR%*%mD
  dmvt(vY,vMu,mSigma,dNu,TRUE)
},dNu,vD=vD,mR = mR, vY=vY, iN=iN,vMu=vMu)

NuScore_mvt(mD, mR, vY, vMu,dNu, iN)

jacobian(function(vPhi,iN){
  mR = MapR_C(vPhi, iN)
  build_vR(mR,iN)
},vPhi,iN=iN)

Jacobian_MapR(vPhi, iN)

## all

vTheta = c(vMu, vD, vR, dNu)

Score_multi(vY, vTheta, iN, "mvt")

grad(function(vTheta,vY,iN){
  dmvt_ThetaParam(vY,vTheta,iN,TRUE)
},vTheta,vY=vY, iN=iN)

iK = NumberParameters("mvt",iN)
vTheta_tilde = UnmapParameters_multi(vTheta,"mvt",iN,iK)


grad(function(vTheta_tilde,vY,iN,iK){
  vTheta = MapParameters_multi(vTheta_tilde,"mvt",iN,iK)
  dmvt_ThetaParam(vY,vTheta,iN,TRUE)
},vTheta_tilde,vY=vY, iN=iN,iK=iK)

vS = Score_multi(vY, vTheta, iN, "mvt")

mJ = MapParametersJacobian_multi(vTheta_tilde, "mvt",iN,iK)

t(mJ)%*%vS

