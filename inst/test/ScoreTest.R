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

