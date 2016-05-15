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
dPhi = 4
dNu  = 5

vTheta = c(dMu,dPhi,dNu)

grad(function(vTheta,dY){
  ddist_univ(dY, vTheta, "std", bLog=TRUE)
},vTheta, dY = dY)

Score_univ(dY, vTheta, "std")


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





