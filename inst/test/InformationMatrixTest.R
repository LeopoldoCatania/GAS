library(GAS)
library(cubature)

#norm
dMu     = 1
dSigma2 = 4

vTheta = c(dMu,dSigma2)

matrix(adaptIntegrate(function(dY,vTheta){
  vScore = Score_univ(dY, vTheta, "norm")
  vScore%*%t(vScore) * ddist_univ(dY, vTheta, "norm", bLog=FALSE)
},lowerLimit = -10, upperLimit = 10, vTheta = vTheta, fDim = length(vTheta)^2)$integral,length(vTheta))

IM_univ(vTheta, "norm")

#std

dMu     = 0
dPhi2   = 1
dNu     = 20

vTheta = c(dMu,dPhi,dNu)

round(matrix(adaptIntegrate(function(dY,vTheta){
  vScore = Score_univ(dY, vTheta, "std")
  vScore%*%t(vScore) * ddist_univ(dY, vTheta, "std", bLog=FALSE)
},lowerLimit = -5, upperLimit = 5, vTheta = vTheta, fDim = length(vTheta)^2)$integral,length(vTheta)),7)

IM_univ(vTheta, "std")


## ast

dMu  = 1
dSigma = 4
dNu1  = 5
dNu2  = 6
dAlpha = 0.4

vTheta = c(dMu,dSigma,dAlpha,dNu1,dNu2)

round(matrix(adaptIntegrate(function(dY,vTheta){
  vScore = Score_univ(dY, vTheta, "ast")
  vScore%*%t(vScore) * ddist_univ(dY, vTheta, "ast", bLog=FALSE)
},lowerLimit = -50, upperLimit = 50, vTheta = vTheta, fDim = length(vTheta)^2)$integral,length(vTheta)),4)

IM_univ(vTheta, "ast")

## ast1

dMu  = 1
dSigma = 4
dNu1  = 5
dAlpha = 0.4

vTheta = c(dMu,dSigma,dAlpha,dNu1)

round(matrix(adaptIntegrate(function(dY,vTheta){
  vScore = Score_univ(dY, vTheta, "ast1")
  vScore%*%t(vScore) * ddist_univ(dY, vTheta, "ast1", bLog=FALSE)
},lowerLimit = -50, upperLimit = 50, vTheta = vTheta, fDim = length(vTheta)^2)$integral,length(vTheta)),4)

IM_univ(vTheta, "ast1")

## snorm

dMu  = 1
dSigma2 = 4
dDelta  = 1.1

vTheta = c(dMu,dSigma2,dDelta)

round(matrix(adaptIntegrate(function(dY,vTheta){
  vScore = Score_univ(dY, vTheta, "snorm")
  vScore%*%t(vScore) * ddist_univ(dY, vTheta, "snorm", bLog=FALSE)
},lowerLimit = -15, upperLimit = 15, vTheta = vTheta, fDim = length(vTheta)^2)$integral,length(vTheta)),4)

IM_univ(vTheta, "snorm")



