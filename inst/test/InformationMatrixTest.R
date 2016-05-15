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

dMu     = 1
dSigma2 = 4
dNu     = 5

vTheta = c(dMu,dSigma2,dNu)

round(matrix(adaptIntegrate(function(dY,vTheta){
  vScore = Score_univ(dY, vTheta, "std")
  vScore%*%t(vScore) * ddist_univ(dY, vTheta, "std", bLog=FALSE)
},lowerLimit = -50, upperLimit = 50, vTheta = vTheta, fDim = length(vTheta)^2)$integral,length(vTheta)),4)

IM_univ(vTheta, "std")   #check

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



