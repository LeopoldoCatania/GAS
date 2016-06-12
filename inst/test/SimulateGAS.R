library(GAS)
## Simulate GAS process

# norm

iT = 1e4
vKappa = c(0,0.1)
mA = matrix(c(0.1,0.0,
              0.0,0.4),2)
mB = matrix(c(0.9,0.0,
              0.0,0.95),2)

lSim = SimulateGAS_univ(iT, vKappa, mA, mB, Dist="norm", ScalingType = "Identity")
plot.ts(lSim$vY)

plot.ts(t(lSim$mTheta))

GASSpec = UnivGASSpec("norm","Identity", GASPar = list(location = F, scale = T, skewness = F, shape = F, shape2 = F))

vY = lSim$vY

Fit     = UniGASFit(GASSpec,vY)

lSim$dLLK

Fit$Estimates$lParList

# std

iT = 1e3
vKappa = c(0,0,log(15))
mA = matrix(c(0.0 , 0.0 , 0.0 ,
              0.0 , 0.0 , 0.0 ,
              0.0 , 0.0 , 0.0),3,byrow = T)

mB = matrix(c(0.7 , 0.0 , 0.0 ,
              0.0 , 0.98, 0.0 ,
              0.0 , 0.0 , 0.9),3,byrow = T)

GASSpec = UnivGASSpec("std","Identity", GASPar = list(location = F, scale = F, skewness = F, shape = F, shape2 = F))

iB = 100
mPar= matrix(,iB,5)

for(i in 1:iB){
  lSim = SimulateGAS_univ(iT, vKappa, mA, mB, Dist="std", ScalingType = "Identity")
  vY = lSim$vY

  Fit = UniGASFit(GASSpec,vY)

  mPar[i,]=c(Fit$Estimates$lParList$vKappa,Fit$Estimates$lParList$mA[2,2],Fit$Estimates$lParList$mB[2,2])
}

## ast

iT = 1e3
vKappa = c(0.0 , 0.1 , 0.0 , 1.0 , 1.5)
mA = matrix(c(0.001 , 0.0 , 0.0 , 0.0 , 0.0,
              0.0 , 0.05, 0.0 , 0.0 , 0.0,
              0.0 , 0.0 , 0.04, 0.0 , 0.0,
              0.0 , 0.0 , 0.0 , 0.1 , 0.0,
              0.0 , 0.0 , 0.0 , 0.0 , 0.2),5,byrow = T)

mB = matrix(c(0.4 , 0.0 , 0.0 , 0.0 , 0.0,
              0.0 , 0.95, 0.0 , 0.0 , 0.0,
              0.0 , 0.0 , 0.97, 0.0 , 0.0,
              0.0 , 0.0 , 0.0 , 0.98, 0.0,
              0.0 , 0.0 , 0.0 , 0.0 , 0.97),5,byrow = T)

lSim = SimulateGAS_univ(iT, vKappa, mA, mB, Dist="ast", ScalingType = "Identity")
plot.ts(lSim$vY)

plot.ts(t(lSim$mTheta))






