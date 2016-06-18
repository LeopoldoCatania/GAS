library(GAS)
## Simulate GAS process

################################
#         UNIVARIATE           #
################################

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

# std

iT = 1e3
vKappa = c(0,0,log(15))
mA = matrix(c(0.001 , 0.0 , 0.0 ,
              0.0 , 0.01 , 0.0 ,
              0.0 , 0.0 , 0.004),3,byrow = T)

mB = matrix(c(0.7 , 0.0 , 0.0 ,
              0.0 , 0.98, 0.0 ,
              0.0 , 0.0 , 0.9),3,byrow = T)

lSim = SimulateGAS_univ(iT, vKappa, mA, mB, Dist="std", ScalingType = "Identity")

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



################################
#        MULTIVARIATE          #
################################

# Simulate 2x2 Normal random vector with time-varying means, variances and correlation

iT = 1e3
iN = 2
vKappa = c(0.0 , 0.1 , 0.0 , 1.0 , 1.5)

mA = matrix(c(0.001 , 0.0 , 0.0 , 0.0 , 0.0,
              0.0 , 0.005, 0.0 , 0.0 , 0.0,
              0.0 , 0.0 , 0.002, 0.0 , 0.0,
              0.0 , 0.0 , 0.0 , 0.001 , 0.0,
              0.0 , 0.0 , 0.0 , 0.0 , 0.2),5,byrow = T)


mB = diag(5)*0.9
lSim = SimulateGAS_multi(iT, iN, vKappa, mA, mB, Dist="mvnorm", ScalingType = "Identity") # only ScalingType = "Identity" supported right now

plot.ts(t(lSim$mY))
plot.ts(t(lSim$mTheta))

# Simulate 2x2 Multi Student random vector with time-varying means, scales, correlation and df

iT = 1e3
iN = 2
vKappa = c(0.0 , 0.1 , 0.0 , 1.0 , 1.5, log(7-2))

mA = matrix(c(0.001 , 0.0 , 0.0 , 0.0 , 0.0, 0.0,
              0.0 , 0.005, 0.0 , 0.0 , 0.0, 0.0,
              0.0 , 0.0 , 0.002, 0.0 , 0.0, 0.0,
              0.0 , 0.0 , 0.0 , 0.001 , 0.0, 0.0,
              0.0 , 0.0 , 0.0 , 0.0 , 0.02,  0.0,
              0.0 , 0.0 , 0.0 , 0.0 , 0.0,  0.1),6,byrow = T)


mB = diag(6)*0.9
lSim = SimulateGAS_multi(iT, iN, vKappa, mA, mB, Dist="mvt", ScalingType = "Identity") # only ScalingType = "Identity" supported right now

plot.ts(t(lSim$mY))
plot.ts(t(lSim$mTheta))
