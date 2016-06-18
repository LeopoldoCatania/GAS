# Simulate 2x2 Normal random vector with time-varying means, variances and correlation

iT = 1e3
iN = 2
vKappa = c(0.0 , 0.1 , 0.0 , 1.0 , 1.5)

mA = matrix(c(0.0 , 0.0 , 0.0 , 0.0 , 0.0,
              0.0 , 0.0 , 0.0 , 0.0 , 0.0,
              0.0 , 0.0 , 0.02, 0.0 , 0.0,
              0.0 , 0.0 , 0.0 , 0.01, 0.0,
              0.0 , 0.0 , 0.0 , 0.0 , 0.2),5,byrow = T)


mB = diag(5)*0.9
lSim = SimulateGAS_multi(iT, iN, vKappa, mA, mB, Dist="mvnorm", ScalingType = "Identity") # only ScalingType = "Identity" supported right now

mY = lSim$mY

## Specification
GASSpec = MultiGASSpec(Dist = "mvnorm", ScalingType = "Identity", GASPar = list(location = F, scale = T, correlation = T, shape = F))

Fit = MultiGASFit(GASSpec,mY)

