library(GAS)
## Simulate GAS process

################################
#         UNIVARIATE           #
################################

# norm

iT = 1e4

help(UnmapParameters_univ)
vKappa = UnmapParameters_univ(c(0,0.1),"norm",2)

mA = matrix(c(0.1,0.0,
              0.0,0.4),2)
mB = matrix(c(0.9,0.0,
              0.0,0.95),2)

help(UniGASSim)
Sim = UniGASSim(iT, vKappa, mA, mB, Dist="norm", ScalingType = "Identity")

plot(Sim)

# std

iT = 1e3
vKappa = UnmapParameters_univ(c(0,0.1,8),"std",3)
mA = matrix(c(0.001 , 0.0 , 0.0 ,
              0.0 , 0.01 , 0.0 ,
              0.0 , 0.0 , 0.04),3,byrow = T)

mB = matrix(c(0.7 , 0.0 , 0.0 ,
              0.0 , 0.98, 0.0 ,
              0.0 , 0.0 , 0.97),3,byrow = T)

Sim = UniGASSim(iT, vKappa, mA, mB, Dist="std", ScalingType = "Identity")

plot(Sim)

## ast

iT = 1e3
vKappa = UnmapParameters_univ(c(0.0 , 0.1 , 0.65 , 5 , 9),"ast",5)
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

Sim = UniGASSim(iT, vKappa, mA, mB, Dist="ast", ScalingType = "Identity")
plot(Sim)

# poi

iT = 1e4

vKappa = log(5)

mA = matrix(c(0.05),1)
mB = matrix(c(0.94),1)

Sim = UniGASSim(iT, vKappa, mA, mB, Dist="poi", ScalingType = "Identity")

plot(Sim)

# beta

iT = 1e4

vKappa = c(log(0.1) , log(0.4))

mA = matrix(c(0.1,0.0,
              0.0,0.4),2)*0.01
mB = matrix(c(0.98,0.0,
              0.0,0.99),2)

Sim = UniGASSim(iT, vKappa, mA, mB, Dist="beta", ScalingType = "Identity")

plot(Sim)


# Bernoulli

iT = 1e4

vKappa = UnmapParameters_univ(0.01,"ber",1)

mA = matrix(c(0.1),1)
mB = matrix(c(0.94),1)

Sim = UniGASSim(iT, vKappa, mA, mB, Dist="ber", ScalingType = "Identity")

plot(Sim)


################################
#        MULTIVARIATE          #
################################

# Simulate from a GAS process with Multivariate Student-t conditional distribution, time-varying locations, scales, correlations and fixed shape parameter.
library(GAS)

set.seed(786)

iT     = 1000 # number of observations to simulate
iN     = 3     # trivariate series
Dist   = "mvt" # conditional Multivariate Studen-t distribution

# build unconditional vector of reparametrised parameters

vMu  = c(0.1,0.2,0.3)   # vector of location parameters (this is not transformed)
vPhi = c(1.0, 1.2, 0.3) # vector of scale parameters for the firs, second and third variables.

vRho = c(0.1,0.2,0.3)   # This represents vec(R), where R is the correlation matrix.
# Note that is up to the user to ensure that vec(R) implies a proper correlation matrix

vTheta = c(vMu, vPhi, vRho, 7) # vector of parameters such that the degrees of freedom are 7.

vKappa = UnmapParameters_multi(vTheta, Dist, iN, iK = length(vTheta))

mA     = matrix(0,length(vKappa),length(vKappa))

diag(mA) = c( 0,0,0 , 0.05,0.01,0.09  , 0.01,0.04,0.07,  0) # update scales and correlations, do not update locations and shape parameters

mB     = matrix(0,length(vKappa),length(vKappa))

diag(mB) = c( 0,0,0 , 0.7,0.7,0.5  , 0.94,0.97,0.92,  0) # update scales and correlations, do not update locations and shape parameters

help(MultiGASSim)
Sim = MultiGASSim(iT,iN, vKappa, mA, mB, Dist, ScalingType = "Identity")

Sim

plot(Sim)

