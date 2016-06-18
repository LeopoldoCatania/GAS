## download/install/load the GAS package
library(devtools)

install_github("LeopoldoCatania/GAS")

library(GAS)

############################
#       UNIVARIATE        #
############################

# Currently available distributions are:
# Gaussian   : Dist = "norm"
# Student--t : Dist = "std"
# Skew Student--t : Dist = "ast1"
# Skew Student--t with two shape parameters : Dist = "ast"
#
# The skew distribution is that of Zhu and Galbraith (2010) JoE
#
# The Order of the distributions parameters is : location, scale, skewness, shape, shape2
#
# ScalingType represents the scaling mechanism for the conditiona score
# possible choices are :
# ScalingType = "Identity" , no scaling
# ScalingType = "Inv"  , invese of the information matrix
# ScalingType = "InvSqrt" , chol of the inverse of the information matrix ( not so stable)
#
#
# Simulate from a Gaussian random variable with time-varying mean and variance

set.seed(7895)

iT = 1e4

vKappa = c(0,0.1) # unconditional reparametrised value of the parameters

mA = matrix(c(0.1,0.0,
              0.0,0.4),2,byrow = T) # matrix that premultiply the score

mB = matrix(c(0.9,0.0,
              0.0,0.95),2,byrow = T) # matrix with autoregressive parameters

Sim = uGASSim(iT, vKappa, mA, mB, Dist = "norm", ScalingType = "Identity")

## summary of the simulation
Sim

## graphical representation
plot(Sim)

## get the simulated data
vY = getObs(Sim)

## build the GAS specification
GASSpec = UniGASSpec(Dist = "norm", ScalingType = "Identity", GASPar = list(location = T, scale = T))

## Estimate the GAS model
Fit = UniGASFit(GASSpec,vY)

## summary of the estimation
Fit

## graphical representation
plot(Fit)

# Simulate from a Student--t random variable with time-varying mean, scale and shape parameters

iT = 1e4
vKappa = c(0,0.1,log(7))

mA = matrix(c(0.1,0.0,0.0,
              0.0,0.4,0.0,
              0.0,0.0,0.01),3,byrow = T)

mB = matrix(c(0.9,0.0,0.0,
              0.0,0.95,0.0,
              0.0,0.0,0.97),3,byrow = T)

Sim = uGASSim(iT, vKappa, mA, mB, Dist = "std", ScalingType = "Identity")

## summary of the simulation
Sim

## graphical representation
plot(Sim)

## get the simulated data
vY = getObs(Sim)

## build the GAS specification
GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = T, scale = T, shape =T))

## Estimate the GAS model
Fit = UniGASFit(GASSpec,vY)

## summary of the estimation
Fit

## graphical representation
plot(Fit)

############################
#       MULTIVARIATE       #
############################
# Currently available distributions are:
# Multivariate Gaussian   : Dist = "mvnorm"
# Multivariate Student--t : Dist = "mvt"
#
# The Order of the distributions parameters is : locations, scales, correlations, shape
#
# only ScalingType = "Identity" supported for multiv dist right now

# Simulate 2x2 Normal random vector with time-varying means, variances and correlation

iT = 1e3
iN = 2
vKappa = c(0.0 , 0.1 , 0.0 , 1.0 , 1.5)

mA = matrix(c(0.0 , 0.0 , 0.0 , 0.0 , 0.0,
              0.0 , 0.0 , 0.0 , 0.0 , 0.0,
              0.0 , 0.0 , 0.02, 0.0 , 0.0,
              0.0 , 0.0 , 0.0 , 0.01, 0.0,
              0.0 , 0.0 , 0.0 , 0.0 , 0.02),5,byrow = T)


mB  = diag(5)*0.9

Sim = mGASSim(iT, iN,vKappa, mA, mB, Dist = "mvnorm", ScalingType = "Identity")

mY = getObs(Sim)
plot(Sim)

## Specification
GASSpec = MultiGASSpec(Dist = "mvnorm", ScalingType = "Identity", GASPar = list(location = F, scale = T, correlation = T, shape = F))

Fit = MultiGASFit(GASSpec,mY)

Fit

plot(Fit)
