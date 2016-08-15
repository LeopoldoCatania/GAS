library(GAS)

# Univariate Example

data("sp500ret")
help(sp500ret)

## build the GAS specification
help(UniGASSpec)
GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = F, scale = T, shape = F))

## Estimate the GAS model
help(UniGASFit)
Fit = UniGASFit(GASSpec,sp500ret)

Fit
plot(Fit)


GASSpec = UniGASSpec(Dist = "sstd", ScalingType = "Identity", GASPar = list(location = F, scale = T))

Fit = UniGASFit(GASSpec,vY)


# Multivariate Example
data("StockIndices")
help(StockIndices)

## Specification mvt
help(MultiGASSpec)
GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity", GASPar = list(location    = FALSE, scale = TRUE,
                                                                             correlation = TRUE,  shape = FALSE),
                       ScalarParameters = FALSE)

help(MultiGASFit)
Fit = MultiGASFit(GASSpec,StockIndices)

Fit

## Specification mvt with scalar parameters
data("StockIndices")
help(StockIndices)

GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity", GASPar = list(location = FALSE, scale = TRUE,
                                                                             correlation = TRUE, shape = FALSE),
                       ScalarParameters = TRUE)

Fit = MultiGASFit(GASSpec,StockIndices)

Fit_FULL = MultiGASFit(GASSpec,StockIndices)


### CPI
data("cpichg")
help(cpichg)

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = T, scale = T, shape = T))

Fit     = UniGASFit(GASSpec,cpichg)

Fit

plot(Fit)


###
library(GAS)
data("StockIndices")
help(StockIndices)

vY = StockIndices["DAX",]

GASSpec = UniGASSpec(Dist = "ald", ScalingType = "Identity", GASPar = list(location = F, scale = T, skewness = T))
Fit     = UniGASFit(GASSpec,vY)

Fit
plot(Fit)


## Integer valued model

data("tqdata")
help(tqdata)

GASSpec = UniGASSpec(Dist = "poi", ScalingType = "Identity", GASPar = list(location = T))

vY = abs(tqdata[,7] - tqdata[,8])  # absolute bid-ask spread

Fit     = UniGASFit(GASSpec,vY)

Fit

plot(Fit)

# Realised volatility

#with gamma distribution

data("sp500rv")
help(sp500rv)

GASSpec = UniGASSpec(Dist = "gamma", ScalingType = "Identity", GASPar = list(shape = T, scale = T))

Fit     = UniGASFit(GASSpec,sp500rv[1:3000])

Fit

plot(Fit)

# with exponential distribution
data("sp500rv")
help(sp500rv)

GASSpec = UniGASSpec(Dist = "exp", ScalingType = "Identity", GASPar = list(location = T))

Fit     = UniGASFit(GASSpec,sp500rv[1:3000])

Fit

plot(Fit)


# unemploiment with BETA - GAS

data(usunp)
help(usunp)

GASSpec = UniGASSpec(Dist = "beta", ScalingType = "Identity", GASPar = list(scale = T, shape = T))

Fit     = UniGASFit(GASSpec,unp)

Fit

plot(Fit)


###

data("dji30ret")

GASSpec = UniGASSpec(Dist = "sstd", ScalingType = "Identity", GASPar = list(scale = T, shape = T, skewness = F))

Fit     = UniGASFit(GASSpec,dji30ret[,1])

Fit

plot(Fit)
