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

# Multivariate Example
data("StockIndex")
help(StockIndex)

## Specification mvt
help(MultiGASSpec)
GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity", GASPar = list(location = F, scale = T, correlation = T, shape = F))

help(MultiGASFit)
Fit = MultiGASFit(GASSpec,StockIndex)

Fit

### CPI
data("cpichg")
help(cpichg)

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = T, scale = T, shape = T))

Fit = UniGASFit(GASSpec,cpichg)

Fit

plot(Fit)

