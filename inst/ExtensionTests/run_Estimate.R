library("GAS")

# Univariate Example

data("sp500ret", package = "GAS")
help(sp500ret)

## build the GAS specification
help(UniGASSpec)
GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = FALSE, scale = TRUE, shape = FALSE))

## Estimate the GAS model
help(UniGASFit)
Fit = UniGASFit(GASSpec, sp500ret)

Fit

# Multivariate Example
data("StockIndices", package = "GAS")
help(StockIndices)

## Specification mvt
help(MultiGASSpec)
GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity", GASPar = list(location    = FALSE, scale = TRUE,
                                                                             correlation = TRUE,  shape = FALSE),
                       ScalarParameters = FALSE)

help(MultiGASFit)
Fit_FULL = MultiGASFit(GASSpec, StockIndices)

Fit_FULL

## Specification mvt with scalar parameters

GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity", GASPar = list(location = FALSE, scale = TRUE,
                                                                             correlation = TRUE, shape = FALSE),
                       ScalarParameters = TRUE)

Fit = MultiGASFit(GASSpec, StockIndices)

Fit

### CPI
data("cpichg", package = "GAS")
help(cpichg)

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = TRUE, scale = TRUE, shape = TRUE))

Fit     = UniGASFit(GASSpec, cpichg)

Fit

###
library(GAS)
data("StockIndices", package = "GAS")
help(StockIndices)

vY = StockIndices[, "DAX"]

GASSpec = UniGASSpec(Dist = "ald", ScalingType = "Identity", GASPar = list(location = FALSE, scale = TRUE, skewness = TRUE))
Fit     = UniGASFit(GASSpec,vY)

Fit

## Integer valued model

data("tqdata", package = "GAS")
help(tqdata)

GASSpec = UniGASSpec(Dist = "poi", ScalingType = "Identity", GASPar = list(location = TRUE))

vY = abs(tqdata[, 6] - tqdata[, 5])  # absolute bid-ask spread

Fit     = UniGASFit(GASSpec,vY)

Fit

# Realised volatility

#with gamma distribution

data("sp500rv", package = "GAS")
help(sp500rv)

GASSpec = UniGASSpec(Dist = "gamma", ScalingType = "Identity", GASPar = list(shape = TRUE, scale = TRUE))

Fit = UniGASFit(GASSpec, sp500rv[1:3000])

Fit

# with exponential distribution
data("sp500rv", package = "GAS")
help(sp500rv)

GASSpec = UniGASSpec(Dist = "exp", ScalingType = "Identity", GASPar = list(location = TRUE))

Fit = UniGASFit(GASSpec, sp500rv[1:3000])

Fit

# unemploiment with BETA - GAS

data("usunp", package = "GAS")

GASSpec = UniGASSpec(Dist = "beta", ScalingType = "Identity", GASPar = list(scale = TRUE, shape = TRUE))

Fit = UniGASFit(GASSpec, usunp)

Fit

