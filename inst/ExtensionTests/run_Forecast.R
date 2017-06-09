library("GAS")

########################
#      UNIVARIATE      #
########################

# Inflation Forecast

data("cpichg")
help(cpichg)

GASSpec   = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = T, scale = T, shape = F))

# Perform H-step ahead forecast with confidence bands

Fit       = UniGASFit(GASSpec,cpichg)
help(UniGASFor)

Forecast  = UniGASFor(Fit, H = 100)

Forecast

# plot(Forecast)

# Perform 1-Step ahead rolling forecast

InsampleData  = cpichg[1:250]
OutSampleData = cpichg[251:276]

Fit       = UniGASFit(GASSpec, InsampleData)

Forecast  = UniGASFor(Fit, Roll = TRUE, out = OutSampleData)

Forecast

# plot(Forecast)

pit(Forecast)
getForecast(Forecast) # see ?UniGASFor and ?uGASFor
LogScore(Forecast)

# Perform 1-step ahead rolling forecast with refit
library(parallel)

cluster = makeCluster(7)

help(UniGASRoll)
Roll = UniGASRoll(as.numeric(cpichg), GASSpec, ForecastLength = 50, RefitEvery = 2, RefitWindow = c("moving"), cluster = cluster)

Roll

########################
#    MULTIVARIATE      #
########################

# data("StockIndices")
data("dji30ret", package = "GAS")

mY = dji30ret[, 1:2]

## Specification mvt
GASSpec = MultiGASSpec(Dist = "mvt", ScalingType = "Identity",
                       GASPar = list(location = FALSE, scale = TRUE, correlation = TRUE, shape = FALSE))

# Perform H-step ahead forecast with confidence bands

# estimation
Fit = MultiGASFit(GASSpec, mY)

#forecast

Forecast  = MultiGASFor(Fit, H = 50)

Forecast

plot(Forecast)

# Perform 1-Step ahead rolling forecast

InSampleData  = mY[1:2521, ]
OutSampleData = mY[2522:5521, ]

# estimation
Fit = MultiGASFit(GASSpec, InSampleData)

Forecast  = MultiGASFor(Fit, Roll = TRUE, out = OutSampleData)

Forecast

plot(Forecast)

# Perform 1-step ahead rolling forecast with refit
library(parallel)

cluster = makeCluster(5)

Roll    = MultiGASRoll(mY, GASSpec, ForecastLength = 500, RefitEvery = 100, RefitWindow = c("moving"), cluster = cluster)

Roll

