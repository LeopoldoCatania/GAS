library(GAS)

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

plot(Forecast)

# Perform 1-Step ahead rolling forecast

InsampleData  = cpichg[1:250]
OutSampleData = cpichg[251:276]

Fit       = UniGASFit(GASSpec,InsampleData)

Forecast  = UniGASFor(Fit, Roll = TRUE, out = OutSampleData)

Forecast

plot(Forecast)

pit(Forecast)
getForecast(Forecast) # see ?UniGASFor and ?uGASFor
LogScore(Forecast)


# Perform 1-step ahead rolling forecast with refit
library(parallel)

cluster = makeCluster(7)


help(UniGASRoll)
Roll = UniGASRoll(cpichg,GASSpec,ForecastLength = 50, RefitEvery = 2, RefitWindow = c("moving"), cluster = cluster)

Roll


########################
#    MULTIVARIATE      #
########################

# data("StockIndices")
data("dji30ret", package = "GAS")

# mY = StockIndices[, 1:2]
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

# InSampleData  = mY[1:1000, ]
# OutSampleData = mY[1001:2404, ]

InSampleData  = mY[1:2521, ]
OutSampleData = mY[2522:5521, ]

# estimation
Fit = MultiGASFit(GASSpec, InSampleData)

Forecast  = MultiGASFor(Fit, Roll = TRUE, out = OutSampleData)

Forecast

plot(Forecast)


DCCRoll = dccroll(spec = mDCCSpec, data = mY, forecast.length = 3000,
                  refit.every = 2998)

aCov = rcov(DCCRoll)
mMu = fitted(DCCRoll)

vShapes = c(sapply(DCCRoll@mforecast, function(x) rep(x@model$pars["mshape", "Level"], ceiling(iH/length(DCCRoll@mforecast)))))

vLLK = numeric(iH)

for (i in 1:iH) {

  dNu = vShapes[i]

  vLLK[i] = dmvt(t(as.numeric(dji30ret_Out[i, vAsset_Multi])),
                     t(as.numeric(mMu[i, ])), aCov[, , i] * (dNu - 2)/dNu, dNu, bLog = TRUE)
}

###



# Perform 1-step ahead rolling forecast with refit
library(parallel)

cluster = makeCluster(7)

Roll    = MultiGASRoll(mY,GASSpec,ForecastLength = 500, RefitEvery = 500, RefitWindow = c("moving"), cluster = cluster)

Roll

## Realised Volatility

data("sp500rv")
help(sp500rv)

GASSpec = UniGASSpec(Dist = "gamma", ScalingType = "Identity", GASPar = list(shape = T, scale = T))

Fit     = UniGASFit(GASSpec,sp500rv[1:3000])

# Forecast  = UniGASFor(Fit, iH = 500, iB = 10000)

Forecast  = UniGASFor(Fit, Roll = TRUE, vOut = sp500rv[3001:4310])

Forecast_Gamma = getMoments(Forecast)[,1]

plot(Forecast)

