library(GAS)

# Inflation Forecast

data("cpichg")
help(cpichg)

GASSpec   = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = T, scale = T, shape = F))

# Perform H-step ahead forecast with confidence bands

Fit       = UniGASFit(GASSpec,cpichg)
help(UniGASFor)

Forecast  = UniGASFor(Fit, iH = 12)

Forecast

# Perform 1-Step ahead rolling forecast

InsampleData  = cpichg[1:250]
OutSampleData = cpichg[251:276]

Fit       = UniGASFit(GASSpec,InsampleData)

Forecast  = UniGASFor(Fit, Roll = TRUE, vOut = OutSampleData)

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



