library(GAS)

data("cpichg")

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = T, scale = T, shape = F))

Fit = UniGASFit(GASSpec,cpichg)

help(ConfidenceBands)
Bands = ConfidenceBands(Fit)

Filtration = getFilteredParameters(Fit)


