## TestThat comparison
## Obtained with R 3.4.1 & GAS 0.2.3
## 20170721

rm(list = ls())
library("GAS")
require("testthat")
options(digits = 10, max.print = 40, prompt = "R> ", warn = 1)
tol <- 1e-5 # tolerance tests

tmp <- sessionInfo()
Plaform_version <- tmp$platform
GAS_Version     <- tmp$otherPkgs$GAS$Version

TestedPlaform_version = c("x86_64-w64-mingw32/x64 (64-bit)",
                          "x86_64-apple-darwin15.6.0 (64-bit)",
                          "x86_64-pc-linux-gnu (64-bit)")

## x86_64-w64-mingw32/x64 (64-bit) results
exp.mean <- c(0.04421838364, 1.84906157709, 1.05189016745, 7.90182294996)
exp.sd   <- c(0.001423688599, 0.551999497058, 0.002433515713, 0.027090204728)
exp.als  <- c(2.0109318694, 1.0537944859, 0.1546565121, 0.3245448897, 0.5248718852, 0.5289226007)

if (!isTRUE(Plaform_version %in% TestedPlaform_version)) {
  mess = paste0("Current platform never tested\n", "Using x86_64-w64-mingw32/x64 (64-bit) results\n")
  warning(mess)
}

if (Plaform_version == "x86_64-apple-darwin15.6.0 (64-bit)") {
  exp.mean <- c(0.04421838386, 1.84906157303, 1.05189016647, 7.90182317122)
  exp.sd   <- c(0.001423688779, 0.551999494038, 0.002433515636, 0.027090303647)
  exp.als  <- c(2.0109318694, 1.0537944857, 0.1546565121, 0.3245448896, 0.5248718851, 0.5289226006)
}

if (Plaform_version == "x86_64-pc-linux-gnu (64-bit)") {
  exp.mean <- c(0.04421838297, 1.84906156745, 1.05189016422, 7.90182306150)
  exp.sd   <- c(0.001423689286, 0.551999496966, 0.002433518239, 0.027090237859)
  exp.als  <- c(2.0109318691, 1.0537944857, 0.1546565121, 0.3245448895, 0.5248718851, 0.5289226005)
}

## asset WMT
data("dji30ret", package = "GAS")

WMTret <- as.vector(dji30ret[, "WMT", drop = TRUE])

## GAS specification
uGASSpec <- GAS::UniGASSpec(Dist = "sstd", ScalingType = "Identity",
                            GASPar = list(location = FALSE, scale = TRUE, shape = FALSE))

## GAS estimation
cat("running estimation (this takes some time...)\n")
uGASRoll <- GAS::UniGASRoll(data = WMTret, GASSpec = uGASSpec, ForecastLength = 300,
                            RefitEvery = 100, RefitWindow = "moving",
                            cluster = NULL, fn.optimizer = GAS::fn.optim)

## testthat
cat("running tests\n")

# GAS parameters mean
est.mean <- apply(GAS::getForecast(uGASRoll), 2, mean)
testthat::expect_true(max(est.mean - exp.mean) < tol)

# GAS parameters std
est.sd <- apply(GAS::getForecast(uGASRoll), 2, sd)
testthat::expect_true(max(est.sd - exp.sd) < tol)

# average log-score
est.als <- GAS::BacktestDensity(uGASRoll, lower = -100, upper = 100)$average
testthat::expect_true(max(est.als - exp.als) < tol)
