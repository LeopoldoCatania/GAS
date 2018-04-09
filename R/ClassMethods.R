setClass("uGASFit", representation(ModelInfo = "list", GASDyn = "list", Estimates = "list", Testing = "list",
                                   Data = "list"))
setClass("mGASFit", representation(ModelInfo = "list", GASDyn = "list", Estimates = "list", Data = "list"))
setClass("uGASSim", representation(ModelInfo = "list", GASDyn = "list", Data = "list"))
setClass("mGASSim", representation(ModelInfo = "list", GASDyn = "list", Data = "list"))
setClass("uGASSpec", representation(Spec = "list"))
setClass("mGASSpec", representation(Spec = "list"))
setClass("uGASFor", representation(Forecast = "list", Bands = "array", Draws = "matrix", Info = "list",
                                   Data = "list"))
setClass("mGASFor", representation(Forecast = "list", Bands = "array", Draws = "array", Info = "list",
                                   Data = "list"))
setClass("uGASRoll", representation(Forecast = "list", Info = "list", Testing = "list", Data = "list"))
setClass("mGASRoll", representation(Forecast = "list", Info = "list", Data = "list"))

setMethod("show", "uGASSpec", function(object) {

  Dist = getDist(object)
  ScalingType = getScalingType(object)
  GASPar = unlist(getGASPar(object))
  GASPar = names(GASPar[GASPar])

  cat("\n-------------------------------------------------------")
  cat("\n-            Univariate GAS Specification             -")
  cat("\n-------------------------------------------------------")
  cat(paste("\nConditional distribution"))
  DistInfo(Dist, FULL = FALSE)
  cat(paste("\nGAS specification"))
  cat("\n-------------------------------------------------------")
  cat(paste("\nScore scaling type: ", ScalingType))
  cat(paste("\nTime varying parameters: ", paste(GASPar, collapse = ", ")))
  #
  cat("\n-------------------------------------------------------")
})

setMethod("show", "mGASSpec", function(object) {

  Dist = getDist(object)

  ScalingType = getScalingType(object)
  GASPar = unlist(getGASPar(object))
  GASPar = names(GASPar[GASPar])
  ScalarParameters = object@Spec$ScalarParameters

  cat("\n-------------------------------------------------------")
  cat("\n-           Multivariate GAS Specification            -")
  cat("\n-------------------------------------------------------")
  cat(paste("\nConditional distribution"))
  DistInfo(Dist, FULL = FALSE)
  cat(paste("\nGAS specification"))
  cat("\n-------------------------------------------------------")
  cat(paste("\nScore scaling type: ", ScalingType))
  cat(paste("\nTime varying parameters: ", paste(GASPar, collapse = ", ")))
  cat(paste("\nScalar Parameters: ", paste(TypeOfParameters(ScalarParameters))))

  #
  cat("\n-------------------------------------------------------")
})

setMethod("show", "uGASFit", function(object) {

  Spec = getSpec(object)
  Dist = getDist(object)
  iT = object@ModelInfo$iT
  iK = NumberParameters(Dist)
  IC = getIC(object)

  ParNames = FullNamesUni(Dist)

  ScalingType = getScalingType(Spec)
  GASPar = unlist(getGASPar(Spec))
  GASPar = names(GASPar[GASPar])

  Inference = object@Estimates$Inference

  vKappa = object@Estimates$lParList$vKappa
  mB = object@Estimates$lParList$mB
  vTheta_Tilde_Unc = solve(diag(iK) - mB) %*% vKappa

  vTheta_Unc = c(MapParameters_univ(vTheta_Tilde_Unc, Dist, iK))
  names(vTheta_Unc) = ParNames

  elapsedTime = object@ModelInfo$elapsedTime

  cat(paste("\n------------------------------------------"))
  cat(paste("\n-          Univariate GAS Fit            -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification:\t")
  cat(paste("\nT = ", iT))
  cat(paste("\nConditional distribution: ", Dist))
  cat(paste("\nScore scaling type: ", ScalingType))
  cat(paste("\nTime varying parameters: ", paste(GASPar, collapse = ", ")))
  #
  cat(paste("\n------------------------------------------"))
  cat(paste("\nEstimates:\n"))
  print(Inference)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nUnconditional Parameters:\n"))
  print(vTheta_Unc)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nInformation Criteria:\n"))
  print(IC)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nConvergence:\t"))
  cat(convergence(object))
  cat(paste("\n------------------------------------------"))
  cat(paste("\n\nElapsed time:", round(as.double(elapsedTime, units = "mins"), 2L), "mins"))
})

setMethod("summary", "uGASFit", function(object) {

  Spec = getSpec(object)
  Dist = getDist(object)
  iT = object@ModelInfo$iT
  iK = NumberParameters(Dist)
  IC = getIC(object)

  ParNames = FullNamesUni(Dist)

  ScalingType = getScalingType(Spec)
  GASPar = unlist(getGASPar(Spec))
  GASPar = names(GASPar[GASPar])

  Inference = object@Estimates$Inference

  vKappa = object@Estimates$lParList$vKappa
  mB = object@Estimates$lParList$mB
  vTheta_Tilde_Unc = solve(diag(iK) - mB) %*% vKappa

  vTheta_Unc = c(MapParameters_univ(vTheta_Tilde_Unc, Dist, iK))
  names(vTheta_Unc) = ParNames

  elapsedTime = object@ModelInfo$elapsedTime

  vRes = residuals(object, standardize = TRUE)
  JB  = JarqueBera(vRes)

  mBoxRes = LjungBox(vRes)
  mBoxRes2 = LjungBox(vRes^2)

  mTest = rbind(JB[c("Statistic", "p-Value")],
                mBoxRes,
                mBoxRes2)

  rownames(mTest) = c(
    " Jarque-Bera Test   R    Chi^2 ",
    " Ljung-Box Test     R    Q(10) ",
    " Ljung-Box Test     R    Q(15) ",
    " Ljung-Box Test     R    Q(20) ",
    " Ljung-Box Test     R^2  Q(10) ",
    " Ljung-Box Test     R^2  Q(15) ",
    " Ljung-Box Test     R^2  Q(20) ")

  cat(paste("\n------------------------------------------"))
  cat(paste("\n-          Univariate GAS Fit            -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification:\t")
  cat(paste("\nT = ", iT))
  cat(paste("\nConditional distribution: ", Dist))
  cat(paste("\nScore scaling type: ", ScalingType))
  cat(paste("\nTime varying parameters: ", paste(GASPar, collapse = ", ")))
  #
  cat(paste("\n------------------------------------------"))
  cat(paste("\nEstimates:\n"))
  print(Inference)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nUnconditional Parameters:\n"))
  print(vTheta_Unc)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nInformation Criteria:\n"))
  print(IC)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nAnalysis of Residuals:\n"))
  print(mTest)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nConvergence:\t"))
  cat(convergence(object))
  cat(paste("\n------------------------------------------"))
  cat(paste("\n\nElapsed time:", round(as.double(elapsedTime, units = "mins"), 2L), "mins"))
})

setMethod("show", "mGASFit", function(object) {

  Spec = getSpec(object)
  iT = object@ModelInfo$iT
  iN = object@ModelInfo$iN
  iK = object@ModelInfo$iK
  IC = getIC(object)

  Dist = getDist(Spec)
  ScalingType = getScalingType(Spec)
  GASPar = unlist(getGASPar(Spec))
  GASPar = names(GASPar[GASPar])

  ParNames = FullNamesMulti(iN, Dist)

  vKappa = object@Estimates$lParList$vKappa
  mB = object@Estimates$lParList$mB
  vTheta_Tilde_Unc = solve(diag(iK) - mB) %*% vKappa

  vTheta_Unc = c(MapParameters_multi(vTheta_Tilde_Unc, Dist, iN, iK))
  names(vTheta_Unc) = ParNames
  Inference = object@Estimates$Inference
  elapsedTime = object@ModelInfo$elapsedTime

  cat(paste("\n------------------------------------------"))
  cat(paste("\n-        Multivariate GAS Fit            -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification:\t")
  cat(paste("\nT = ", iT))
  cat(paste("\nConditional distribution: ", Dist))
  cat(paste("\nScore scaling type: ", ScalingType))
  cat(paste("\nTime varying parameters: ", paste(GASPar, collapse = ", ")))
  #
  cat(paste("\n------------------------------------------"))
  cat(paste("\nEstimates:\n"))
  print(Inference)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nUnconditional Parameters:\n"))
  print(vTheta_Unc)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nInformation Criteria:\n"))
  print(IC)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nConvergence:\t"))
  cat(convergence(object))
  cat(paste("\n------------------------------------------"))
  cat(paste("\n\nElapsed time:", round(as.double(elapsedTime, units = "mins"), 2L), "mins"))
})

setMethod("summary", "mGASFit", function(object) {

  Spec = getSpec(object)
  iT = object@ModelInfo$iT
  iN = object@ModelInfo$iN
  iK = object@ModelInfo$iK
  IC = getIC(object)

  Dist = getDist(Spec)
  ScalingType = getScalingType(Spec)
  GASPar = unlist(getGASPar(Spec))
  GASPar = names(GASPar[GASPar])

  ParNames = FullNamesMulti(iN, Dist)

  vKappa = object@Estimates$lParList$vKappa
  mB = object@Estimates$lParList$mB
  vTheta_Tilde_Unc = solve(diag(iK) - mB) %*% vKappa

  vTheta_Unc = c(MapParameters_multi(vTheta_Tilde_Unc, Dist, iN, iK))
  names(vTheta_Unc) = ParNames
  Inference = object@Estimates$Inference
  elapsedTime = object@ModelInfo$elapsedTime

  cat(paste("\n------------------------------------------"))
  cat(paste("\n-        Multivariate GAS Fit            -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification:\t")
  cat(paste("\nT = ", iT))
  cat(paste("\nConditional distribution: ", Dist))
  cat(paste("\nScore scaling type: ", ScalingType))
  cat(paste("\nTime varying parameters: ", paste(GASPar, collapse = ", ")))
  #
  cat(paste("\n------------------------------------------"))
  cat(paste("\nEstimates:\n"))
  print(Inference)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nUnconditional Parameters:\n"))
  print(vTheta_Unc)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nInformation Criteria:\n"))
  print(IC)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nConvergence:\t"))
  cat(convergence(object))
  cat(paste("\n------------------------------------------"))
  cat(paste("\n\nElapsed time:", round(as.double(elapsedTime, units = "mins"), 2L), "mins"))
})

setMethod("show", "uGASSim", function(object) {

  iT = object@ModelInfo$iT

  Dist = getDist(object)
  ScalingType = getScalingType(object)
  ParNames = FullNamesUni(Dist)
  iK = NumberParameters(Dist)

  mA = object@ModelInfo$mA
  mB = object@ModelInfo$mB
  vKappa = object@ModelInfo$vKappa
  names(vKappa) = ParNames

  vTheta_Tilde_Unc = solve(diag(iK) - mB) %*% vKappa

  vTheta_Unc = c(MapParameters_univ(vTheta_Tilde_Unc, Dist, iK))
  names(vTheta_Unc) = ParNames

  cat(paste("\n------------------------------------------"))
  cat(paste("\n-          Univariate GAS Sim            -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification:\t")
  cat(paste("\nT = ", iT))
  cat(paste("\nConditional distribution: ", Dist))
  cat(paste("\nScore scaling type: ", ScalingType))
  #
  cat(paste("\n------------------------------------------"))
  cat(paste("\nParameters:\n"))
  cat("vKappa:\n")
  print(vKappa)
  cat("mA:\n")
  print(mA)
  cat("mB:\n")
  print(mB)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nUnconditional Parameters:\n"))
  print(vTheta_Unc)
})

setMethod("show", "mGASSim", function(object) {

  iT = object@ModelInfo$iT
  iN = object@ModelInfo$iN

  Dist = getDist(object)
  ScalingType = getScalingType(object)
  ParNames = FullNamesMulti(iN, Dist)
  iK = NumberParameters(Dist, iN)

  mA = object@ModelInfo$mA
  mB = object@ModelInfo$mB
  vKappa = object@ModelInfo$vKappa
  names(vKappa) = ParNames

  vTheta_Tilde_Unc = solve(diag(iK) - mB) %*% vKappa

  vTheta_Unc = c(MapParameters_multi(vTheta_Tilde_Unc, Dist, iN, iK))
  names(vTheta_Unc) = ParNames

  cat(paste("\n------------------------------------------"))
  cat(paste("\n-          Univariate GAS Sim            -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification:\t")
  cat(paste("\nT = ", iT))
  cat(paste("\nN = ", iN))
  cat(paste("\nConditional distribution: ", Dist))
  cat(paste("\nScore scaling type: ", ScalingType))
  #
  cat(paste("\n------------------------------------------"))
  cat(paste("\nParameters:\n"))
  cat("vKappa:\n")
  print(vKappa)
  cat("mA:\n")
  print(mA)
  cat("mB:\n")
  print(mB)
  cat(paste("\n------------------------------------------"))
  cat(paste("\nUnconditional Parameters:\n"))
  print(vTheta_Unc)
})

setMethod("show", "uGASFor", function(object) {

  iH = object@Info$iH

  Dist = getDist(object)
  ScalingType = getScalingType(object)
  Roll = object@Info$Roll
  PointForecast = getForecast(object)

  if (Roll) {
    PointForecast = cbind(PointForecast, realized = object@Data$vOut)
  }

  cat(paste("\n------------------------------------------"))
  cat(paste("\n-        Univariate GAS Forecast         -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification")
  cat(paste("\nConditional distribution: ", Dist))
  cat(paste("\nScore scaling type: ", ScalingType))
  cat(paste("\nHorizon: ", iH))
  cat(paste("\nRolling forecast: ", Roll))
  #
  cat(paste("\n------------------------------------------"))
  cat(paste("\nParameters forecast:\n"))
  if (nrow(PointForecast) > 10L ) {
    print(head(PointForecast, 5L))
    cat(paste("\n....................\n"))
    print(tail(PointForecast, 5L))
  } else {
    print(PointForecast)
  }
})

setMethod("show", "mGASFor", function(object) {

  iH = object@Info$iH
  iN = object@Info$iN

  Dist = getDist(object)
  ScalingType = getScalingType(object)
  Roll = object@Info$Roll
  PointForecast = getForecast(object)

  if (Roll) {
    PointForecast = cbind(PointForecast, realized = object@Data$vOut)
  }

  cat(paste("\n------------------------------------------"))
  cat(paste("\n-      Multivariate GAS Forecast         -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification")
  cat(paste("\nConditional distribution: ", Dist))
  cat(paste("\nScore scaling type: ", ScalingType))
  cat(paste("\nHorizon: ", iH))
  cat(paste("\nNumber of series: ", iN))
  cat(paste("\nRolling forecast: ", Roll))
  #
  cat(paste("\n------------------------------------------"))
  cat(paste("\nParameters forecast:\n"))
  if (nrow(PointForecast) > 10L) {
    print(head(PointForecast, 5L))
    cat(paste("\n....................\n"))
    print(tail(PointForecast, 5L))
  } else {
    print(PointForecast)
  }
})

setMethod("show", "uGASRoll", function(object) {

  Dist = getDist(object)
  ScalingType = getScalingType(object)
  PointForecast = getForecast(object)
  elapsedTime = object@Info$elapsedTime

  cat(paste("\n------------------------------------------"))
  cat(paste("\n-    Univariate GAS Rolling Forecast     -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification")
  cat(paste("\nConditional distribution: ", Dist))
  cat(paste("\nScore scaling type: ", ScalingType))
  #
  cat(paste("\n------------------------------------------"))
  cat(paste("\nParameters forecast:\n"))
  if (nrow(PointForecast) > 10L) {
    print(head(PointForecast, 5L))
    cat(paste("\n....................\n"))
    print(tail(PointForecast, 5L))
  } else {
    print(PointForecast)
  }
  cat(paste("\n------------------------------------------"))
  cat(paste("\n\nElapsed time:", round(as.double(elapsedTime, units = "mins"), 2L), "mins"))
})

setMethod("show", "mGASRoll", function(object) {

  Dist = getDist(object)
  ScalingType = getScalingType(object)
  PointForecast = getForecast(object)
  elapsedTime = object@Info$elapsedTime

  cat(paste("\n------------------------------------------"))
  cat(paste("\n-   Multivariate GAS Rolling Forecast    -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification")
  cat(paste("\nConditional distribution: ", Dist))
  cat(paste("\nScore scaling type: ", ScalingType))
  #
  cat(paste("\n------------------------------------------"))
  cat(paste("\nParameters forecast:\n"))
  if (nrow(PointForecast) > 10L) {
    print(head(PointForecast, 5L))
    cat(paste("\n....................\n"))
    print(tail(PointForecast, 5L))
  } else {
    print(PointForecast)
  }
  cat(paste("\n------------------------------------------"))
  cat(paste("\n\nElapsed time:", round(as.double(elapsedTime, units = "mins"), 2L), "mins"))
})

setMethod("plot", signature(x = "uGASFit", y = "missing"), function(x, which = NULL) {
  iK = x@ModelInfo$iK
  iT = x@ModelInfo$iT
  vY = x@Data$vY

  FilteredParameters = getFilteredParameters(x)[1:iT, , drop = FALSE]
  Moments = getMoments(x)[1:iT, , drop = FALSE]
  vU = pit(x)

  if (is(vY, "xts")) {
    vDates = as.Date(index(vY))
  } else {
    vDates = 1:length(vY)
  }
  PlotType = 1L
  while (PlotType > 0L) {
    if (is.null(which)) {
      vMenu = PlotMenu(x)

      cat(paste("Print 1-",length(vMenu)," or 0 to exit", sep = ""))
      PlotType = menu(vMenu)
    } else {
      PlotType = which
    }

    if (PlotType == 1L) {
      PlotMultipleSeries(FilteredParameters, iK, iT, vDates)
    }

    if (PlotType == 2L) {
      PlotMultipleSeries(Moments, 4L, iT, vDates)
    }
    if (PlotType == 3L) {
      PlotPit(vU, x@Testing$PitTest$Hist)
    }
    if (PlotType == 4L) {
      PlotSingleSeries(vY, iT, vDates)
    }
    if (PlotType == 5L) {
      mRealVsFiltered = cbind(Moments[, 1L], vY)
      PlotForecastVsRealized_Univ(mRealVsFiltered, vDates, x)
    }
    if (!is.null(which)) {
      PlotType = 0L
    }
  }
})

setMethod("plot", signature(x = "mGASFit", y = "missing"), function(x, which = NULL) {
  iK = x@ModelInfo$iK
  iN = x@ModelInfo$iN
  iT = x@ModelInfo$iT
  mY = t(x@Data$mY)

  vSeriesName = colnames(mY)
  if (is.null(vSeriesName)) {
    vSeriesName = paste("series", 1:iN, sep = "")
  }

  if (is(mY, "xts")) {
    vDates = as.Date(index(mY))
  } else {
    vDates = rownames(mY)
    vDates = try(as.Date(vDates), silent = TRUE)
    if (is.null(vDates) || is(vDates, "try-error")) {
      vDates = 1:nrow(mY)
    }
  }
  PlotType = 1L
  while (PlotType > 0L) {
    if (is.null(which)) {
      cat(paste("Print 1-3 or 0 to exit"))
      PlotType = menu(PlotMenu(x))
    } else {
      PlotType = which
    }
    if (PlotType == 1L) {
      series2plot = getFilteredParameters(x)[1:iT, , drop = FALSE]
    }
    if (PlotType == 2L) {
      lMoments = getMoments(x)
      mMean = lMoments[["mean"]][1:iT, , drop = FALSE]
      dimnames(mMean) = list(vDates, paste(vSeriesName, "mean", sep = "."))

      aCov  = lMoments[["cov"]][,, 1:iT, drop = FALSE]
      dimnames(aCov) = list(vSeriesName, vSeriesName, vDates)
      mCov = Array2Matrix(aCov, type = 2L)

    }
    if (PlotType == 3L) {
      series2plot = mY
    }

    if (PlotType == 1L) {
      PlotMultipleSeries(series2plot, iK, iT, vDates)
    }
    if (PlotType == 2L) {

      PlotMultipleSeries(mMean, iN, iT, vDates)
      foo = readline("Print enter to plot covariances\n:")
      PlotMultipleSeries(mCov, ncol(mCov), iT, vDates)

    }
    if (PlotType == 3L) {
      PlotMultipleSeries(series2plot, iN, iT, vDates)
    }
    if (!is.null(which)) {
      PlotType = 0L
    }
  }
})

setMethod("plot", signature(x = "uGASSim", y = "missing"), function(x, which = NULL) {
  iK = x@ModelInfo$iK
  iT = x@ModelInfo$iT
  vY = x@Data$vY

  vDates = 1:iT

  PlotType = 1L
  while (PlotType > 0L) {
    if (is.null(which)) {
      cat(paste("Print 1-3 or 0 to exit"))
      PlotType = menu(PlotMenu(x))
    } else {
      PlotType = which
    }
    if (PlotType == 1L) {
      series2plot = getFilteredParameters(x)[1:iT, , drop = FALSE]
    }
    if (PlotType == 2L) {
      series2plot = getMoments(x)[1:iT, , drop = FALSE]
    }
    if (PlotType == 3L) {
      series2plot = vY
    }

    if (PlotType == 1L) {
      PlotMultipleSeries(series2plot, iK, iT, vDates)
    }
    if (PlotType == 2L) {
      PlotMultipleSeries(series2plot, iK, iT, vDates)
    }
    if (PlotType == 3L) {
      PlotSingleSeries(series2plot, iT, vDates)
    }
    if (!is.null(which)) {
      PlotType = 0L
    }
  }
})

setMethod("plot", signature(x = "mGASSim", y = "missing"), function(x, which = NULL) {
  iK = x@ModelInfo$iK
  iT = x@ModelInfo$iT
  iN = x@ModelInfo$iN
  mY = t(x@Data$mY)

  vDates = 1:iT
  vSeriesName = colnames(mY)

  PlotType = 1L
  while (PlotType > 0L) {
    if (is.null(which)) {
      cat(paste("Print 1-3 or 0 to exit"))
      PlotType = menu(PlotMenu(x))
    } else {
      PlotType = which
    }
    if (PlotType == 1L) {
      series2plot = getFilteredParameters(x)[1:iT, , drop = FALSE]
    }
    if (PlotType == 2L) {
      lMoments = getMoments(x)
      mMean = lMoments[["mean"]][1:iT, , drop = FALSE]
      dimnames(mMean) = list(vDates, paste(vSeriesName, "mean", sep = "."))

      aCov  = lMoments[["cov"]][,, 1:iT, drop = FALSE]
      dimnames(aCov) = list(vSeriesName, vSeriesName, vDates)
      mCov = Array2Matrix(aCov, type = 2L)

    }
    if (PlotType == 3L) {
      series2plot = mY
    }

    if (PlotType == 1L) {
      PlotMultipleSeries(series2plot, iK, iT, vDates)
    }
    if (PlotType == 2L) {

      PlotMultipleSeries(mMean, iN, iT, vDates)
      foo = readline("Print enter to plot covariances\n:")
      PlotMultipleSeries(mCov, ncol(mCov), iT, vDates)

    }
    if (PlotType == 3L) {
      PlotMultipleSeries(series2plot, iN, iT, vDates)
    }
    if (!is.null(which)) {
      PlotType = 0L
    }
  }
})

setMethod("plot", signature(x = "uGASFor", y = "missing"), function(x, which = NULL) {
  iK = x@Info$iK
  vY = x@Data$vY
  iH = x@Info$iH
  iT = length(vY)

  Roll = x@Info$Roll
  vOut = x@Data$vOut

  FilteredParameters = x@Data$FilteredParameters
  FilteredParameters = FilteredParameters[-nrow(FilteredParameters), ]  #remove one step ahead forecast
  ParametersForecast = getForecast(x)
  cBands = x@Bands

  Dist = getDist(x)

  vLS = LogScore(x)

  if (is(vY, "xts")) {
    vDates_is = as.Date(index(vY))
    if (Roll) {
      vDates_os = as.Date(index(vOut))
    } else {
      DiffTime = vDates_is[2L] - vDates_is[1L]
      vDates_os = seq(tail(vDates_is, 1L) + DiffTime, tail(vDates_is, 1L) + DiffTime * iH, by = DiffTime)
    }
    ParametersForecast = xts(ParametersForecast, vDates_os)
    FilteredParameters = xts(FilteredParameters, vDates_is)
  } else {
    vDates_is = 1:iT
    vDates_os = (iT + 1L):(iT + iH)
  }

  PlotType = 1L
  while (PlotType > 0L) {
    if (!Roll) {
      if (is.null(which)) {
        cat(paste("Print 1-6 or 0 to exit"))
        PlotType = menu(PlotMenu(x))
      } else {
        PlotType = which
      }
      if (PlotType == 1L) {
        PlotMultipleSeries(ParametersForecast, iK, iH, vDates_os)
      }
      if (PlotType == 2L) {
        PlotMultipleSeries_Bands(ParametersForecast, iK, iH, vDates_os, cBands)
      }
      if (PlotType == 3L) {
        PlotMultipleSeries_wis(FilteredParameters, ParametersForecast, iK, iH, vDates_os, vDates_is)
      }
      if (PlotType == 4L) {
        PlotMultipleSeries_Bands_wis(FilteredParameters, ParametersForecast, iK, iH, vDates_os,
                                     vDates_is, cBands)
      }
      if (PlotType == 5L) {
        Moments_is = getMoments(x)
        PlotMultipleSeries(Moments_is, 4L, iH, vDates_os)
      }
      if (PlotType == 6L) {
        Moments_os = getMoments(x)
        Moments_is = EvalMoments_univ(t(FilteredParameters), Dist)
        colnames(Moments_is) = paste("M", 1:4, sep = "")
        PlotMultipleSeries_wis(Moments_is, Moments_os, 4L, iH, vDates_os, vDates_is)
      }
    } else {
      if (is.null(which)) {
        cat(paste("Print 1-4 or 0 to exit"))
        PlotType = menu(PlotMenu(x))
      } else {
        PlotType = which
      }
      if (PlotType == 1L) {
        PlotMultipleSeries(ParametersForecast, iK, iH, vDates_os)
      }
      if (PlotType == 2L) {
        Moments_os = getMoments(x)
        Mu = Moments_os[, 1L]
        mRealVsForecast = cbind(Mu, vOut)
        PlotForecastVsRealized_Univ(mRealVsForecast, vDates_os, x)
      }
      if (PlotType == 3L) {
        Moments_os = getMoments(x)
        PlotMultipleSeries(Moments_os, 4L, iH, vDates_os)
      }
      if (PlotType == 4L) {
        PlotSingleSeries(vLS, iH, vDates_os)
      }
    }
    if (!is.null(which)) {
      PlotType = 0L
    }
  }
})

setMethod("plot", signature(x = "mGASFor", y = "missing"), function(x, which = NULL) {
  iK = x@Info$iK
  iN = x@Info$iN
  mY = x@Data$mY
  mOut = x@Data$mOut
  iH = x@Info$iH
  iT = ncol(mY)

  Roll = x@Info$Roll
  vLS = LogScore(x)

  FilteredParameters = x@Data$FilteredParameters
  FilteredParameters = FilteredParameters[-nrow(FilteredParameters), ]  #remove one step ahead forecast
  ParametersForecast = getForecast(x)
  cBands = x@Bands

  if (is(mY, "xts")) {
    vDates_is = as.Date(index(mY))
    if (Roll) {
      vDates_os = as.Date(index(mOut))
    } else {
      DiffTime = vDates_is[2L] - vDates_is[1L]
      vDates_os = seq(tail(vDates_is, 1L) + DiffTime, tail(vDates_is, 1L) + DiffTime * iH, by = DiffTime)
    }
    ParametersForecast = xts(ParametersForecast, vDates_os)
    FilteredParameters = xts(FilteredParameters, vDates_is)
  } else {
    vDates_is = 1:iT
    vDates_os = (iT + 1L):(iT + iH)
  }

  PlotType = 1L
  while (PlotType > 0L) {
    if (!Roll) {
      if (is.null(which)) {
        cat(paste("Print 1-4 or 0 to exit"))
        PlotType = menu(PlotMenu(x))
      } else {
        PlotType = which
      }
      if (PlotType == 1L) {
        PlotMultipleSeries(ParametersForecast, iK, iH, vDates_os)
      }
      if (PlotType == 2L) {
        PlotMultipleSeries_Bands(ParametersForecast, iK, iH, vDates_os, cBands)
      }
      if (PlotType == 3L) {
        PlotMultipleSeries_wis(FilteredParameters, ParametersForecast, iK, iH, vDates_os, vDates_is)
      }
      if (PlotType == 4L) {
        PlotMultipleSeries_Bands_wis(FilteredParameters, ParametersForecast, iK, iH, vDates_os,
                                     vDates_is, cBands)
      }
    } else {
      if (is.null(which)) {
        cat(paste("Print 1-4 or 0 to exit"))
        PlotType = menu(PlotMenu(x))
      } else {
        PlotType = which
      }
      if (PlotType == 1L) {
        PlotMultipleSeries(ParametersForecast, iK, iT, vDates_os)
      }
      if (PlotType == 2L) {
        Moments_os = getMoments(x)
        mForcasted = Moments_os[["mean"]]
        PlotForecastVsRealized_Multi(t(mOut), mForcasted, iN, vDates_os, x)
      }
      if (PlotType == 3L) {
        Moments_os = getMoments(x)
        mMean = Moments_os[["mean"]]
        colnames(mMean) = colnames(t(mOut))
        cCov = Moments_os[["cov"]]

        PlotMultipleSeries(mMean, iN, iH, vDates_os)
        foo = readline("Print enter to plot covariances\n:")
        PlotCovariances(cCov, iN, iH, vDates_os, colnames(t(mOut)))
      }
      if (PlotType == 4L) {
        PlotSingleSeries(vLS, iH, vDates_os)
      }
    }
    if (!is.null(which)) {
      PlotType = 0L
    }
  }
})

setMethod("plot", signature(x = "uGASRoll", y = "missing"), function(x, which = NULL) {
  iK = x@Info$iK
  vY = x@Data$vY
  iH = x@Info$ForecastLength
  vOut = tail(vY, iH)
  iT = length(vY)

  ParametersForecast = getForecast(x)
  vU = pit(x)

  if (is(vY, "xts")) {
    vDates_os = tail(as.Date(index(vY)), iH)
    ParametersForecast = xts(ParametersForecast, vDates_os)
  } else {
    vDates_os = 1:iH
  }

  PlotType = 1L
  while (PlotType > 0L) {
    if (is.null(which)) {
      cat(paste("Print 1-4 or 0 to exit"))
      PlotType = menu(PlotMenu(x))
    } else {
      PlotType = which
    }
    if (PlotType == 1L) {
      PlotMultipleSeries(ParametersForecast, iK, iT, vDates_os)
    }
    if (PlotType == 2L) {
      Moments_os = getMoments(x)
      Mu = Moments_os[, 1L]
      mRealVsForecast = cbind(Mu, vOut)
      PlotForecastVsRealized_Univ(mRealVsForecast, vDates_os, x)
    }
    if (PlotType == 3L) {
      Moments_os = getMoments(x)
      PlotMultipleSeries(Moments_os, 4L, iH, vDates_os)
    }
    if (PlotType == 4L) {
      PlotPit(vU, x@Testing$PitTest$Hist)
    }
    if (!is.null(which)) {
      PlotType = 0L
    }
  }
})

setMethod("plot", signature(x = "mGASRoll", y = "missing"), function(x, which = NULL) {
  iN = x@Info$iN
  iK = x@Info$iK
  mY = x@Data$mY
  iH = x@Info$ForecastLength
  mOut = t(tail(t(mY), iH))
  iT = ncol(mY)

  ParametersForecast = getForecast(x)

  if (is(mY, "xts")) {
    vDates_os = tail(as.Date(index(mY)), iH)
    ParametersForecast = xts(ParametersForecast, vDates_os)
  } else {
    vDates_os = 1:iH
  }

  PlotType = 1L
  while (PlotType > 0L) {
    if (is.null(which)) {
      cat(paste("Print 1-4 or 0 to exit"))
      PlotType = menu(PlotMenu(x))
    } else {
      PlotType = which
    }
    if (PlotType == 1L) {
      PlotMultipleSeries(ParametersForecast, iK, iT, vDates_os)
    }
    if (PlotType == 2L) {
      Moments_os = getMoments(x)
      mForcasted = Moments_os[["mean"]]
      PlotForecastVsRealized_Multi(t(mOut), mForcasted, iN, vDates_os, x)
    }
    if (PlotType == 3L) {
      Moments_os = getMoments(x)
      mMean = Moments_os[["mean"]]
      colnames(mMean) = colnames(t(mOut))
      cCov = Moments_os[["cov"]]

      PlotMultipleSeries(mMean, iN, iH, vDates_os)
      foo = readline("Print enter to plot covariances\n:")
      PlotCovariances(cCov, iN, iH, vDates_os, colnames(t(mOut)))
    }
    if (!is.null(which)) {
      PlotType = 0L
    }
  }
})


getFilteredParameters = function(object) {
  UseMethod("getFilteredParameters")
}
.getFilteredParameters <- function(object) {
  if (is(object, "uGASFit") | is(object, "mGASFit"))
    mTheta = object@GASDyn$mTheta
  if (is(object, "uGASSim") | is(object, "mGASSim"))
    mTheta = object@GASDyn$mTheta

  mTheta = t(mTheta)
  parNames = getParNames(object)
  colnames(mTheta) = parNames

  return(mTheta)
}
setMethod("getFilteredParameters", signature(object = "uGASFit"), .getFilteredParameters)
setMethod("getFilteredParameters", signature(object = "mGASFit"), .getFilteredParameters)
setMethod("getFilteredParameters", signature(object = "uGASSim"), .getFilteredParameters)
setMethod("getFilteredParameters", signature(object = "mGASSim"), .getFilteredParameters)

getObs = function(object) {
  UseMethod("getObs")
}
.getObs <- function(object) {
  if (is(object, "uGASFit"))
    Data = object@Data$vY
  if (is(object, "mGASFit"))
    Data = object@Data$mY

  if (is(object, "uGASSim"))
    Data = object@Data$vY
  if (is(object, "mGASSim"))
    Data = object@Data$mY

  if (is(object, "uGASFor"))
    Data = object@Data$vY

  if (is(object, "uGASRoll"))
    Data = object@Data$vY

  return(Data)
}
setMethod("getObs", signature(object = "uGASFit"), .getObs)
setMethod("getObs", signature(object = "mGASFit"), .getObs)
setMethod("getObs", signature(object = "uGASSim"), .getObs)
setMethod("getObs", signature(object = "mGASSim"), .getObs)
setMethod("getObs", signature(object = "uGASFor"), .getObs)
setMethod("getObs", signature(object = "uGASRoll"), .getObs)


getMoments = function(object) {
  UseMethod("getMoments")
}
.getMoments <- function(object) {
  if (is(object, "uGASFit")) {
    Moments = object@Estimates$Moments
  }
  if (is(object, "mGASFit")) {
    Moments = object@Estimates$Moments
  }
  if (is(object, "uGASSim")) {
    Moments = object@Data$Moments
  }
  if (is(object, "mGASSim")) {
    Moments = object@Data$Moments
  }

  if (is(object, "uGASFor")) {
    Moments = object@Forecast$Moments
  }
  if (is(object, "mGASFor")) {
    Moments = object@Forecast$Moments
  }
  if (is(object, "uGASRoll")) {
    Moments = object@Forecast$Moments
  }
  if (is(object, "mGASRoll")) {
    Moments = object@Forecast$Moments
  }

  Dist = getDist(object)

  if (!is.null(Moments) & DistType(Dist) != "multivariate") {
    colnames(Moments) = paste("M", 1:4, sep = "")
  }

  return(Moments)
}
setMethod("getMoments", signature(object = "uGASFit"), .getMoments)
setMethod("getMoments", signature(object = "mGASFit"), .getMoments)
setMethod("getMoments", signature(object = "uGASSim"), .getMoments)
setMethod("getMoments", signature(object = "mGASSim"), .getMoments)
setMethod("getMoments", signature(object = "uGASFor"), .getMoments)
setMethod("getMoments", signature(object = "mGASFor"), .getMoments)
setMethod("getMoments", signature(object = "uGASRoll"), .getMoments)
setMethod("getMoments", signature(object = "mGASRoll"), .getMoments)

.getCoef <- function(object, do.list = FALSE) {
  if (is(object, "uGASFit") | is(object, "mGASFit")) {
    if (do.list) {
      ans = list(lCoef = object@Estimates$lParList, mCoef = object@Estimates$Inference)
    } else {
      ans = object@Estimates$Inference[, "Estimate"]
    }
  }
  if (is(object, "uGASSim") | is(object, "mGASSim")) {
    ans = list(vKappa = object@ModelInfo$vKappa, mA = object@ModelInfo$mA, mB = object@ModelInfo$mB)
  }
  if (is(object, "uGASRoll") | is(object, "mGASRoll")) {
    ans = object@Forecast$Coef
    if (!do.list) {
      ans = as.matrix(do.call(rbind, lapply(ans, function(x) x$mCoef[, "Estimate"])))
    }
    rownames(ans) = paste("refit", 1:nrow(ans))
  }

  return(ans)

}

setMethod("coef", signature(object = "uGASFit"), .getCoef)
setMethod("coef", signature(object = "mGASFit"), .getCoef)
setMethod("coef", signature(object = "uGASSim"), .getCoef)
setMethod("coef", signature(object = "mGASSim"), .getCoef)

.getQuantile <- function(x, probs = c(0.01, 0.05)) {

  if (is(x, "uGASFit"))
    mTheta = getFilteredParameters(x)
  if (is(x, "uGASSim"))
    mTheta = getFilteredParameters(x)
  if (is(x, "uGASRoll"))
    mTheta = getForecast(x)

  Dist = getDist(x)

  mQuantile = Quantiles(t(mTheta), Dist, probs)
  colnames(mQuantile) = paste("q.", probs, sep = "")

  return(mQuantile)

}

.getQuantile_Sim <- function(x, probs = c(0.01, 0.05)) {

  bRoll = x@Info$Roll

  iH = x@Info$iH

  if (bRoll) {

    mTheta_tph = getForecast(x)

    mQuantile = matrix(NA, iH, length(probs), dimnames = list(paste("T+", 1:iH, sep = ""),
                                                              paste("q.", probs, sep = "")))

    for (h in 1:iH) {
      mQuantile[h, ] = Quantiles(t(mTheta_tph[h, , drop = FALSE]), Dist = getDist(x), probs)
    }

  } else {

    mDraws = x@Draws

    if (is.null(mDraws) & iH > 1) {
      stop("ReturnDraws = TRUE needs to be selected in the
            UniGASFor function for multistep ahead quantile evaluation.")
    }

    mQuantile = matrix(NA, iH, length(probs), dimnames = list(paste("T+", 1:iH, sep = ""),
                                                              paste("q.", probs, sep = "")))

    ## one step ahead
    vTheta_tp1 = getForecast(x)[1, ,drop = FALSE]
    mQuantile[1, ] = Quantiles(t(vTheta_tp1), Dist = getDist(x), probs)

    ## multi step ahead
    if (iH > 1) {
      for (h in 2:iH) {
        mQuantile[h, ] = quantile(mDraws[, h], probs)
      }
    }

  }

  return(mQuantile)

}

.getES <- function(object, probs = c(0.01, 0.05)) {

  if (is(object, "uGASFit"))
    mTheta = getFilteredParameters(object)
  if (is(object, "uGASSim"))
    mTheta = getFilteredParameters(object)
  if (is(object, "uGASFor"))
    mTheta = getForecast(object)
  if (is(object, "uGASRoll"))
    mTheta = getForecast(object)

  Dist = getDist(object)

  mQuantile = Quantiles(t(mTheta), Dist, probs)
  colnames(mQuantile) = paste("q.", probs, sep = "")

  mES = mQuantile

  for (i in 1:nrow(mES)) {
    for (j in 1:ncol(mES)) {
      mES[i, j] = adaptIntegrate(Quantiles, lower = 1e-7, upper = probs[j],
                            mTheta = t(mTheta[i, , drop = FALSE]), Dist = Dist)$integral

    }
  }

  mES = t(t(mES) / probs)

  return(mES)

}

.getES_Sim <- function(object, probs = c(0.01, 0.05)) {

  iH = object@Info$iH

  bRoll = x@Info$Roll

  iH = x@Info$iH

  Dist = getDist(object)

  if (bRoll) {

    mTheta_tph = getForecast(x)

    mQuantile = quantile(object, probs)

    mES = mQuantile

    vTheta_tp1 = getForecast(object)[1, ,drop = FALSE]

    for (h in 1:iH) {
      for (j in 1:ncol(mES)) {
        mES[h, j] = integrate(Quantiles, lower = 1e-7, upper = probs[j],
                              mTheta = t(mTheta_tph[h, , drop = FALSE]), Dist = Dist)$value/probs[j]
      }
    }

  } else {

    mDraws = object@Draws

    if (is.null(mDraws) & iH > 1) {
      stop("ReturnDraws = TRUE needs to be selected in the
         UniGASFor function for multistep ahead quantile evaluation.")
    }

    mQuantile = quantile(object, probs)
    mES = mQuantile

    ## one step ahead
    vTheta_tp1 = getForecast(object)[1, ,drop = FALSE]
    for (j in 1:ncol(mES)) {
      mES[1, j] = integrate(Quantiles, lower = 1e-7, upper = probs[j],
                            mTheta = t(vTheta_tp1), Dist = Dist)$value/probs[j]
    }

    ## multi step ahead
    if (iH > 1) {
      for (h in 2:iH) {
        vDraws = mDraws[, h]
        for (j in 1:ncol(mES)) {
          mES[h, j] = mean(vDraws[vDraws < mQuantile[h, j]])
        }
      }
    }
  }

  return(mES)

}

setMethod("quantile", signature(x = "uGASFit"), .getQuantile)
setMethod("quantile", signature(x = "uGASSim"), .getQuantile)
setMethod("quantile", signature(x = "uGASFor"), .getQuantile_Sim)
setMethod("quantile", signature(x = "uGASRoll"), .getQuantile)

ES = function(object, ...) {
  UseMethod("ES")
}

setMethod("ES", signature(object = "uGASFit"), .getES)
setMethod("ES", signature(object = "uGASSim"), .getES)
setMethod("ES", signature(object = "uGASFor"), .getES_Sim)
setMethod("ES", signature(object = "uGASRoll"), .getES)

pit = function(object) {
  UseMethod("pit")
}

setMethod("pit", signature(object = "uGASFit"), function(object) object@Estimates$vU)
setMethod("pit", signature(object = "uGASFor"), function(object) object@Forecast$vU)
setMethod("pit", signature(object = "uGASRoll"), function(object) object@Forecast$vU)


LogScore = function(object) {
  UseMethod("LogScore")
}

setMethod("LogScore", signature(object = "uGASFor"), function(object) object@Forecast$vLS)
setMethod("LogScore", signature(object = "mGASFor"), function(object) object@Forecast$vLS)
setMethod("LogScore", signature(object = "uGASRoll"), function(object) object@Forecast$vLS)
setMethod("LogScore", signature(object = "mGASRoll"), function(object) object@Forecast$vLS)


getForecast = function(object) {
  UseMethod("getForecast")
}
setMethod("getForecast", signature(object = "uGASFor"), function(object) return(object@Forecast$PointForecast))
setMethod("getForecast", signature(object = "mGASFor"), function(object) return(object@Forecast$PointForecast))
setMethod("getForecast", signature(object = "uGASRoll"), function(object) return(object@Forecast$PointForecast))
setMethod("getForecast", signature(object = "mGASRoll"), function(object) return(object@Forecast$PointForecast))


getPwNames = function(object) {
  UseMethod("getPwNames")
}
setMethod("getPwNames", signature(object = "uGASSpec"), function(object) return(object@Spec$PwNames))
setMethod("getPwNames", signature(object = "mGASSpec"), function(object) return(object@Spec$PwNames))


residuals = function(object, ...) {
  UseMethod("residuals")
}

setMethod("residuals", signature(object = "uGASFit"), function(object, standardize = TRUE) {

  vY = object@Data$vY
  iT = length(vY)

  mMoments = getMoments(object)[1:iT, ]
  vRes     = vY - mMoments[, 1L]

  if (standardize) {
    vRes = vRes/sqrt(mMoments[, 2L])
  }

  return(vRes)

})
setMethod("residuals", signature(object = "mGASFit"), function(object, standardize = TRUE) {

  mY = object@Data$mY
  iT = ncol(mY)

  lMoments = getMoments(object)
  mMean = lMoments$mean[1:iT, ]

  mRes     = t(mY) - mMean

  if (standardize) {
    aCov = lMoments$cov
    for (i in 1:iT) {
      mRes[i, ] = solve(chol(aCov[,, i])) %*% t(mRes[i, ,drop = FALSE])
    }
  }

  return(mRes)

})
setMethod("residuals", signature(object = "uGASRoll"), function(object, standardize = TRUE) {

  iH = object@Info$ForecastLength
  vY = tail(object@Data$vY, iH)

  mMoments = getMoments(object)
  vRes     = vY - mMoments[, 1L]

  if (standardize) {
    vRes = vRes/sqrt(mMoments[, 2L])
  }

  return(vRes)

})
setMethod("residuals", signature(object = "mGASRoll"), function(object, standardize = TRUE) {

  iH = object@Info$ForecastLength
  mY = tail(t(object@Data$mY), iH)

  lMoments = getMoments(object)
  mMean = lMoments$mean

  mRes     = mY - mMean

  if (standardize) {
    aCov = lMoments$cov
    for (i in 1:iH) {
      mRes[i, ] = solve(chol(aCov[,, i])) %*% t(mRes[i, ,drop = FALSE])
    }
  }

  return(mRes)

})

convergence = function(object) {
  UseMethod("convergence")
}
setMethod("convergence", signature(object = "uGASFit"), function(object) return(object@ModelInfo$convergence))
setMethod("convergence", signature(object = "mGASFit"), function(object) return(object@ModelInfo$convergence))

