UniGASFor <- function(uGASFit, H = NULL, Roll = FALSE, out = NULL, B = 10000, Bands = c(0.1, 0.15, 0.85, 0.9), ReturnDraws = FALSE) {

    vOut = out
    if (Roll) {
        H = length(vOut)
    }

    iH = H
    iB = B

    vBands = Bands
    bReturnDraws = ReturnDraws

    iK = uGASFit@ModelInfo$iK
    ScalingType = getScalingType(uGASFit)
    lParList = coef(uGASFit, do.list = TRUE)$lCoef
    Dist = getDist(uGASFit)
    GASPar = getGASPar(uGASFit)
    vY = getObs(uGASFit)

    FilteredParameters = getFilteredParameters(uGASFit)

    if (Roll) {

        iT = length(vY)
        vYf = c(vY, vOut)
        iT = iT + iH

        GASDyn = GASFilter_univ(vYf, lParList$vKappa, lParList$mA, lParList$mB, iT, iK, Dist, ScalingType)

        PointForecast = t(GASDyn$mTheta[, (iT - iH + 1):(iT)])
        cBands = array(0, dim = c(1, 1, 1))
        mY = matrix(0, 1, 1)
        vU = EvaluatePit_Univ(t(PointForecast), vOut, Dist, iH)
        vLS = GASDyn$vLLK[(iT - iH + 1):(iT)]

    } else {

        vTheta_tp1 = tail(getFilteredParameters(uGASFit), 1)
        Forc = uGASMultiForcast(vTheta_tp1, lParList$vKappa, lParList$mA, lParList$mB, iH, iB, iK, Dist,
            ScalingType, bReturnDraws)

        PointForecast = matrix(0, iH, iK)
        cBands = array(0, dim = c(iH, length(vBands), iK), dimnames = list(1:iH, paste("q.", vBands,
            sep = ""), colnames(vTheta_tp1)))

        if (iH > 1) {
            for (k in 1:iK) {
                PointForecast[, k] = apply(Forc$cTheta[k, , ], 1, function(x) median(na.omit(x)))
                for (q in vBands) {
                  cBands[, paste("q.", q, sep = ""), k] = apply(Forc$cTheta[k, , ], 1, function(x, q) quantile(na.omit(x),
                    probs = q), q = q)
                }
            }
        } else {
            PointForecast[1, ] = vTheta_tp1
            for (k in 1:iK) cBands[, , k] = vTheta_tp1[k]
        }
        if (bReturnDraws) {
            mY = Forc$mY
        } else {
            mY = matrix(0, 1, 1)
        }

        vU = vLS = NULL
    }

    mMoments = EvalMoments_univ(t(PointForecast), Dist)

    colnames(PointForecast) = names(GASPar)
    rownames(PointForecast) = paste("T+", 1:iH, sep = "")
    rownames(mMoments)      = paste("T+", 1:iH, sep = "")

    if ( Roll &  any(class(vOut)[1] == c("zoo", "ts", "xts") | !is.null(names(vOut)))) {
      if ( any(class(vOut)[1] == c("zoo", "ts", "xts"))) {
        PointForecast = xts(PointForecast, order.by = index(vOut))
        mMoments      = xts(mMoments, order.by = index(vOut))
      } else {
        rownames(PointForecast) = names(vOut)
        rownames(mMoments)      = names(vOut)
      }
    }

    Out <- new("uGASFor",
               Forecast = list(PointForecast = PointForecast,
                               Moments = mMoments,
                               vU = vU,
                               vLS = vLS),
               Bands = cBands,
               Draws = mY,
               Info = list(iH = iH,
                           Roll = Roll,
                           iB = iB,
                           vBands = vBands,
                           bReturnDraws = bReturnDraws,
                           GASPar = GASPar,
                           Dist = Dist,
                           ScalingType = ScalingType,
                           iK = iK),
               Data = list(vY = vY,
                           FilteredParameters = FilteredParameters,
                           vOut = vOut))

    return(Out)

}

MultiGASFor <- function(mGASFit, H = NULL, Roll = FALSE, out = NULL, B = 10000, Bands = c(0.1, 0.15, 0.85, 0.9),
    ReturnDraws = FALSE) {

    if (Roll) {

      if (is.null(out)) stop("When Roll = TRUE, 'out' must be submitted")

      mOut = t(out)
      H = ncol(mOut)

      if (any(class(out)[1] == c("ts", "zoo", "xts"))) {
        vOutDate = index(out)
      } else {
        vOutDate = paste("T+", 1:H, sep = "")
      }

    } else {
        mOut = NULL
        vOutDate = paste("T+", 1:H, sep = "")
    }

    iH = H
    iB = B

    vBands = Bands
    bReturnDraws = ReturnDraws

    iK = mGASFit@ModelInfo$iK
    iN = mGASFit@ModelInfo$iN
    ScalingType = getScalingType(mGASFit)
    lParList = coef(mGASFit, do.list = TRUE)$lCoef
    Dist = getDist(mGASFit)
    GASPar = getGASPar(mGASFit)
    mY = getObs(mGASFit)

    FilteredParameters = getFilteredParameters(mGASFit)

    if (Roll) {

        iH = ncol(mOut)
        iT = ncol(mY)
        mYf = cbind(mY, mOut)
        iT = iT + iH

        GASDyn = GASFilter_multi(mYf, lParList$vKappa, lParList$mA, lParList$mB, iT, iN, iK, Dist, ScalingType)

        PointForecast = t(GASDyn$mTheta[, (iT - iH + 1):(iT)])
        cBands = array(0, dim = c(1, 1, 1))
        cY = array(0, dim = c(1, 1, 1))
        vLS = GASDyn$vLLK[(iT - iH + 1):(iT)]

    } else {

        vTheta_tp1 = tail(getFilteredParameters(mGASFit), 1)
        Forc = mGASMultiForcast(vTheta_tp1, lParList$vKappa, lParList$mA, lParList$mB, iH, iB, iK, iN,
            Dist, ScalingType, bReturnDraws)

        PointForecast = matrix(0, iH, iK)
        cBands = array(0, dim = c(iH, length(vBands), iK), dimnames = list(1:iH, paste("q.", vBands,
            sep = ""), colnames(vTheta_tp1)))

        if (iH > 1) {
            for (k in 1:iK) {
                PointForecast[, k] = apply(Forc$cTheta[k, , ], 1, function(x) median(na.omit(x)))
                for (q in vBands) {
                  cBands[, paste("q.", q, sep = ""), k] = apply(Forc$cTheta[k, , ], 1, function(x, q) quantile(na.omit(x),
                    probs = q), q = q)
                }
            }
        } else {
            PointForecast[1, ] = vTheta_tp1
            for (k in 1:iK) cBands[, , k] = vTheta_tp1[k]
        }

        if (bReturnDraws) {
            cY = Forc$cY
        } else {
            cY = array(0, dim = c(1, 1, 1))
        }

        vLS = NULL
    }

    colnames(PointForecast) = FullNamesMulti(iN, Dist)
    rownames(PointForecast) = vOutDate

    if (any(class(out)[1] == c("xts", "ts", "zoo"))) {
      PointForecast = xts(PointForecast, order.by = index(out))
    }

    lMoments = EvalMoments_multi(t(PointForecast), Dist, iN)

    Out <- new("mGASFor",
               Forecast = list(PointForecast = PointForecast,
                               Moments = lMoments,
                               vLS = vLS),
               Bands = cBands,
               Draws = cY,
               Info = list(iH = iH,
                           iN = iN,
                           Roll = Roll,
                           iB = iB,
                           vBands = vBands,
                           bReturnDraws = bReturnDraws,
                           GASPar = GASPar,
                           Dist = Dist,
                           ScalingType = ScalingType,
                           iK = iK),
               Data = list(mY = mY,
                           FilteredParameters = FilteredParameters,
                           mOut = mOut))

    return(Out)

}

UniGASRoll <- function(data, GASSpec, ForecastLength = 500, Nstart = NULL, RefitEvery = 23, RefitWindow = c("moving",
    "recursive"), cluster = NULL, Compute.SE = FALSE, ...) {

    StartTime = Sys.time()

    vY = data

    iT = length(vY)

    Dist = getDist(GASSpec)
    iK = NumberParameters(Dist)

    if (!is.null(cluster)) {
        clusterEvalQ(cluster, {
            library(GAS)
        })
    }
    if (!is.null(Nstart)) {
        iStart = Nstart
    } else {
        iStart = iT - ForecastLength
    }

    FitIndex = seq(iStart, iT, RefitEvery)
    if (tail(FitIndex, 1) == iT) {
        FitIndex = FitIndex[-length(FitIndex)]
    }

    lFits = list()
    lForecasts = list()
    lData = list()
    lOut = list()

    if (RefitWindow[1] == "recursive") {
        for (i in 1:length(FitIndex)) {
            lData[[i]] = vY[1:FitIndex[i]]
        }
    }
    if (RefitWindow[1] == "moving") {
        for (i in 1:length(FitIndex)) {
            lData[[i]] = vY[(FitIndex[i] - iStart + 1):FitIndex[i]]
        }
    }
    # fits
    if (is.null(cluster)) {
        lFits = lapply(lData, UniGASFit, GASSpec = GASSpec, ... = ..., Compute.SE = Compute.SE)
    }
    if (!is.null(cluster)) {
        lFits = parLapply(cluster, lData, UniGASFit, GASSpec = GASSpec, ... = ..., Compute.SE = Compute.SE)
    }

    # coef
    lCoef = lapply(lFits, coef)

    if (RefitEvery == 1) {

        mForc = do.call(rbind, lapply(lFits, function(x) tail(getFilteredParameters(x), 1)))
        vU = EvaluatePit_Univ(t(mForc), tail(vY, ForecastLength), Dist, ForecastLength)
        vLS = EvaluateLogScore_Univ(t(mForc), tail(vY, ForecastLength), Dist, ForecastLength)
        Moments = EvalMoments_univ(t(mForc), Dist)

    } else {

        for (i in 1:length(FitIndex)) {
            if (i != length(FitIndex)) {
                lOut[[i]] = vY[(FitIndex[i] + 1):(FitIndex[i + 1])]
            } else {
                lOut[[i]] = vY[(FitIndex[i] + 1):iT]
            }
        }

        lForecasts = lapply(1:length(lOut), function(i, lFits, lOut) {

            UniGASFor(uGASFit = lFits[[i]], out = lOut[[i]], Roll = TRUE)

        }, lFits = lFits, lOut = lOut)

        mForc   = do.call(rbind, lapply(lForecasts, getForecast))
        vU      = do.call(c, lapply(lForecasts, pit))
        vLS     = do.call(c, lapply(lForecasts, LogScore))
        Moments = do.call(rbind, lapply(lForecasts, getMoments))
    }

    if (all(class(vY)[1] != c("zoo", "ts", "xts"))) {
      rownames(mForc)  = paste("T+", 1:ForecastLength, sep = "")
      rownames(Moments) = paste("T+", 1:ForecastLength, sep = "")
    }

    PitTest = PIT_test(vU, G = 20, alpha = 0.05, plot = FALSE)

    elapsedTime = Sys.time() - StartTime

    Out <- new("uGASRoll",
               Forecast = list(PointForecast = mForc,
                               vU = vU,
                               vLS = vLS,
                               Moments = Moments,
                               Coef = lCoef),
               Info = list(GASSpec = GASSpec,
                           ForecastLength = ForecastLength,
                           RefitEvery = RefitEvery,
                           RefitWindow = RefitWindow[1],
                           iT = iT,
                           iK = iK,
                           elapsedTime = elapsedTime),
               Testing = list(PitTest = PitTest),
               Data = list(vY = vY))

    return(Out)

}

MultiGASRoll <- function(data, GASSpec, ForecastLength = 500, Nstart = NULL, RefitEvery = 23, RefitWindow = c("moving",
    "recursive"), cluster = NULL, Compute.SE = FALSE, ...) {

    StartTime = Sys.time()

    mY = t(data)

    iT = ncol(mY)
    iN = nrow(mY)

    Dist = getDist(GASSpec)
    iK = NumberParameters(Dist, iN)

    if (!is.null(cluster))
        clusterEvalQ(cluster, {
            library(GAS)
        })
    if (!is.null(Nstart)) {
        iStart = Nstart
    } else {
        iStart = iT - ForecastLength
    }

    FitIndex = seq(iStart, iT, RefitEvery)
    if (tail(FitIndex, 1) == iT)
        FitIndex = FitIndex[-length(FitIndex)]

    lFits = list()
    lForecasts = list()
    lData = list()
    lOut = list()

    if (RefitWindow[1] == "recursive") {
        for (i in 1:length(FitIndex)) {
            lData[[i]] = t(mY[, 1:FitIndex[i]])
        }
    }
    if (RefitWindow[1] == "moving") {
        for (i in 1:length(FitIndex)) {
            lData[[i]] = t(mY[, (FitIndex[i] - iStart + 1):FitIndex[i]])
        }
    }
    # fits
    if (is.null(cluster)) {
      lFits = lapply(lData, MultiGASFit, GASSpec = GASSpec, ... = ..., Compute.SE = Compute.SE)
    }
    if (!is.null(cluster)) {
        lFits = parLapply(cluster, lData, MultiGASFit, GASSpec = GASSpec, ... = ..., Compute.SE = Compute.SE)
    }

    # coef
    lCoef = lapply(lFits, coef)

    if (RefitEvery == 1) {

        mForc = do.call(rbind, lapply(lFits, function(x) tail(getFilteredParameters(x), 1)))
        vU = NULL
        vLS = EvaluateLogScore_Multi(t(mForc), t(tail(t(mY), ForecastLength)), Dist, ForecastLength,
            iN)
        Moments = EvalMoments_multi(t(mForc), Dist, iN)
    } else {

        for (i in 1:length(FitIndex)) {
            if (i != length(FitIndex)) {
                lOut[[i]] = t(mY[, (FitIndex[i] + 1):(FitIndex[i + 1])])
            } else {
                lOut[[i]] = t(mY[, (FitIndex[i] + 1):iT])
            }
        }

        lForecasts = lapply(1:length(lOut), function(i, lFits, lOut) {

            MultiGASFor(mGASFit = lFits[[i]], out = lOut[[i]], Roll = TRUE)

        }, lFits = lFits, lOut = lOut)

        mForc = do.call(rbind, lapply(lForecasts, getForecast))
        vU = NULL
        vLS = do.call(c, lapply(lForecasts, LogScore))
        Moments = EvalMoments_multi(t(mForc), Dist, iN)
    }

    if (all(class(mY)[1] != c("zoo", "ts", "xts"))) {
      rownames(mForc)  = paste("T+", 1:ForecastLength, sep = "")
      rownames(Moments) = paste("T+", 1:ForecastLength, sep = "")
    }

    elapsedTime = Sys.time() - StartTime

    Out <- new("mGASRoll",
               Forecast = list(PointForecast = mForc,
                               vU = vU,
                               vLS = vLS,
                               Moments = Moments,
                               Coef = lCoef),
               Info = list(GASSpec = GASSpec,
                           ForecastLength = ForecastLength,
                           RefitEvery = RefitEvery,
                           RefitWindow = RefitWindow[1],
                           iT = iT,
                           iK = iK,
                           iN = iN,
                           elapsedTime = elapsedTime),
               Data = list(mY = mY))

    return(Out)
}

