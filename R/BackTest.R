BacktestVaR <- function(data, VaR, alpha, Lags = 4) {

    vY = data
    vVaR = VaR
    dTau = alpha

    vY = as.numeric(vY)
    vVaR = as.numeric(vVaR)

    Hit = HitSequence(vY, vVaR)
    LRuc = Kupiec(Hit, dTau)
    LRcc = Christoffersen(Hit, dTau)
    DQ = DQOOStest(vY, vVaR, dTau, Lags)
    AE = ActualOverExpected(Hit, dTau)
    AD = AbsoluteDeviation(Hit, vY, vVaR)
    Loss = QLoss(vY, vVaR, dTau)

    lOut = list(LRuc = LRuc, LRcc = LRcc, AE = AE, AD = AD, DQ = DQ, Loss = Loss)

    return(lOut)
}

BacktestDensity <- function(Roll, lower, upper, K = 1000, a = 0, b = 1) {

    dLower = lower
    dUpper = upper
    iK = K
    dA = a
    dB = b

    iH = Roll@Info$ForecastLength
    vY = tail(getObs(Roll), iH)
    Dist = getDist(Roll)
    mTheta = getForecast(Roll)

    iT = length(vY)
    vLS = EvaluateLogScore_Univ(t(mTheta), vY, Dist, iT)

    mWCRPS = mWCRPS_backtest(vY, t(mTheta), Dist, dLower, dUpper, iK, dA, dB)
    colnames(mWCRPS) = c("uniform", "center", "tails", "tail_r", "tail_l")

    vAvg = c(NLS = -mean(vLS), apply(mWCRPS, 2, mean))

    lOut = list()
    lOut[["series"]] = list(LS = vLS, WCRPS = mWCRPS)
    lOut[["average"]] = vAvg

    return(lOut)
}

