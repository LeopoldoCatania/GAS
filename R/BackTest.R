BacktestVaR <- function(data, VaR, alpha, Lags = 4L) {

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

BacktestDensity <- function(Roll, lower, upper, K = 1000L, a = NULL, b = NULL) {

  dLower = lower
  dUpper = upper
  iK = K

  iH = Roll@Info$ForecastLength

  vY     = getObs(Roll)
  vY_oos = tail(vY, iH)
  vY_is  = vY[1:(length(vY) - iH)]

  Dist = getDist(Roll)
  mTheta = getForecast(Roll)

  if (is.null(a)) {
    a = mean(vY_is)
  }
  if (is.null(b)) {
    b = sd(vY_is)
  }

  dA = a
  dB = b

  iT = length(vY_oos)
  vLS = EvaluateLogScore_Univ(t(mTheta), vY_oos, Dist, iT)

  mWCRPS = mWCRPS_backtest(vY_oos, t(mTheta), Dist, dLower, dUpper, iK, dA, dB)
  colnames(mWCRPS) = c("uniform", "center", "tails", "tail_l", "tail_r")

  vAvg = c(NLS = -mean(vLS), apply(mWCRPS, 2L, mean))

  lOut = list()
  lOut[["series"]] = list(LS = vLS, WCRPS = mWCRPS)
  lOut[["average"]] = vAvg

  return(lOut)
}

FZLoss <- function(data, VaR, ES, alpha) {

  vY = data
  vVaR = VaR
  dTau = alpha
  vES  = ES

  vY = as.numeric(vY)
  vVaR = as.numeric(vVaR)
  vES  = as.numeric(vES)

  vHit = HitSequence(vY, vVaR)

  vLoss = -vHit/(dTau * vES) * (vVaR - vY) + vVaR/vES + log(-vES) - 1.0

  return(vLoss)

}

