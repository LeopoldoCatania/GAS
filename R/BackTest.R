BacktestVaR <- function(vY, vVaR, dTau, alphaTest = 0.95, cLags = 4) {
  vY = as.numeric(vY)
  vVaR = as.numeric(vVaR)

  Hit = HitSequence(vY, vVaR)
  LRuc = Kupiec(Hit, dTau, alphaTest = alphaTest)
  LRcc = Christoffersen(Hit, dTau, alphaTest = alphaTest)
  DQ = DQOOStest(vY, vVaR, dTau, cLags)
  AE = ActualOverExpected(Hit, dTau)
  AD = AbsoluteDeviation(Hit, vY, vVaR)
  Loss = QLoss(vY, vVaR, dTau)

  lOut = list(LRuc = LRuc, LRcc = LRcc, AE = AE, AD = AD, DQ = DQ, Loss = Loss)

  return(lOut)
}

BacktestDensity <- function(Roll, dLower, dUpper, iK = 1000L, dA = 0.0, dB = 1.0){

  iH = Roll@Info$ForecastLength
  vY = tail(getObs(Roll), iH)
  Dist = getDist(Roll)
  mTheta = getForecast(Roll)

  iT  = length(vY)
  vLS = EvaluateLogScore_Univ(t(mTheta), vY, Dist, iT)

  mWCRPS = mWCRPS_backtest(vY, t(mTheta), Dist, dLower, dUpper, iK, dA, dB)
  colnames(mWCRPS) = c("uniform", "center", "tails", "tail_r", "tail_l")

  vAvg = c(LS = -mean(vLS), apply(mWCRPS, 2, mean))

  lOut = list()
  lOut[["series"]]  = list(LS = vLS, WCRPS = mWCRPS )
  lOut[["average"]] = vAvg

  return(lOut);
}

