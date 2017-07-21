
################################ DENSITY #

############ PIT#################

BinTest <- function(pit, g = 20L, alpha = 0.05, plot = FALSE) {

    h = hist(pit, nclass = g, plot = FALSE)
    n_i = h$counts
    test = sum((n_i - mean(n_i))^2/mean(n_i))
    crit = qchisq(1.0 - alpha, g - 1L)
    pvalue = 1 - pchisq(test, g - 1L)
    confidence = mean(n_i) + c(-qnorm(1 - alpha) * sqrt(mean(n_i)), +qnorm(1 - alpha) * sqrt(mean(n_i)))

    if (plot) {
        plot(h, col = "blue", ylim = c(0, max(confidence * 1.2, n_i * 1.2)))
        abline(h = confidence, col = "red", lwd = 2L, xlim = c(0, 1))
    }

    out = list(test = test, crit = crit, pvalue = pvalue, hist = h, confidence = confidence)
    return(out)
}

iidTest <- function(pit, alpha = 0.05) {

    N = length(pit)

    m1 = as.numeric(pit - mean(pit))
    m2 = as.numeric((pit - mean(pit))^2)
    m3 = as.numeric((pit - mean(pit))^3)
    m4 = as.numeric((pit - mean(pit))^4)

    data1 = do.call(cbind, lapply(1:20, function(i) c(m1[-(1:i)], rep(NA, i))))
    data1 = data.frame(head(data1, N - 20L))
    data2 = do.call(cbind, lapply(1:20, function(i) c(m2[-(1:i)], rep(NA, i))))
    data2 = data.frame(head(data1, N - 20L))
    data3 = do.call(cbind, lapply(1:20, function(i) c(m3[-(1:i)], rep(NA, i))))
    data3 = data.frame(head(data1, N - 20L))
    data4 = do.call(cbind, lapply(1:20, function(i) c(m4[-(1:i)], rep(NA, i))))
    data4 = data.frame(head(data1, N - 20L))

    m1 = head(m1, N - 20L)
    m2 = head(m2, N - 20L)
    m3 = head(m3, N - 20L)
    m4 = head(m4, N - 20L)

    fit1 = lm(m1 ~ ., data = data1)
    fit2 = lm(m2 ~ ., data = data2)
    fit3 = lm(m3 ~ ., data = data3)
    fit4 = lm(m4 ~ ., data = data4)

    test1 = (N - 20) * summary(fit1)$r.squared
    test2 = (N - 20) * summary(fit2)$r.squared
    test3 = (N - 20) * summary(fit3)$r.squared
    test4 = (N - 20) * summary(fit4)$r.squared

    crit = qchisq(1 - alpha, 20L)

    pvalue1 = 1 - pchisq(test1, 20L)
    pvalue2 = 1 - pchisq(test2, 20L)
    pvalue3 = 1 - pchisq(test3, 20L)
    pvalue4 = 1 - pchisq(test4, 20L)


    out = list(test = c(test1 = test1, test2 = test2, test3 = test3, test4 = test4), crit = crit, pvalue = c(pvalue1,
        pvalue2, pvalue3, pvalue4))
    return(out)
}

PIT_test <- function(U, G = 20L, alpha = 0.05, plot = FALSE) {

    iG = G
    dAlpha = alpha
    vU = U

    if (length(vU) < 100L)
        iG = 5L

    Hist = BinTest(pit = vU, g = iG, alpha = dAlpha, plot = plot)
    IID = iidTest(pit = vU, alpha = dAlpha)

    return(list(Hist = Hist, IID = IID))
}

############################### QUANTILE #

DQOOStest <- function(y, VaR, tau, cLags) {

    cT = length(y)
    vHit = numeric(cT)
    vHit[y < VaR] = 1 - tau
    vHit[y > VaR] = -tau

    vConstant = rep(1, (cT - cLags))
    vHIT = vHit[(cLags + 1):cT]
    vVaRforecast = VaR[(cLags + 1):cT]
    mZ = matrix(0, cT - cLags, cLags)


    for (st in 1:cLags) {

        mZ[, st] = vHit[st:(cT - (cLags + 1L - st))]

    }

    mX = cbind(vConstant, vVaRforecast, mZ)

    dDQstatOut = (t(vHIT) %*% mX %*% ginv(t(mX) %*% mX) %*% t(mX) %*% (vHIT))/(tau * (1 - tau))

    dDQpvalueOut = 1 - pchisq(dDQstatOut, ncol(mX))

    out = list(stat = dDQstatOut, pvalue = dDQpvalueOut)
}

HitSequence <- function(returns_X, VaR_X) {
    N = length(returns_X)
    Hit_X = numeric(N)
    Hit_X[which(returns_X <= VaR_X)] = 1L
    return(Hit_X)
}

Kupiec <- function(Hit, tau) {
    N = length(Hit)
    x = sum(Hit)
    rate = x/N
    test = -2 * log(((1 - tau)^(N - x) * tau^x)/((1 - rate)^(N - x) * rate^x))
    if (is.nan(test))
        test = -2 * ((N - x) * log(1 - tau) + x * log(tau) - (N - x) * log(1 - rate) - x * log(rate))
    # threshold = qchisq(alphaTest, df = 1)
    pvalue = 1 - pchisq(test, df = 1)

    LRpof = c(test, pvalue)
    names(LRpof) = c("Test", "Pvalue")
    return(LRpof)
}

Christoffersen <- function(Hit, tau) {
    n00 = n01 = n10 = n11 = 0
    N = length(Hit)
    for (i in 2:N) {
        if (Hit[i] == 0L & Hit[i - 1L] == 0L)
            n00 = n00 + 1
        if (Hit[i] == 0L & Hit[i - 1L] == 1L)
            n01 = n01 + 1
        if (Hit[i] == 1L & Hit[i - 1L] == 0L)
            n10 = n10 + 1
        if (Hit[i] == 1L & Hit[i - 1L] == 1L)
            n11 = n11 + 1
    }
    pi0 = n01/(n00 + n01)
    pi1 = n11/(n10 + n11)
    pi = (n01 + n11)/(n00 + n01 + n10 + n11)
    LRind = -2 * log(((1 - pi)^(n00 + n10) * pi^(n01 + n11))/((1 - pi0)^n00 * pi0^n01 * (1 - pi1)^n10 *
        pi1^n11))
    if (is.nan(LRind))
        LRind = -2 * ((n00 + n10) * log(1 - pi) + (n01 + n11) * log(pi) - n00 * log(1 - pi0) - n01 *
            log(pi0) - n10 * log(1 - pi1) - n11 * log(pi1))
    LRpof = Kupiec(Hit, tau)["Test"]
    LRcc = LRpof + LRind
    pvalue = 1 - pchisq(LRcc, df = 2L)
    LRcc = c(LRcc, pvalue)
    names(LRcc) = c("Test", "Pvalue")
    return(LRcc)
}

ActualOverExpected <- function(Hit, tau) {
    N = length(Hit)
    x = sum(Hit)
    Actual = x
    AovE = Actual/(tau * N)
}

AbsoluteDeviation <- function(Hit, returns_X, VaR_X) {
    series = abs(VaR_X - returns_X)
    series = series[which(Hit == 1L)]
    ADmean = mean(series)
    ADmax = max(series)

    out = c(ADmean, ADmax)
    names(out) = c("ADmean", "ADmax")
    return(out)
}

QLoss <- function(vY, vVaR, dTau) {
    vHit = HitSequence(vY, vVaR)
    vLoss = (vY - vVaR) * (dTau - vHit)
    dLoss = mean(vLoss)

    return(list(Loss = dLoss, LossSeries = vLoss))
}

############ Gaussianity

JarqueBera <- function(vRes, dAlpha = 0.05) {

  iT = length(vRes)

  m1 = mean(vRes)
  m2 = mean((vRes - m1)^2.0)
  m3 = mean((vRes - m1)^3.0)
  m4 = mean((vRes - m1)^4.0)

  dS = m3 / m2^(3.0 / 2.0)
  dK = m4 / m2^2.0

  dStat = iT / 6.0 * (dS^2.0 + (dK - 3.0)^2.0 / 4.0)

  dPval = 1.0 - pchisq(dStat, 2)

  dCritical = qchisq(1.0 - dAlpha, 2L)

  return(c("Statistic" = dStat, "p-Value" = dPval, "critical" = dCritical))

}

LjungBox <- function(vRes, vLag = c(10, 15, 20)) {

  iP    = length(vLag)
  mTest = matrix(data = NA, iP, 2L, dimnames =
                   list(vLag, c("Statistic", "p-Value")))

  for (p in 1:iP) {
    Test = Box.test(vRes, vLag[p], type = "Ljung-Box")
    mTest[p, ] = c(Test$statistic, Test$p.value)
  }

  return(mTest)

}




