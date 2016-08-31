
################################ DENSITY #

############ PIT#################

BinTest <- function(pit, g = 20, alpha = 0.05, plot = F) {

    h = hist(pit, nclass = g, plot = F)
    n_i = h$counts
    test = sum((n_i - mean(n_i))^2/mean(n_i))
    crit = qchisq(1 - alpha, g - 1)
    pvalue = 1 - pchisq(test, g - 1)
    confidence = mean(n_i) + c(-qnorm(1 - alpha) * sqrt(mean(n_i)), +qnorm(1 - alpha) * sqrt(mean(n_i)))

    if (plot) {
        plot(h, col = "blue", ylim = c(0, max(confidence * 1.2, n_i * 1.2)))
        abline(h = confidence, col = "red", lwd = 2, xlim = c(0, 1))
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
    data1 = data.frame(head(data1, N - 20))
    data2 = do.call(cbind, lapply(1:20, function(i) c(m2[-(1:i)], rep(NA, i))))
    data2 = data.frame(head(data1, N - 20))
    data3 = do.call(cbind, lapply(1:20, function(i) c(m3[-(1:i)], rep(NA, i))))
    data3 = data.frame(head(data1, N - 20))
    data4 = do.call(cbind, lapply(1:20, function(i) c(m4[-(1:i)], rep(NA, i))))
    data4 = data.frame(head(data1, N - 20))

    m1 = head(m1, N - 20)
    m2 = head(m2, N - 20)
    m3 = head(m3, N - 20)
    m4 = head(m4, N - 20)

    fit1 = lm(m1 ~ ., data = data1)
    fit2 = lm(m2 ~ ., data = data2)
    fit3 = lm(m3 ~ ., data = data3)
    fit4 = lm(m4 ~ ., data = data4)

    test1 = (N - 20) * summary(fit1)$r.squared
    test2 = (N - 20) * summary(fit2)$r.squared
    test3 = (N - 20) * summary(fit3)$r.squared
    test4 = (N - 20) * summary(fit4)$r.squared

    crit = qchisq(1 - alpha, 20)

    pvalue1 = 1 - pchisq(test1, 20)
    pvalue2 = 1 - pchisq(test2, 20)
    pvalue3 = 1 - pchisq(test3, 20)
    pvalue4 = 1 - pchisq(test4, 20)


    out = list(test = c(test1 = test1, test2 = test2, test3 = test3, test4 = test4), crit = crit, pvalue = c(pvalue1,
        pvalue2, pvalue3, pvalue4))
    return(out)
}

PIT_test <- function(U, G = 20, alpha = 0.05, plot = FALSE) {

    iG = G
    dAlpha = alpha
    vU = U

    if (length(vU) < 100)
        iG = 5

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

        mZ[, st] = vHit[st:(cT - (cLags + 1 - st))]

    }

    mX = cbind(vConstant, vVaRforecast, mZ)

    dDQstatOut = (t(vHIT) %*% mX %*% ginv(t(mX) %*% mX) %*% t(mX) %*% (vHIT))/(tau * (1 - tau))

    dDQpvalueOut = 1 - pchisq(dDQstatOut, ncol(mX))

    out = list(stat = dDQstatOut, pvalue = dDQpvalueOut)
}

HitSequence <- function(returns_X, VaR_X) {
    N = length(returns_X)
    Hit_X = numeric(N)
    Hit_X[which(returns_X <= VaR_X)] = 1
    return(Hit_X)
}
Kupiec <- function(Hit, tau, alphaTest = 0.95) {
    N = length(Hit)
    x = sum(Hit)
    rate = x/N
    test = -2 * log(((1 - tau)^(N - x) * tau^x)/((1 - rate)^(N - x) * rate^x))
    if (is.nan(test))
        test = -2 * ((N - x) * log(1 - tau) + x * log(tau) - (N - x) * log(1 - rate) - x * log(rate))
    threshold = qchisq(alphaTest, df = 1)
    pvalue = 1 - pchisq(test, df = 1)

    LRpof = c(test, pvalue)
    names(LRpof) = c("Test", "Pvalue")
    return(LRpof)
}

Christoffersen <- function(Hit, tau, alphaTest = 0.95) {
    n00 = n01 = n10 = n11 = 0
    N = length(Hit)
    for (i in 2:N) {
        if (Hit[i] == 0 & Hit[i - 1] == 0)
            n00 = n00 + 1
        if (Hit[i] == 0 & Hit[i - 1] == 1)
            n01 = n01 + 1
        if (Hit[i] == 1 & Hit[i - 1] == 0)
            n10 = n10 + 1
        if (Hit[i] == 1 & Hit[i - 1] == 1)
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
    LRpof = Kupiec(Hit, tau, alphaTest = alphaTest)["Test"]
    LRcc = LRpof + LRind
    threshold = qchisq(alphaTest, df = 2)
    pvalue = 1 - pchisq(LRcc, df = 2)
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
    series = series[which(Hit == 1)]
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


