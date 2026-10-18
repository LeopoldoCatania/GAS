
################################ DENSITY #

############ PIT#################


.validatePIT <- function(pit, alpha) {
    if (!is.numeric(pit) || !length(pit) || any(!is.finite(pit)) ||
        any(pit < 0 | pit > 1))
        stop("PIT values must be a nonempty finite numeric vector in [0,1].")
    if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
        stop("alpha must lie strictly between zero and one.")
}

BinTest <- function(pit, g = 20L, alpha = 0.05, plot = FALSE) {
    .validatePIT(pit, alpha)
    if (length(g) != 1L || !is.finite(g) || g < 2 || g != floor(g) ||
        g >= .Machine$integer.max) stop("g must be an integer of at least two.")
    h <- hist(pit, breaks = seq(0, 1, length.out = g + 1L),
              include.lowest = TRUE, right = TRUE, plot = FALSE)
    expected <- length(pit) / g
    test <- sum((h$counts - expected)^2 / expected)
    adequate <- expected >= 5
    crit <- if (adequate) qchisq(1 - alpha, g - 1L) else NA_real_
    pvalue <- if (adequate) pchisq(test, g - 1L, lower.tail = FALSE) else NA_real_
    # Marginal bin counts have binomial, not Poisson, variance.
    confidence <- expected + c(-1, 1) * qnorm(1 - alpha/2) *
        sqrt(length(pit) / g * (1 - 1/g))
    confidence <- pmin(length(pit), pmax(0, confidence))
    if (plot) {
        plot(h, freq = TRUE, col = "blue",
             ylim = c(0, max(1, confidence, h$counts) * 1.2))
        abline(h = confidence, col = "red", lwd = 2L)
    }
    list(test = test, crit = crit, pvalue = pvalue, hist = h,
         confidence = confidence,
         reason = if (adequate) NULL else
             "Chi-squared calibration unavailable: expected bin count is below five.")
}

iidTest <- function(pit, alpha = 0.05) {
    .validatePIT(pit, alpha)
    L <- 20L
    tests <- setNames(rep(NA_real_, 4L), paste0("test", 1:4))
    pvalues <- setNames(rep(NA_real_, 4L), paste0("pvalue", 1:4))
    reasons <- rep(NA_character_, 4L)
    # The lagged design has an intercept plus L regressors.
    if (length(pit) - L <= L + 1L)
        return(list(test = tests, crit = NA_real_, pvalue = pvalues,
                    reason = "Insufficient observations for the 20-lag regressions."))
    z <- as.numeric(pit - mean(pit))
    for (k in 1:4) {
        dat <- embed(z^k, L + 1L)
        y <- dat[, 1L]
        X <- cbind(1, dat[, -1L, drop = FALSE])
        fit <- lm.fit(X, y)
        tss <- sum((y - mean(y))^2)
        if (fit$rank != ncol(X) || !is.finite(tss) || tss <= 0) {
            reasons[k] <- "Rank-deficient lag regression or constant dependent variable."
            next
        }
        r2 <- max(0, min(1, 1 - sum(fit$residuals^2)/tss))
        tests[k] <- nrow(X) * r2
        pvalues[k] <- pchisq(tests[k], L, lower.tail = FALSE)
    }
    list(test = tests, crit = qchisq(1 - alpha, L), pvalue = pvalues,
         reason = reasons)
}

PIT_test <- function(U, G = 20L, alpha = 0.05, plot = FALSE) {
    list(Hist = BinTest(U, g = G, alpha = alpha, plot = plot),
         IID = iidTest(U, alpha = alpha))
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
    vY2_lag = y[cLags:(cT - 1)]^2

    for (st in 1:cLags) {

        mZ[, st] = vHit[st:(cT - (cLags + 1L - st))]

    }

    mX = cbind(vConstant, vVaRforecast, mZ, vY2_lag)

    dDQstatOut = (t(vHIT) %*% mX %*% MASS::ginv(t(mX) %*% mX) %*% t(mX) %*% (vHIT))/(tau * (1 - tau))

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




