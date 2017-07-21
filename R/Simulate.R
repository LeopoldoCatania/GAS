
UniGASSim <- function(fit = NULL, T.sim = 1000L, kappa = NULL, A = NULL, B = NULL, Dist = NULL, ScalingType = NULL) {

  if (!is.null(fit)) {
    lCoef = coef(fit, do.list = TRUE)
    kappa = lCoef$lCoef$vKappa
    A = lCoef$lCoef$mA
    B = lCoef$lCoef$mB
    Dist = getDist(fit)
    ScalingType = getScalingType(fit)
  } else {
    if (any(is.null(kappa), is.null(A), is.null(B), is.null(Dist), is.null(ScalingType))) {
      stop("If an uGASFit object is not provided, arguments kappa, A, B,
           Dist and ScalingType, have to be provided")
    }
  }

    iT = T.sim

    vKappa = kappa
    mA = A
    mB = B

    lSim = SimulateGAS_univ(iT, vKappa, mA, mB, Dist, ScalingType)

    iK = NumberParameters(Dist)

    mMoments = EvalMoments_univ(lSim$mTheta, Dist)

    Sim <- new("uGASSim",
               ModelInfo = list(iT = iT,
                                iK = iK,
                                vKappa = vKappa,
                                mA = mA,
                                mB = mB,
                                Dist = Dist,
                                ScalingType = ScalingType),
               GASDyn = lSim,
               Data = list(vY = lSim$vY,
                           Moments = mMoments))

    return(Sim)
}

MultiGASSim <- function(fit = NULL, T.sim = 1000L, N = NULL, kappa = NULL, A = NULL, B = NULL, Dist = NULL, ScalingType = NULL) {

  if (!is.null(fit)) {
    lCoef = coef(fit, do.list = TRUE)
    kappa = lCoef$lCoef$vKappa
    A = lCoef$lCoef$mA
    B = lCoef$lCoef$mB
    Dist = getDist(fit)
    ScalingType = getScalingType(fit)
    N = fit@ModelInfo$iN
  } else {
    if (any(is.null(kappa), is.null(N), is.null(A), is.null(B), is.null(Dist), is.null(ScalingType))) {
      stop("If an uGASFit object is not provided, arguments kappa, A, B,
           Dist and ScalingType, have to be provided")
    }
  }

    iT = T.sim

    vKappa = kappa
    mA = A
    mB = B
    iN = N

    lSim = SimulateGAS_multi(iT, iN, vKappa, mA, mB, Dist, ScalingType)

    iK = NumberParameters(Dist, iN)
    mY = lSim$mY
    rownames(mY) = paste("Series", 1:iN)

    lMoments = EvalMoments_multi(lSim$mTheta, Dist, iN)

    Sim <- new("mGASSim",
               ModelInfo = list(iT = iT,
                                iN = iN,
                                iK = iK,
                                vKappa = vKappa,
                                mA = mA,
                                mB = mB,
                                Dist = Dist,
                                ScalingType = ScalingType),
               GASDyn = lSim,
               Data = list(mY = mY,
                           Moments = lMoments))

    return(Sim)
}




