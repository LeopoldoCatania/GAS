
StaticMLFIT <- function(vY, Dist, fn.optimizer) {

  iT = length(vY)
  iK = NumberParameters(Dist)

  vTheta_tilde = StaticStarting_Uni(vY, Dist, iK)

  lArguments = list(Dist    = Dist,
                    iT      = iT,
                    iK      = iK)

  # optimiser = fn.optimizer(par0 = vTheta_tilde, data = vY, GASSpec = lArguments, FUN = wrapper_StaticLLKoptimizer_Uni)

  optimiser = StaticOptimizationLink_Univ(vTheta_tilde, vY, lArguments, fn.optimizer)

  vTheta_tilde = optimiser$pars

  vTheta = as.numeric(MapParameters_univ(vTheta_tilde, Dist, iK))
  names(vTheta) = FullNamesUni(Dist)

  out = list(vTheta = vTheta, dLLK = -optimiser$value, optimiser = optimiser)

  return(out)
}

StaticMLFIT_Multiv <- function(mY, Dist, fn.optimizer) {

  iT = ncol(mY)
  iN = nrow(mY)
  iK = NumberParameters(Dist, iN)

  vTheta_tilde = StaticStarting_Multi(mY, Dist, iN)

  lArguments = list(Dist    = Dist,
                    iT      = iT,
                    iK      = iK,
                    iN      = iN)

  optimiser = fn.optimizer(par0 = vTheta_tilde, data = mY, GASSpec = lArguments, FUN = wrapper_StaticLLKoptimizer_Multi)

  vTheta_tilde = optimiser$pars

  vTheta = as.numeric(MapParameters_multi(vTheta_tilde, Dist, iN, iK))
  names(vTheta) = FullNamesMulti(iN, Dist)

  out = list(vTheta = vTheta, dLLK = -optimiser$value, optimiser = optimiser)

  return(out)
}

UniGASFit <- function(GASSpec, data, fn.optimizer = fn.optim, Compute.SE = TRUE) {

  vY = data

  Start = Sys.time()

  iT = length(vY)

  Dist = getDist(GASSpec)
  ScalingType = getScalingType(GASSpec)
  GASPar = getGASPar(GASSpec)
  iK = NumberParameters(Dist)

  # starting par
  lStarting = UniGAS_Starting(vY, iT, iK, Dist, ScalingType, GASPar, fn.optimizer)
  vPw = lStarting$vPw
  StaticFit = lStarting$StaticFit

  # fixed par
  FixedPar = GetFixedPar_Uni(Dist, GASPar)
  vPw = RemoveFixedPar(vPw, FixedPar)

  GASSpec@Spec$PwNames = names(vPw)

  # optimise
  optimiser = fn.optimizer(par0 = vPw, data = vY, GASSpec = GASSpec, FUN = UniGASOptimiser)

  vPw = optimiser$pars

  if (is.null(names(vPw))) {
    names(vPw) = getPwNames(GASSpec)
  }

  lParList = vPw2lPn_Uni(vPw, iK)
  lParList = AddFixedPar(lParList)

  mHessian = optimiser$hessian

  if (Compute.SE) {
  if (is.null(mHessian)) {
    mHessian = numDeriv::hessian(UniGASOptimiser, vPw, data = vY, GASSpec = GASSpec)
    if (any(!is.finite(mHessian)) || min(eigen(mHessian)$values) < 0) {
      mHessian = optim(vPw, UniGASOptimiser, data = vY, GASSpec = GASSpec, hessian = TRUE, control = list(maxit = 1))$hessian
    }
    if (min(eigen(mHessian)$values) < 0) {
      warning("Hessian matrix is not positive definite. Standard errors are not reliable")
    }
  }
  }

  Inference = InferenceFun_Uni(mHessian, vPw, iK, Compute.SE)

  GASDyn = GASFilter_univ(vY, lParList$vKappa, lParList$mA, lParList$mB, iT, iK, Dist, ScalingType)

  IC = ICfun(-optimiser$value, length(optimiser$pars), iT)

  vU = EvaluatePit_Univ(GASDyn$mTheta, vY, Dist, iT)
  PitTest = PIT_test(vU, G = 20, alpha = 0.05, plot = FALSE)

  mMoments = EvalMoments_univ(GASDyn$mTheta, Dist)

  elapsedTime = Sys.time() - Start

  Out <- new("uGASFit",
             ModelInfo = list(Spec = GASSpec,
                              iT = iT,
                              iK = iK,
                              elapsedTime = elapsedTime,
                              Date = Start,
                              convergence = optimiser$convergence),
             GASDyn = GASDyn,
             Estimates = list(lParList = lParList,
                              optimiser = optimiser,
                              StaticFit = StaticFit,
                              Inference = Inference,
                              IC = IC,
                              vU = vU,
                              Moments = mMoments),
             Testing = list(PitTest = PitTest),
             Data = list(vY = vY))

  return(Out)
}

# mY is NxT
MultiGASFit <- function(GASSpec, data, fn.optimizer = fn.optim, Compute.SE = TRUE) {

  mY = t(data)

  Start = Sys.time()
  # getInfo
  Dist = getDist(GASSpec)
  ScalingType = getScalingType(GASSpec)
  GASPar = getGASPar(GASSpec)
  ScalarParameters = GASSpec@Spec$ScalarParameters

  # dimension par
  iT = ncol(mY)
  iN = nrow(mY)
  iK = NumberParameters(Dist, iN)

  if (is.null(rownames(mY))) {
    rownames(mY) = paste("Series", 1:iN)
  }

  # starting par
  vPw = MultiGAS_Starting(mY, iT, iN, iK, Dist, GASPar, ScalingType, ScalarParameters, fn.optimizer)

  # fixed par
  FixedPar = GetFixedPar_Multi(Dist, GASPar, iN, ScalarParameters)
  vPw = RemoveFixedPar(vPw, FixedPar)

  GASSpec@Spec$PwNames = names(vPw)

  # optimise
  optimiser = fn.optimizer(par0 = vPw, data = mY, GASSpec = GASSpec, FUN = MultiGASOptimiser)

  vPw = optimiser$pars

  if (is.null(names(vPw))) {
    names(vPw) = getPwNames(GASSpec)
  }

  lParList = vPw2lPn_Multi(vPw, Dist, iK, iN, ScalarParameters)
  lParList = AddFixedPar(lParList)

  mHessian = optimiser$hessian

  if (Compute.SE) {
  if (is.null(mHessian)) {
    mHessian = numDeriv::hessian(MultiGASOptimiser, vPw, data = mY, GASSpec = GASSpec)
    if (any(!is.finite(mHessian)) || min(eigen(mHessian)$values < 0)) {
      mHessian = optim(vPw, MultiGASOptimiser, data = mY, GASSpec = GASSpec, hessian = TRUE, control = list(maxit = 1))$hessian
    }
    if (min(eigen(mHessian)$values) < 0) {
      warning("Hessian matrix is not positive definite. Standard errors are not reliable")
    }
  }
  }

  Inference = InferenceFun_Multi(mHessian, Dist, vPw, iK, iN, ScalarParameters, Compute.SE)

  GASDyn = GASFilter_multi(mY, lParList$vKappa, lParList$mA, lParList$mB, iT, iN, iK, Dist, ScalingType)

  IC = ICfun(-optimiser$value, length(optimiser$pars), iT)

  ## Moments
  mMoments = EvalMoments_multi(GASDyn$mTheta, Dist, iN)

  elapsedTime = Sys.time() - Start

  Out <- new("mGASFit",
             ModelInfo = list(Spec = GASSpec,
                              iT = iT,
                              iN = iN,
                              iK = iK,
                              elapsedTime = elapsedTime,
                              Date = Start,
                              convergence = optimiser$convergence),
             GASDyn = GASDyn,
             Estimates = list(lParList = lParList,
                              optimiser = optimiser,
                              Inference = Inference,
                              IC = IC,
                              Moments = mMoments),
             Data = list(mY = mY))
  return(Out)
}
