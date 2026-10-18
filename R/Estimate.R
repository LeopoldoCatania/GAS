
StaticMLFIT <- function(vY, Dist, fn.optimizer) {

  iT = length(vY)
  iK = NumberParameters(Dist)

  vTheta_tilde = StaticStarting_Uni(vY, Dist, iK)

  lArguments = list(Dist    = Dist,
                    iT      = iT,
                    iK      = iK)

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

  optimiser <- GAS_finalize_optimiser(optimiser, vPw, vY, GASSpec,
                                      UniGASOptimiser, Compute.SE)

  vPw = optimiser$pars

  if (is.null(names(vPw))) {
    names(vPw) = getPwNames(GASSpec)
  }

  lParList = vPw2lPn_Uni(vPw, iK)
  lParList = AddFixedPar(lParList)

  mHessian = optimiser$hessian
  valid_inference <- Compute.SE && isTRUE(optimiser$diagnostics$inference_ok)

  Inference = InferenceFun_Uni(mHessian, vPw, iK, valid_inference)

  GASDyn = GASFilter_univ(vY, lParList$vKappa, lParList$mA, lParList$mB, iT, iK, Dist, ScalingType)

  IC = ICfun(-optimiser$value, length(optimiser$pars), iT)

  vU = EvaluatePit_Univ(GASDyn$mTheta, vY, Dist, iT)
  PitTest = PIT_test(vU, G = 20L, alpha = 0.05, plot = FALSE)

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

  analytic_gaussian <- Dist == "mvnorm" && !any(unlist(GASPar))
  if (analytic_gaussian) {
    mu <- rowMeans(mY)
    centered <- sweep(mY, 1L, mu)
    S <- tcrossprod(centered)/iT
    if (any(!is.finite(S)) || is.null(tryCatch(chol(S), error = function(e) NULL)))
      stop("A static multivariate Gaussian MLE requires a full-rank sample covariance.")
    Theta <- c(mu, sqrt(diag(S)), as.numeric(build_vR(stats::cov2cor(S), iN)))
    vPw <- as.numeric(MultiUnmapParameters(Theta, Dist, iN))
    names(vPw) <- paste0("kappa.", FullNamesMulti(iN, Dist))
  } else {
    vPw = MultiGAS_Starting(mY, iT, iN, iK, Dist, GASPar, ScalingType, ScalarParameters, fn.optimizer)
    FixedPar = GetFixedPar_Multi(Dist, GASPar, iN, ScalarParameters)
    vPw = RemoveFixedPar(vPw, FixedPar)
  }

  GASSpec@Spec$PwNames = names(vPw)

  # optimise
  if (analytic_gaussian) {
    optimiser <- list(pars = vPw,
      value = MultiGASOptimiser(vPw, data = mY, GASSpec = GASSpec),
      convergence = 0L, hessian = NULL, method = "analytical Gaussian MLE")
  } else {
    optimiser = fn.optimizer(par0 = vPw, data = mY, GASSpec = GASSpec, FUN = MultiGASOptimiser)
  }

  optimiser <- GAS_finalize_optimiser(optimiser, vPw, mY, GASSpec,
                                      MultiGASOptimiser, Compute.SE)

  vPw = optimiser$pars

  if (is.null(names(vPw))) {
    names(vPw) = getPwNames(GASSpec)
  }

  lParList = vPw2lPn_Multi(vPw, Dist, iK, iN, ScalarParameters)
  lParList = AddFixedPar(lParList)

  mHessian = optimiser$hessian
  valid_inference <- Compute.SE && isTRUE(optimiser$diagnostics$inference_ok)

  Inference = InferenceFun_Multi(mHessian, Dist, vPw, iK, iN, ScalarParameters, valid_inference)

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
