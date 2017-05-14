StaticLLKoptimizer_Uni <- function(vTheta_tilde, vY, Dist, iT, iK) {
    vTheta = MapParameters_univ(vTheta_tilde, Dist, iK)
    dLLK = StaticLLK_Univ(vY, vTheta, iT, Dist)

    if (is.na(dLLK)) {
        dLLK = -1e+50
    }
    return(-dLLK)
}

# Note that in wrapper_StaticLLKoptimizer_Uni() the argument GASSpec
# is used as an artificial variable which contains additional parameters
# needed by StaticLLKoptimizer_Uni(). This is introduced in order to use
# the same optimzer for the choice of starting values and model estimation.
wrapper_StaticLLKoptimizer_Uni <- function(vTheta_tilde, data, GASSpec) {

  Dist = GASSpec$Dist
  iT   = GASSpec$iT
  iK   = GASSpec$iK

  dmLLK = StaticLLKoptimizer_Uni(vTheta_tilde = vTheta_tilde,
                                 vY   = data,
                                 Dist = Dist,
                                 iT   = iT,
                                 iK   = iK)

  return(dmLLK)

}


StaticLLKoptimizer_Multi <- function(vTheta_tilde, mY, Dist, iT, iN, iK) {
    vTheta = MapParameters_multi(vTheta_tilde, Dist, iN, iK)
    dLLK = StaticLLK_Multi(mY, vTheta, iT, iN, Dist)

    if (is.na(dLLK)) {
        dLLK = -1e+50
    }
    return(-dLLK)
}

# Note that in wrapper_StaticLLKoptimizer_Multi() the argument GASSpec
# is used as an artificial variable which contains additional parameters
# needed by StaticLLKoptimizer_Uni(). This is introduced in order to use
# the same optimzer for the choice of starting values and model estimation.
wrapper_StaticLLKoptimizer_Multi <- function(vTheta_tilde, data, GASSpec) {

  Dist = GASSpec$Dist
  iT   = GASSpec$iT
  iK   = GASSpec$iK
  iN   = GASSpec$iN

  dmLLK = StaticLLKoptimizer_Multi(vTheta_tilde = vTheta_tilde,
                                 mY   = data,
                                 Dist = Dist,
                                 iT   = iT,
                                 iK   = iK,
                                 iN   = iN)

  return(dmLLK)

}

UniGASOptimiser <- function(vPw, data, GASSpec) {

  if (is.null(names(vPw))) {
    names(vPw) = GAS:::getPwNames(GASSpec)
  }

  Dist = getDist(GASSpec)
  ScalingType = getScalingType(GASSpec)
  iK = NumberParameters(Dist)
  iT = length(data)

  lParList = vPw2lPn_Uni(vPw, iK)

  lParList = AddFixedPar(lParList)

  dLLK = try(GASFilter_univ(data, lParList$vKappa, lParList$mA, lParList$mB, iT, iK, Dist, ScalingType)$dLLK,
             silent = TRUE)

  if (!is(dLLK, "try-error")) {
    dMLLK = -dLLK
  } else {
    dMLLK = 1e+50
  }

  if (!is.finite(dMLLK)) {
    dMLLK = 1e+50
  }

  return(dMLLK)

}

MultiGASOptimiser <- function(vPw, data, GASSpec){

  if (is.null(names(vPw))) {
    names(vPw) = GAS:::getPwNames(GASSpec)
  }

  Dist = getDist(GASSpec)
  ScalingType = getScalingType(GASSpec)
  iT = ncol(data)
  iN = nrow(data)
  iK = NumberParameters(Dist, iN)
  ScalarParameters = GASSpec@Spec$ScalarParameters

  lParList = vPw2lPn_Multi(vPw, Dist, iK, iN, ScalarParameters)
  lParList = AddFixedPar(lParList)

  dLLK = try(GASFilter_multi(data, lParList$vKappa, lParList$mA, lParList$mB, iT, iN, iK, Dist, ScalingType)$dLLK,
             silent = TRUE)

  if (!is(dLLK, "try-error")) {
    dMLLK = -dLLK
  } else {
    dMLLK = 1e+50
  }
  if (!is.finite(dMLLK)) {
    dMLLK = 1e+50
  }
  return(dMLLK)
}
