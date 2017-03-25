StaticLLKoptimizer_Uni <- function(vTheta_tilde, vY, Dist, iT, iK) {
    vTheta = MapParameters_univ(vTheta_tilde, Dist, iK)
    dLLK = StaticLLK_Univ(vY, vTheta, iT, Dist)

    if (is.na(dLLK)) {
        dLLK = -1e+50
    }
    return(-dLLK)
}
StaticLLKoptimizer_Multi <- function(vTheta_tilde, mY, Dist, iT, iN, iK) {
    vTheta = MapParameters_multi(vTheta_tilde, Dist, iN, iK)
    dLLK = StaticLLK_Multi(mY, vTheta, iT, iN, Dist)

    if (is.na(dLLK)) {
        dLLK = -1e+50
    }
    return(-dLLK)
}

UniGASOptimiser <- function(vPw, data, GASSpec) {

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
