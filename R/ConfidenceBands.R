ConfidenceBands <- function(object, B = 10000L, probs = c(0.01, 0.1, 0.9, 0.99), ...) {

    if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1 ||
        B != floor(B) || B > .Machine$integer.max)
      stop("B must be a positive integer.")
    if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
        any(probs < 0 | probs > 1)) stop("probs must contain finite probabilities in [0,1].")
    if (!isTRUE(convergence(object) == 0))
      stop("Confidence bands require a fit that passed convergence checks.")
    iB = B

    vPw = object@Estimates$optimiser$pars
    iT = object@ModelInfo$iT
    iK = object@ModelInfo$iK
    Data = getObs(object)
    Dist = getDist(object)
    DType = DistType(Dist)
    ScalingType = getScalingType(object)

    objective <- if (DType == "univariate") UniGASOptimiser else MultiGASOptimiser
    gradient <- GAS_objective_gradient(objective, vPw, Data, object@ModelInfo$Spec)
    if (any(!is.finite(gradient)) || max(abs(gradient)) > 1e-5)
      stop("Confidence bands require an adequately optimized likelihood.")
    vCov <- GAS_checked_covariance(object@Estimates$optimiser$opg, length(vPw))
    if (is.null(vCov)) stop("Confidence bands require a finite, nonsingular OPG information matrix. Refit with Compute.SE = TRUE.")
    # Use the checked covariance directly; do not repair it with chol_safe().
    mPw <- sweep(matrix(stats::rnorm(iB*length(vPw)), nrow = iB) %*%
                   chol(vCov), 2L, vPw, "+")

    mTheta = getFilteredParameters(object)

    cTheta = array(0, dim = c(iT + 1L, iB + 1L, iK), dimnames = list(1:(iT + 1L), 1:(iB + 1L), colnames(mTheta)))
    cTheta[, 1L, ] = mTheta

    for (b in 2:(iB + 1L)) {
        vPw_foo = mPw[b - 1L, ]
        names(vPw_foo) = names(vPw)

        if (DType == "univariate") {

          lParList <- vPw2lPn_Uni(vPw_foo, iK)

        } else if (DType == "multivariate") {

          lParList <- vPw2lPn_Multi(
            vPw = vPw_foo,
            Dist = Dist,
            iK = iK,
            iN = object@ModelInfo$iN,
            ScalarParameters =
              object@ModelInfo$Spec@Spec$ScalarParameters
          )

        } else {
          stop("Unsupported distribution type: ", DType)
        }

        lParList <- AddFixedPar(lParList)

        if (DType == "univariate")
            cTheta[, b, ] = t(GASFilter_univ(Data, lParList$vKappa, lParList$mA, lParList$mB, iT, iK,
                Dist, ScalingType)$mTheta)
        if (DType == "multivariate")
            cTheta[, b, ] = t(GASFilter_multi(Data, lParList$vKappa, lParList$mA, lParList$mB, iT, nrow(Data),
                iK, Dist, ScalingType)$mTheta)

    }

    cQuantile = array(0, dim = c(iT + 1L, length(probs), iK), dimnames = list(1:(iT + 1L), paste("q",
        probs, sep = ""), colnames(mTheta)))

    for (k in seq_len(iK)) {
      simulated_paths <- matrix(
        cTheta[, -1L, k, drop = FALSE],
        nrow = iT + 1L,
        ncol = iB
      )

      cQuantile[, , k] <- t(
        apply(simulated_paths, 1L, quantile, probs = probs)
      )
    }

    return(cQuantile)
}




