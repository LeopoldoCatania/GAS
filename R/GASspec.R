
UniGASSpec <- function(Dist = "norm", ScalingType = "Identity", GASPar = list(location = FALSE, scale = TRUE,
    skewness = FALSE, shape = FALSE, shape2 = FALSE)) {

    iK = NumberParameters(Dist)

    if (is.null(GASPar$location)) {
        GASPar$location = FALSE
    }
    if (is.null(GASPar$scale)) {
        GASPar$scale = FALSE
    }
    if (is.null(GASPar$skewness)) {
        GASPar$skewness = FALSE
    }
    if (is.null(GASPar$shape)) {
        GASPar$shape = FALSE
    }
    if (is.null(GASPar$shape2)) {
        GASPar$shape2 = FALSE
    }

    DistPar = DistParameters(Dist)

    GASPar = GASPar[DistPar]

    Spec <- new("uGASSpec", Spec = list(Dist = Dist, ScalingType = ScalingType, iK = iK, GASPar = GASPar))
    return(Spec)
}

GAS_validate_multi_scaling <- function(ScalingType) {
    if (!is.character(ScalingType) || length(ScalingType) != 1L ||
        is.na(ScalingType) || !ScalingType %in% c("Identity", "Inv", "InvSqrt"))
        stop('ScalingType must be "Identity", "Inv", or "InvSqrt".',
             call. = FALSE)
    invisible(ScalingType)
}

MultiGASSpec <- function(Dist = "mvnorm", ScalingType = "Identity", GASPar = list(location = FALSE,
    scale = TRUE, correlation = FALSE, shape = FALSE), ScalarParameters = TRUE) {

    GAS_validate_multi_scaling(ScalingType)

    if (is.null(GASPar$location)) {
        GASPar$location = FALSE
    }
    if (is.null(GASPar$scale)) {
        GASPar$scale = FALSE
    }
    if (is.null(GASPar$correlation)) {
        GASPar$correlation = FALSE
    }
    if (is.null(GASPar$shape)) {
        GASPar$shape = FALSE
    }

    DistPar = DistParameters(Dist)
    GASPar = GASPar[DistPar]

    Spec <- new("mGASSpec", Spec = list(Dist = Dist, ScalingType = ScalingType, GASPar = GASPar, ScalarParameters = ScalarParameters))
    return(Spec)
}


