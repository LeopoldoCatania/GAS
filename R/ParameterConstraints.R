
GetFixedPar_Uni <- function(Dist, GASPar) {
    FixedPar = NULL
    if (Dist == "norm") {
        if (!GASPar$location)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
        if (!GASPar$scale)
            FixedPar = c(FixedPar, a2 = 0, b2 = 0)
    }
    if (Dist == "snorm") {
        if (!GASPar$location)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
        if (!GASPar$scale)
            FixedPar = c(FixedPar, a2 = 0, b2 = 0)
        if (!GASPar$skewness)
            FixedPar = c(FixedPar, a3 = 0, b3 = 0)
    }
    if (Dist == "std") {
        if (!GASPar$location)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
        if (!GASPar$scale)
            FixedPar = c(FixedPar, a2 = 0, b2 = 0)
        if (!GASPar$shape)
            FixedPar = c(FixedPar, a3 = 0, b3 = 0)
    }
    if (Dist == "sstd") {
        if (!GASPar$location)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
        if (!GASPar$scale)
            FixedPar = c(FixedPar, a2 = 0, b2 = 0)
        if (!GASPar$skewness)
            FixedPar = c(FixedPar, a3 = 0, b3 = 0)
        if (!GASPar$shape)
            FixedPar = c(FixedPar, a4 = 0, b4 = 0)
    }
    if (Dist == "ast") {
        if (!GASPar$location)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
        if (!GASPar$scale)
            FixedPar = c(FixedPar, a2 = 0, b2 = 0)
        if (!GASPar$skewness)
            FixedPar = c(FixedPar, a3 = 0, b3 = 0)
        if (!GASPar$shape)
            FixedPar = c(FixedPar, a4 = 0, b4 = 0)
        if (!GASPar$shape2)
            FixedPar = c(FixedPar, a5 = 0, b5 = 0)
    }
    if (Dist == "ast1") {
        if (!GASPar$location)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
        if (!GASPar$scale)
            FixedPar = c(FixedPar, a2 = 0, b2 = 0)
        if (!GASPar$skewness)
            FixedPar = c(FixedPar, a3 = 0, b3 = 0)
        if (!GASPar$shape)
            FixedPar = c(FixedPar, a4 = 0, b4 = 0)
    }
    if (Dist == "ald") {
        if (!GASPar$location)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
        if (!GASPar$scale)
            FixedPar = c(FixedPar, a2 = 0, b2 = 0)
        if (!GASPar$skewness)
            FixedPar = c(FixedPar, a3 = 0, b3 = 0)
    }
    if (Dist == "ghskt") {
      if (!GASPar$location)
        FixedPar = c(FixedPar, a1 = 0, b1 = 0)
      if (!GASPar$scale)
        FixedPar = c(FixedPar, a2 = 0, b2 = 0)
      if (!GASPar$skewness)
        FixedPar = c(FixedPar, a3 = 0, b3 = 0)
      if (!GASPar$shape)
        FixedPar = c(FixedPar, a4 = 0, b4 = 0)
    }
    if (Dist == "poi") {
        if (!GASPar$location)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
    }
    if (Dist == "ber") {
        if (!GASPar$location)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
    }
    if (Dist == "gamma") {
        if (!GASPar$scale)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
        if (!GASPar$shape)
            FixedPar = c(FixedPar, a2 = 0, b2 = 0)
    }
    if (Dist == "exp") {
        if (!GASPar$location)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
    }
    if (Dist == "beta") {
        if (!GASPar$scale)
            FixedPar = c(FixedPar, a1 = 0, b1 = 0)
        if (!GASPar$shape)
            FixedPar = c(FixedPar, a2 = 0, b2 = 0)

    }
    if (Dist == "negbin") {
      if (!GASPar$location)
        FixedPar = c(FixedPar, a1 = 0, b1 = 0)
      if (!GASPar$scale)
        FixedPar = c(FixedPar, a2 = 0, b2 = 0)
    }
    if (Dist == "skellam") {
      if (!GASPar$location)
        FixedPar = c(FixedPar, a1 = 0, b1 = 0)
      if (!GASPar$scale)
        FixedPar = c(FixedPar, a2 = 0, b2 = 0)
    }
    return(FixedPar)
}

MultiFixedScale <- function(iN, Dist, ScalarParameters) {

    if (ScalarParameters) {

        FixedPar = c(0, 0)
        if (Dist == "mvnorm")
            names(FixedPar) = c("a.sigma", "b.sigma")
        if (Dist == "mvt")
            names(FixedPar) = c("a.phi", "b.phi")

    } else {

        FixedPar = rep(0, iN * 2L)
        if (Dist == "mvnorm")
            names(FixedPar) = c(paste("a.sigma", 1:iN, sep = ""), paste("b.sigma", 1:iN, sep = ""))
        if (Dist == "mvt")
            names(FixedPar) = c(paste("a.phi", 1:iN, sep = ""), paste("b.phi", 1:iN, sep = ""))

    }
    return(FixedPar)
}

MultiFixedLocation <- function(iN, ScalarParameters) {
    if (ScalarParameters) {

        FixedPar = c(0, 0)
        names(FixedPar) = c("a.mu", "b.mu")

    } else {

        FixedPar = rep(0, iN * 2L)
        names(FixedPar) = c(paste("a.mu", 1:iN, sep = ""), paste("b.mu", 1:iN, sep = ""))

    }

    return(FixedPar)
}
MultiFixedCorrelation <- function(iN, ScalarParameters) {
    if (ScalarParameters) {

        FixedPar = c(0, 0)
        names(FixedPar) = c("a.rho", "b.rho")

    } else {

        FixedPar = rep(0, iN * (iN - 1L))
        vRhoNames = RhoNames(iN)
        names(FixedPar) = c(paste("a.", vRhoNames, sep = ""), paste("b.", vRhoNames, sep = ""))

    }

    return(FixedPar)
}

GetFixedPar_Multi <- function(Dist, GASPar, iN, ScalarParameters) {
    FixedPar = NULL

    if (!GASPar$location) {
        FixedPar = c(FixedPar, MultiFixedLocation(iN, ScalarParameters))
    }
    if (!GASPar$scale) {
        FixedPar = c(FixedPar, MultiFixedScale(iN, Dist, ScalarParameters))
    }
    if (!GASPar$correlation) {
        FixedPar = c(FixedPar, MultiFixedCorrelation(iN, ScalarParameters))
    }
    if (Dist == "mvt") {
      if (!GASPar$shape) {
        FixedPar = c(FixedPar, a.nu = 0, b.nu = 0)
      }
    }

    return(FixedPar)
}
RemoveFixedPar <- function(vPw, FixedPar) {
    if (!is.null(FixedPar))
        vPw = vPw[-which(names(vPw) %in% names(FixedPar))]
    return(vPw)
}
AddFixedPar <- function(lParList) {
    for (i in 1:length(lParList)) {
        lParList[[i]][is.na(lParList[[i]])] = 0
    }
    return(lParList)
}


FixedDynamicPar_Multi <- function(Dist, iN, GASPar) {

    ParNames = FullNamesMulti(iN, Dist)
    vBool = rep(FALSE, length(ParNames))
    names(vBool) = ParNames

    if (GASPar$location) {
        vBool[1:iN] = TRUE
    }
    if (GASPar$scale) {
        vBool[(iN + 1L):(2L * iN)] = TRUE
    }
    if (GASPar$correlation) {
        vBool[(2L * iN + 1L):(2L * iN + iN * (iN - 1L)/2L)] = TRUE
    }
    if (Dist == "mvt") {
      if (GASPar$shape) {
        vBool[2L * iN + iN * (iN - 1L)/2L + 1L] = TRUE
      }
    }

    return(vBool)
}


MatrixCoefficientStructure_Multi <- function(Dist, iN, iK, GASPar) {

    lStructure = list()

    for (i in 1:length(GASPar)) lStructure[[names(GASPar)[i]]] = matrix(FALSE, iK, iK)

    if (GASPar$location) {
        diag(lStructure[["location"]])[1:iN] = TRUE
    }
    if (GASPar$scale) {
        diag(lStructure[["scale"]])[(iN + 1L):(2L * iN)] = TRUE
    }
    if (GASPar$correlation) {
        diag(lStructure[["correlation"]])[(2L * iN + 1L):(2L * iN + iN * (iN - 1L)/2L)] = TRUE
    }
    if (Dist == "mvt") {
      if (GASPar$shape) {
        diag(lStructure[["shape"]])[(1L + 2L * iN + iN * (iN - 1L)/2L)] = TRUE
      }
    }

    return(lStructure)

}

