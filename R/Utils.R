LowerA <- function() {
    return(0)
}
UpperA <- function() {
    return(10)
}

LowerB <- function() {
    return(0)
}
UpperB <- function() {
    return(0.9999)
}

LowerNu <- function() return(4.0)
UpperNu <- function() return(50)

Array2Matrix <- function(aArray, type) {

  iN = dim(aArray)[1L]
  iT = dim(aArray)[3L]

  vSeriesName = dimnames(aArray)[[1L]]
  vDates = dimnames(aArray)[[3L]]

  if (is.null(vSeriesName)) {
    vSeriesName = paste("series", 1:iN, sep = "")
  }

  if (is.null(vDates)) {
    vDates = 1:iT
  }

  vColNames = NULL

  # keep all the lements
  if (type == 1L) {
    iL = iN * iN
    mMat = matrix(data = NA, iT, iL)
    iC = 1L
    for (i in 1:iN) {
      for (j in 1:iN) {
        mMat[, iC] = aArray[i, j, ]
        iC = iC + 1L
        vColNames = c(vColNames, paste(vSeriesName[i], vSeriesName[j], sep = "."))
      }
    }
  }
  # keep the lower triangular elements
  if (type == 2L) {
    iL = iN * (iN + 1L) / 2L
    mMat = matrix(data = NA, iT, iL)
    iC = 1L
    for (i in 1:iN) {
      for (j in 1:i) {
        mMat[, iC] = aArray[i, j, ]
        iC = iC + 1L
        vColNames = c(vColNames, paste(vSeriesName[i], vSeriesName[j], sep = "."))
      }
    }
  }
  # keep the lower triangular elements
  # without the main diagonal elements
  if (type == 3L) {
    iL = iN * (iN + 1L) / 2L
    mMat = matrix(data = NA, iT, iL)
    iC = 1L
    for (i in 1:iN) {
      for (j in 1:i) {
        if (i != j) {
          mMat[, iC] = aArray[i, j, ]
          iC = iC + 1L
          vColNames = c(vColNames, paste(vSeriesName[i], vSeriesName[j], sep = "."))
        }
      }
    }
  }
  colnames(mMat) = vColNames
  rownames(mMat) = vDates

  return(mMat)
}
