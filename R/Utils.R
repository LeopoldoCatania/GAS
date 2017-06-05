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
    return(0.9999999)
}

LowerNu <- function() return(2.5)
UpperNu <- function() return(50)

Array2Matrix <- function(aArray, type) {

  iN = dim(aArray)[1]
  iT = dim(aArray)[3]

  vSeriesName = dimnames(aArray)[[1]]
  vDates = dimnames(aArray)[[3]]

  if (is.null(vSeriesName)) {
    vSeriesName = paste("series", 1:iN, sep = "")
  }

  if (is.null(vDates)) {
    vDates = 1:iT
  }

  vColNames = NULL

  # keep all the lements
  if (type == 1) {
    iL = iN * iN
    mMat = matrix(data = NA, iT, iL)
    iC = 1
    for (i in 1:iN) {
      for (j in 1:iN) {
        mMat[, iC] = aArray[i, j, ]
        iC = iC + 1
        vColNames = c(vColNames, paste(vSeriesName[i], vSeriesName[j], sep = "."))
      }
    }
  }
  # keep the lower triangular elements
  if (type == 2) {
    iL = iN * (iN + 1) / 2
    mMat = matrix(data = NA, iT, iL)
    iC = 1
    for (i in 1:iN) {
      for (j in 1:i) {
        mMat[, iC] = aArray[i, j, ]
        iC = iC + 1
        vColNames = c(vColNames, paste(vSeriesName[i], vSeriesName[j], sep = "."))
      }
    }
  }
  # keep the lower triangular elements
  # without the main diagonal elements
  if (type == 3) {
    iL = iN * (iN + 1) / 2
    mMat = matrix(data = NA, iT, iL)
    iC = 1
    for (i in 1:iN) {
      for (j in 1:i) {
        if (i != j) {
          mMat[, iC] = aArray[i, j, ]
          iC = iC + 1
          vColNames = c(vColNames, paste(vSeriesName[i], vSeriesName[j], sep = "."))
        }
      }
    }
  }
  colnames(mMat) = vColNames
  rownames(mMat) = vDates

  return(mMat)
}
