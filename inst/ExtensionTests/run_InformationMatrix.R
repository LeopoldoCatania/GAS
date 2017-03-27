require("cubature")
require("GAS")
require("testthat")

tol = 1e-4

# norm

dMu     = 1
dSigma2 = 4

vTheta = c(dMu, dSigma2)

mApprox = round(matrix(cubature::adaptIntegrate(function(dY, vTheta) {
  vScore = GAS::Score_Uni(dY, vTheta, "norm")
  vScore %*% t(vScore) * GAS::ddist_Uni(dY, vTheta, "norm", log = FALSE)
}, lowerLimit = -10, upperLimit = 10, vTheta = vTheta, fDim = length(vTheta)^2)$integral,
length(vTheta)), 7)

mGAS = GAS:::IM_Uni(vTheta, "norm")

expect_equal(dim(mApprox), c(2, 2))
expect_equal(dim(mGAS), c(2, 2))
expect_true(max(mApprox - mGAS) < tol)

# std

dMu   = 0
dPhi2 = 1
dNu   = 7

vTheta = c(dMu, dPhi2, dNu)

mApprox = round(matrix(cubature::adaptIntegrate(function(dY, vTheta) {
  vScore = GAS::Score_Uni(dY, vTheta, "std")
  vScore %*% t(vScore) * GAS::ddist_Uni(dY, vTheta, "std", log = FALSE)
}, lowerLimit = -15, upperLimit = 15, vTheta = vTheta, fDim = length(vTheta)^2, maxEval = 2e4)$integral,
length(vTheta)), 7)

mGAS = GAS:::IM_Uni(vTheta, "std")

expect_equal(dim(mApprox), c(3, 3))
expect_equal(dim(mGAS), c(3, 3))
expect_true(max(mApprox - mGAS) < tol)

# ast

dMu    = 1
dSigma = 4
dNu1   = 5
dNu2   = 6
dAlpha = 0.4

vTheta = c(dMu, dSigma, dAlpha, dNu1, dNu2)

mApprox = round(matrix(cubature::adaptIntegrate(function(dY, vTheta) {
  vScore = GAS::Score_Uni(dY, vTheta, "ast")
  vScore %*% t(vScore) * GAS::ddist_Uni(dY, vTheta, "ast", log = FALSE)
}, lowerLimit = -50, upperLimit = 50, vTheta = vTheta, fDim = length(vTheta)^2)$integral,
length(vTheta)), 4)

mGAS = GAS:::IM_Uni(vTheta, "ast")

expect_equal(dim(mApprox), c(5, 5))
expect_equal(dim(mGAS), c(5, 5))
expect_true(max(mApprox - mGAS) < tol)

# ast1

dMu    = 1
dSigma = 4
dNu1   = 5
dAlpha = 0.4

vTheta = c(dMu, dSigma, dAlpha, dNu1)

mApprox = round(matrix(cubature::adaptIntegrate(function(dY, vTheta) {
  vScore = GAS::Score_Uni(dY, vTheta, "ast1")
  vScore %*% t(vScore) * GAS::ddist_Uni(dY, vTheta, "ast1", log = FALSE)
}, lowerLimit = -50, upperLimit = 50, vTheta = vTheta, fDim = length(vTheta)^2)$integral,
length(vTheta)), 4)

mGAS = GAS:::IM_Uni(vTheta, "ast1")

expect_equal(dim(mApprox), c(4, 4))
expect_equal(dim(mGAS), c(4, 4))
expect_true(max(mApprox - mGAS) < tol)

