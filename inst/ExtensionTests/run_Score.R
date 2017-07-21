require("numDeriv")
require("GAS")
require("testthat")


tol = 1e-04

## univariate

# norm

dY = 0
dMu = 1
dSigma2 = 4

vTheta = c(dMu, dSigma2)

vApprox = matrix(numDeriv::grad(function(vTheta, dY) {
  GAS::ddist_Uni(dY, vTheta, "norm", log = TRUE)
}, vTheta, dY = dY), ncol = 1)

vGAS = GAS:::Score_Uni(dY, vTheta, "norm")

expect_equal(dim(vApprox), c(2, 1))
expect_equal(dim(vGAS), c(2, 1))
expect_true(max(vApprox - vGAS) < tol)

# std

dY = 0
dMu = 1
dPhi2 = 5
dNu = 5

vTheta = c(dMu, dPhi2, dNu)

vApprox = matrix(numDeriv::grad(function(vTheta, dY) {
  GAS::ddist_Uni(dY, vTheta, "std", log = TRUE)
}, vTheta, dY = dY), ncol = 1)

vGAS = GAS:::Score_Uni(dY, vTheta, "std")

expect_equal(dim(vApprox), c(3, 1))
expect_equal(dim(vGAS), c(3, 1))
expect_true(max(vApprox - vGAS) < tol)

# ast

dY = 0
dMu = 1
dSigma = 4
dNu1 = 5
dNu2 = 6
dAlpha = 0.4

vTheta = c(dMu, dSigma, dAlpha, dNu1, dNu2)

vApprox = matrix(numDeriv::grad(function(vTheta, dY) {
  GAS::ddist_Uni(dY, vTheta, "ast", log = TRUE)
}, vTheta, dY = dY), ncol = 1)

vGAS = GAS:::Score_Uni(dY, vTheta, "ast")

expect_equal(dim(vApprox), c(5, 1))
expect_equal(dim(vGAS), c(5, 1))
expect_true(max(vApprox - vGAS) < tol)

# ast1

dY = 0
dMu = 1
dSigma = 4
dNu = 5
dAlpha = 0.4

vTheta = c(dMu, dSigma, dAlpha, dNu)

vApprox = matrix(numDeriv::grad(function(vTheta, dY) {
  GAS::ddist_Uni(dY, vTheta, "ast1", log = TRUE)
}, vTheta, dY = dY), ncol = 1)

vGAS = GAS:::Score_Uni(dY, vTheta, "ast1")

expect_equal(dim(vApprox), c(4, 1))
expect_equal(dim(vGAS), c(4, 1))
expect_true(max(vApprox - vGAS) < tol)

# snorm

dY = 0.9
dMu = 1
dSigma = 4
dXi = 1

vTheta = c(dMu, dSigma, dXi)

vApprox = matrix(numDeriv::grad(function(vTheta, dY) {
  GAS::ddist_Uni(dY, vTheta, "snorm", log = TRUE)
}, vTheta, dY = dY, method = "Richardson"), ncol = 1)

vGAS = GAS:::Score_Uni(dY, vTheta, "snorm")

expect_equal(dim(vApprox), c(3, 1))
expect_equal(dim(vGAS), c(3, 1))
expect_true(max(vApprox - vGAS) < tol)

# sstd

dY = 0.8
dMu = 1
dSigma = 4
dXi = 1.2
dNu = 7

vTheta = c(dMu, dSigma, dXi, dNu)

vApprox = matrix(numDeriv::grad(function(vTheta, dY) {
  GAS::ddist_Uni(dY, vTheta, "sstd", log = TRUE)
}, vTheta, dY = dY), ncol = 1)

vGAS = GAS:::Score_Uni(dY, vTheta, "sstd")

abs(vApprox[2] - vGAS[2])

expect_equal(dim(vApprox), c(4, 1))
expect_equal(dim(vGAS), c(4, 1))
expect_true(abs(vApprox[2] - vGAS[2]) < tol)

# Skellam

dY = 4
dMu = 1
dSigma2 = 4

vTheta = c(dMu, dSigma2)

vApprox = matrix(numDeriv::grad(function(vTheta, dY) {
  GAS::ddist_Uni(dY, vTheta, "skellam", log = TRUE)
}, vTheta, dY = dY, method = "simple"), ncol = 1)


vGAS = GAS:::Score_Uni(dY, vTheta, "skellam")


# ghskt

dY = 1.5
dMu = 9.0
dSigma = 1.7
dBetaBar = 1.5
dNu = 5.0

vTheta = c(dMu, dSigma, dBetaBar, dNu)

vApprox = matrix(numDeriv::grad(function(vTheta, dY) {
  GAS::ddist_Uni(dY, vTheta, "ghskt", log = TRUE)
}, vTheta, dY = dY), ncol = 1)

vGAS = GAS:::Score_Uni(dY, vTheta, "ghskt")


## multivariate

# mvt

iN = 3

vPhi = rnorm(iN * (iN - 1)/2)
mR = GAS:::MapR_C(vPhi, iN)

vD = exp(rnorm(iN))/10
mD = diag(vD)
vMu = rnorm(iN)
dNu = 4

mSigma = mD %*% mR %*% mD

vY = rnorm(iN)

vR = GAS:::build_vR(mR, iN)

vTheta = c(vMu, vD, vR, dNu)

vGAS = GAS:::Score_Multi(vY, vTheta, "mvt")

vApprox = matrix(numDeriv::grad(function(vTheta, vY, iN) {
  GAS::ddist_Multi(vY, vTheta, "mvt", TRUE)
}, vTheta, vY = vY, iN = iN), ncol = 1)

expect_equal(dim(vApprox), c(10, 1))
expect_equal(dim(vGAS), c(10, 1))
expect_true(max(vApprox - vGAS) < tol)

