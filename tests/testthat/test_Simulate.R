context("Test Simulate")

tol = 1e-4

test_that("Simulation", {

  ## UNIVARIATE

  # norm

  iT = 1e4

  A = matrix(c(0.1,0.0,
               0.0,0.4),2)
  B = matrix(c(0.9,0.0,
               0.0,0.95),2)

  Kappa = (diag(2) - B) %*% GAS::UniUnmapParameters(c(0, 0.1), "norm")

  set.seed(123)
  Sim = GAS::UniGASSim(fit = NULL, iT, Kappa, A, B, Dist = "norm", ScalingType = "Identity")

  tmp = getObs(Sim)
  expect_equal(dim(tmp), c(iT, 1))
  expect_true(abs(sum(tmp) - -129.6511) < tol)

  tmp = max(abs(Sim@GASDyn$mTheta[,100] - c(2.63350183, 0.07058801)))
  expect_true(tmp < tol)

  # snorm

  iT = 1e3

  A = matrix(c(0.0 , 0.0 , 0.0 ,
               0.0 , 0.1 , 0.0 ,
               0.0 , 0.0 , 0.0), 3, byrow = TRUE)

  B = matrix(c(0.7 , 0.0 , 0.0 ,
               0.0 , 0.98, 0.0 ,
               0.0 , 0.0 , 0.97), 3, byrow = TRUE)

  Kappa = (diag(3) - B) %*% GAS::UniUnmapParameters(c(0, 0.1, 1.1), "snorm")

  set.seed(123)
  Sim = GAS::UniGASSim(fit = NULL, iT, Kappa, A, B, Dist = "snorm", ScalingType = "Identity")

  tmp = getObs(Sim)
  expect_equal(dim(tmp), c(iT, 1))
  expect_true(abs(sum(tmp) - 5.506904) < tol)

  tmp = max(abs(Sim@GASDyn$mTheta[,100] - c(0.00000000, 0.1052736, 1.10000000)))
  expect_true(tmp < tol)

  # std

  iT = 1e3
  A = matrix(c(0.001 , 0.0 , 0.0 ,
               0.0 , 0.01 , 0.0 ,
               0.0 , 0.0 , 0.04),3,byrow = T)

  B = matrix(c(0.7 , 0.0 , 0.0 ,
               0.0 , 0.98, 0.0 ,
               0.0 , 0.0 , 0.97),3,byrow = T)

  Kappa = (diag(3) - B) %*% GAS::UniUnmapParameters(c(0, 0.1, 8), "std")

  set.seed(123)
  Sim = GAS::UniGASSim(fit = NULL, iT, Kappa, A, B, Dist = "std", ScalingType = "Identity")

  tmp = getObs(Sim)
  expect_equal(dim(tmp), c(iT, 1))
  expect_true(abs(sum(tmp) - -8.063126) < tol)

  tmp = max(abs(Sim@GASDyn$mTheta[,100] - c(0.0001108389, 0.0991014929, 7.9800370961)))
  expect_true(tmp < tol)

  # sstd

  iT = 1e3

  A = matrix(c(0.0 , 0.0 , 0.0 , 0.0,
               0.0 , 0.1 , 0.0 , 0.0,
               0.0 , 0.0 , 0.00,  0.0,
               0.0 , 0.0 , 0.00,  0.0),4,byrow = T)

  B = matrix(c(0.0 , 0.0 , 0.0 , 0.0,
               0.0 , 0.98, 0.0 , 0.0,
               0.0 , 0.0 , 0.00, 0.0,
               0.0 , 0.0 , 0.00, 0.0),4,byrow = T)

  Kappa = (diag(4) - B) %*% GAS::UniUnmapParameters(c(0, 0.1, 1.1, 8), "sstd")

  set.seed(123)
  Sim = GAS::UniGASSim(fit = NULL, iT, Kappa, A, B, Dist = "sstd", ScalingType = "Identity")

  tmp = getObs(Sim)
  expect_equal(dim(tmp), c(iT, 1))
  expect_true(abs(sum(tmp) - 4.379524) < tol)

  tmp = max(abs(Sim@GASDyn$mTheta[,100] - c(0.0000000, 0.1473415, 1.1000000, 8.0000000)))
  expect_true(tmp < tol)

  ## MULTIVARIATE

  iT   = 1000 # number of observations to simulate
  N    = 3     # trivariate series
  Dist = "mvt" # conditional Multivariate Studen-t distribution

  # build unconditional vector of reparametrised parameters

  Mu  = c(0.1,0.2,0.3)   # vector of location parameters (this is not transformed)
  Phi = c(1.0, 1.2, 0.3) # vector of scale parameters for the firs, second and third variables.

  Rho = c(0.1,0.2,0.3)   # This represents vec(R), where R is the correlation matrix.
  # Note that is up to the user to ensure that vec(R) implies a proper correlation matrix

  Theta = c(Mu, Phi, Rho, 7) # vector of parameters such that the degrees of freedom are 7.

  A = matrix(0, length(Theta), length(Theta))

  diag(A) = c(0, 0, 0, 0.05, 0.01, 0.09, 0.01, 0.04, 0.07, 0) # update scales and correlations, do not update locations and shape parameters

  B = matrix(0, length(Theta), length(Theta))

  diag(B) = c(0, 0, 0, 0.7, 0.7, 0.5, 0.94, 0.97, 0.92, 0) # update scales and correlations, do not update locations and shape parameters

  Kappa = (diag(length(Theta)) - B) %*% GAS::MultiUnmapParameters(Theta, Dist, N)

  set.seed(123)
  Sim = GAS::MultiGASSim(fit = NULL, iT, N, Kappa, A, B, Dist, ScalingType = "Identity")

  tmp = getObs(Sim)
  expect_equal(dim(tmp), c(N, iT))

  tmp = max(abs(Sim@GASDyn$mTheta[, 100] - c(0.10000, 0.20000, 0.30000, 0.92775,
                                             1.21015, 0.27259, 0.08726, 0.09672, 0.50185, 7.00000)))
  expect_true(tmp < tol)

})

