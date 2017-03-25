context("Test Estimate")

tol = 1e-4

test_that("Estimation", {

  # Univariate Example

  data("sp500ret", package = "GAS")

  GASSpec = GAS::UniGASSpec(Dist = "std", ScalingType = "Identity",
                            GASPar = list(location = FALSE, scale = TRUE, shape = FALSE))

  Fit = GAS::UniGASFit(GASSpec, sp500ret)

  tmp = max(abs(Fit@Estimates$optimiser$pars - c(0.056934913, -0.007766315, -2.634285043, -4.028353135, 4.579973587)))
  expect_true(tmp < tol)

  # Multivariate Example

  data("StockIndices", package = "GAS")
  GASSpec = GAS::MultiGASSpec(Dist = "mvt", ScalingType = "Identity",
                              GASPar = list(location    = FALSE, scale = TRUE,
                                            correlation = TRUE,  shape = FALSE),
                              ScalarParameters = FALSE)

  Fit = GAS::MultiGASFit(GASSpec, StockIndices)

  tmp = max(abs(Fit@Estimates$optimiser$pars - c(6.675485e-02,  3.642856e-02,  1.320821e-02,  1.353661e-03,
                                                 9.434927e-06,  4.982442e-03,  1.154864e-02,  2.973417e-02,
                                                 9.272084e-03, -2.696278e+00, -7.175445e+00, -7.757717e+00,
                                                 -6.790044e+00, -8.268374e+00, -7.038423e+00, -8.112272e+00,
                                                 4.600537e+00,  7.655267e+00,  4.089359e+00,  3.410833e+00,
                                                 2.834004e+00,  4.721346e+00 )))
  expect_true(tmp < tol)


  # Specification mvt with scalar parameters

  GASSpec = GAS::MultiGASSpec(Dist = "mvt", ScalingType = "Identity",
                              GASPar = list(location = FALSE, scale = TRUE,
                                            correlation = TRUE, shape = FALSE),
                              ScalarParameters = TRUE)

  Fit = GAS::MultiGASFit(GASSpec, StockIndices)

  tmp = max(abs(Fit@Estimates$optimiser$pars - c(0.0689883072,  0.0373772602,  0.0173364589,  0.0009332931,
                                                 0.0012413557,  0.0019411287,  0.1871876929,  0.2694115828,
                                                 0.5252761292, -2.7707074972, -7.0965094385, -7.2230917292,
                                                 4.8487286766, -0.0137374433)))
  expect_true(tmp < tol)

  # CPI

  data("cpichg", package = "GAS")

  GASSpec = GAS::UniGASSpec(Dist = "std", ScalingType = "Identity",
                            GASPar = list(location = TRUE, scale = TRUE, shape = TRUE))

  Fit = GAS::UniGASFit(GASSpec, cpichg)

  tmp = max(abs(Fit@Estimates$optimiser$pars - c(0.04955147, -0.23864318, -0.34812790, -4.86529113,
                                                 -3.12422661, -1.86188241,  2.68124270,  1.94666403,
                                                 3.89099313)))
  expect_true(tmp < tol)
})

