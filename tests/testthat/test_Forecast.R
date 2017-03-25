context("Test Forecast") 

tol = 1e-4

test_that("Forecasting", { 
  
  data("cpichg")
  
  GASSpec = GAS::UniGASSpec(Dist = "std", ScalingType = "Identity", 
                       GASPar = list(location = TRUE, scale = TRUE, shape = FALSE))
  
  Fit = GAS::UniGASFit(GASSpec, cpichg)
  
  set.seed(1234)
  Forecast = GAS::UniGASFor(Fit, H = 100)
  
  tmp = max(abs(Forecast@Forecast$Moments[100,] - c(0.04033691, 1.06576392, 0.00000000, 5.37511686)))
  expect_true(tmp < tol)
  
  tmp = max(abs(Forecast@Forecast$PointForecast[100,] - c(0.04033691, 0.73915263, 6.52619149)))
  expect_true(tmp < tol)
  
  # Perform 1-Step ahead rolling forecast
  
  InsampleData  = cpichg[1:250]
  OutSampleData = cpichg[251:276]
  
  Fit = GAS::UniGASFit(GASSpec, InsampleData)
  
  set.seed(1234)
  Forecast = GAS::UniGASFor(Fit, Roll = TRUE, out = OutSampleData)
  
  tmp = max(abs(Forecast@Forecast$Moments[10,] - c(0.6129428, 0.229406, 0, 5.526927)))
  expect_true(tmp < tol)
  
  tmp = max(abs(Forecast@Forecast$PointForecast[10,] - c(0.6129428, 0.157429, 6.374426)))
  expect_true(tmp < tol)

}) 