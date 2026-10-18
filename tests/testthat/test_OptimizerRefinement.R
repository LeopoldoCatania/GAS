test_that("small intercepts receive an accurate final optimization", {
  # A skewed, sharply curved likelihood with a known optimum. Default
  # optim differences (1e-3) bias its intercept gradient even near the optimum.
  objective <- function(p, data, GASSpec) {
    z <- 100 * (p["intercept"] - .001)
    GASSpec$iT * (3 + exp(z) - z - 1 + (p["other"] - .3)^2)
  }
  start <- c(intercept = .00101, other = .30001)
  opt <- fn.optim(start, NULL, list(iT = 50L), objective)
  # Check against analytic derivatives, independently of numDeriv.
  analytic_gradient <- c(100 * expm1(100 * (opt$pars["intercept"] - .001)),
                         2 * (opt$pars["other"] - .3))
  expect_identical(names(opt$pars), names(start))
  expect_equal(opt$convergence, 0)
  expect_lt(max(abs(analytic_gradient)), 1e-5)
  expect_equal(unname(opt$value), 150, tolerance = 1e-9)
  expect_equal(unname(opt$pars), c(.001, .3), tolerance = 1e-5)
})
