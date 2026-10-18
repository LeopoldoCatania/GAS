test_that("Gaussian OPG matches analytic scores and sample-size scaling", {
  y <- qnorm((seq_len(300)-.5)/300) + .3
  spec <- UniGASSpec(Dist="norm", GASPar=list(location=FALSE, scale=FALSE))
  fit <- UniGASFit(spec, y)
  opt <- fit@Estimates$optimiser
  mu <- opt$pars["kappa1"]
  v <- exp(opt$pars["kappa2"])
  S <- cbind((y-mu)/v, .5*((y-mu)^2/v-1))
  expect_equal(unname(opt$scores), unname(S), tolerance=1e-6)
  expect_equal(unname(opt$opg), unname(crossprod(S)), tolerance=1e-6)
  V <- solve(crossprod(S))
  expect_equal(unname(opt$covariance), unname(V), tolerance=1e-6)
  expect_equal(unname(fit@Estimates$Inference[,"Std. Error"]),
               unname(sqrt(diag(V))), tolerance=1e-6)
  twice <- GAS:::GAS_opg(opt$pars, rep(y,2), fit@ModelInfo$Spec)
  expect_equal(unname(solve(twice$opg)), unname(V/2), tolerance=1e-6)
  expect_equal(-colMeans(opt$scores),
               setNames(opt$diagnostics$gradient_per_observation, names(opt$pars)),
               tolerance=1e-6)
  expect_true(opt$diagnostics$inference_ok)
  expect_identical(opt$covariance_method, "OPG")
})

test_that("dynamic scores differentiate initialization and the complete recursion", {
  spec <- UniGASSpec(Dist="norm", GASPar=list(location=FALSE,scale=TRUE))
  p <- c(kappa1=.2,kappa2=-.03,
         a2=qlogis((.15-LowerA())/(UpperA()-LowerA())),
         b2=qlogis((.8-LowerB())/(UpperB()-LowerB())))
  spec@Spec$PwNames <- names(p)
  y <- c(-.8,.1,1.5,-2,.2,.7,-.3,.8,.4,-1.2)
  reference <- function(p) {
    a <- LowerA()+(UpperA()-LowerA())*plogis(p["a2"])
    b <- LowerB()+(UpperB()-LowerB())*plogis(p["b2"])
    f <- p["kappa2"]/(1-b)
    out <- numeric(length(y))
    for(t in seq_along(y)) {
      v <- exp(f); z <- y[t]-p["kappa1"]
      out[t] <- dnorm(y[t], p["kappa1"], sqrt(v), log=TRUE)
      f <- p["kappa2"] + a*.5*(z^2/v-1) + b*f
    }
    out
  }
  opg <- GAS:::GAS_opg(p,y,spec)
  expect_equal(GAS:::GAS_loglik_contributions(p,y,spec), reference(p),tolerance=1e-10)
  # Independent central differences of the R recursion.
  h <- 1e-5
  S <- sapply(seq_along(p),function(j) {
    plus <- minus <- p; plus[j] <- p[j]+h; minus[j] <- p[j]-h
    (reference(plus)-reference(minus))/(2*h)
  })
  expect_equal(unname(opg$scores), S, tolerance=1e-6)
  expect_equal(sum(reference(p)), -GAS:::UniGASOptimiser(p,y,spec),tolerance=1e-10)
  # Check the delta method against the analytic logistic derivatives.
  a <- .15; b <- .8
  J <- diag(c(1,1,(a-LowerA())*(UpperA()-a)/(UpperA()-LowerA()),
              (b-LowerB())*(UpperB()-b)/(UpperB()-LowerB())))
  tab <- GAS:::InferenceFun_Uni(opg$opg,p,2L)
  expect_equal(unname(tab[,"Std. Error"]),
               sqrt(diag(J%*%solve(opg$opg)%*%t(J))),tolerance=1e-6)
})

test_that("multivariate OPG uses joint observation contributions", {
  set.seed(313)
  y <- matrix(rnorm(800),ncol=2)
  y[,2] <- .4*y[,1]+.8*y[,2]
  spec <- MultiGASSpec(Dist="mvnorm",GASPar=list(location=FALSE,scale=FALSE,correlation=FALSE))
  fit <- MultiGASFit(spec,y)
  opt <- fit@Estimates$optimiser
  p <- opt$pars
  independent <- function(p) {
    theta <- as.numeric(MultiMapParameters(p,"mvnorm",2L))
    mu <- theta[1:2]; sd <- theta[3:4]; rho <- theta[5]
    cov <- diag(sd)%*%matrix(c(1,rho,rho,1),2)%*%diag(sd)
    z <- sweep(y,2,mu)
    -log(2*pi)-.5*as.numeric(determinant(cov,logarithm=TRUE)$modulus)-
      .5*rowSums((z%*%solve(cov))*z)
  }
  expect_equal(GAS:::GAS_loglik_contributions(p,t(y),fit@ModelInfo$Spec), independent(p), tolerance=1e-8)
  expect_equal(unname(opt$scores), numDeriv::jacobian(independent,p), tolerance=1e-6)
  expect_equal(-unname(colMeans(opt$scores)),opt$diagnostics$gradient_per_observation,tolerance=1e-6)
  expect_true(opt$diagnostics$inference_ok)
  expect_true(all(is.finite(fit@Estimates$Inference[,"Std. Error"])))
  expect_identical(dim(opt$scores),c(nrow(y),length(p)))
  expect_true(all(is.finite(ConfidenceBands(fit,B=3,probs=c(.1,.9)))))
})

test_that("inference ignores optimizer Hessians and respects Compute.SE", {
  y <- qnorm((seq_len(200)-.5)/200)
  spec <- UniGASSpec(Dist="norm",GASPar=list(location=FALSE,scale=FALSE))
  bad_hessian <- function(par0,data,GASSpec,FUN) {
    out <- fn.optim(par0,data,GASSpec,FUN)
    out$hessian <- diag(-1,length(par0))
    out
  }
  fit <- UniGASFit(spec,y,fn.optimizer=bad_hessian)
  expect_true(fit@Estimates$optimiser$diagnostics$inference_ok)
  expect_true(all(is.finite(fit@Estimates$Inference[,"Std. Error"])))
  expect_true(all(diag(fit@Estimates$optimiser$hessian)<0))
  expect_true(all(is.finite(ConfidenceBands(fit,B=3,probs=c(.1,.9)))))
  skip <- UniGASFit(spec,y,Compute.SE=FALSE)
  expect_null(skip@Estimates$optimiser$scores)
  expect_null(skip@Estimates$optimiser$opg)
  expect_false(skip@Estimates$optimiser$diagnostics$inference_ok)
  expect_true(all(is.na(skip@Estimates$Inference[,"Std. Error"])))
  expect_error(ConfidenceBands(skip,B=2),"Refit with Compute.SE")
  expect_false(GAS:::GAS_information_info(matrix(1,2,2))$ok)
  expect_null(GAS:::GAS_checked_covariance(matrix(1,2,2),2L))
})
