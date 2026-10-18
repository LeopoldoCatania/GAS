conditional_fixture <- function() {
  spec <- UniGASSpec(Dist="norm", GASPar=list(location=TRUE,scale=FALSE))
  b <- .7; mu <- .3; v <- 1.2
  p <- c(kappa1=mu*(1-b), kappa2=log(v), a1=-1000,
         b1=qlogis((b-LowerB())/(UpperB()-LowerB())))
  spec@Spec$PwNames <- names(p)
  set.seed(481)
  y <- sample(qnorm((seq_len(400)-.5)/400))*sqrt(v)+mu
  list(spec=spec,p=p,y=y,b=b,mu=mu,v=v,
       opg=GAS:::GAS_opg(p,y,spec)$opg)
}

test_that("inactive dynamics yield the analytic conditional Gaussian covariance", {
  f <- conditional_fixture()
  z <- f$y-f$mu
  S <- cbind(z/(f$v*(1-f$b)), .5*(z^2/f$v-1))
  expected <- solve(crossprod(S))
  result <- GAS:::GAS_opg_inference(f$opg,f$p,f$spec)
  expect_false(result$full_info$ok)
  expect_true(result$info$ok)
  expect_true(result$conditional)
  expect_identical(result$excluded,c("a1","b1"))
  expect_identical(result$included,c("kappa1","kappa2"))
  expect_equal(unname(result$covariance[1:2,1:2]),expected,tolerance=1e-6)
  expect_true(all(is.na(result$covariance[c("a1","b1"),])))
  table <- GAS:::InferenceFun_Uni(f$opg,f$p,2L,TRUE,result$included)
  expect_equal(unname(table[1:2,"Std. Error"]),sqrt(diag(expected)),tolerance=1e-6)
  expect_true(all(is.na(table[c("a1","b1"),-1])))
  expect_equal(unname(table[c("a1","b1"),"Estimate"]),c(0,f$b),tolerance=1e-10)
})

test_that("fallback is limited to near-zero pairs and a valid retained block", {
  f <- conditional_fixture()
  disabled <- GAS:::GAS_opg_inference(f$opg,f$p,f$spec,SE.A.tol=0)
  expect_false(disabled$info$ok)
  expect_null(disabled$covariance)
  expect_length(disabled$excluded,0)
  # A near-zero loading alone must not cause exclusions if the full OPG passes.
  joint <- GAS:::GAS_opg_inference(diag(4),f$p,f$spec)
  expect_true(joint$info$ok)
  expect_false(joint$conditional)
  expect_length(joint$excluded,0)
  p <- f$p; p["a1"] <- qlogis(.1/10)
  unrelated <- GAS:::GAS_opg_inference(f$opg,p,f$spec)
  expect_false(unrelated$info$ok)
  expect_length(unrelated$excluded,0)
  bad <- f$opg; bad[2,] <- bad[1,]; bad[,2] <- bad[,1]
  reduced_bad <- GAS:::GAS_opg_inference(bad,f$p,f$spec)
  expect_false(reduced_bad$info$ok)
  expect_false(reduced_bad$conditional)
  expect_null(reduced_bad$covariance)
  for(tol in list(-1,Inf,NA_real_,c(1,2),"small"))
    expect_error(GAS:::GAS_validate_se_a_tol(tol),"SE.A.tol")
})

test_that("multivariate exclusions preserve scalar and diagonal coefficient grouping", {
  for(scalar in c(TRUE,FALSE)) {
    spec <- MultiGASSpec(Dist="mvnorm",ScalarParameters=scalar,
                        GASPar=list(location=TRUE,scale=TRUE,correlation=FALSE))
    suffix <- if(scalar) c("mu","sigma") else c("mu1","mu2")
    nm <- c("kappa.mu1","kappa.mu2",paste0("a.",suffix),paste0("b.",suffix))
    p <- setNames(c(.1,.2,-1000,qlogis(.1/10),0,0),nm)
    opg <- diag(6); opg[3,3] <- 0
    result <- GAS:::GAS_opg_inference(opg,p,spec)
    expect_identical(result$excluded,nm[c(3,5)])
    expect_identical(result$included,nm[c(1,2,4,6)])
    expect_true(result$conditional)
    expect_true(result$info$ok)
  }
})

test_that("conditional bands hold excluded coefficients fixed in every draw", {
  f <- conditional_fixture()
  y <- f$y
  p <- f$p
  p["kappa1"] <- mean(y)*(1-f$b)
  p["kappa2"] <- log(mean((y-mean(y))^2))
  static <- UniGASSpec(Dist="norm",GASPar=list(location=FALSE,scale=FALSE))
  fit <- UniGASFit(static,y)
  fit@ModelInfo$Spec <- f$spec
  expect_warning(opt <- GAS:::GAS_finalize_optimiser(
    list(pars=p,convergence=0L),p,y,f$spec,GAS:::UniGASOptimiser,TRUE),
    "Conditional OPG standard errors")
  fit@Estimates$optimiser <- opt
  par <- GAS:::AddFixedPar(GAS:::vPw2lPn_Uni(p,2L))
  fit@GASDyn <- GAS:::GASFilter_univ(y,par$vKappa,par$mA,par$mB,
                                    length(y),2L,"norm","Identity")
  B <- 20L; probs <- c(.1,.9)
  set.seed(1234)
  expect_warning(bands <- ConfidenceBands(fit,B=B,probs=probs),
                 "Conditional confidence bands")
  set.seed(1234)
  draws <- sweep(matrix(rnorm(B*2),nrow=B)%*%chol(opt$covariance[1:2,1:2]),2,p[1:2],"+")
  expected <- cbind(quantile(draws[,1]/(1-f$b),probs),
                    quantile(exp(draws[,2]),probs))
  expect_equal(unname(bands[1,,]),unname(expected),tolerance=1e-8)
  expect_equal(unname(bands[length(y)+1,,]),unname(expected),tolerance=1e-8)
  expect_true(all(is.finite(bands)))
})
