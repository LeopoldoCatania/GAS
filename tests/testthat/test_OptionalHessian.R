# Count actual stats::optimHess calls, including initialization. A NULL
# returned Hessian alone would not establish that its calculation was skipped.
with_hessian_counter <- function(code) {
  counter <- new.env(parent=emptyenv())
  counter$n <- 0L
  suppressMessages(trace('optimHess', where=asNamespace('stats'), print=FALSE,
    tracer=function() { counter$n <- counter$n+1L }))
  on.exit(suppressMessages(untrace('optimHess',where=asNamespace('stats'))))
  code(counter)
}

test_that('optional Hessians have the correct scale and are truly skipped', {
  objective <- function(p,data,GASSpec) GASSpec$iT*sum((p-c(.2,-.3))^2)
  with_hessian_counter(function(counter) {
    off <- fn.optim(c(a=1,b=2),NULL,list(iT=40L),objective)
    expect_identical(counter$n,0L)
    expect_null(off$hessian)
    on <- fn.optim(c(a=1,b=2),NULL,list(iT=40L),objective,Compute.Hessian=TRUE)
    expect_identical(counter$n,2L)
    expect_equal(unname(on$hessian),diag(80,2),tolerance=1e-6)
    expect_equal(off$pars,on$pars,tolerance=1e-9)
    expect_equal(off$value,on$value,tolerance=1e-9)
  })
})

test_that('the option reaches static and dynamic fits independently of OPG', {
  set.seed(91)
  y <- rgamma(200,shape=4,rate=2)
  spec <- UniGASSpec('gamma',GASPar=list(location=FALSE,scale=FALSE))
  with_hessian_counter(function(counter) {
    off <- UniGASFit(spec,y,Compute.SE=TRUE)
    expect_identical(counter$n,0L)
    expect_null(off@Estimates$optimiser$hessian)
    expect_null(off@Estimates$StaticFit$optimiser$hessian)
    expect_true(off@Estimates$optimiser$diagnostics$inference_ok)
    expect_true(all(is.finite(off@Estimates$Inference[,'Std. Error'])))
    on <- UniGASFit(spec,y,Compute.SE=TRUE,Compute.Hessian=TRUE)
    expect_gte(counter$n,4L)
    expect_true(is.matrix(on@Estimates$optimiser$hessian))
    expect_true(is.matrix(on@Estimates$StaticFit$optimiser$hessian))
    expect_equal(off@Estimates$optimiser$pars,on@Estimates$optimiser$pars,tolerance=1e-8)
    expect_equal(off@Estimates$optimiser$opg,on@Estimates$optimiser$opg,tolerance=1e-8)
    no_se <- UniGASFit(spec,y,Compute.SE=FALSE,Compute.Hessian=TRUE)
    expect_true(is.matrix(no_se@Estimates$optimiser$hessian))
    expect_null(no_se@Estimates$optimiser$opg)
  })
})

test_that('analytical multivariate Gaussian fits respect the Hessian option', {
  set.seed(163)
  y <- matrix(rnorm(400),200,2)
  spec <- MultiGASSpec('mvnorm',GASPar=list(location=FALSE,scale=FALSE,correlation=FALSE))
  with_hessian_counter(function(counter) {
    off <- MultiGASFit(spec,y)
    expect_identical(counter$n,0L)
    expect_null(off@Estimates$optimiser$hessian)
    expect_true(off@Estimates$optimiser$diagnostics$inference_ok)
    on <- MultiGASFit(spec,y,Compute.Hessian=TRUE)
    expect_identical(counter$n,1L)
    expect_true(all(is.finite(on@Estimates$optimiser$hessian)))
    expect_equal(off@Estimates$optimiser$opg,on@Estimates$optimiser$opg)
  })
})

test_that('legacy custom optimizers keep their contract and opt-in ones receive the flag', {
  legacy <- function(par0,data,GASSpec,FUN) list(pars=par0)
  expect_identical(GAS:::GAS_optimizer_hessian(legacy,FALSE),legacy)
  expect_identical(GAS:::GAS_optimizer_hessian(legacy,TRUE),legacy)
  for(flag in c(FALSE,TRUE)) {
    seen <- logical()
    custom <- function(par0,data,GASSpec,FUN,Compute.Hessian) {
      seen <<- c(seen,Compute.Hessian)
      fn.optim(par0,data,GASSpec,FUN,Compute.Hessian=Compute.Hessian)
    }
    set.seed(18)
    fit <- UniGASFit(UniGASSpec('gamma',GASPar=list(location=FALSE,scale=FALSE)),
                     rgamma(100,shape=4,rate=2),fn.optimizer=custom,Compute.SE=FALSE,Compute.Hessian=flag)
    expect_identical(seen,c(flag,flag))
    expect_identical(is.null(fit@Estimates$optimiser$hessian),!flag)
  }
})

test_that('invalid Hessian flags fail before optimization', {
  for(flag in list(NA,NULL,1,'FALSE',c(TRUE,FALSE))) {
    expect_error(UniGASFit(UniGASSpec(),1:20,Compute.Hessian=flag),'Compute.Hessian must')
    expect_error(MultiGASFit(MultiGASSpec(),matrix(1:40,20),Compute.Hessian=flag),'Compute.Hessian must')
    expect_error(fn.optim(NULL,NULL,NULL,NULL,Compute.Hessian=flag),'Compute.Hessian must')
    expect_error(fn.solnp(NULL,NULL,NULL,NULL,Compute.Hessian=flag),'Compute.Hessian must')
  }
})

test_that('SOLNP retains its solution while making Hessian output optional', {
  objective <- function(p,data,GASSpec) sum((p-c(.2,-.3))^2)
  off <- fn.solnp(c(a=1,b=2),NULL,list(),objective)
  on <- fn.solnp(c(a=1,b=2),NULL,list(),objective,Compute.Hessian=TRUE)
  expect_null(off$hessian)
  expect_true(is.matrix(on$hessian))
  expect_true(all(is.finite(on$hessian)))
  expect_equal(off$pars,on$pars)
})
