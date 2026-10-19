# Independent standardized two-piece density, for numerical score checks.
skew_reference_parts <- function(theta, Dist) {
  xi <- theta[3]
  student <- Dist == 'sstd'
  nu <- if(student) theta[4] else Inf
  m1 <- if(student) exp(log(2)+.5*log(nu-2)-log(nu-1)-lbeta(.5,nu/2)) else sqrt(2/pi)
  m <- m1*(xi-1/xi)
  c <- sqrt(1+(1-m1^2)*(xi-1/xi)^2)
  b <- theta[2]/c * if(student) sqrt((nu-2)/nu) else 1
  list(a=theta[1]-theta[2]*m/c, b=b, nu=nu)
}

skew_reference_log_density <- function(theta, y, Dist) {
  p <- skew_reference_parts(theta,Dist)
  z <- (y-p$a)/p$b
  x <- ifelse(z < 0,z*theta[3],z/theta[3])
  base <- if(Dist=='sstd') dt(x,p$nu,log=TRUE) else dnorm(x,log=TRUE)
  log(2)-log(theta[3]+1/theta[3])-log(p$b)+base
}

test_that('skew-family scores agree with independent density derivatives', {
  for(Dist in c('snorm','sstd')) for(xi in c(.55,1,1.45)) {
    theta <- c(.2,1.3,xi,if(Dist=='sstd') 3.2)
    for(y in c(-5,-.7,.6,4)) {
      expect_equal(ddist_Uni(y,theta,Dist,log=TRUE),
                   skew_reference_log_density(theta,y,Dist),tolerance=1e-10)
      numeric_score <- numDeriv::grad(function(p) skew_reference_log_density(p,y,Dist),theta)
      expect_equal(as.numeric(Score_Uni(y,theta,Dist)),numeric_score,tolerance=2e-7)
    }
  }
})

test_that('closed-form skew information equals integrated score outer products', {
  for(nu in c(Inf,2.1,7,49.9)) for(xi in c(.55,1,1.45)) {
    Dist <- if(is.infinite(nu)) 'snorm' else 'sstd'
    theta <- c(.2,1.3,xi,if(Dist=='sstd') nu)
    p <- skew_reference_parts(theta,Dist)
    K <- length(theta)
    expected <- matrix(0,K,K)
    mean_score <- numeric(K)
    # Conditional on each side, X has the absolute base distribution.
    for(sign in c(-1,1)) {
      lambda <- xi^sign
      weight <- lambda/(xi+1/xi)
      integrand <- function(x,i,j=NULL) vapply(x,function(xx) {
        s <- as.numeric(Score_Uni(p$a+sign*p$b*lambda*xx,theta,Dist))
        density <- if(Dist=='sstd') dt(xx,nu) else dnorm(xx)
        2*weight*density*s[i]*if(is.null(j)) 1 else s[j]
      },numeric(1))
      for(i in seq_len(K)) {
        mean_score[i] <- mean_score[i]+integrate(integrand,0,Inf,i=i,
          rel.tol=1e-8,abs.tol=1e-9,subdivisions=500L)$value
        for(j in seq_len(K)) expected[i,j] <- expected[i,j]+integrate(integrand,0,Inf,i=i,j=j,
          rel.tol=1e-8,abs.tol=1e-9,subdivisions=500L)$value
      }
    }
    information <- IM_Uni(theta,Dist)
    unit <- sqrt(outer(diag(information),diag(information)))
    expect_lt(max(abs(information-expected)/unit),2e-6)
    expect_lt(max(abs(mean_score)/sqrt(diag(information))),2e-7)
    expect_gt(min(eigen(information,symmetric=TRUE,only.values=TRUE)$values),0)
  }
})

test_that('skew information remains positive definite near parameter bounds', {
  for(nu in c(Inf,2.0001,2.01,49.999)) for(xi in c(.5001,1,1.4999))
    for(sigma in c(.01,1,100)) {
      Dist <- if(is.infinite(nu)) 'snorm' else 'sstd'
      theta <- c(0,sigma,xi,if(Dist=='sstd') nu)
      information <- IM_Uni(theta,Dist)
      normalized <- information/sqrt(outer(diag(information),diag(information)))
      expect_true(all(is.finite(information)))
      expect_gt(min(eigen(normalized,symmetric=TRUE,only.values=TRUE)$values),0)
    }
  # At symmetry, the standardized normal's mean and variance scores are
  # orthogonal to skewness; skewness still has positive information.
  expect_equal(IM_Uni(c(0,1,1),'snorm'),diag(c(1,2,3-8/pi)),tolerance=1e-12)
  expect_error(IM_Uni(c(0,1,1,2),'sstd'),'Invalid parameters')
})

test_that('all skew scalings follow the full transformed-state recursion', {
  y <- c(-2,.3,1.7,-.2,2.5)
  for(Dist in c('snorm','sstd')) for(scaling in c('Identity','Inv','InvSqrt')) {
    theta <- c(.2,1.3,.8,if(Dist=='sstd') 7)
    f <- as.numeric(UniUnmapParameters(theta,Dist))
    K <- length(f); A <- diag(.001,K); B <- diag(.8,K); kappa <- .2*f
    filtered <- GAS:::GASFilter_univ(y,kappa,A,B,length(y),K,Dist,scaling)
    path <- matrix(NA_real_,K,length(y)+1L); path[,1] <- theta
    for(t in seq_along(y)) {
      theta <- as.numeric(UniMapParameters(f,Dist))
      score <- numDeriv::grad(function(p) skew_reference_log_density(p,y[t],Dist),theta)
      J <- numDeriv::jacobian(function(p) as.numeric(UniMapParameters(p,Dist)),f)
      score <- as.numeric(t(J)%*%score)
      information <- t(J)%*%IM_Uni(theta,Dist)%*%J
      innovation <- switch(scaling,Identity=score,Inv=solve(information,score),
                           InvSqrt=solve(t(chol(information)),score))
      expect_equal(as.numeric(filtered$mInnovations[,t]),as.numeric(innovation),tolerance=2e-6)
      f <- as.numeric(kappa+A%*%innovation+B%*%f)
      path[,t+1L] <- UniMapParameters(f,Dist)
    }
    expect_equal(filtered$mTheta,path,tolerance=2e-6)
    expect_match(DistScalingType(Dist),scaling,fixed=TRUE)
    expect_s4_class(UniGASSpec(Dist,scaling),'uGASSpec')
  }
})

test_that('skew simulation and forecasts use the same scaled recursion', {
  set.seed(428)
  for(Dist in c('snorm','sstd')) for(scaling in c('Inv','InvSqrt')) {
    theta <- c(.2,1.3,.8,if(Dist=='sstd') 7)
    f <- as.numeric(UniUnmapParameters(theta,Dist)); K <- length(f)
    A <- diag(.001,K); B <- diag(.8,K); kappa <- .2*f
    sim <- UniGASSim(T.sim=10L,kappa=kappa,A=A,B=B,Dist=Dist,ScalingType=scaling)
    filtered <- GAS:::GASFilter_univ(sim@GASDyn$vY,kappa,A,B,10L,K,Dist,scaling)
    expect_equal(sim@GASDyn$mTheta,filtered$mTheta,tolerance=1e-12)
    forecast <- uGASMultiForcast(theta,kappa,A,B,5L,2L,K,Dist,scaling,TRUE)
    for(b in 1:2) {
      path <- GAS:::GASFilter_univ(forecast$mY[b,],kappa,A,B,5L,K,Dist,scaling)
      expect_equal(forecast$cTheta[,,b],path$mTheta[,1:5],tolerance=1e-10)
    }
  }
})
