# Independent R log density, used for numerical score and information checks.
multi_reference_log_density <- function(theta, y, Dist) {
  n <- ncol(y)
  R <- diag(n)
  k <- 2L*n
  for (i in seq_len(n-1L)) for (j in (i+1L):n) {
    k <- k+1L
    R[i,j] <- R[j,i] <- theta[k]
  }
  D <- diag(theta[n+seq_len(n)])
  Sigma <- D %*% R %*% D
  z <- sweep(y, 2L, theta[seq_len(n)])
  q <- rowSums((z %*% solve(Sigma))*z)
  logdet <- as.numeric(determinant(Sigma, logarithm=TRUE)$modulus)
  if (Dist == "mvnorm") return(-.5*(n*log(2*pi)+logdet+q))
  nu <- tail(theta,1L)
  lgamma((nu+n)/2)-lgamma(nu/2)-.5*n*log(nu*pi)-.5*logdet-
    .5*(nu+n)*log1p(q/nu)
}

test_that("multivariate information agrees with independent density scores", {
  set.seed(481)
  N <- 30000L
  theta0 <- c(.2,-.1,.3,1.1,.8,1.3,.3,-.15,.2)
  R <- matrix(c(1,.3,-.15,.3,1,.2,-.15,.2,1),3L)
  D <- diag(theta0[4:6])
  for (nu in c(Inf,1.5,7,30)) {
    Dist <- if(is.infinite(nu)) "mvnorm" else "mvt"
    theta <- if(Dist=="mvt") c(theta0,nu) else theta0
    z <- matrix(rnorm(3L*N),N,3L) %*% chol(D%*%R%*%D)
    if (Dist=="mvt") z <- z/sqrt(rchisq(N,nu)/nu)
    y <- sweep(z,2L,theta[1:3],"+")
    S <- numDeriv::jacobian(function(p) multi_reference_log_density(p,y,Dist),theta)
    H <- GAS:::IM_multi(theta,Dist,3L)
    empirical <- crossprod(S)/N
    # Monte Carlo uncertainty, including the off-diagonal and shape blocks.
    se <- sqrt(pmax(0,crossprod(S^2)/N-empirical^2)/N)
    expect_lt(max(abs(H-empirical)/(se+1e-10)),6.5)
    expect_gt(min(eigen(H,symmetric=TRUE,only.values=TRUE)$values),0)
  }
})

test_that("one-dimensional information reduces to existing univariate formulas", {
  for(nu in c(.0011,.1,1,2,4,10,49.9)) {
    # std uses scale squared; mvt uses scale.
    J <- diag(c(1,2*1.3,1))
    expected <- t(J)%*%IM_Uni(c(.2,1.3^2,nu),"std")%*%J
    expect_equal(GAS:::IM_multi(c(.2,1.3,nu),"mvt",1L),expected,tolerance=1e-10)
  }
  expect_equal(GAS:::IM_multi(c(.2,1.3),"mvnorm",1L),diag(c(1,2)/1.3^2))
})

test_that("all multivariate scaling mechanisms follow the full state recursion", {
  y <- rbind(c(.4,-.6,.2),c(-.1,.8,-.7),c(.6,.3,-.2))
  theta0 <- c(.1,-.2,.05,1.1,.8,1.3,.3,-.15,.2)
  for (Dist in c("mvnorm","mvt")) for (scaling in c("Identity","Inv","InvSqrt")) {
    theta <- if(Dist=="mvt") c(theta0,7) else theta0
    f <- as.numeric(MultiUnmapParameters(theta,Dist,3L))
    K <- length(f)
    A <- diag(.001,K); B <- diag(.8,K); kappa <- .2*f
    filtered <- GAS:::GASFilter_multi(t(y),kappa,A,B,nrow(y),3L,K,Dist,scaling)
    path <- matrix(NA_real_,K,nrow(y)+1L)
    path[,1L] <- theta
    for(t in seq_len(nrow(y))) {
      theta <- as.numeric(MultiMapParameters(f,Dist,3L))
      score <- numDeriv::grad(function(p) multi_reference_log_density(p,y[t,,drop=FALSE],Dist),theta)
      J <- numDeriv::jacobian(function(p) as.numeric(MultiMapParameters(p,Dist,3L)),f)
      score <- as.numeric(t(J)%*%score)
      H <- t(J)%*%GAS:::IM_multi(theta,Dist,3L)%*%J
      innovation <- switch(scaling,Identity=score,Inv=solve(H,score),
                           InvSqrt=solve(t(chol(H)),score))
      expect_equal(as.numeric(filtered$mInnovations[,t]),as.numeric(innovation),tolerance=2e-6)
      f <- as.numeric(kappa+A%*%innovation+B%*%f)
      path[,t+1L] <- MultiMapParameters(f,Dist,3L)
    }
    expect_equal(filtered$mTheta,path,tolerance=2e-6)
    expect_equal(as.numeric(filtered$vLLK),vapply(seq_len(nrow(y)),function(t)
      multi_reference_log_density(path[,t],y[t,,drop=FALSE],Dist),numeric(1)),tolerance=2e-6)
  }
})

test_that("multivariate scaling validation and distribution metadata agree", {
  for (Dist in c("mvnorm","mvt")) for (scaling in c("Identity","Inv","InvSqrt")) {
    expect_s4_class(MultiGASSpec(Dist,scaling),"mGASSpec")
    expect_match(DistScalingType(Dist),scaling,fixed=TRUE)
  }
  for(scaling in list("bad",NA_character_,NULL,c("Inv","InvSqrt")))
    expect_error(MultiGASSpec(ScalingType=scaling),"ScalingType must")
})

test_that('multivariate simulation and forecasts use the same scaled recursion', {
  set.seed(512)
  for(Dist in c('mvnorm','mvt')) for(scaling in c('Inv','InvSqrt')) {
    theta <- c(.1,-.2,1.1,.8,.3)
    if(Dist=='mvt') theta <- c(theta,7)
    f <- as.numeric(MultiUnmapParameters(theta,Dist,2L)); K <- length(f)
    A <- diag(seq(.0001,.0005,length.out=K)); B <- diag(.8,K); kappa <- .2*f
    sim <- MultiGASSim(T.sim=8L,N=2L,kappa=kappa,A=A,B=B,Dist=Dist,ScalingType=scaling)
    filtered <- GAS:::GASFilter_multi(sim@GASDyn$mY,kappa,A,B,8L,2L,K,Dist,scaling)
    expect_equal(sim@GASDyn$mTheta,filtered$mTheta,tolerance=1e-12)
    # The exported C++ forecast wrapper starts from the supplied natural state.
    forecast <- mGASMultiForcast(theta,kappa,A,B,3L,2L,K,2L,Dist,scaling,TRUE)
    for(b in 1:2) {
      path <- GAS:::GASFilter_multi(forecast$cY[,,b],kappa,A,B,3L,2L,K,Dist,scaling)
      expect_equal(forecast$cTheta[,,b],path$mTheta[,1:3],tolerance=1e-10)
    }
  }
})

test_that('small dispersion matrices are evaluated without likelihood repairs', {
  y <- matrix(c(.4,-.6,.2),1L)
  for(s in c(1e-7,1e-12)) for(nu in c(.0011,.5,7)) {
    theta <- c(0,0,0,rep(s,3),.3,-.15,.2,nu)
    expect_equal(ddist_Multi(as.numeric(y),theta,'mvt',log=TRUE),
                 as.numeric(multi_reference_log_density(theta,y,'mvt')),tolerance=1e-10)
  }
  # Non-positive-definite matrices are not silently changed into valid ones.
  bad <- c(0,0,1,1,1.1)
  expect_error(ddist_Multi(c(.2,-.1),bad,'mvnorm',log=TRUE),'positive definite')
  expect_error(ddist_Multi(c(.2,-.1),c(bad,7),'mvt',log=TRUE),'positive definite')
})
