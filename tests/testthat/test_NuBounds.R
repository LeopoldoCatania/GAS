test_that("R bounds delegate to the C++ definitions", {
  expected <- c(std=.001, ast=.001, ast1=.001, mvt=.001, sstd=2, ghskt=4)
  for (d in names(expected)) {
    expect_equal(LowerNu(d), unname(expected[d]))
    expect_equal(UpperNu(d), 50)
    expect_identical(LowerNu(d), GAS:::LowerNu_cpp(d))
    expect_identical(UpperNu(d), GAS:::UpperNu_cpp(d))
  }
  expect_equal(LowerNu(), .001)
  expect_equal(UpperNu(), 50)
  expect_error(LowerNu("norm"), "not defined")
  expect_error(UpperNu("unknown"), "not defined")
})

nu_fixtures <- list(std=c(.2,1.3,8), sstd=c(.2,1.3,1.1,8),
                    ast=c(.2,1.3,.4,8,12), ast1=c(.2,1.3,.4,8),
                    ghskt=c(.2,1.3,.2,9))
nu_indices <- list(std=3L,sstd=4L,ast=4:5,ast1=4L,ghskt=4L)

test_that("each shape map and inverse uses the correct open interval", {
  for(d in names(nu_fixtures)) {
    ix<-nu_indices[[d]];lo<-LowerNu(d);hi<-UpperNu(d)
    vals<-unique(c(lo+.01, c(.5,1.5,2.5,3.5,4.005,8,49)))
    vals<-vals[vals>lo & vals<hi]
    for(v in vals) {
      th<-nu_fixtures[[d]];th[ix]<-v
      f<-as.numeric(UniUnmapParameters(th,d))
      expect_equal(as.numeric(UniMapParameters(f,d)),th,tolerance=1e-10)
      expect_equal(f[ix],rep(log((v-lo)/(hi-v)),length(ix)),tolerance=1e-10)
    }
    for(x in c(-1000,-100,-10,0,10,100,1000)) {
      f<-as.numeric(UniUnmapParameters(nu_fixtures[[d]],d));f[ix]<-x
      nu<-as.numeric(UniMapParameters(f,d))[ix]
      expect_true(all(is.finite(nu) & nu>lo & nu<hi))
    }
    for(v in c(lo,lo-1,hi,hi+1,NA_real_,Inf)) {
      th<-nu_fixtures[[d]];th[ix[1]]<-v
      expect_error(UniUnmapParameters(th,d),"strictly between")
    }
  }
  for(v in c(.01,.5,1.5,2,8,49)) {
    th<-c(.2,-.1,1,1.3,.2,v)
    f<-MultiUnmapParameters(th,"mvt",2L)
    expect_equal(as.numeric(MultiMapParameters(f,"mvt",2L)),th,tolerance=1e-10)
    expect_equal(tail(as.numeric(f),1),log((v-LowerNu("mvt"))/(50-v)))
  }
  for(x in c(-1000,1000)) {
    f<-as.numeric(MultiUnmapParameters(c(.2,-.1,1,1.3,.2,8),"mvt",2L));f[6]<-x
    v<-as.numeric(MultiMapParameters(f,"mvt",2L))[6]
    expect_true(v>LowerNu("mvt") && v<UpperNu("mvt"))
  }
})

test_that("filter shape Jacobians agree with numerical map derivatives", {
  for(d in names(nu_fixtures)) {
    th<-nu_fixtures[[d]];ix<-nu_indices[[d]]
    th[ix]<-LowerNu(d)+1.5
    f<-as.numeric(UniUnmapParameters(th,d));k<-length(f)
    J<-numDeriv::jacobian(function(z) as.numeric(UniMapParameters(z,d)),f)
    expected<-as.numeric(crossprod(J,Score_Uni(.7,th,d)))
    out<-GAS:::GASFilter_univ(.7,f,diag(.01,k),matrix(0,k,k),1L,k,d,"Identity")
    expect_equal(as.numeric(out$mInnovations[ix,1]),expected[ix],tolerance=1e-6)
    expect_true(is.finite(out$dLLK))
  }
  th<-c(.2,-.1,1,1.3,.2,1.5);f<-as.numeric(MultiUnmapParameters(th,"mvt",2L))
  J<-numDeriv::jacobian(function(z)as.numeric(MultiMapParameters(z,"mvt",2L)),f)
  y<-c(.7,-.3);expected<-as.numeric(crossprod(J,Score_Multi(y,th,"mvt")))
  out<-GAS:::GASFilter_multi(matrix(y,2),f,diag(.01,6),matrix(0,6,6),1L,2L,6L,"mvt","Identity")
  expect_equal(as.numeric(out$mInnovations[6,1]),expected[6],tolerance=1e-6)
  expect_true(is.finite(out$dLLK))
})

test_that("mvt initialization remains finite below nu = 2", {
  q<-qt(seq(.005,.995,length.out=200),df=1)
  y<-rbind(q,q[c(seq(2,200,2),seq(1,199,2))])
  nu<-GAS:::StartingNu(c(y),"mvt")
  expect_true(nu>LowerNu("mvt") && nu<2)
  f<-GAS:::StaticStarting_Multi(y,"mvt",2L)
  expect_true(all(is.finite(f)))
  th<-as.numeric(MultiMapParameters(f,"mvt",2L))
  expect_equal(th[6],nu,tolerance=1e-7)
  expect_true(all(th[3:4]>0))
  expect_true(is.finite(ddist_Multi(c(0,0),th,"mvt",log=TRUE)))
})

test_that("Student-t dynamics match an independent R recursion", {
  th0<-c(0,.1,8);A<-diag(c(.001,.01,.04));B<-diag(c(.7,.98,.97))
  f0<-as.numeric(UniUnmapParameters(th0,"std"));kap<-as.numeric((diag(3)-B)%*%f0)
  set.seed(123)
  sim<-UniGASSim(fit=NULL,T.sim=1000,kappa=kap,A=A,B=B,Dist="std",ScalingType="Identity")
  y<-as.numeric(getObs(sim));ref<-matrix(NA_real_,3,1001);f<-f0
  lo<-.001;hi<-50
  for(t in seq_len(1001)) {
    mu<-f[1];v<-exp(f[2]);nu<-lo+(hi-lo)*plogis(f[3]);ref[,t]<-c(mu,v,nu)
    if(t<=1000) {
      z<-y[t]-mu;r<-z^2/(nu*v)
      score<-c((nu+1)*z/(nu*v+z^2),
               -1/(2*v)+(nu+1)*z^2/(2*v^2*nu*(1+r)),
               .5*digamma((nu+1)/2)-.5*digamma(nu/2)-.5/nu-
                 .5*log1p(r)+.5*(nu+1)*r/(nu*(1+r)))
      jac<-c(1,v,(nu-lo)*(hi-nu)/(hi-lo))
      f<-as.numeric(kap+A%*%(jac*score)+B%*%f)
    }
  }
  expect_equal(unname(sim@GASDyn$mTheta),ref,tolerance=1e-8)
})
