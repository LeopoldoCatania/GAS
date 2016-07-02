#'@export
BinTest<-function(pit,g=20,alpha=0.05,plot=F){

  h=hist(pit,nclass=g,plot=F)
  n_i=h$counts
  test=sum((n_i-mean(n_i))^2/mean(n_i))
  crit=qchisq(1-alpha,g-1)
  pvalue=1-pchisq(test,g-1)
  confidence=mean(n_i)+c(-qnorm(1-alpha)*sqrt(mean(n_i)),+qnorm(1-alpha)*sqrt(mean(n_i)))

  if(plot) {
    plot(h,col="blue",ylim=c(0,max(confidence*1.2,n_i*1.2)))
    abline(h=confidence,col="red",lwd=2,xlim=c(0,1))
  }

  out=list(test=test,crit=crit,pvalue=pvalue,hist=h, confidence=confidence)
  return(out)
}
#'@export
iidTest<-function(pit,alpha=0.05){

  N=length(pit)

  m1=as.numeric(pit-mean(pit))
  m2=as.numeric((pit-mean(pit))^2)
  m3=as.numeric((pit-mean(pit))^3)
  m4=as.numeric((pit-mean(pit))^4)

  data1=do.call(cbind,lapply(1:20, function(i) c( m1[-(1:i)]  ,rep(NA,i) )  ))
  data1=data.frame(head(data1,N-20))
  data2=do.call(cbind,lapply(1:20, function(i) c( m2[-(1:i)]  ,rep(NA,i) )  ))
  data2=data.frame(head(data1,N-20))
  data3=do.call(cbind,lapply(1:20, function(i) c( m3[-(1:i)]  ,rep(NA,i) )  ))
  data3=data.frame(head(data1,N-20))
  data4=do.call(cbind,lapply(1:20, function(i) c( m4[-(1:i)]  ,rep(NA,i) )  ))
  data4=data.frame(head(data1,N-20))

  m1=head(m1,N-20);m2=head(m2,N-20);m3=head(m3,N-20);m4=head(m4,N-20)

  fit1=lm(m1~.,data=data1)
  fit2=lm(m2~.,data=data2)
  fit3=lm(m3~.,data=data3)
  fit4=lm(m4~.,data=data4)

  test1=(N-20)*summary(fit1)$r.squared
  test2=(N-20)*summary(fit2)$r.squared
  test3=(N-20)*summary(fit3)$r.squared
  test4=(N-20)*summary(fit4)$r.squared

  crit=qchisq(1-alpha,20)

  pvalue1=1-pchisq(test1,20)
  pvalue2=1-pchisq(test2,20)
  pvalue3=1-pchisq(test3,20)
  pvalue4=1-pchisq(test4,20)


  out=list(test=c(test1=test1,test2=test2,test3=test3,test4=test4),crit=crit,pvalue=c(pvalue1,pvalue2,pvalue3,pvalue4))
  return(out)
}

PIT_test<-function(vU, iG=20, dAlpha=0.05, dBeta = 0.05, plot=F){

  if(length(vU)<100) iG = 5

  Hist = BinTest(pit = vU,g = iG, alpha = dAlpha, plot=plot)
  IID  = iidTest(pit = vU, alpha = dBeta)

  return(list(Hist = Hist, IID = IID))
}
