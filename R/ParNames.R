mvnormParNames<-function(iN){
  foo = c(paste("mu",1:iN,sep=""),paste("sigma",1:iN,sep=""))
  baz = RhoNames(iN)
  foo = c(foo,baz)
  return(foo)
}

mvtParNames<-function(iN){
  foo = c(paste("mu",1:iN,sep=""),paste("phi",1:iN,sep=""))
  baz = RhoNames(iN)
  foo = c(foo,baz,"nu")
  return(foo)
}

RhoNames<-function(iN){
  baz = numeric(iN*(iN-1)/2)
  iC  = 1
  for(i in 1:iN){
    for(j in i:iN){
      if(i!=j){
        baz[iC] = paste("rho",i,j,sep="")
        iC = iC + 1
      }
    }
  }
  return(baz)
}

