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

FullNamesMulti<-function(iN,Dist){
  vRhoNames = RhoNames(iN)
  if(Dist == "mvnorm")  vNames = c(paste("mu",1:iN,sep=""),paste("sigma",1:iN,sep=""),vRhoNames)
  if(Dist == "mvt")     vNames = c(paste("mu",1:iN,sep=""),paste("phi",1:iN,sep=""),vRhoNames,"nu")
  return(vNames)
}

FullNamesUni<-function(Dist){
  vNames = c("location","scale","skewness","shape","shape2")

  if(Dist == "norm")  vNames = vNames[c(1,2)]
  if(Dist == "std")   vNames = vNames[c(1,2,4)]
  if(Dist == "ast")   vNames = vNames
  if(Dist == "ast1")  vNames = vNames[c(1,2,3,4)]
  if(Dist == "poi")   vNames = vNames[1]
  if(Dist == "gamma") vNames = vNames[c(2, 4)]
  if(Dist == "exp")   vNames = vNames[1]
  if(Dist == "beta") vNames = vNames[c(2, 4)]
  return(vNames)
}

getParNamesUni<-function(object){
  Dist = getDist(object)
  parNames = FullNamesUni(Dist)
  return(parNames)
}
getParNamesMulti<-function(object){
  Dist = getDist(object)
  iN   = object@ModelInfo$iN
  parNames = FullNamesMulti(iN,Dist)
  return(parNames)
}

getParNames<-function(object){
  if(is(object,"uGASFit") | is(object,"uGASSim")) parNames = getParNamesUni(object)
  if(is(object,"mGASFit") | is(object,"mGASSim")) parNames = getParNamesMulti(object)
  return(parNames)
}


