
GetFixedPar_Univ<-function(Dist,GASPar){
  FixedPar = NULL
  if(Dist=="norm"){
    if(!GASPar$location)    FixedPar = c(FixedPar, a1=0,b1=0)
    if(!GASPar$scale) FixedPar = c(FixedPar, a2=0,b2=0)
  }
  if(Dist=="std"){
    if(!GASPar$location)    FixedPar = c(FixedPar, a1=0,b1=0)
    if(!GASPar$scale) FixedPar = c(FixedPar, a2=0,b2=0)
    if(!GASPar$shape)    FixedPar = c(FixedPar, a3=0,b3=0)
  }
  if(Dist=="ast"){
    if(!GASPar$location) FixedPar = c(FixedPar, a1=0,b1=0)
    if(!GASPar$scale)    FixedPar = c(FixedPar, a2=0,b2=0)
    if(!GASPar$skewness) FixedPar = c(FixedPar, a3=0,b3=0)
    if(!GASPar$shape)    FixedPar = c(FixedPar, a4=0,b4=0)
    if(!GASPar$shape2)   FixedPar = c(FixedPar, a5=0,b5=0)
  }
  if(Dist=="ast1"){
    if(!GASPar$location) FixedPar = c(FixedPar, a1=0,b1=0)
    if(!GASPar$scale)    FixedPar = c(FixedPar, a2=0,b2=0)
    if(!GASPar$skewness) FixedPar = c(FixedPar, a3=0,b3=0)
    if(!GASPar$shape)    FixedPar = c(FixedPar, a4=0,b4=0)
  }
  return(FixedPar)
}

MultiFixedScale<-function(iN,Dist){
  FixedPar = rep(0,iN*2);
  if(Dist == "mvnorm") names(FixedPar) = c(paste("a.sigma",1:iN,sep=""),paste("b.sigma",1:iN,sep=""))
  if(Dist == "mvt")    names(FixedPar) = c(paste("a.phi",1:iN,sep=""),paste("b.phi",1:iN,sep=""))
  return(FixedPar)
}

MultiFixedLocation<-function(iN){
  FixedPar = rep(0,iN*2);
  names(FixedPar) = c(paste("a.mu",1:iN,sep=""),paste("b.mu",1:iN,sep=""))
  return(FixedPar)
}
MultiFixedCorrelation<-function(iN){
  FixedPar = rep(0,iN*(iN-1));
  vRhoNames = RhoNames(iN)
  names(FixedPar) = c(paste("a.",vRhoNames,sep=""),paste("b.",vRhoNames,sep=""))
  return(FixedPar)
}
GetFixedPar_Multi<-function(Dist,GASPar,iN){
  FixedPar = NULL
  if(Dist=="mvnorm"){
    if(!GASPar$location)    FixedPar = c(FixedPar, MultiFixedLocation(iN))
    if(!GASPar$scale)       FixedPar = c(FixedPar, MultiFixedScale(iN,"mvnorm"))
    if(!GASPar$correlation) FixedPar = c(FixedPar, MultiFixedCorrelation(iN))

  }
  if(Dist=="mvt"){
    if(!GASPar$location)    FixedPar = c(FixedPar, MultiFixedLocation(iN))
    if(!GASPar$scale)       FixedPar = c(FixedPar, MultiFixedScale(iN,"mvt"))
    if(!GASPar$correlation) FixedPar = c(FixedPar, MultiFixedCorrelation(iN))
    if(!GASPar$shape)       FixedPar = c(FixedPar, "a.nu"=0,"b.nu"=0)
  }
  return(FixedPar)
}
RemoveFixedPar<-function(vPw, FixedPar){
  if(!is.null(FixedPar)) vPw = vPw[-which(names(vPw)%in%names(FixedPar))]
  return(vPw)
}
AddFixedPar<-function(lParList){
  for(i in 1:length(lParList)){
    lParList[[i]][is.na(lParList[[i]])]=0
  }
  return(lParList)
}
