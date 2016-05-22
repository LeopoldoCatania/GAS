
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
