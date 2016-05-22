getParNames<-function(Dist){
  if(Dist=="norm") out=c("mu","sigma2")
  # if(Dist=="jsu") out=c("mu","sigma","tau","df")
  # if(Dist=="ghsstd" | dist=="FSsstd") out=c("mu","delta","beta","df")
  if(Dist=="ast") out=c("mu","sigma","alpha","nu1","nu2")
  if(Dist=="ast1") out=c("mu","sigma","alpha","nu1")
  if(Dist=="std") out=c("mu","phi","df")
  # if(dist=="std2") out=c("mu","lambda","df")
  return(out)
}


LowerA<-function(){
  return(0)
}
UpperA<-function(){
  return(5)
}

LowerB<-function(){
  return(0)
}
UpperB<-function(){
  return(0.9999)
}


