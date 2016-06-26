DistLabels<-function(){return(c("norm","std","ast","ast1","mvnorm","mvt"))}

DistName<-function(DistLabel){
  if(DistLabel=="norm") return("Gaussian")
  if(DistLabel=="std")  return("Student-t")
  if(DistLabel=="ast")  return("Asymmetric Student-t with two tail decay parameters")
  if(DistLabel=="ast1") return("Asymmetric Student-t with two one decay parameter")
  if(DistLabel=="mvnorm") return("Multivariate Gaussian")
  if(DistLabel=="mvt") return("Multivariate Student-t")
}

DistReference<-function(DistLabel){
  if(DistLabel=="norm") return("")
  if(DistLabel=="std")  return("")
  if(DistLabel=="ast")  return("Zhu, D., & Galbraith, J. W. (2010). A generalized asymmetric Student-t distribution with application to financial econometrics. Journal of Econometrics, 157(2), 297-305.")
  if(DistLabel=="ast1")  return("Zhu, D., & Galbraith, J. W. (2010). A generalized asymmetric Student-t distribution with application to financial econometrics. Journal of Econometrics, 157(2), 297-305.")
  if(DistLabel=="mvnorm") return("")
  if(DistLabel=="mvt") return("")
}

DistParameters<-function(DistLabel){
  if(DistLabel=="norm") return("location, scale")
  if(DistLabel=="std")  return("location, scale, shape")
  if(DistLabel=="ast")  return("location, scale, skewness, shape, shape2")
  if(DistLabel=="ast1") return("location, scale, skewness, shape")
  if(DistLabel=="mvnorm") return("locations, scales, correlations")
  if(DistLabel=="mvt") return("locations, scales, correlations, shape")
}

DistType<-function(DistLabel){
  if(DistLabel=="norm") return("univariate")
  if(DistLabel=="std")  return("univariate")
  if(DistLabel=="ast")  return("univariate")
  if(DistLabel=="ast1") return("univariate")
  if(DistLabel=="mvnorm") return("multivariate")
  if(DistLabel=="mvt") return("multivariate")
}

DistScalingType<-function(DistLabel){
    if(DistLabel=="norm") return("Identity, Inv, Inv.Sqrt")
    if(DistLabel=="std")  return("Identity, Inv, Inv.Sqrt")
    if(DistLabel=="ast")  return("Identity, Inv, Inv.Sqrt")
    if(DistLabel=="ast1") return("Identity, Inv, Inv.Sqrt")
    if(DistLabel=="mvnorm") return("Identity")
    if(DistLabel=="mvt") return("Identity")
  }

DistInfo<-function(DistLabel = NULL, N = 2){
  if(is.null(DistLabel)) DistLabel = DistLabels()
  for(i in 1:length(DistLabel)){
    cat("\n-------------------------------------------------------")
    cat(paste("\nName:\t",DistName(DistLabel[i]),sep=""))
    cat(paste("\nLabel:\t",DistLabel[i],sep=""))
    cat(paste("\nType:\t",DistType(DistLabel[i]),sep=""))
    cat(paste("\nParameters:\t",DistParameters(DistLabel[i]),sep=""))
    if(DistType(DistLabel[i])=="univariate") cat(paste("\nNumber of Parameters:\t",NumberParameters(DistLabel[i]),sep=""))
    if(DistType(DistLabel[i])=="multivariate") cat(paste("\nNumber of Parameters:\t",NumberParameters(DistLabel[i],N)," with N = ",N,sep=""))
    cat(paste("\nScaling Type:\t",DistScalingType(DistLabel[i])))
    cat(paste("\nReferences:\t",DistReference(DistLabel[i]),sep = ""))
    cat("\n-------------------------------------------------------")
  }
}






