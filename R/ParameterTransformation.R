##############################
#         UNIVARIATE         #
##############################
vPw2lPn_Univ<-function(vPw,iK){

  vA_tilde = vPw[paste("a",1:iK,sep = "")]
  vB_tilde = vPw[paste("b",1:iK,sep = "")]

  mA = diag(c(Map_Vec(vA_tilde, LowerA(), UpperA())))
  mB = diag(c(Map_Vec(vB_tilde, LowerB(), UpperB())))

  lParList = list(vKappa = vPw[paste("kappa",1:iK,sep = "")], mA=mA, mB=mB)

  return(lParList)
}

vPw2vPn_Univ<-function(vPw,iK){

  vA = c(Map_Vec(vPw[paste("a",1:iK,sep = "")], LowerA(), UpperA())); names(vA) = paste("a",1:iK,sep="")
  vB = c(Map_Vec(vPw[paste("b",1:iK,sep = "")], LowerB(), UpperB())); names(vB) = paste("b",1:iK,sep="")

  vParList = c(vPw[paste("kappa",1:iK,sep = "")], vA, vB)

  vParList = vParList[!is.na(vParList)] # nas are fixed parameters. the order is preserved

  return(vParList)
}
##############################
#       MULTIVARIATE         #
##############################


vPw2lPn_Multi<-function(vPw,Dist,iK,iN){

  vFullNames = FullNamesMulti(iN,Dist)

  vA_tilde = vPw[paste("a.",vFullNames,sep = "")]
  vB_tilde = vPw[paste("b.",vFullNames,sep = "")]

  mA = diag(c(Map_Vec(vA_tilde, LowerA(), UpperA())))
  mB = diag(c(Map_Vec(vB_tilde, LowerB(), UpperB())))

  lParList = list(vKappa = vPw[paste("kappa.",vFullNames,sep = "")], mA=mA, mB=mB)

  return(lParList)
}

vPw2vPn_Multi<-function(vPw,Dist,iK,iN){

  vFullNames = FullNamesMulti(iN,Dist)

  vA_tilde = vPw[paste("a.",vFullNames,sep = "")]
  vB_tilde = vPw[paste("b.",vFullNames,sep = "")]

  vA = c(Map_Vec(vA_tilde, LowerA(), UpperA())) ; names(vA) = paste("a.",vFullNames,sep = "")
  vB = c(Map_Vec(vB_tilde, LowerB(), UpperB())) ; names(vB) = paste("b.",vFullNames,sep = "")

  vKappa   = vPw[paste("kappa.",vFullNames,sep = "")]

  vParList = c(vKappa, vA, vB)

  vParList = vParList[!is.na(vParList)] # nas are fixed parameters. the order is preserved

  return(vParList)
}
