##############################
#         UNIVARIATE         #
##############################
vPw2lPn_Uni<-function(vPw,iK){

  vA_tilde = vPw[paste("a",1:iK,sep = "")]
  vB_tilde = vPw[paste("b",1:iK,sep = "")]

  if(iK>1){
    mA = diag(c(Map_Vec(vA_tilde, LowerA(), UpperA())))
    mB = diag(c(Map_Vec(vB_tilde, LowerB(), UpperB())))
  }else{
    mA = matrix(c(Map_Vec(vA_tilde, LowerA(), UpperA())), iK , iK)
    mB = matrix(c(Map_Vec(vB_tilde, LowerB(), UpperB())), iK , iK)
  }
  lParList = list(vKappa = vPw[paste("kappa",1:iK,sep = "")], mA=mA, mB=mB)

  return(lParList)
}

vPw2vPn_Uni<-function(vPw,iK){

  vA = c(Map_Vec(vPw[paste("a",1:iK,sep = "")], LowerA(), UpperA())); names(vA) = paste("a",1:iK,sep="")
  vB = c(Map_Vec(vPw[paste("b",1:iK,sep = "")], LowerB(), UpperB())); names(vB) = paste("b",1:iK,sep="")

  vParList = c(vPw[paste("kappa",1:iK,sep = "")], vA, vB)

  vParList = vParList[!is.na(vParList)] # nas are fixed parameters. the order is preserved

  return(vParList)
}
##############################
#       MULTIVARIATE         #
##############################


vPw2lPn_Multi<-function(vPw,Dist,iK,iN, ScalarParameters){

  vA_tilde = vPw[FullNamesCoefMulti(iN,Dist, "a", ScalarParameters)]
  vB_tilde = vPw[FullNamesCoefMulti(iN,Dist, "b", ScalarParameters)]

  if(ScalarParameters){

    if(Dist == "mvt"){

      dA_tilde_nu = vA_tilde[4]
      dB_tilde_nu = vB_tilde[4]

    }else{

      dA_tilde_nu = NULL
      dB_tilde_nu = NULL

    }

    vA_tilde = c(rep(vA_tilde[1],iN), rep(vA_tilde[2],iN), rep(vA_tilde[3], iN), dA_tilde_nu)
    vB_tilde = c(rep(vB_tilde[1],iN), rep(vB_tilde[2],iN), rep(vB_tilde[3], iN), dB_tilde_nu)

  }

  mA = diag(c(Map_Vec(vA_tilde, LowerA(), UpperA())))
  mB = diag(c(Map_Vec(vB_tilde, LowerB(), UpperB())))

  lParList = list(vKappa = vPw[paste("kappa.",FullNamesMulti(iN, Dist),sep = "")], mA=mA, mB=mB)

  return(lParList)
}

vPw2vPn_Multi<-function(vPw,Dist,iK,iN, ScalarParameters){

  vA_tilde = vPw[FullNamesCoefMulti(iN,Dist, "a", ScalarParameters)]
  vB_tilde = vPw[FullNamesCoefMulti(iN,Dist, "b", ScalarParameters)]

  vA = c(Map_Vec(vA_tilde, LowerA(), UpperA())) ; names(vA) = FullNamesCoefMulti(iN,Dist, "a", ScalarParameters)
  vB = c(Map_Vec(vB_tilde, LowerB(), UpperB())) ; names(vB) = FullNamesCoefMulti(iN,Dist, "b", ScalarParameters)

  vKappa   = vPw[paste("kappa.",FullNamesMulti(iN, Dist),sep = "")]

  vParList = c(vKappa, vA, vB)

  vParList = vParList[!is.na(vParList)] # nas are fixed parameters. the order is preserved

  return(vParList)
}
