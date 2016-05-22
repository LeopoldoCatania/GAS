vPw2lPn_Univ<-function(vPw,iK){

  vA_tilde = vPw[paste("a",1:iK,sep = "")]
  vB_tilde = vPw[paste("b",1:iK,sep = "")]

  mA = diag(c(Map_Vec(vA_tilde, LowerA(), UpperA())))
  mB = diag(c(Map_Vec(vB_tilde, LowerB(), UpperB())))

  lParList = list(vKappa = vPw[paste("kappa",1:iK,sep = "")], mA=mA, mB=mB)

  return(lParList)
}
