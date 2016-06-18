#'@export
InferenceFun_Univ<-function(mHessian,vPw, iK){

  vPn  = vPw2vPn_Univ(vPw,iK)
  iK_s = length(vPn)

  out =matrix(NA,iK_s,4,dimnames=list(names(vPn),c("Estimates","Standard Errors","Test","p-value")))

  mJacob      = jacobian(vPw2vPn_Univ,vPw,iK=iK)
  mInvHessian = ginv(mHessian)
  mSandwitch  = t(mJacob)%*%mInvHessian%*%mJacob

  vSE      = sqrt(diag(mSandwitch))
  vTest    = vPn/vSE
  vPvalues = 1.0 - pnorm(abs(vTest))

  out[,"Estimates"]       = vPn
  out[,"Standard Errors"] = vSE
  out[,"Test"]            = vTest
  out[,"p-value"]         = vPvalues

  return(out)
}

#'@export
InferenceFun_Multi<-function(mHessian,Dist,vPw, iK, iN){

  vPn  = vPw2vPn_Multi(vPw,Dist,iK,iN)
  iK_s = length(vPn)

  out =matrix(NA,iK_s,4,dimnames=list(names(vPn),c("Estimates","Standard Errors","Test","p-value")))

  mJacob      = jacobian(vPw2vPn_Multi,vPw,iK=iK,iN=iN,Dist=Dist)
  mInvHessian = ginv(mHessian)
  mSandwitch  = t(mJacob)%*%mInvHessian%*%mJacob

  vSE      = sqrt(diag(mSandwitch))
  vTest    = vPn/vSE
  vPvalues = 1.0 - pnorm(abs(vTest))

  out[,"Estimates"]       = vPn
  out[,"Standard Errors"] = vSE
  out[,"Test"]            = vTest
  out[,"p-value"]         = vPvalues

  return(out)
}
