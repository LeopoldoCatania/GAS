#'@export
InferenceFun<-function(mHessian,vPw, iK){

  vPn  = vPw2vPn_Univ(vPw,iK)
  iK_s = length(vPn)

  out =matrix(NA,iK_s,4,dimnames=list(names(vPn),c("Estimates","Standard Errors","Test","p-value")))

  mJacob      = jacobian(vPw2vPn_Univ,vPw,iK=iK)
  mInvHessian = ginv(mHessian)
  mSandwitch  = t(mJacob)%*%mInvHessian%*%mJacob

  vTest    = vPn/vSE
  vSE      = sqrt(diag(mSandwitch))
  vPvalues = 1.0 - pnorm(abs(vTest))

  out[,"Estimates"]       = vPn
  out[,"Standard Errors"] = vSE
  out[,"Test"]            = vTest
  out[,"p-value"]         = vPvalues

  return(out)
}




