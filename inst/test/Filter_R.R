GASFilter_multi_R<-function(mY, vKappa, mA, mB, iT, iN, iK, Dist, ScalingType){

  vLLK = numeric(iT);
  dLLK = 0;

  #//initialise parameter
  mTheta_tilde = matrix(0,iK,iT+1);
  mTheta = matrix(0,iK,iT+1);
  mInnovations = matrix(0,iK,iT);

  #//initialise Dynamics
  vIntercept = ( diag(iK) - mB) %*% vKappa;
  mTheta_tilde[,1] = vKappa;

  mTheta[,1] = MapParameters_multi(mTheta_tilde[,1],Dist, iN, iK);

  # //initialise Likelihood
  vLLK[1] = ddist_multi(mY[,1], mTheta[,1],iN, Dist, TRUE);
  dLLK    = dLLK + vLLK[1];

  # // Dynamics
  for(i in 2:(iT+1)){
    mInnovations[,i-1] = GASInnovation_multi(mY[,i-1], mTheta[,i-1], mTheta_tilde[,i-1],iN, iK, Dist, ScalingType);
    mTheta_tilde[,i]   = vIntercept + mA %*% mInnovations[,i-1] + mB %*% mTheta_tilde[,i-1];
    mTheta[,i]         = MapParameters_multi(mTheta_tilde[,i],Dist,iN, iK);
    if(i<iT+1){
      vLLK[i] = ddist_multi(mY[,i], mTheta[,i], iN,Dist, TRUE);
      dLLK    = dLLK + vLLK[i];
    }
  }

  out = list();

  out["mTheta"]       = mTheta;
  out["mInnovations"] = mInnovations;
  out["mTheta_tilde"] = mTheta_tilde;

  out["vLLK"] = vLLK;
  out["dLLK"] = dLLK;

  return(out);
}





