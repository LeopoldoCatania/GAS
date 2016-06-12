#include <RcppArmadillo.h>
#include "GASInnovation.h"
#include "DistWrap.h"
#include "Mapping.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List GASFilter_univ(arma::vec vY, arma::vec vKappa, arma::mat mA, arma::mat mB, int iT, int iK, std::string Dist, std::string ScalingType){

  int i;
  arma::vec vLLK(iT);
  double dLLK = 0;

  //initialise parameter
  arma::mat mTheta_tilde(iK,iT+1);
  arma::mat mTheta(iK,iT+1);
  arma::mat mInnovations(iK,iT);

  //initialise Dynamics
  arma::vec vIntercept = ( eye(iK,iK) - mB) * vKappa;
  mTheta_tilde.col(0) = vKappa;

  mTheta.col(0) = MapParameters_univ(mTheta_tilde.col(0),Dist, iK);

  //initialise Likelihood
  vLLK(0) = ddist_univ(vY(0), mTheta.col(0), Dist, true);
  dLLK   += vLLK(0);

  // Dynamics
  for(i=1;i<iT+1;i++){
    mInnovations.col(i-1) = GASInnovation_univ(vY(i-1), mTheta.col(i-1), mTheta_tilde.col(i-1), iK, Dist, ScalingType);
    mTheta_tilde.col(i)   = vIntercept + mA * mInnovations.col(i-1) + mB *  mTheta_tilde.col(i-1);
    mTheta.col(i)         = MapParameters_univ(mTheta_tilde.col(i),Dist, iK);
    if(i<iT){
      vLLK(i) = ddist_univ(vY(i), mTheta.col(i), Dist, true);
      dLLK   += vLLK(i);
    }
  }

  List out;

  out["mTheta"]       = mTheta;
  out["mInnovations"] = mInnovations;
  out["mTheta_tilde"] = mTheta_tilde;

  out["vLLK"] = vLLK;
  out["dLLK"] = dLLK;

  return out;
}

//[[Rcpp::export]]
List GASFilter_multi(arma::mat mY, arma::vec vKappa, arma::mat mA, arma::mat mB, int iT, int iN, int iK, std::string Dist, std::string ScalingType){

  int i;
  arma::vec vLLK(iT);
  double dLLK = 0;

  //initialise parameter
  arma::mat mTheta_tilde(iK,iT+1);
  arma::mat mTheta(iK,iT+1);
  arma::mat mInnovations(iK,iT);

  //initialise Dynamics
  arma::vec vIntercept = ( eye(iK,iK) - mB) * vKappa;
  mTheta_tilde.col(0) = vKappa;

  mTheta.col(0) = MapParameters_multi(mTheta_tilde.col(0),Dist, iN, iK);

  //initialise Likelihood
  vLLK(0) = ddist_multi(mY.col(0), mTheta.col(0),iN, Dist, true);
  dLLK   += vLLK(0);

  // Dynamics
  for(i=1;i<iT+1;i++){
    mInnovations.col(i-1) = GASInnovation_multi(mY.col(i-1), mTheta.col(i-1), mTheta_tilde.col(i-1),iN, iK, Dist, ScalingType);
    mTheta_tilde.col(i)   = vIntercept + mA * mInnovations.col(i-1) + mB *  mTheta_tilde.col(i-1);
    mTheta.col(i)         = MapParameters_multi(mTheta_tilde.col(i),Dist,iN, iK);
    if(i<iT){
      vLLK(i) = ddist_multi(mY.col(i), mTheta.col(i), iN,Dist, true);
      dLLK   += vLLK(i);
    }
  }

  List out;

  out["mTheta"]       = mTheta;
  out["mInnovations"] = mInnovations;
  out["mTheta_tilde"] = mTheta_tilde;

  out["vLLK"] = vLLK;
  out["dLLK"] = dLLK;

  return out;
}

List FFBS(arma::mat allprobs, arma::vec delta, arma::mat mGamma, int iJ, int iT){

  arma::mat lalpha=zeros(iJ,iT);
  arma::mat lbeta=zeros(iJ,iT);

  arma::vec foo(iJ);
  double sumfoo,lscale;
  int i;

  foo    = delta % allprobs.row(0).t();
  sumfoo = sum(foo);
  lscale = log(sumfoo);
  foo    = foo/sumfoo ;

  lalpha.col(0) = log(foo)+lscale;
  for(i=1;i<iT;i++){
    foo           = (foo.t() * mGamma).t() % allprobs.row(i).t();
    sumfoo        = sum(foo);
    lscale        = lscale+log(sumfoo);
    foo           = foo/sumfoo;
    lalpha.col(i) = log(foo)+lscale;
  }
  for(i=0;i<iJ;i++) {
    foo(i)=1.0/iJ;
  }
  lscale = log(iJ);
  for(i=iT-2;i>=0;i--){
    foo          = mGamma * (allprobs.row(i+1).t() % foo);
    lbeta.col(i) = log(foo)+lscale;
    sumfoo       = sum(foo);
    foo          = foo/sumfoo;
    lscale       = lscale+log(sumfoo);
  }

  List FS;
  FS["lalpha"]=lalpha;
  FS["lbeta"]=lbeta;

  return FS;
}

// List HMMlalphabeta(arma::vec vY, arma::mat mGamma, arma::vec vMu, arma::vec vSigma2, int T, int K){
//
//   arma::vec vDelta=getDelta( mGamma, K);
//
//   arma::mat allprobs = GaussianLk(vY, vMu, vSigma2, K, T, 0);
//
//   List FB=FFBS(allprobs, vDelta, mGamma, K, T);
//
//   FB["allprobs"]=allprobs;
//
//   return FB;
// }

