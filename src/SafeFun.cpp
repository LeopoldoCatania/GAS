#include <RcppArmadillo.h>
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

double lgammafn_safe(double dX){
  dX = InfRemover(dX, 1e50);
  double dOut = Rf_lgammafn(dX);
  return dOut;
}

arma::mat chol_safe(arma::mat mX){

  arma::vec vEigval = eig_sym( mX );

  double dMinEigen = min(vEigval);

  if(dMinEigen<1e-15){
    mX.diag() = mX.diag() + abs3(dMinEigen) + 1e-10;
  }

  return chol(mX) ;
}
