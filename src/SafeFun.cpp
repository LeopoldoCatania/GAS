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
  arma::mat U;
  if (!mX.is_finite() || !arma::chol(U, mX))
    Rcpp::stop("The dispersion matrix must be finite and positive definite.");
  return U;
}
