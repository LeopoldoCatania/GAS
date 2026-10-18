#include <RcppArmadillo.h>
#include "norm.h"
#include "snorm.h"
#include "std.h"
#include "sstd.h"
#include "ast.h"
#include "ald.h"
#include "poi.h"
#include "ber.h"
#include "gamma.h"
#include "exp.h"
#include "beta.h"
#include "negbin.h"
#include "skellam.h"
#include "Utils.h"


using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat IM_univ(arma::vec vTheta,std::string Dist){
  arma::mat mIM;
  if(Dist == "norm") mIM = norm_IM(vTheta);
  if(Dist == "snorm") mIM = snorm_IM(vTheta);
  if(Dist == "std")  mIM = std_IM(vTheta);
  if(Dist == "sstd")  mIM = sstd_IM(vTheta);
  if(Dist == "ast")  mIM = ast_IM(vTheta);
  if(Dist == "ast1") mIM = ast1_IM(vTheta);
  if(Dist == "ald") mIM = ald_IM(vTheta);
  if(Dist == "poi")  mIM = poi_IM(vTheta(0));
  if(Dist == "ber")  mIM = ber_IM(vTheta(0));
  if(Dist == "gamma") mIM = gamma_IM(vTheta);
  if(Dist == "exp")  mIM = exp_IM(vTheta(0));
  if(Dist == "beta") mIM = beta_IM(vTheta);
  if(Dist == "negbin") mIM = negbin_IM(vTheta);
  if(Dist == "skellam") mIM = skellam_IM(vTheta);
  return mIM;
}

// Expected information for location, marginal scales, correlations and nu.
// Sigma is the t dispersion matrix, not its covariance matrix. These formulas
// hold for every nu > 0 and include the scale/correlation--nu cross block.
// See Lange, Little and Taylor (1989), Appendix B.
//[[Rcpp::export]]
arma::mat IM_multi(arma::vec vTheta, std::string Dist, int iN) {
  if (Dist != "mvnorm" && Dist != "mvt")
    Rcpp::stop("Multivariate Fisher information is available for mvnorm and mvt.");
  if (iN < 1) Rcpp::stop("The dimension must be positive.");
  const bool student = Dist == "mvt";
  const int nCov = iN + iN*(iN-1)/2;
  const int iK = iN + nCov + (student ? 1 : 0);
  if (vTheta.n_elem != static_cast<unsigned>(iK) || !vTheta.is_finite())
    Rcpp::stop("Invalid multivariate distribution parameters.");
  arma::vec scales = vTheta.subvec(iN, 2*iN-1);
  if (arma::any(scales <= 0)) Rcpp::stop("Scales must be positive.");
  arma::mat R = arma::eye(iN, iN);
  int k = 2*iN;
  for (int i = 0; i < iN; ++i) for (int j = i+1; j < iN; ++j)
    R(i,j) = R(j,i) = vTheta(k++);
  arma::mat Sigma = arma::diagmat(scales) * R * arma::diagmat(scales);
  arma::mat precision;
  if (!arma::inv_sympd(precision, Sigma))
    Rcpp::stop("Fisher information requires positive-definite dispersion.");
  const double nu = student ? vTheta(iK-1) : 1.0;
  if (student && nu <= 0) Rcpp::stop("Degrees of freedom must be positive.");
  const double factor = student ? (nu+iN)/(nu+iN+2.0) : 1.0;
  arma::mat information(iK, iK, arma::fill::zeros);
  information.submat(0,0,iN-1,iN-1) = factor * precision;

  // W_j = Sigma^{-1} (d Sigma / d theta_j).
  arma::cube W(iN, iN, nCov, arma::fill::zeros);
  arma::vec traces(nCov);
  for (int j = 0; j < iN; ++j) {
    arma::mat derivative(iN, iN, arma::fill::zeros);
    derivative.row(j) += Sigma.row(j)/scales(j);
    derivative.col(j) += Sigma.col(j)/scales(j);
    W.slice(j) = precision * derivative;
  }
  k = iN;
  for (int i = 0; i < iN; ++i) for (int j = i+1; j < iN; ++j) {
    arma::mat derivative(iN, iN, arma::fill::zeros);
    derivative(i,j) = derivative(j,i) = scales(i)*scales(j);
    W.slice(k++) = precision * derivative;
  }
  for (int j = 0; j < nCov; ++j) traces(j) = arma::trace(W.slice(j));
  for (int j = 0; j < nCov; ++j) {
    for (int k = j; k < nCov; ++k) {
      double entry = 0.5 * factor * arma::accu(W.slice(j) % W.slice(k).t());
      if (student) entry -= traces(j)*traces(k)/(2.0*(nu+iN+2.0));
      information(iN+j,iN+k) = information(iN+k,iN+j) = entry;
    }
    if (student) {
      double cross = -traces(j)/((nu+iN)*(nu+iN+2.0));
      information(iN+j,iK-1) = information(iK-1,iN+j) = cross;
    }
  }
  if (student)
    information(iK-1,iK-1) = 0.25*(R::trigamma(nu/2.0)-R::trigamma((nu+iN)/2.0))
      - iN*(nu+iN+4.0)/(2.0*nu*(nu+iN)*(nu+iN+2.0));
  return information;
}

