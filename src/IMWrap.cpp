#include <RcppArmadillo.h>
#include "norm.h"
#include "snorm.h"
#include "std.h"
#include "ast.h"
#include "ald.h"
#include "poi.h"
#include "ber.h"
#include "gamma.h"
#include "exp.h"
#include "beta.h"



using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat IM_univ(arma::vec vTheta,std::string Dist){
  arma::mat mIM;
  if(Dist == "norm") mIM = norm_IM(vTheta);
  if(Dist == "snorm") mIM = snorm_IM(vTheta);
  if(Dist == "std")  mIM = std_IM(vTheta);
  if(Dist == "ast")  mIM = ast_IM(vTheta);
  if(Dist == "ast1") mIM = ast1_IM(vTheta);
  if(Dist == "ald") mIM = ald_IM(vTheta);
  if(Dist == "poi")  mIM = poi_IM(vTheta(0));
  if(Dist == "ber")  mIM = ber_IM(vTheta(0));
  if(Dist == "gamma") mIM = gamma_IM(vTheta);
  if(Dist == "exp")  mIM = exp_IM(vTheta(0));
  if(Dist == "beta") mIM = beta_IM(vTheta);
  return mIM;
}

