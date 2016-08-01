#include <RcppArmadillo.h>
#include "norm.h"
#include "std.h"
#include "sstd.h"
#include "ast.h"
#include "ald.h"
#include "poi.h"
#include "ber.h"
#include "gamma.h"
#include "exp.h"
#include "beta.h"
#include "mvnorm.h"
#include "mvt.h"


using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec Score_univ(double dY, arma::vec vTheta,std::string Dist){
  arma::vec vScore;
  if(Dist == "norm") vScore = norm_Score(dY,vTheta);
  if(Dist == "std")  vScore = std_Score(dY,vTheta);
  if(Dist == "sstd") vScore = sstd_Score(dY,vTheta);
  if(Dist == "ast")  vScore = ast_Score(dY,vTheta);
  if(Dist == "ald")  vScore = ald_Score(dY,vTheta);
  if(Dist == "ast1") vScore = ast1_Score(dY,vTheta);
  if(Dist == "poi")  vScore = poi_Score(dY,vTheta(0));
  if(Dist == "ber")  vScore = ber_Score(dY,vTheta(0));
  if(Dist == "gamma") vScore = gamma_Score(dY,vTheta);
  if(Dist == "exp")  vScore = exp_Score(dY,vTheta(0));
  if(Dist == "beta") vScore = beta_Score(dY,vTheta);

  return vScore;

}
//[[Rcpp::export]]
arma::vec Score_multi(arma::vec vY, arma::vec vTheta, int iN,std::string Dist){
  arma::vec vScore;
  if(Dist == "mvnorm") vScore = mvnorm_Score(vY, vTheta, iN);
  if(Dist == "mvt")    vScore = mvt_Score(vY, vTheta, iN);
  return vScore;
}
