#include <RcppArmadillo.h>
#include "Mapping.h"
#include "IMWrap.h"
#include "Utils.h"
#include "ScoreWrap.h"

using namespace Rcpp;
using namespace arma;

arma::vec GASInnovation_univ(double dY, arma::vec vTheta, arma::vec vTheta_tilde, int iK, std::string Dist, std::string ScalingType){

  arma::vec vS_tilde(iK);

  arma::vec vScore       = Score_univ(dY,vTheta,Dist);
  arma::mat mJ           = MapParametersJacobian_univ(vTheta_tilde, Dist, iK);
  arma::vec Score_tilde  = mJ.t() * vScore ;

  if (ScalingType == "Identity") {

    vS_tilde = Score_tilde;

  } else if (ScalingType == "Inv" || ScalingType == "InvSqrt") {

    // Fisher information in the original parameterization.
    arma::mat mIM = IM_univ(vTheta, Dist);

    // Fisher information in the unconstrained state parameterization.
    arma::mat mIM_tilde = mJ.t() * mIM * mJ;

    // Upper-triangular factor U satisfying I_f = U.t() * U.
    arma::mat U;

    if (!mIM_tilde.is_finite() || !arma::chol(U, mIM_tilde)) {
      Rcpp::stop(
        "Score scaling requires a finite, positive-definite "
        "transformed Fisher information matrix."
      );
    }

    // Whiten the transformed score: z = inverse(U.t()) * s_f.
    arma::vec z = arma::solve(
      arma::trimatl(U.t()),
      Score_tilde,
      arma::solve_opts::no_approx
    );

    if (ScalingType == "Inv") {

      // inverse(I_f) * s_f = inverse(U) * inverse(U.t()) * s_f.
      vS_tilde = arma::solve(
        arma::trimatu(U),
        z,
        arma::solve_opts::no_approx
      );

    } else {

      // Inverse-square-root scaling using Cholesky whitening.
      vS_tilde = z;

    }

  } else {

    Rcpp::stop("Unknown ScalingType.");

  }

  //
  vS_tilde = Thresholding_vec(vS_tilde, 1e5);
  vS_tilde = NaN2Zero(vS_tilde);
  vS_tilde = InfRemover_vec(vS_tilde, 1e5);

  return vS_tilde;
}

arma::vec GASInnovation_multi(arma::vec vY, arma::vec vTheta, arma::vec vTheta_tilde, int iN, int iK,  std::string Dist, std::string ScalingType){

  arma::vec vScore   = Score_multi(vY, vTheta, iN, Dist);
  arma::mat mJ       = MapParametersJacobian_multi(vTheta_tilde, Dist, iN, iK);
  arma::vec vS_tilde = mJ.t() * vScore ;

  vS_tilde = Thresholding_vec(vS_tilde, 1e5);
  vS_tilde = NaN2Zero(vS_tilde);
  vS_tilde = InfRemover_vec(vS_tilde, 1e5);

  return vS_tilde;
}
