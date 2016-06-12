#ifndef GASINNOVATION_H
#define GASINNOVATION_H

arma::vec GASInnovation_univ(double dY, arma::vec vTheta, arma::vec vTheta_tilde, int iK, std::string Dist, std::string ScalingType);
arma::vec GASInnovation_multi(arma::vec vY, arma::vec vTheta, arma::vec vTheta_tilde, int iN, int iK,  std::string Dist, std::string ScalingType);
#endif
