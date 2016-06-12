#ifndef MAPPING_H
#define MAPPING_H

arma::vec MapParameters_univ(arma::vec vTheta_tilde, std::string Dist, int iK);
arma::vec UnmapParameters_univ(arma::vec vTheta, std::string Dist, int iK);
arma::vec MapParameters_multi(arma::vec vTheta_tilde, std::string Dist,int iN, int iK);
arma::vec MapParametersJacobian_univ(arma::vec vTheta_tilde, std::string Dist, int iK);
arma::mat MapParametersJacobian_multi(arma::vec vTheta_tilde, std::string Dist, int iN, int iK);
#endif
