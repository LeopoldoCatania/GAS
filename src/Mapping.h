#ifndef MAPPING_H
#define MAPPING_H

arma::vec MapParameters(arma::vec vTheta_tilde, std::string Dist, int iK);
arma::vec UnmapParameters(arma::vec vTheta, std::string Dist, int iK);
arma::vec MapParametersJacobian(arma::vec vTheta_tilde, std::string Dist, int iK);
#endif
