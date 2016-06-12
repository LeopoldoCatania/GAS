#ifndef SCOREWRAP_H
#define SCOREWRAP_H

arma::vec Score_univ(double dY, arma::vec vTheta,std::string Dist);
arma::vec Score_multi(arma::vec vY, arma::vec vTheta, int iN,std::string Dist);
#endif
