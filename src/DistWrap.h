#ifndef DISTWRAP_H
#define DISTWRAP_H

double ddist_univ(double dY, arma::vec vTheta, std::string Dist, bool bLog);
double ddist_multi(arma::vec vY, arma::vec vTheta, int iN,std::string Dist, bool bLog);
double rdist_univ(arma::vec vTheta, std::string Dist);
arma::vec rdist_multi(arma::vec vTheta, int iN,std::string Dist);
#endif
