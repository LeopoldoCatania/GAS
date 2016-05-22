#ifndef DISTWRAP_H
#define DISTWRAP_H

double ddist_univ(double dY, arma::vec vTheta, std::string Dist, bool bLog);
double rdist_univ(arma::vec vTheta, std::string Dist);
#endif
