#ifndef DISTWRAP_H
#define DISTWRAP_H

double ddist_univ(double dY, arma::vec vTheta, std::string Dist, bool bLog);
double ddist_multi(arma::vec vY, arma::vec vTheta, int iN,std::string Dist, bool bLog);
double rdist_univ(arma::vec vTheta, std::string Dist);
arma::vec rdist_multi(arma::vec vTheta, int iN,std::string Dist);
double pdist_univ(double dQ, arma::vec vTheta, std::string Dist);
double qdist_univ(double dP, arma::vec vTheta, std::string Dist);
arma::vec mdist_univ(arma::vec vTheta, std::string Dist);
arma::vec mdist_multi_mean(arma::vec vTheta, std::string Dist, int iN);
arma::mat mdist_multi_cov(arma::vec vTheta, std::string Dist, int iN);
// Silent paths for aggregating warnings over an entire time series.
arma::vec mdist_univ_quiet(arma::vec vTheta, std::string Dist);
arma::vec mdist_multi_mean_quiet(arma::vec vTheta, std::string Dist, int iN);
arma::mat mdist_multi_cov_quiet(arma::vec vTheta, std::string Dist, int iN);
#endif
