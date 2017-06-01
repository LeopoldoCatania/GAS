#ifndef SKELLAM_H
#define SKELLAM_H

double dSKELLAM(double dY, double dMu1, double dMu2, bool bLog = false);
double pSKELLAM(double dY, double dMu1, double dMu2) ;
double qSKELLAM(double dP, double dMu1, double dMu2) ;
double rSKELLAM(double dMu1, double dMu2);
arma::vec mSKELLAM(double dMu1, double dMu2);
arma::vec skellam_Score(double dY, arma::vec vTheta);
arma::mat skellam_IM(arma::vec vTheta);
#endif
