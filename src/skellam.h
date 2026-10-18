#ifndef SKELLAM_H
#define SKELLAM_H

double dSKELLAM(double dY, double dMu, double dSigma2, bool bLog = false);
double pSKELLAM(double dY, double dMu, double dSigma2) ;
double qSKELLAM(double dP, double dMu, double dSigma2) ;
double rSKELLAM(double dMu, double dSigma2);
arma::vec mSKELLAM(double dMu, double dSigma2);
arma::vec skellam_Score(double dY, arma::vec vTheta);
arma::mat skellam_IM(arma::vec vTheta);
#endif
