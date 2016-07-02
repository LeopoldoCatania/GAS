#ifndef GAMMA_H
#define GAMMA_H

double dGAMMA(double dY, double dAlpha, double dBeta, bool bLog=false);
double pGAMMA(double dY, double dAlpha, double dBeta);
double qGAMMA(double dP, double dAlpha, double dBeta);
double rGAMMA(double dAlpha, double dBeta);
arma::vec mGAMMA(double dAlpha, double dBeta);
arma::vec gamma_Score(double dY, arma::vec vTheta);
arma::mat gamma_IM(arma::vec vTheta);

#endif
