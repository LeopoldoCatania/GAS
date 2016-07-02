#ifndef BETA_H
#define BETA_H

double dBETA(double dY, double dAlpha,  double dBeta, bool bLog=false);
double pBETA(double dY, double dAlpha,  double dBeta) ;
double qBETA(double dP, double dAlpha,  double dBeta);
double rBETA(double dAlpha,  double dBeta);
arma::vec mBETA(double dAlpha,  double dBeta);
arma::vec beta_Score(double dY, arma::vec vTheta);
arma::mat beta_IM(arma::vec vTheta);
#endif
