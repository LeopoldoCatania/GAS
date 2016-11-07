#ifndef SNORM_H
#define SNORM_H
double dSNORM(double dY, double dMu, double dSigma, double dXi, bool bLog = false);
double pSNORM(double dY,double dMu, double dSigma, double dXi);
double rSNORM(double dMu, double dSigma, double dXi);
double qSNORM(double dP, double dMu, double dSigma, double dXi);
arma::vec mSNORM(double dMu, double dSigma, double dXi);
arma::vec snorm_Score(double dY, arma::vec vTheta);
arma::mat snorm_IM(arma::vec vTheta);
#endif
