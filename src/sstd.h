#ifndef SSTD_H
#define SSTD_H
double dSSTD(double dY, double dMu, double dSigma, double dXi, double dNu, bool bLog = false);
double pSSTD(double dY, double dMu, double dSigma, double dXi, double dNu);
double rSSTD(double dMu, double dSigma, double dXi, double dNu);
double qSSTD(double dP, double dMu, double dSigma, double dXi, double dNu);
arma::vec mSSTD(double dMu, double dSigma, double dXi, double dNu);
arma::vec sstd_Score(double dY, arma::vec vTheta);
arma::mat sstd_IM(arma::vec vTheta);
#endif
