#ifndef NORM_H
#define NORM_H

arma::vec norm_Score(double dY, arma::vec vTheta);
arma::mat norm_IM(arma::vec vTheta);
double dNORM(double dY, double dMu, double dSigma2, bool bLog = false);
#endif
