#ifndef SNORM_H
#define SNORM_H
double dSNORM(double dY, double dMu, double dSigma2, double dDelta, bool bLog = false);
double pSNORM(double dY,double dMu, double dSigma2, double dDelta);
double rSNORM(double dMu, double dSigma2, double dDelta);
double qSNORM(double dP, double dMu, double dSigma2, double dDelta,
            double lower=-150, double upper=150, int maxiter=1e4, double eps=1e-7);
arma::vec mSNORM(double dMu, double dSigma2, double dDelta);
arma::vec snorm_Score(double dY, arma::vec vTheta);
arma::mat snorm_IM(arma::vec vTheta);
#endif
