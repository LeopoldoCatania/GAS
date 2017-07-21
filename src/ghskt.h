#ifndef GHSKT_H
#define GHSKT_H

double dGHSKT(double dY, double dMuBar, double dSigma, double dBetaBar, double dNu, bool bLog = false);
double pGHSKT(double dY, double dMuBar, double dSigma, double dBetaBar, double dNu);
double rGHSKT(double dMuBar, double dSigma, double dBetaBar, double dNu);
double qGHSKT(double dP, double dMuBar, double dSigma, double dBetaBar, double dNu,
              int maxiter = 1e4, double eps = 1e-7);
arma::vec mGHSKT(double dMuBar, double dSigma, double dBetaBar, double dNu);
arma::vec ghskt_Score(double dY, arma::vec vTheta);
arma::mat ghskt_IM(arma::vec vTheta);
#endif
