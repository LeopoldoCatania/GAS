#ifndef NEGBIN_H
#define NEGBIN_H

double dNEGBIN(double dY, double dPi, double dNu, bool bLog = false);
double pNEGBIN(double dY, double dPi, double dNu)  ;
double qNEGBIN(double dP, double dPi, double dNu);
double rNEGBIN(double dPi, double dNu);
arma::vec mNEGBIN(double dPi, double dNu);
arma::vec negbin_Score(double dY, arma::vec vTheta);
arma::mat negbin_IM(arma::vec vTheta);
#endif
