#ifndef STD_H
#define STD_H

arma::vec std_Score(double dY, arma::vec vTheta);
arma::mat std_IM( arma::vec vTheta);
double dSTD(double dY, double dMu, double dPhi , double dNu, bool bLog=false);
#endif
