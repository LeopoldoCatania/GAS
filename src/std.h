#ifndef STD_H
#define STD_H

arma::vec std_Score(double dY, arma::vec vTheta);
arma::mat std_IM( arma::vec vTheta);
double dSTD(double dY, double dMu, double dPhi , double dNu, bool bLog=false);
double pSTD(double dY, double dMu, double dPhi2 , double dNu);
double qSTD(double dP, double dMu, double dPhi2 , double dNu);
double rSTD(double dMu, double dPhi2 , double dNu);
arma::vec mSTD(double dMu, double dPhi2, double dNu);
#endif
