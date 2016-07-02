#ifndef EXP_H
#define EXP_H

double dEXP(double dY, double dMu, bool bLog=false);
double pEXP(double dY, double dMu) ;
double qEXP(double dP, double dMu);
double rEXP(double dMu);
arma::vec mEXP(double dMu);
arma::vec exp_Score(double dY, double dMu);
arma::mat exp_IM(double dMu);
#endif
