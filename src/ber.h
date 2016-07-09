#ifndef BER_H
#define BER_H

double dBER(double dY, double dPi, bool bLog=false);
double pBER(double dY, double dPi) ;
double qBER(double dP, double dPi);
double rBER(double dPi);
arma::vec mBER(double dPi);
arma::vec ber_Score(double dY, double dPi);
arma::mat ber_IM(double dPi);
#endif
