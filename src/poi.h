#ifndef POI_H
#define POI_H

double dPOI(double dY, double dMu, bool bLog=false);
double pPOI(double dY, double dMu) ;
double qPOI(double dP, double dMu);
double rPOI(double dMu);
arma::vec mPOI(double dMu);
arma::vec poi_Score(double dY, double dMu);
arma::mat poi_IM(double dMu);
#endif
