#ifndef UTILS_H
#define UTILS_H

double sign_C(double dX);
double maxDouble_C(double a, double b);
double minDouble_C(double a, double b);
double abs3(double x);
int NumberParameters(std::string Dist, int iN = 1);
arma::vec InfRemover_vec(arma::vec vX, double dTol=1e50);
double InfRemover(double dX, double dTol=1e50);
arma::mat build_mR(arma::vec vR, int iN);
arma::vec build_vR(arma::mat mR, int iN);
arma::mat FillUpperTriangular(arma::vec vX,int iN);
arma::mat cumprodMat_removeLastRow(arma::mat Mat);
arma::mat Up_rbind_C(arma::mat Mat, arma::vec Vec);
arma::vec ZeroRemover_v(arma::vec vX);
arma::vec NaN2Zero(arma::vec vX, double To=0);
#endif
