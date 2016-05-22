#ifndef UTILS_H
#define UTILS_H

double sign_C(double dX);
double maxDouble_C(double a, double b);
double minDouble_C(double a, double b);
double abs3(double x);
int NumberParameters(std::string Dist);
arma::vec InfRemover_vec(arma::vec vX, double dTol=1e50);
double InfRemover(double dX, double dTol=1e50);
#endif
