#ifndef AST_H
#define AST_H

double dAST(double dY, double dMu, double dSigma, double dAlpha, double dNu1,double dNu2 , bool bLog=false);
double rAST(double dMu, double dSigma,double dAlpha, double dNu1, double dNu2);
arma::vec ast_Score(double dY, arma::vec vTheta);
arma::vec ast1_Score(double dY, arma::vec vTheta);
arma::mat ast_IM(arma::vec vTheta);
arma::mat ast1_IM(arma::vec vTheta);
#endif
