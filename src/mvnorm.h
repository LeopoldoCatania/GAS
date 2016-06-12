#ifndef MVNORM_H
#define MVNORM_H

arma::mat rmvnorm_mat(int iN, arma::vec vMu, arma::mat mSigma);
arma::vec mvnorm_Score(arma::vec vTheta, arma::vec vY, int iN);
#endif
