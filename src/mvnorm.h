#ifndef MVNORM_H
#define MVNORM_H

arma::mat rmvnorm_mat(int iN, arma::vec vMu, arma::mat mSigma);
arma::vec mvnorm_Score(arma::vec vY, arma::vec vTheta, int iN);
double dmvnorm(arma::vec vY,
               arma::vec vMu,
               arma::mat mSigma,
               bool bLog = false);
double dmvnorm_ThetaParam(arma::vec vY,
                          arma::vec vTheta,
                          int iN,
                          bool bLog = false);
arma::mat rmvnorm_ThetaParam(arma::vec vTheta,int iN, int iJ) ;
arma::vec mMVNORM_mean(arma::vec vTheta, int iN);
arma::mat mMVNORM_cov(arma::vec vTheta, int iN);
#endif
