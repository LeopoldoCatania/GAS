#ifndef MVT_H
#define MVT_H
arma::vec mvt_Score(arma::vec vY,arma::vec vTheta,  int iN);
double dmvt(arma::vec vY,
            arma::vec vMu,
            arma::mat mSigma,
            double dNu,
            bool bLog= false);
double dmvt_ThetaParam(arma::vec vY,
                       arma::vec vTheta,
                       int iN,
                       bool bLog = false);
arma::mat rmvt_mat(int iN, arma::vec vMu, arma::mat mSigma, double dNu);
arma::mat rmvt_ThetaParam(arma::vec vTheta, int iN, int iJ);
#endif
