#ifndef ALD_H
#define ALD_H
double dALD(double dY, double dTheta, double dSigma, double dKappa, bool bLog = false);
double pALD(double dY, double dTheta, double dSigma, double dKappa);
double rALD(double dTheta, double dSigma, double dKappa);
double qALD(double dP, double dTheta, double dSigma, double dKappa,
            double lower=-150, double upper=150, int maxiter=1e4, double eps=1e-7);
arma::vec mALD(double dTheta, double dSigma, double dKappa);
arma::vec ald_Score(double dY, arma::vec vTheta);
arma::mat ald_IM(arma::vec vTheta);
#endif
