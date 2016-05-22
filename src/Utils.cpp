#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
int NumberParameters(std::string Dist){
  int iK = 0;
  if(Dist == "norm") iK = 2;
  if(Dist == "std")  iK = 3;
  if(Dist == "ast")  iK = 5;
  if(Dist == "ast1") iK = 4;
  return iK;
}

double sign_C(double dX){
  double dSign = 0.0;
  if(dX<0) dSign = -1;
  if(dX>0) dSign = 1;
  return dSign;
}
double maxDouble_C(double a, double b){
  double out=a;
  if(b>a){
    out=b;
  }
  return out;
}
double minDouble_C(double a, double b){
  double out=a;
  if(b<a){
    out=b;
  }
  return out;
}
double abs3(double x){
  if(x<0) x=-1.0*x;
  return x;
}
double InfRemover(double dX, double dTol=1e50){
  double INF= arma::datum::inf ;
  if(dX==INF)    dX=dTol;
  if(dX==-1*INF) dX=-dTol;
  return dX;
}
arma::vec InfRemover_vec(arma::vec vX, double dTol=1e50){
  double INF= arma::datum::inf ;
  int i,iZ = vX.size();
  for(i=0;i<iZ;i++){
    if(vX(i)==INF)    vX(i)=dTol;
    if(vX(i)==-1*INF) vX(i)=-dTol;
  }
  return vX;
}

