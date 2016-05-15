#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

int NumberParameters(std::string Dist){
  int iK = 0;
  if(Dist == "norm") iK = 2;
  if(Dist == "std")  iK = 3;
  if(Dist == "ast")  iK = 5;
  if(Dist == "ast1") iK = 4;
  return iK;
}

// NumericVector sign_C(NumericVector x){
//   int N = x.size();
//   NumericVector out(N);
//
//   for(int i=0;i<N;i++){
//     if(x[i]<0) out[i] = -1;
//     if(x[i]>0) out[i] = 1;
//   }
//   return out;
// }
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
