#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
int NumberParameters(std::string Dist, int iN = 1){
  int iK = 0;
  if(Dist == "norm") iK = 2;
  if(Dist == "std")  iK = 3;
  if(Dist == "ast")  iK = 5;
  if(Dist == "ast1") iK = 4;
  if(Dist == "poi")  iK = 1;
  if(Dist == "gamma")  iK = 2;
  if(Dist == "exp")  iK = 1;
  if(Dist == "beta")  iK = 2;
  if(Dist == "mvnorm") iK = 2*iN + iN*(iN-1)/2;
  if(Dist == "mvt")    iK = 2*iN + iN*(iN-1)/2 + 1;

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
double ZeroRemover(double dX){
  if(dX<1e-50){
    dX = 1e-50;
  }
  return dX;
}

arma::vec ZeroRemover_v(arma::vec vX){
  int iK = vX.size();
  for(int i=0;i<iK;i++){
    vX(i) = ZeroRemover(vX(i));
  }
  return vX;
}

arma::mat build_mR(arma::vec vR, int iN){

  arma::mat mR = eye(iN,iN);
  int i,j, iC = 0;

  for(i = 0;i<iN;i++){
    for(j = i;j<iN;j++){
      if(i!=j){
        mR(i,j) = vR(iC);
        mR(j,i) = vR(iC);
        iC     += 1;
      }
    }
  }

  return mR;
}
//[[Rcpp::export]]
arma::vec build_vR(arma::mat mR, int iN){
  arma::vec vR(iN*(iN-1)/2);
  int i,j,iC = 0;
  for(i = 0;i<iN;i++){
    for(j = i;j<iN;j++){
      if(i!=j){
        vR(iC) = mR(i,j);
        iC    += 1;
      }
    }
  }
  return vR;
}
arma::mat FillUpperTriangular(arma::vec vX,int iN){
  arma::mat Mat = zeros(iN,iN);

  int i,j;

  int iC=0;
  for(i=0;i<iN-1;i++){
    for(j=i+1;j<iN;j++){
      Mat(i,j)=vX(iC);
      iC+=1 ;
    }
  }
  return Mat;
}
arma::vec cumprod_removeLast(arma::vec mvec) {
  int nElem = mvec.n_elem;
  double cprod = mvec(0);
  arma::vec out(nElem-1);
  out(0)=cprod;
  for (int i = 1; i < nElem-1; ++i) {
    cprod *= mvec(i);
    out(i) = cprod;
  }
  return out;
}
arma::mat cumprodMat_removeLastRow(arma::mat Mat){

  int k = Mat.n_cols;

  arma::mat cprodMat=zeros(k-1,k);

  for(int i = 0;i<k;i++){
    cprodMat.col(i) = cumprod_removeLast(Mat.col(i));
  }
  return(cprodMat);
}
arma::mat Up_rbind_C(arma::mat Mat, arma::vec Vec){
  int n=Mat.n_rows;
  int k=Mat.n_cols;
  arma::mat out = zeros(n+1,k);
  for(int i =1;i<n+1;i++){
    out.row(i)=Mat.row(i-1);
  }
  out.row(0)=Vec.t();
  return(out);
}
arma::vec NaN2Zero(arma::vec vX, double To=0){
  int iN = vX.size();
  for(int i=0;i<iN;i++){
      if(R_IsNaN(vX(i)) ) vX(i) = To ;
  }
  return vX;
}
