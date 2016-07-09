library(numDeriv)


IndexFinder<-function(p,q, iN){

  if(q==1){
    iC = (p-q)
  }else{
    iC = (q-1)*iN - sum(1:(q-1)) + (p-q)
  }

  return(iC)
}

IndexFinder(4,1,4)

Jacobian_element<-function(mS,mC,mPhi,j,i,h,k){

  dJ=0

  if(i == 1){
    if(h == j & k == i) {
      dJ = -mS[j,i]
    }else{
      dJ = 0
    }
  }

  if(i == 2){
    if( h == j & k == 2){
      dJ = -mS[j,2]*mS[2,1]*mS[j,1]
    }else if( h == j & k == 1 ){
      dJ = -mS[j,1]*mC[2,1] + mC[j,2]*mS[2,1]*mC[j,1]
    }else{
      dJ = 0
    }
  }

  dFooProd1 = 1.0
  dFooProd2 = 1.0
  dFooSum  = 0.0

  if(i > 2){

    if(h == j & k == i){

      for(l in 1:(i-1)) dFooProd1 = dFooProd1*mS[i,l]*mS[j,l]
      dJ = -mS[j,i] * dFooProd1

    }else if(h == j & k == 1 & i == 3){
      for(m in 2:(i-1)){
        dFooSum  = dFooSum + mC[i,m]*mC[j,m]
        dFooProd1 = dFooProd1 * mS[i,m]*mS[j,m]
      }
      dJ = -mS[j,1]*mC[i,1] + mC[j,1]*mS[i,1] * dFooSum + mC[j,i]*mC[j,1]*mS[i,1] * dFooProd1
    }else if(h == j & k == 1 & i > 3){
      for(m in 3:(i-1)){
        dFooProd1 = 1.0
        for(l in 2:(m-1)){
          dFooProd1 = dFooProd1*mS[i,l]*mS[j,l]
        }
        dFooSum   = dFooSum + mC[i,m]*mC[j,m] * dFooProd2
        dFooProd2 = dFooProd2 * mS[i,m]*mS[j,m]
      }
      dJ = -mS[j,1]*mC[i,1] + mS[j,1]*mC[j,1]*(mC[j,2]*mS[i,2] + dFooSum)  + mC[j,i]*mC[j,1]*mS[i,1] * dFooProd2
    }else{
      dJ = 0.0
    }

  }


  return(dJ)
}


Jacobian_MapR2<-function(vPhi, iN){

  iL = iN*(iN - 1)/2

  mJ = matrix(0,iL,iL)

  mPhi = matrix(0,iN,iN)
  mS   = matrix(0,iN,iN)
  mC   = matrix(0,iN,iN)

  iC = 1
  for(i in 1:iN){
    for(j in 1:iN){
      if(j>i){
      mPhi[j,i] = vPhi[iC]
      mS[j,i]   = sin(vPhi[iC])
      mC[j,i]   = cos(vPhi[iC])

      iC = iC + 1
      }
    }
  }

  iRow = 0
  iCol = 0

  for(j in 1:iN){
    for(i in 1:iN){
      for(h in 1:iN){
        for(k in 1:iN){
          if( j > i & h > k & h <= j){
            iRow = IndexFinder(j,i,iN)
            iCol = IndexFinder(h,k,iN)

            mJ[iRow, iCol] = Jacobian_element(mS,mC,mPhi,j,i,h,k)

          }
        }
      }
    }
  }

  return(mJ)

}



### Correlations

iN = 3

mR    = cor(matrix(rnorm(iN*100),100,iN))

vRho  = build_vR( mR, iN)

vPhi = UnMapR_C(vRho, iN)

MapR_C(vPhi, iN) - mR

jacobian(function(x,iN) {
  mR = MapR_C(x,iN)
  build_vR(mR,iN)
}, vPhi, iN = iN)

Jacobian_MapR2(vPhi, iN)


