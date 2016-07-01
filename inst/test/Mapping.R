library(numDeriv)


### Correlations

iN = 4

mR    = cor(matrix(rnorm(iN*100),100,iN))

vRho  = build_vR( mR, iN)

vPhi = UnMapR2_C(mR, iN)

MapR_C(vPhi, iN) - mR

jacobian(function(x,iN) {
  mR = MapR_C(x,iN)
  build_vR(mR,iN)
}, vPhi, iN = iN)

-sin(vPhi)

Jacobian_MapR(vPhi, iN)
Jacobian_MapR2(vPhi, iN)


IndexesFinder(1, 2)

