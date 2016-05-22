
# ast

iT = 1e4
vTheta = c(0.2 , 1.0 , 0.6 , 5 , 7)

vY = numeric(iT)
for(i in 1:iT) vY[i] = rdist_univ(vTheta, Dist = "ast")

Fit = StaticMLFIT(vY,Dist = "ast")
Fit$vTheta
# std

iT = 5e4
vTheta = c(0.2 , 1.0 , 15)

vY = numeric(iT)
for(i in 1:iT) vY[i] = rdist_univ(vTheta, Dist = "std")

Fit = StaticMLFIT(vY,Dist = "std")
Fit$vTheta
