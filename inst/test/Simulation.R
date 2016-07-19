library(GAS)
library(parallel)

mA = matrix(c(0.01 , 0.0 , 0.0 ,
              0.0 , 0.4 , 0.0 ,
              0.0 , 0.0 , 0.00),3,byrow = T)

mB = matrix(c(0.9 , 0.0 , 0.0 ,
              0.0 , 0.98, 0.0 ,
              0.0 , 0.0 , 0.0),3,byrow = T)

vKappa = (diag(3) - mB) %*% UnmapParameters_univ(c(0,0.1,8),"std",3)

vPar = c(vKappa, mA[1, 1], mA[2, 2], mB[1, 1], mB[2, 2]) ;
names(vPar) = c("kappa1", "kappa2", "kappa3", "a1", "a2", "b1", "b2")

### MC design
iB = 1e2 # number of resamples
iT = 1e3 # number of observations

### GAS specification
GASSpec = UniGASSpec(Dist = "std", GASPar = list(location = TRUE, scale = TRUE))

### parallelisation
cluster = makeCluster(7)
clusterEvalQ(cluster, library(GAS))
clusterExport(cluster, c("vKappa", "mA", "mB", "iT", "GASSpec"))

## Experiment

lCoef = parLapply(cluster, 1:iB, function(b){

  Sim = UniGASSim(iT, vKappa, mA, mB, Dist="std", ScalingType = "Identity")

  vY = getObs(Sim)

  Fit = UniGASFit(GASSpec, vY)

  coef(Fit)$mCoef[,1]
})

mCoef = do.call(rbind, lCoef)

### Figures

## Graph parameters

layout(matrix(1:9,3,3)
       ,heights=rep(1.5, 9)
)

par(mar = c(1.5, 1.5, 1.5, 1.5))

for(i in 1:length(vPar)){

  dens_y = density(mCoef[,i], kernel = "gaussian")$y
  dens_x = density(mCoef[,i], kernel = "gaussian")$x
  Title =  names(vPar)[i]

  plot(dens_x,dens_y,main = Title,type="l",lwd = 2,#ylim = c(0,max_y),
       #xlim=x_lim,
       xlab = "", ylab = "", tcl = -0.2, mgp = c(3, 0.1, 0))

  lines(dens_x, dens_y,  type="l", col = "purple", lwd = 2, xlab = "", ylab = "")

  abline(v=vPar[i],col="red",lty = 2)

  grid(nx = 10, ny = 10, col = "gray", lty = "dotted")

}

stopCluster(cluster)
