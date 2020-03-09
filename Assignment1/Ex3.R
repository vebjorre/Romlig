library(geoR)
library(akima)
library(MASS)
library(fields)

# The rdist function is used to find the distance between all the grid points.
n=30
kombi = expand.grid(seq(1,n, by=1), seq(1,n, by=1)) #Make n x n grid
dist = rdist(kombi, kombi)

#Correlation matrix, exponential covariance function
phi <- 3
rho=matrix(data=0,nrow=(n*n),ncol=(n*n))
for (i in 1:(n*n))     
{
  for (j in 1:(n*n))
  {
    rho[i,j] = exp(-dist[i,j]/ohi)
  }
}

#The covariance matrix is made from the correlation matrix.
sigma1=2 #Variance of GRF.
sigma=sigma1*rho
mu  = rep(0,n*n) #Expected value 0.

Realisation = mvrnorm(n = 1, mu, sigma) #Sample one realisation.
Result=interp(x=kombi$Var1,y=kombi$Var2,z=Realisation,xo=seq(1,30),yo=seq(1,30)) #Match element in Realisation to its position.

# Display the realisation
image.plot(Result, asp=1)

#Variogram
variogram.full <- variog(coords = kombi, data=Realisation)
x=seq(1,41)
correct.variogram <- sigma1*(1-(exp(-x/epsilon)))
plot(variogram.full,type="l")
lines(correct.variogram,type="l")

#d
variogram <- function(m){
  location.x <- round(runif(m,1,30))
  location.y <- round(runif(m,1,30))
  coords <- cbind(location.x,location.y)
  data=diag(Result$z[location.x,location.y])
  variogram.obs <- variog(coords = coords, data=data)
  plot(variogram.obs,type="l")
  lines(correct.variogram,type="l")
  lik.obs <- likfit(coords = coords, data=data, ini.cov.pars = cbind(2,3))
  lik.full <- likfit(coords = kombi, data=Realisation, ini.cov.pars = cbind(2,3))
  sigma.obs <- lik.obs$sigmasq
  sigma.full <-lik.full$sigmasq
  phi.obs <- lik.obs$phi
  phi.full <- lik.full$phi
  ml.variogram.obs <- sigma.obs*(1-(exp(-x/phi.obs)))
  ml.variogram.full <- sigma.full*(1-(exp(-x/phi.full)))
  plot(ml.variogram.obs, type="l")
  lines(ml.variogram.full, type="l")
  lines(correct.variogram,type="l")
}
variogram(100)
#e


