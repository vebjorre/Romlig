library(geoR)
library(akima)
library(MASS)
library(fields)

par(mar = c(5.1, 5.1, 4.1, 2.1))

#a.
n=30
kombi = expand.grid(seq(1,n, by=1), seq(1,n, by=1)) #Make n x n grid
dist = rdist(kombi, kombi) #find the distance between all the grid points

#Correlation matrix, exponential covariance function
phi <- 3
rho=matrix(data=0,nrow=(n*n),ncol=(n*n))
for (i in 1:(n*n))     
{
  for (j in 1:(n*n))
  {
    rho[i,j] = exp(-dist[i,j]/phi)
  }
}

#The covariance matrix is made from the correlation matrix.
sigma1=2 #Variance of GRF.
sigma=sigma1*rho
mu  = rep(0,n*n) #Expected value 0.

Realisation = mvrnorm(n = 1, mu, sigma) #Sample one realisation.
Result=interp(x=kombi$Var1,y=kombi$Var2,z=Realisation,xo=seq(1,30),yo=seq(1,30)) #Match element in Realisation to its position.

# Display the realisation
image.plot(Result,xlab=bquote(x),ylab=bquote(y))

#b - Variogram
x=seq(1,41)
variogram.full <- variog(coords = kombi, data=Realisation)
correct.variogram <- sigma1*(1-(exp(-x/phi)))
plot(variogram.full,type="l",xlab=bquote(tau),ylab=bquote(gamma[r](tau)))
lines(correct.variogram,type="l")

#c - repeat a) and b) three times

#d
generate.locations <- function(m){
  length <- 1
  coords <- matrix(0,m,2)
  while (length<=m){
    logical=TRUE
    x=round(runif(1,1,30))
    y=round(runif(1,1,30))
    coord=cbind(x,y)
    check=which(coords[,1]==x)
    if (length(check)>0){
      for (i in check){
        in.coord=which(coords[i,]==y)
        if (length(in.coord)>0){
          logical=FALSE
        }
      }
    }
    if (logical==TRUE){
      coords[length,]=cbind(x,y)
      length = length+1
    }
  }
  return (coords)
}

variogram <- function(m){
  coords <- generate.locations(m)
  data=Result$z[coords]
  variogram.obs <- variog(coords = coords, data=data)
  plot(variogram.obs,type="l",ylim=c(0,4))
  lines(correct.variogram,type="l",col = "red")
  lik.obs <- likfit(coords = coords, data=data, ini.cov.pars = cbind(2,3))
  lik.full <- likfit(coords = kombi, data=Realisation, ini.cov.pars = cbind(2,3))
  sigma.obs <- lik.obs$sigmasq
  sigma.full <-lik.full$sigmasq
  phi.obs <- lik.obs$phi
  phi.full <- lik.full$phi
  ml.variogram.obs <- sigma.obs*(1-(exp(-x/phi.obs)))
  ml.variogram.full <- sigma.full*(1-(exp(-x/phi.full)))
  plot(ml.variogram.obs, type="l",ylim=c(0,3))
  lines(ml.variogram.full, lty=2)
  lines(correct.variogram,type="l",col = "red")
}
variogram(36)

#e
variogram(9)
variogram(64)
variogram(100)


