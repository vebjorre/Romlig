library(MASS)
library(spatial)

#cells <- ppinit("cells.dat") #repulsive
#redwood <- ppinit("redwood.dat") #clustered
#pines <- ppinit("pines.dat") #Stationary

cells <- as.list(read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/cells.dat",col.names= c('x', 'y')))
redwood <- as.list(read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/redwood.dat",col.names= c('x', 'y')))
pines <- ppinit("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/pines.dat")

cells_df <- data.frame(x=cells$x,y=cells$y)
redwood_df <- data.frame(x=redwood$x,y=redwood$y)
pines_df <- data.frame(x=pines$x,y=pines$y)

par(mfrow=c(1,1))
plot(cells_df, xlim=c(0,1), ylim=c(0,1))
plot(redwood_df, xlim=c(0,1), ylim=c(0,1))
plot(pines_df, xlim=c(0,1), ylim=c(0,1))

#L-function
L_cells <- Kfn(cells,1)
L_redwood <- Kfn(redwood,1)
L_pines <- Kfn(pines,1)

#Plot L-function
par(mfrow=c(1,1))
plot(L_cells, type="l", xlab="t", ylab="L(t)", ylim=c(0,.7))
lines(L_cells$x,L_cells$x,lty=2,col="red")
plot(L_redwood, type="l", xlab="t", ylab="L(t)", ylim=c(0,.7))
lines(L_redwood$x,L_redwood$x,lty=2,col="red")
plot(L_pines, type="l", xlab="t", ylab="L(t)", ylim=c(0,.7))
lines(L_pines$x,L_pines$x,lty=2,col="red")

#Plot J-function
par(mfrow=c(2,2))
plot(L_cells$y/L_cells$x)
plot(L_redwood$y/L_redwood$x)
plot(L_pines$y/L_pines$x)

#MCMC-test
n_cells <- length(cells_df$x)
n_redwood <- length(redwood_df$x)
n_pines <- length(pines_df$x)

MCMC_test <- function(n, L_hat){
  L <- matrix(0,100,70)
  for (i in 1:100){
    x <- runif(n)
    y <- runif(n)
    xy <- list(x=x,y=y)
    L[i,] <- Kfn(xy,1)$y
  }
  q <- matrix(0,2,70)
  for (j in 1:70){
    q[,j] <- quantile(L[,j], c(0.05,0.95))
  }
  plot(L_hat, type="l", xlab="t", ylab="L(t)", ylim=c(0,0.7), )
  lines(L_hat$x,q[1,], lty=2, col="red")
  lines(L_hat$x,q[2,], lty=2, col="red")
}

MCMC_test(n_cells,L_cells)
MCMC_test(n_redwood, L_redwood)
MCMC_test(n_pines, L_pines)
