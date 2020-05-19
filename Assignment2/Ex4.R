library(spatstat)
library(spatial)
cells <- as.list(read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/cells.dat",col.names= c('x', 'y')))
cells_df <- data.frame(x=cells$x,y=cells$y)
ppregion(0,1,0,1)

#Plot cells
plot(cells_df,xlab="x", ylab="y",xlim=c(0,1), ylim=c(0,1))
L_cells <- Kfn(cells_df,1)

#Strauss Event RF
Strauss <- function(k,tau0, phi0, phi1,n=1000){
  X <- matrix(runif(k*2,0,1),k,2)
  for (i in 1:n){
    u <- sample(1:k,1)
    x_p <- runif(2,0,1)
    sum <- 0
    for (j in 1:k){
      if (j!=u){
        sum <- sum + (phi(x_p-X[j,], tau0, phi0, phi1)-phi(X[u,]-X[j,], tau0, phi0, phi1))
      }
    }
    alpha <- min(1,exp(-sum))
    if (alpha >= runif(1,0,1)){
      X[u,] <- x_p
    }
  }
  return(X)
}

#Interaction function 
phi <- function(tau_vector, tau0, phi0, phi1){
  tau <- sqrt(sum(tau_vector^2))
  if (tau0 >= tau){ #if events are located close to each other 
    return (phi0)
  }
  else{
    return(phi0*exp(-phi1*(tau-tau0)))
  }
}

#MCMC-test
MCMC_test_Strauss <- function(L_hat, tau_0, phi_0, phi_1){
  L_S <- matrix(0,100,70)
  for (i in 1:100){
    S <- Strauss(k, tau_0, phi_0, phi_1)
    S_df <- data.frame(x=S[,1],y=S[,2])
    L_S[i,] <- Kfn(S_df,1)$y
  }
  q <- matrix(0,2,70)
  for (j in 1:70){
    q[,j] <- quantile(L_S[,j], c(0.05,0.95))
  }
  plot(L_hat, type="l", xlim=c(0,0.7),ylim=c(0,0.7), xlab="t", ylab="L(t)")
  lines(L_hat$x,q[1,], lty=2, col="red")
  lines(L_hat$x,q[2,], lty=2, col="red")
}
  
#Plot one realisation of Strauss event RF with guestimated parameters
k <- 42 #observations
tau0 <- min(nndist(cells))
phi0 <- 10
phi1 <- 100
S <- Strauss(k, tau0, phi0, phi1, n=10000)
plot(S, xlab="x", ylab="y",xlim=c(0,1), ylim=c(0,1))

#L-function 
S_df <- data.frame(x=S[,1],y=S[,2])
L_S <- Kfn(S_df,1)
plot(L_S, xlab="t", ylab="L(t)")
lines(L_S$x,L_S$x)

#Evaluate the parameter values by MCMC-test on the L-interaction function
MCMC_test_Strauss(L_cells, tau_0, phi_0, phi1)

#Iterate our gestimate procedure to improve the fit
tau0_new <- 0.095
phi0_new <- 7
phi1_new <- 100 
S_new <- Strauss(k, tau0_new, phi0_new, phi1_new, n=10000)
plot(S_new, xlab="x", ylab="y",xlim=c(0,1), ylim=c(0,1))
MCMC_test_Strauss(L_cells,tau0_new, phi0_new, phi1_new)