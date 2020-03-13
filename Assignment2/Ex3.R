library(dae)
library(spatial)

redwood <- as.list(read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/redwood.dat",col.names= c('x', 'y')))
redwood_df <- data.frame(x=redwood$x,y=redwood$y)
ppregion(0,1,0,1)

#MCMC-test
MCMC_test_NS <- function(L_hat, lambda_M, sigma_c, lambda_c){
  L_NS <- matrix(0,100,70)
  for (i in 1:100){
    NS <- Neuman_Scott(lambda_M, sigma_c, lambda_c)
    NS_df <- data.frame(x=NS[,1],y=NS[,2])
    L_NS[i,] <- Kfn(NS_df,1)$y
  }
  q <- matrix(0,2,70)
  for (j in 1:70){
    q[,j] <- quantile(L_NS[,j], c(0.05,0.95))
  }
  plot(L_hat, type="l", ylim=c(0,0.7), xlab="t", ylab="L(t)")
  lines(L_hat$x,q[1,], lty=2, col="red")
  lines(L_hat$x,q[2,], lty=2, col="red")
}

Neuman_Scott <- function(lambda_M,sigma_c,lambda_c){
  k <- 0 #total number of events
  D <- 1 #domain
  lambda_MD <- lambda_M*D 
  k_M <- rpois(1,lambda_MD) #Mother count
  x <- matrix(0,0,2)
  for (j in 1:k_M){ #iterate through the mothers
    x_M <- runif(2,0,1) #Sample uniformly from D to get mother location
    k_c <- rpois(1,lambda_c) #Child count
    for (i in 1:k_c){ #iterate through the children in mother[j]
      x_p <- rmvnorm(x_M,sigma_c*diag(2)) #children location with mother[j] as senter
      if (x_p[1]>0 && x_p[1]<1 && x_p[2]>0 && x_p[2]<1 ){ #Check if children is inside D
        x <- rbind(x,x_p)
        k <- k+1
      }
    }
  }
  k_D <- k-1 
  return(x)
}

#Plot redwood
plot(redwood_df,xlab="x", ylab="y",xlim=c(0,1), ylim=c(0,1))
L_redwood <- Kfn(redwood_df,1)

#Plot one realisation of Neuman Scott event RF with guestimated parameters
lambda_M <- 8 #number of clusters in redwood
sigma_c <- 0.002 #trial and error - spred in the cluster 
lambda_c <- 7.75 #number of observations divided by number of clusters  - 62/8
NS <- Neuman_Scott(lambda_M, sigma_c,lambda_c)
plot(NS, xlab="x", ylab="y", xlim=c(0,1), ylim=c(0,1))


#L-function 
NS_df <- data.frame(x=NS[,1],y=NS[,2])
L_NS <- Kfn(NS_df,1)
plot(L_NS, xlab="t", ylab="L(t)")
lines(L_NS$x,L_NS$x)

#Evaluate the parameter values by MCMC-test on the L-interaction function
MCMC_test_NS(L_redwood, lambda_M, sigma_c, lambda_c)

#Iterate our gestimate procedure to improve the fit
lambda_M_new <- 16
lambda_c_new <- 4
sigma_c_new <- 0.002
NS_new <- Neuman_Scott(lambda_M_new, sigma_c_new,lambda_c_new)
plot(NS_new, xlab="x", ylab="y")
NS_df_new <- data.frame(x=NS_new[,1],y=NS_new[,2])
L_NS_new <- Kfn(NS_df_new,1)
MCMC_test_NS(L_redwood,lambda_M_new, sigma_c_new, lambda_c_new)