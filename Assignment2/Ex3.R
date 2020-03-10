library(dae)
redwood <- as.list(read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/redwood.dat",col.names= c('x', 'y')))
redwood_df <- data.frame(x=redwood$x,y=redwood$y)

#MCMC-test
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
  plot(L_hat, type="l", ylim=c(0,0.7))
  lines(L_hat$x,q[1,], lty=2, col="red")
  lines(L_hat$x,q[2,], lty=2, col="red")
}

Neuman_Scott <- function(lambda_M,sigma_c,lambda_c){
  set.seed(222)
  k <- 0 #total number of events
  D <- 1 #domain
  lambda_MD <- lambda_M*D 
  k_M <- rpois(1,lambda_MD) #Mother count
  x <- matrix(0,0,2)
  for (j in 1:k_M){ #iterate through the mothers
    x_M <- runif(2,0,1) #Sample uniformly from D to get mother location
    k_c <- rpois(1,lambda_c) #sample to get child count
    x_j <- matrix(0,k_c,2) #Child events after torus border condition
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
plot(redwood_df)

#Plot one realisation of Neuman Scott event RF with guestimated parameters
lambda_M <- 9 #number of clusters in redwood
sigma_c <- 0.002 #try and fail - spredning
lambda_c <- 6.9 #number of observations divided by number of clusters  - 62/9
NS <- Neuman_Scott(lambda_M, sigma_c,lambda_c)
plot(NS)


#L-function 
NS_df <- data.frame(x=NS[,1],y=NS[,2])
L_NS <- Kfn(NS_df,1)
plot(L_NS)
lines(L_NS$x,L_NS$x)

#Evaluate the parameter values by MCMC-test on the L-interaction function
MCMC_test(length(NS_df$x),L_NS)

#Iterate our gestimate procedure to improve the fit
lambda_M <- 10
lambda_c <- 16
sigma_c <- 0.0015
NS <- Neuman_Scott(lambda_M, sigma_c,lambda_c)
plot(NS)
NS_df <- data.frame(x=NS[,1],y=NS[,2])
L_NS <- Kfn(NS_df,1)
MCMC_test(length(NS_df$x),L_NS)


