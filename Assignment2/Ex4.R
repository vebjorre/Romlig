cells <- as.list(read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/cells.dat",col.names= c('x', 'y')))
cells_df <- data.frame(x=cells$x,y=cells$y)

#Plot cells
plot(cells_df,xlab="x", ylab="y")
L_cells <- Kfn(cells_df,1)

#MCMC-test
MCMC_test_Strauss <- function(L_hat, tau_0, phi_0, phi_1){
  L_S <- matrix(0,100,70)
  for (i in 1:100){
    S <- Strauss(tau_0, phi_0, phi_1)
    S_df <- data.frame(x=S[,1],y=S[,2])
    L_S[i,] <- Kfn(S_df,1)$y
  }
  q <- matrix(0,2,70)
  for (j in 1:70){
    q[,j] <- quantile(L_S[,j], c(0.05,0.95))
  }
  plot(L_hat, type="l", ylim=c(0,0.7), xlab="t", ylab="L(t)")
  lines(L_hat$x,q[1,], lty=2, col="red")
  lines(L_hat$x,q[2,], lty=2, col="red")
}

#Strauss Event RF
Strauss <- function(tau_0, phi_0, phi_1){
  x <- matrix(0,k,2)
}
  
  
  
  