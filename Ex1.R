library(geoR)
library(tidyverse)
library(lazyeval)
library(MASS)

set.seed(0)

mu = 0
x <- seq(1,50)
sigma <- c(1,5)
kappa_m <- c(1,3)
kappa_p <- c(1,1.9)

H <- abs(x%*%t(rep(1,50))-rep(1,50) %*% t(x))

#correlation matrix
rho_matern1 <- cov.spatial(H,cov.model="matern",cov.pars=c(1,10),kappa=kappa_m[1])
rho_matern2 <- cov.spatial(H,cov.model="matern",cov.pars=c(1,10),kappa=kappa_m[2])
rho_pow.exp1 <- cov.spatial(H,cov.model="powered.exponential",cov.pars=c(1,10),kappa=kappa_p[1])
rho_pow.exp2 <- cov.spatial(H,cov.model="powered.exponential",cov.pars=c(1,10),kappa=kappa_p[2])

plot(rho_matern1[,1], type="l")
plot(rho_matern2[,1], type="l")
plot(rho_pow.exp1[,1], type="l")
plot(rho_pow.exp2[,1], type="l")

#variogram vectors
gamma_matern1 <- sigma[1]*(1-rho_matern1[,1])
gamma_matern2 <- sigma[1]*(1-rho_matern2[,1])
gamma_pow.exp1 <- sigma[1]*(1-rho_pow.exp1[,1])
gamma_pow.exp2 <- sigma[1]*(1-rho_pow.exp2[,1])

plot(gamma_matern1, type="l")
plot(gamma_matern2, type="l")
plot(gamma_pow.exp1, type="l")
plot(gamma_pow.exp2, type="l")

Sigma1_m1 <- sigma[1]*rho_matern1
Sigma1_m2 <- sigma[1]*rho_matern2
Sigma1_p1 <- sigma[1]*rho_pow.exp1
Sigma1_p2 <- sigma[1]*rho_pow.exp2
Sigma2_m1 <- sigma[2]*rho_matern1
Sigma2_m2 <- sigma[2]*rho_matern2
Sigma2_p1 <- sigma[2]*rho_pow.exp1
Sigma2_p2 <- sigma[2]*rho_pow.exp2

r1 <- matrix(0,50,4)
L <- t(chol(Sigma1_m1))
r1[,1] <- L %*% rnorm(50) + mu
plot(r1[,1], type="l", ylim=c(-2,2))
for (i in 2:4){
  r1[,i] <- L %*% rnorm(50) + mu
  lines(r1[,i])
}

r2 <- matrix(0,50,4)
L <- t(chol(Sigma1_m2))
r2[,1] <- L %*% rnorm(50) + mu
plot(r2[,1], type="l", ylim=c(-2,2))
for (i in 2:4){
  r2[,i] <- L %*% rnorm(50) + mu
  lines(r2[,i])
}

r3 <- matrix(0,50,4)
L <- t(chol(Sigma1_p1))
r3[,1] <- L %*% rnorm(50) + mu
plot(r3[,1], type="l", ylim=c(-2,2))
for (i in 2:4){
  r3[,i] <- L %*% rnorm(50) + mu
  lines(r3[,i])
}

r4 <- matrix(0,50,4)
L <- t(chol(Sigma1_p2))
r4[,1] <- L %*% rnorm(50) + mu
plot(r4[,1], type="l", ylim=c(-2,2))
for (i in 2:4){
  r4[,i] <- L %*% rnorm(50) + mu
  lines(r4[,i])
}

r5 <- matrix(0,50,4)
L <- t(chol(Sigma2_m1))
r5[,1] <- L %*% rnorm(50) + mu
plot(r5[,1], type="l", ylim=c(-2,2))
for (i in 2:4){
  r5[,i] <- L %*% rnorm(50) + mu
  lines(r5[,i])
}

r6 <- matrix(0,50,4)
L <- t(chol(Sigma2_m2))
r6[,1] <- L %*% rnorm(50) + mu
plot(r6[,1], type="l", ylim=c(-4,4))
for (i in 2:4){
  r6[,i] <- L %*% rnorm(50) + mu
  lines(r6[,i])
}

r7 <- matrix(0,50,4)
L <- t(chol(Sigma2_p1))
r7[,1] <- L %*% rnorm(50) + mu
plot(r7[,1], type="l", ylim=c(-5,5))
for (i in 2:4){
  r7[,i] <- L %*% rnorm(50) + mu
  lines(r7[,i])
}

r8 <- matrix(0,50,4)
L <- t(chol(Sigma2_p2))
r8[,1] <- L %*% rnorm(50) + mu
plot(r8[,1], type="l", ylim=c(-5,5))
for (i in 2:4){
  r8[,i] <- L %*% rnorm(50) + mu
  lines(r8[,i])
}

##d)

F <- matrix(0,3,50)
F[1,10] <- 1
F[2,25] <- 1
F[3,30] <- 1

sigma_obs <- c(0,0.25)
y1 <- F %*% r6[,1] + sigma_obs[1] * rnorm(3)
C1=F %*%Sigma2_m2%*%t(F) # + 0
rhat1=Sigma2_m2 %*% t(F) %*% solve(C1,y1)

Vhat1=Sigma2_m2-Sigma2_m2 %*% t(F)%*% solve(C1,(F %*%Sigma2_m2))
var_r1=diag(Vhat1)

rlow1=rhat1-1.28*sqrt(var_r1)
rupp1=rhat1+1.28*sqrt(var_r1)
df=data.frame(r6[,1],rhat1,rlow1,rupp1,x)

ggplot(data=df, aes(x=x, y=rhat1)) +
  geom_line(aes(x=x, y=rhat1)) +
  geom_line(aes(x=x, y=r6[,1]),linetype=2,color=2) +
  geom_line(aes(x=x, y=rlow1),linetype=3,colour=4) +
  geom_line(aes(x=x, y=rupp1),linetype=3, colour=4) 

#### New sigma_obs

y2 <- F %*% r6[,1] + sigma_obs[2] * rnorm(3)
C2=F %*%Sigma2_m2%*%t(F)+diag(sigma_obs[2],nrow=3,ncol=3)
rhat2=Sigma2_m2 %*% t(F) %*% solve(C2,y2)

Vhat2=Sigma2_m2-Sigma2_m2 %*% t(F)%*% solve(C2,(F %*%Sigma2_m2))
var_r2=diag(Vhat2)

rlow2=rhat2-1.28*sqrt(var_r2)
rupp2=rhat2+1.28*sqrt(var_r2)
df=data.frame(r6[,1],rhat2,rlow2,rupp2,x)

ggplot(data=df, aes(x=x, y=rhat2)) +
  geom_line(aes(x=x, y=rhat2)) +
  geom_line(aes(x=x, y=r6[,1]),linetype=2,color=2) +
  geom_line(aes(x=x, y=rlow2),linetype=3,colour=4) +
  geom_line(aes(x=x, y=rupp2),linetype=3, colour=4) 

##e)

mu_post1 <- rhat1
Sigma_post1 <- Vhat1
r_post1 <- matrix(0,50,100)
for (i in 1:100){
  r_post1[,i] <- mu_post1 + t(chol(Sigma_post1)) %*% rnorm(50)
}

mu_post2 <- rhat2
Sigma_post2 <- Vhat2
r_post2 <- matrix(0,50,100)
for (i in 1:100){
  r_post2[,i] <- mu_post2 + t(chol(Sigma_post2)) %*% rnorm(50)
}

mean.emp2 <- rowMeans(r_post2)
