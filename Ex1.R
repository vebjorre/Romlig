library(geoR)
library(tidyverse)
library(MASS)

par(mar = c(5.1, 5.1, 4.1, 2.1))

mu = 0
x <- seq(1,50)
sigma <- c(1,5)
kappa_m <- c(1,3)
kappa_p <- c(1,1.9)

H <- abs(x%*%t(rep(1,50))-rep(1,50) %*% t(x))

#correlation matrix
rho_matern1 <- cov.spatial(H,cov.model="matern",cov.pars=c(sigma[1],10),kappa=kappa_m[1])
rho_matern2 <- cov.spatial(H,cov.model="matern",cov.pars=c(sigma[1],10),kappa=kappa_m[2])
rho_pow.exp1 <- cov.spatial(H,cov.model="powered.exponential",cov.pars=c(sigma[1],10),kappa=kappa_p[1])
rho_pow.exp2 <- cov.spatial(H,cov.model="powered.exponential",cov.pars=c(sigma[1],10),kappa=kappa_p[2])

# rho_matern3 <- cov.spatial(H,cov.model="matern",cov.pars=c(sigma[2],10),kappa=kappa_m[1])
# rho_matern4 <- cov.spatial(H,cov.model="matern",cov.pars=c(sigma[2],10),kappa=kappa_m[2])
# rho_pow.exp3 <- cov.spatial(H,cov.model="powered.exponential",cov.pars=c(sigma[2],10),kappa=kappa_p[1])
# rho_pow.exp4 <- cov.spatial(H,cov.model="powered.exponential",cov.pars=c(sigma[2],10),kappa=kappa_p[2])

plot(rho_matern1[,1], type="l",xlab=bquote(tau),ylab=bquote(rho[r](tau)),cex.lab=1.4)
lines(rho_matern2[,1], lty=2)

plot(rho_pow.exp1[,1], type="l",xlab=bquote(tau),ylab=bquote(rho[r](tau)),cex.lab=1.4)
lines(rho_pow.exp2[,1], lty=2)

# plot(rho_matern3[,1], type="l")
# plot(rho_matern4[,1], type="l")
# plot(rho_pow.exp3[,1], type="l")
# plot(rho_pow.exp4[,1], type="l")


#variogram vectors
gamma_matern1 <- sigma[1]*(1-rho_matern1[,1])
gamma_matern2 <- sigma[1]*(1-rho_matern2[,1])
gamma_pow.exp1 <- sigma[1]*(1-rho_pow.exp1[,1])
gamma_pow.exp2 <- sigma[1]*(1-rho_pow.exp2[,1])

plot(gamma_matern1, type="l",xlab=bquote(tau),ylab=bquote(gamma[r](tau)),cex.lab=1.4)
lines(gamma_matern2, lty=2)
plot(gamma_pow.exp1, type="l",xlab=bquote(tau),ylab=bquote(gamma[r](tau)),cex.lab=1.4)
lines(gamma_pow.exp2, lty=2)


##b)

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
set.seed(1)
r1[,1] <- L %*% rnorm(50) + mu
plot(r1[,1], type="l", ylim=c(-3,3), xlab=bquote(italic(x)), ylab=bquote(italic(r(x))), cex.lab=1.4)
for (i in 2:4){
  set.seed(i)
  r1[,i] <- L %*% rnorm(50) + mu
  lines(r1[,i])
}

r2 <- matrix(0,50,4)
L <- t(chol(Sigma1_m2))
set.seed(1)
r2[,1] <- L %*% rnorm(50) + mu
plot(r2[,1], type="l", ylim=c(-3,3), xlab=bquote(italic(x)), ylab=bquote(italic(r(x))), cex.lab=1.4)
for (i in 2:4){
  set.seed(i)
  r2[,i] <- L %*% rnorm(50) + mu
  lines(r2[,i])
}

r3 <- matrix(0,50,4)
L <- t(chol(Sigma1_p1))
set.seed(1)
r3[,1] <- L %*% rnorm(50) + mu
plot(r3[,1], type="l", ylim=c(-3,3), xlab=bquote(italic(x)), ylab=bquote(italic(r(x))), cex.lab=1.4)
for (i in 2:4){
  set.seed(i)
  r3[,i] <- L %*% rnorm(50) + mu
  lines(r3[,i])
}

r4 <- matrix(0,50,4)
L <- t(chol(Sigma1_p2))
set.seed(1)
r4[,1] <- L %*% rnorm(50) + mu
plot(r4[,1], type="l", ylim=c(-3,3), xlab=bquote(italic(x)), ylab=bquote(italic(r(x))), cex.lab=1.4)
for (i in 2:4){
  set.seed(i)
  r4[,i] <- L %*% rnorm(50) + mu
  lines(r4[,i])
}

r5 <- matrix(0,50,4)
L <- t(chol(Sigma2_m1))
set.seed(1)
r5[,1] <- L %*% rnorm(50) + mu
plot(r5[,1], type="l", ylim=c(-5,6), xlab=bquote(italic(x)), ylab=bquote(italic(r(x))), cex.lab=1.4)
for (i in 2:4){
  set.seed(i)
  r5[,i] <- L %*% rnorm(50) + mu
  lines(r5[,i])
}

r6 <- matrix(0,50,4)
L <- t(chol(Sigma2_m2))
set.seed(1)
r6[,1] <- L %*% rnorm(50) + mu
plot(r6[,1], type="l", ylim=c(-5,6), xlab=bquote(italic(x)), ylab=bquote(italic(r(x))), cex.lab=1.4)
for (i in 2:4){
  set.seed(i)
  r6[,i] <- L %*% rnorm(50) + mu
  lines(r6[,i])
}

r7 <- matrix(0,50,4)
L <- t(chol(Sigma2_p1))
set.seed(1)
r7[,1] <- L %*% rnorm(50) + mu
plot(r7[,1], type="l", ylim=c(-5,6), xlab=bquote(italic(x)), ylab=bquote(italic(r(x))), cex.lab=1.4)
for (i in 2:4){
  set.seed(i)
  r7[,i] <- L %*% rnorm(50) + mu
  lines(r7[,i])
}

r8 <- matrix(0,50,4)
L <- t(chol(Sigma2_p2))
set.seed(1)
r8[,1] <- L %*% rnorm(50) + mu
plot(r8[,1], type="l", ylim=c(-5,6), xlab=bquote(italic(x)), ylab=bquote(italic(r(x))), cex.lab=1.4)
for (i in 2:4){
  set.seed(i)
  r8[,i] <- L %*% rnorm(50) + mu
  lines(r8[,i])
}

##d)
#Choose first realization with powwered exponential covariance function and kappa=1
#I.e r7[,1]

r <- r7[,1]

F <- matrix(0,3,50)
F[1,10] <- 1
F[2,25] <- 1
F[3,30] <- 1

sigma_obs <- c(0,0.25)
y1 <- F %*% r + sigma_obs[1] * rnorm(3)
C1=F %*%Sigma2_p1%*%t(F) # + 0
rhat1=Sigma2_p1 %*% t(F) %*% solve(C1,y1)

Vhat1=Sigma2_p1-Sigma2_p1 %*% t(F)%*% solve(C1,(F %*%Sigma2_p1))
var_r1=diag(Vhat1)

rlow1=rhat1-1.28*sqrt(var_r1)
rupp1=rhat1+1.28*sqrt(var_r1)
df=data.frame(r,rhat1,rlow1,rupp1,x)

ggplot(data=df, aes(x=x, y=rhat1)) +
  geom_line(aes(y=rhat1), size=.8) +
  geom_line(aes(y=r),linetype=2, color=2, size=.8) +
  geom_line(aes(y=rlow1),linetype=3, color=4, size=.8) +
  geom_line(aes(y=rupp1),linetype=3, color=4, size=.8) +
  xlab(bquote(italic(x))) +
  ylab(bquote(italic(r(x)))) +
  theme(axis.title=element_text(size=15))

#### New sigma_obs

y2 <- F %*% r + sigma_obs[2] * rnorm(3)
C2=F %*%Sigma2_p1%*%t(F)+diag(sigma_obs[2],nrow=3,ncol=3)
rhat2=Sigma2_p1 %*% t(F) %*% solve(C2,y2)

Vhat2=Sigma2_p1-Sigma2_p1 %*% t(F)%*% solve(C2,(F %*%Sigma2_p1))
var_r2=diag(Vhat2)

rlow2=rhat2-1.28*sqrt(var_r2)
rupp2=rhat2+1.28*sqrt(var_r2)
df=data.frame(r,rhat2,rlow2,rupp2,x)

ggplot(data=df, aes(x=x, y=rhat2)) +
  geom_line(aes(y=rhat2), size=.8) +
  geom_line(aes(y=r),linetype=2, color=2, size=.8) +
  geom_line(aes(y=rlow2),linetype=3, color=4, size=.8) +
  geom_line(aes(y=rupp2),linetype=3, color=4, size=.8) +
  xlab(bquote(italic(x))) +
  ylab(bquote(italic(r(x)))) +
  theme(axis.title=element_text(size=15))

##e)

mu_post1 <- rhat1
Sigma_post1 <- Vhat1
Sigma_post1[10,10] <- 0.000001
Sigma_post1[25,25] <- 0.000001
Sigma_post1[30,30] <- 0.000001
r_post1 <- matrix(0,100,50)
for (i in 1:100){
  r_post1[i,] <- mu_post1 + t(chol(Sigma_post1)) %*% rnorm(50)
}

mean.emp1 <- colMeans(r_post1)
var.emp1 <- diag(cov(r_post1))
mean.emp1_low <- mean.emp1 - 1.28 * sqrt(var.emp1)
mean.emp1_upp <- mean.emp1 + 1.28 * sqrt(var.emp1)
real <- sort(rep(seq(1,100),50))
post1.vec <- as.vector(t(r_post1))

post1.df <- data.frame(real,x=rep(x,100),r,post1.vec,mean.emp1,mean.emp1_low,mean.emp1_upp)
ggplot(data=post1.df, aes(x=x, y=mean.emp1)) +
  geom_path(aes(y=post1.vec, group=real), colour=8) +
  geom_path(aes(y=mean.emp1, group=real), colour=1, size=.8) +
  geom_path(aes(y=r, group=real), colour=2, linetype=2, size=.8) +
  geom_path(aes(y=mean.emp1_low, group=real), colour=4, linetype=3, size=.8) +
  geom_path(aes(y=mean.emp1_upp, group=real), colour=4, linetype=3, size=.8) +
  xlab(bquote(italic(x))) +
  ylab(bquote(italic(r(x)))) +
  theme(axis.title=element_text(size=15))

mu_post2 <- rhat2
Sigma_post2 <- Vhat2
r_post2 <- matrix(0,100,50)
for (i in 1:100){
  r_post2[i,] <- mu_post2 + t(chol(Sigma_post2)) %*% rnorm(50)
}

mean.emp2 <- colMeans(r_post2)
var.emp2 <- diag(cov(r_post2))
mean.emp2_low <- mean.emp2 - 1.28 * sqrt(var.emp2)
mean.emp2_upp <- mean.emp2 + 1.28 * sqrt(var.emp2)
post2.vec <- as.vector(t(r_post2))

post2.df <- data.frame(real,x=rep(x,100),r,post2.vec,mean.emp2,mean.emp2_low,mean.emp2_upp)
ggplot(data=post2.df, aes(x=x, y=mean.emp2)) +
  geom_path(aes(x=x, y=post2.vec, group=real), color=8) +
  geom_path(aes(x=x, y=mean.emp2, group=real), color=1, size=.8) +
  geom_path(aes(x=x, y=r, group=real), linetype=2, color=2, size=.8) +
  geom_path(aes(x=x, y=mean.emp2_low, group=real),linetype=3,colour=4, size=.8) +
  geom_path(aes(x=x, y=mean.emp2_upp, group=real),linetype=3,colour=4, size=.8) +
  xlab(bquote(italic(x))) +
  ylab(bquote(italic(r(x)))) +
  theme(axis.title=element_text(size=15)) 

integral1 <- 0
for (i in 1:50){
  for (j in 1:100){
    if (r_post2[j,i] > 2){
      integral1 = integral1 + 1
    }
  }
}
integral1 = integral1/100

integral2 <- 0
for (i in 1:50){
    if (rhat1[i] > 2){
      integral2 = integral2 + 1
    }
}

integral3 <- 0
for (i in 1:50){
  if (r[i] > 2){
    integral3 = integral3 + 1
  }
}