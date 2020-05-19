library(spatial)
library(tidyverse)
library(reshape)
library(scales)
library(matrixStats)
#Plots: 600x500

### a)
seismic <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/seismic.dat")

#Make grid
L <- vapply(seq(1,75), rep, rep(1,75), times=75)
seismic.df <- melt(L)
seismic.df$value <- seismic$V1

#Plot data
p <- ggplot(seismic.df, aes(x=X1,y=X2))
fill.c <- scale_fill_gradient2(midpoint = 0.04)
theme <- theme_minimal()
p+geom_raster(aes(fill=value))+fill.c+theme+xlab("x")+ylab("y")

### b)
n <- length(seismic$V1)
nsims <- 6
#p(d|l)
likelihood0 <- dnorm(seismic$V1,0.02,0.06)
likelihood1 <- dnorm(seismic$V1,0.08,0.06)
#p(d)
p_di <- .5*(likelihood0 + likelihood1)
#p(l|d) = E[l|d]
posterior <- 0.5 * likelihood1/p_di
#Var[l|d] = p(l|d)(1-p(l|d))
post.var <- posterior*(1-posterior)

#Generate realizations
posterior.realizations <- replicate(nsims, as.integer(rbernoulli(n,posterior)))

#Plot realizations
reals.df <-
  as_tibble(posterior.realizations) %>%
  bind_cols(dplyr::select(seismic.df, X1, X2)) %>%
  gather(sim, value, -X1, -X2)

ggplot(reals.df, aes(x=X1, y=X2, fill=factor(value))) +
  facet_wrap(~sim) +
  geom_raster() +
  scale_fill_brewer(palette=4,type="qual") +
  theme_minimal() +
  labs(fill="Value") +
  xlab("x") +
  ylab("y")

#Plot Mean, variance and MMAP
seismic.df["Mean"] <- posterior
seismic.df["Variance"] <- post.var
seismic.df["MMAP"] <- as.factor(as.integer(posterior>0.5))
p <- ggplot(seismic.df, aes(x=X1,y=X2))
fill.m <- scale_fill_gradient2(midpoint=0.5)
p + geom_raster(aes(fill=Mean)) + fill.m + theme + xlab("x") + ylab("y")
fill.v <- scale_fill_gradient()
p + geom_raster(aes(fill=Variance)) + fill.v + theme + xlab("x") + ylab("y")
fill.d <- scale_fill_brewer(palette=4,type="qual")
p + geom_raster(aes(fill=MMAP)) + fill.d + theme + xlab("x") + ylab("y")

### c)
complit <- as.matrix(read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/complit.dat"))
#Make grid
Lc <- vapply(seq(1,66), rep, rep(1,66), times=66)
#Plot data
complit.df <- melt(Lc)
complit.df$value <- as.factor(complit)
p <- ggplot(complit.df, aes(x=X1,y=X2))
p+geom_raster(aes(fill=value))+fill.d+theme+xlab("x")+ylab("y")

#Returns index of nearest neighbours on each side of point i
get_neighbours <- function(l, dim, i){
  n <- length(l)
  neighbours <- rep(NA, 4)
  #c(up,right,down,left)
  moves <- c(-dim,1,dim,-1)
  valid <- ifelse((i + moves) < 1 | (i + moves) > n, FALSE,TRUE)
  return(i+moves[valid])
}

#Returns index of nearest neighbours on each side of point i
#with torus border conditions
get_neighbours_bc <- function(l, dim, i){
  n <- length(l)
  neighbours <- rep(NA, 4)
  #c(up,right,down,left)
  moves <- c(-dim,1,dim,-1)
  if (i <= dim){ neighbours[1] <- n-dim + i }
  else { neighbours[1] <- i + moves[1] }
  if (i %% dim == 0){ neighbours[2] <- (i %/% dim - 1)*dim + 1 }
  else { neighbours[2] <- i + moves[2] }
  if (i > dim*(dim-1) && i %% dim != 0) { neighbours[3] <- i %% dim }
  else if (i > dim*(dim-1) && i %% dim == 0) { neighbours[3] <- dim  } 
  else { neighbours[3] <- i + moves[3] }
  if (i %% dim == 1){ neighbours[4] <- (i %/% dim + 1)*dim - 1 }
  else { neighbours[4] <- i + moves[4]}
  return (neighbours)
}

#loglikelihood of observation i
loglik_i <- function(beta,l,dim,i){
  neighb <- get_neighbours_bc(l,dim,i)
  part1 <- sum(l[i]==l[neighb])*log(beta)
  part2 <- log(beta^(sum(l[neighb]==0))+beta^(sum(l[neighb]==1)))
  return (part1-part2)
}

#total loglikelihood to be maximised wrt beta
loglik <- function(beta,l,dim){
  n <- length(l)
  ll <- rep(NA,n)
  for (i in 1:n){
    ll[i] <- loglik_i(beta,l,dim,i)
  }
  return (sum(ll))
}

#optimisation of loglikelihood wrt beta
betahat <- optimize(loglik,c(0,5),complit.df$value,66, maximum=TRUE)[[1]]

###Plotting loglikelihood
# beta <- seq(0,5,0.1)
# logliklist <- rep(NA,length(beta))
# for (i in 1:length(beta)){
#   logliklist[i] <- loglik(beta[i],complit,66)
# }
# plot(beta,logliklist)

#One gibbs "sweep"
gibbs_iteration <- function(d,l,beta,dim){
  n = length(d)
  #Sample n locations in l (with replacement)
  idx <- sample(1:n,n,replace=TRUE)
  for (i in idx){
    #p(di|li)
    likelihood0 <- dnorm(d[i], mean=0.02, sd=0.06)
    likelihood1 <- dnorm(d[i], mean=0.08, sd=0.06)
    
    #get neighbours (with torus border conditions)
    neighb <- get_neighbours_bc(l,dim,i)
    
    #p(li|l-i)
    prior0 <- beta^(sum(l[neighb]==0))
    prior1 <- beta^(sum(l[neighb]==1))
    
    #p(li|di)
    posterior0 <- prior0*likelihood0
    posterior1 <- prior1*likelihood1
    
    #Probability of sand
    p0 <- posterior0/(posterior0+posterior1)
    if (runif(1)<p0){
      l[i] <- 0
    }
    else{
      l[i] <- 1
    }
  }
  return(l)
}

dim <- 75
d <- seismic$V1
#Set initial l
l <- round(runif(n))
num_it <- 2000
MCMC_res <- matrix(0,ncol=num_it,nrow=n)
for (i in 1:num_it){
  l <- gibbs_iteration(d,l,betahat,dim)
  MCMC_res[,i] <- l
}
#Compute sand proportion
shalerate = colSums(MCMC_res)/n
sandrate = 1-shalerate

#Set burn-in period
burnin <- c(1:500)
#Choose which realizations to keep
MCMC_reals <- MCMC_res[,c(500,800,1100,1400,1700,2000)]

#Plot realizations
MCMC.df <-
  as_tibble(MCMC_reals) %>%
  bind_cols(dplyr::select(seismic.df, X1, X2)) %>%
  gather(sim, value, -X1, -X2)

ggplot(MCMC.df, aes(x=X1, y=X2, fill=factor(value))) +
  facet_wrap(~sim) +
  geom_raster() +
  scale_fill_brewer(palette=4,type="qual") +
  theme_minimal() +
  labs(x="x", y="y", fill="Value")

#Trace plot for sand proportion
ratio.df <- data.frame("Iteration"=1:num_it, sandrate)
ggplot(data=ratio.df, aes(Iteration,sandrate)) + geom_line() + ylab("Sand proportion")

#Estimate moments of MCMC sample
MCMC_mean <- rowMeans(MCMC_res[,-burnin])
MCMC_var <- rowVars(MCMC_res[,-burnin])
MCMC_MMAP <- as.integer(rowSums(MCMC_res[,-burnin])/(num_it-length(burnin)) > 0.5)
seismic.df["MCMC_m"] <- MCMC_mean
seismic.df["MCMC_v"] <- MCMC_var
seismic.df["MCMC_p"] <- as.factor(MCMC_MMAP)

#Plot moments of MCMC sample
p <- ggplot(data=seismic.df, aes(X1,X2))
fill.l <- scale_fill_gradient2(midpoint=0.5)
p + geom_raster(aes(fill=MCMC_m)) + fill.l + theme + labs(x="x", y="y", fill="Mean")
p + geom_raster(aes(fill=MCMC_v)) + fill.v + theme + labs(x="x", y="y", fill="Variance")
p + geom_raster(aes(fill=MCMC_p)) + fill.d + theme + labs(x="x", y="y", fill="MMAP")
