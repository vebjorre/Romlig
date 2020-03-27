library(MASS)
library(spatial)
library(tidyverse)
library(gridExtra)

# Load data
data <- as_tibble(read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/obsprob.txt", header=TRUE))
data["pines"] <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/obspines.txt", header=TRUE)[3]

# Plot number of observed pines and observation probabilities
p <- ggplot(data, aes(x, y))
fill <- scale_fill_gradient(low="burlywood", high="darkgreen")
theme <- theme_minimal()
p + geom_raster(aes(fill=pines)) + fill + theme
fill <- scale_fill_gradient(low="black", high="lightgrey")
p + geom_raster(aes(fill=alpha)) + fill + theme

#Estimate lambda
node_area <- 100
n <- nrow(data)
pihat <- n/sum(data$alpha)
lambdahat <- pihat / node_area

#Simulate count model from prior distribution
nsims <- 6
prior.sims <- replicate(nsims, rpois(n,pihat))

#Make plots of count models (Not included in report)
prior.plots <- vector("list", nsims)
fill <- scale_fill_gradient(low="burlywood", high="darkgreen")
for (i in 1:nsims){
  d <- data
  d["pines"] <- prior.sims[,i]
  p <- ggplot(d, aes(x,y))
  prior.plots[[i]] <- p + geom_raster(aes(fill=pines)) + fill + theme
}
prior.plots

#Function for simulating location model given count model
counts2loc <- function(sim){
  x.sim <- matrix(NA,nrow=sum(sim),ncol=2)
  nk <- 0
  for (i in 1:n){
    k <- sim[i]
    if (k>0){
      x <- runif(k,(i%%30)*10, ((i%%30)*10+10))
      y <- runif(k,((i-1)%/%30)*10, (((i-1)%/%30)*10+10))
      x.sim[(nk+1):(nk+k),] <- cbind(x,y)
      nk = nk + k
    }
  }
  return (x.sim)
}

#Simulate prior and plot location models
prior.loc <- vector("list", nsims)
for (i in 1:nsims){
  prior.loc[[i]] <- counts2loc(prior.sims[,i])
  plot(prior.loc[[i]], col="darkgreen", xlab="x", ylab="y")
}

#Posterior intensity
intensity.post <- pihat*(1 - data$alpha)

#Simulate posterior count models
### CHANGED:
post.counts <- replicate(nsims, rpois(n, intensity.post)+data$pines)

#Plot posterior count models (not included in report) 
post.plots <- vector("list",nsims)
for (i in 1:nsims){
  d <- data
  d["pines"] <- post.counts[,i]
  p <- ggplot(d, aes(x,y))
  post.plots[[i]] <- p + geom_raster(aes(fill=pines)) + fill + theme
}
post.plots

#Simulate and plot posterior location models
### DIFFERENT BECAUSE OF CHANGE
post.loc <- vector("list", nsims)
for (i in 1:nsims){
  post.loc[[i]] <- counts2loc(post.counts[,i])
  plot(post.loc[[i]], col="darkgreen", xlab="x", ylab="y")
}

#Simulate 100 count models from prior and posterior
nsims <- 100
sample.prior <- replicate(nsims, rpois(n,pihat))
mean.prior <- rowMeans(sample.prior)
### CHANGED:
sample.post <- replicate(nsims, rpois(n,intensity.post)+data$pines)
mean.post <- rowMeans(sample.post)

#Mean of 100 realizations
### DIFFERENT BECAUSE OF CHANGE
mean.plots <- vector("list", 2)
d <- data
d["pines"] <- mean.prior
p <- ggplot(d, aes(x,y))
mean.plots[[1]] <- p + geom_raster(aes(fill=pines)) + fill + theme
d["pines"] <- mean.post
p <- ggplot(d, aes(x,y))
mean.plots[[2]] <- p + geom_raster(aes(fill=pines)) + fill + theme

mean.plots
