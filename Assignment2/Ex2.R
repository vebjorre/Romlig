library(MASS)
library(spatial)
library(tidyverse)
library(gridExtra)

# Load data
data <- as_tibble(read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/obsprob.txt", header=TRUE))
data["pines"] <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise2/obspines.txt", header=TRUE)[3]
tot_area <- 300^2
node_area <- 100
n <- nrow(data)

# Plot number of observed pines and observation probabilities
p <- ggplot(data, aes(x, y))
fill <- scale_fill_gradient(low="burlywood", high="darkgreen")
theme <- theme_minimal()
p + geom_raster(aes(fill=pines)) + fill + theme
fill <- scale_fill_gradient(low="darkgreen", high="burlywood")
p + geom_raster(aes(fill=alpha)) + fill + theme

pihat <- n/sum(data$alpha)
lambdahat <- pihat / node_area

nsims <- 6
prior.sims <- replicate(nsims, rpois(n,pihat))
prior.plots <- vector("list", nsims)
fill <- scale_fill_gradient(low="burlywood", high="darkgreen")

for (i in 1:nsims){
  d <- data
  d["pines"] <- prior.sims[,i]
  p <- ggplot(d, aes(x,y))
  prior.plots[[i]] <- p + geom_raster(aes(fill=pines)) + fill + theme
}
prior.plots

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
prior.loc <- vector("list", nsims)
for (i in 1:nsims){
  prior.loc[[i]] <- counts2loc(prior.sims[,i])
  plot(prior.loc[[i]], col="darkgreen")
}

intensity.post <- pihat*(1 - data$alpha)
post.counts <- replicate(nsims, rpois(n, intensity.post))
post.plots <- vector("list",nsims)

for (i in 1:nsims){
  d <- data
  d["pines"] <- post.counts[,i]
  p <- ggplot(d, aes(x,y))
  post.plots[[i]] <- p + geom_raster(aes(fill=pines)) + fill + theme
}
post.plots

post.loc <- vector("list", nsims)
for (i in 1:nsims){
  post.loc[[i]] <- counts2loc(post.counts[,i])
  plot(post.loc[[i]], col="darkgreen")
}

nsims <- 100
sample.prior <- replicate(nsims, rpois(n,pihat))
mean.prior <- rowMeans(sample.prior)
sample.post <- replicate(nsims, rpois(n,intensity.post))
mean.post <- rowMeans(sample.post)

mean.plots <- vector("list", 2)
d <- data
d["pines"] <- mean.prior
p <- ggplot(d, aes(x,y))
mean.plots[[1]] <- p + geom_raster(aes(fill=pines)) + fill + theme
d["pines"] <- mean.post
p <- ggplot(d, aes(x,y))
mean.plots[[2]] <- p + geom_raster(aes(fill=pines)) + fill + theme

mean.plots
