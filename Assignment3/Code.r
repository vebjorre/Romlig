library(spatial)
library(tidyverse)
library(reshape)
library(scales)
#Plots: 600x500

### a)
seismic <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/seismic.dat")

L <- vapply(seq(1,75), rep, rep(1,75), times=75)
seismic.df <- melt(L)
seismic.df$value <- seismic$V1

p <- ggplot(seismic.df, aes(x=X1,y=X2))
fill.c <- scale_fill_gradient2(midpoint = 0.04)
theme <- theme_minimal()

p+geom_raster(aes(fill=value))+fill.c+theme+xlab("x")+ylab("y")

### b)
n <- length(seismic$V1)
nsims <- 6
likelihood0 <- dnorm(seismic$V1,0.02,0.06)
likelihood1 <- dnorm(seismic$V1,0.08,0.06)
p_di <- .5*(likelihood0 + likelihood1)
posterior <- 0.5 * likelihood1/p_di
post.var <- posterior*(1-posterior)

posterior.realizations <- replicate(nsims, as.integer(rbernoulli(n,posterior)))
posterior.plots <- vector("list",nsims)
fill.d <- scale_fill_brewer(palette=4,type="qual")
for (i in 1:nsims){
  df <- seismic.df
  df["value"] <- as.factor(posterior.realizations[,i])
  p <- ggplot(df, aes(x=X1,y=X2))
  posterior.plots[[i]] <- p + geom_raster(aes(fill=value)) + fill.d + theme + xlab("x") + ylab("y")
}
posterior.plots

b_data <-
  as_tibble(posterior.realizations) %>%
  bind_cols(dplyr::select(seismic.df, X1, X2)) %>%
  gather(sim, value, -X1, -X2)

ggplot(b_data, aes(x=X1, y=X2, fill=factor(value))) +
  facet_wrap(~sim) +
  geom_raster() +
  scale_fill_brewer(palette=4,type="qual") +
  theme_minimal() +
  labs(fill="Value") +
  xlab("x") +
  ylab("y")

seismic.df["Mean"] <- posterior
seismic.df["Variance"] <- post.var
seismic.df["MMAP"] <- as.factor(as.integer(posterior>0.5))
p <- ggplot(seismic.df, aes(x=X1,y=X2))
fill.m <- scale_fill_gradient2(midpoint=0.5)
p + geom_raster(aes(fill=Mean)) + fill.m + theme + xlab("x") + ylab("y")
fill.v <- scale_fill_gradient()
p + geom_raster(aes(fill=Variance)) + fill.v + theme + xlab("x") + ylab("y")
p + geom_raster(aes(fill=MMAP)) + fill.d + theme + xlab("x") + ylab("y")

### c)
complit <- as.matrix(read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/complit.dat"))
Lc <- vapply(seq(1,66), rep, rep(1,66), times=66)
complit.df <- melt(Lc)

### Rotate?
# rotate <- function(x) t(apply(x, 2, rev))
# complit <- complit[,c(66:1),drop = FALSE]
# complit.df$value <- as.factor(rotate(complit))

complit.df$value <- as.factor(complit)

p <- ggplot(complit.df, aes(x=X1,y=X2))
p+geom_raster(aes(fill=value))+fill.d+theme+xlab("x")+ylab("y")







