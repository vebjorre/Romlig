library(akima)
library(fields)

data <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise1/topo.dat")

field <- interp(data$x,data$y,data$z)
image.plot(field)
contour(field)
image(field)

x <- seq(1,315)
L <- expand.grid(x=x,y=x)
coords <- matrix(data=cbind(data$x,data$y),ncol=2)


control <- krige.control(type.krige="OK",cov.model="powered.exponential",cov.pars=c(2500,100),kappa=1.5)
krig <- krige.conv(coords=coords, data=data$z, krige=control, locations=L)
krig.predict <- matrix(krig$predict, nrow=315, ncol=315)
krig.var <- matrix(krig$krige.var, nrow=315, ncol=315)
image(krig.predict)
image.plot(krig.predict)
contour(krig.predict)

image(krig.var)
image.plot(krig.var)


##d)
H <- abs(x%*%t(rep(1,315))-rep(1,315) %*% t(x))
Sigma.r <- exp(-(0.01*H)^1.5) #ikke gange med sigma?
Sigma.d <- Sigma.r[data$x,data$y] #??

g <- function(coords){
  return(cbind(rep(1,52),coords[,1],coords[,2],coords[,1]*coords[,2],coords[,1]^2,coords[,2]^2))
}

G <- g(coords)

control2 <- krige.control(type.krige="OK",cov.model="powered.exponential",cov.pars=c(2500,100),kappa=1.5, trend.d="2nd", trend.l="2nd")
krig2 <- krige.conv(coords=coords, data=data$z, krige=control2, locations=L)

krig2.predict <- matrix(krig2$predict, nrow=315, ncol=315)
krig2.var <- matrix(krig2$krige.var, nrow=315, ncol=315)
image(krig2.predict)
image.plot(krig2.predict)
contour(krig2.predict)

image(krig2.var)
image.plot(krig2.var)

##e)
x_0 <- 850
mu_0 <- krig.predict[100,100]
sigma_0 <- sqrt(krig.var[100,100])
pnorm(x_0, mu_0, sigma_0, lower.tail=FALSE)
