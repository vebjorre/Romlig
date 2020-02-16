library(geoR)
library(akima)
library(fields)

par(mar = c(5.1, 5.1, 4.1, 2.1))

data <- read.table("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise1/topo.dat")

field <- interp(data$x,data$y,data$z)
image.plot(field, xlab="x", ylab="y", col=terrain.colors(30), cex.lab=1.4)
contour(field, xlab="x", ylab="y", cex.lab=1.4)

x <- seq(1,315)
L <- expand.grid(x=x,y=x)
coords <- matrix(data=cbind(data$x,data$y),ncol=2)


control <- krige.control(type.krige="OK",cov.model="powered.exponential",cov.pars=c(2500,100),kappa=1.5, trend.d="cte", trend.l="cte")
krig <- krige.conv(coords=coords, data=data$z, krige=control, locations=L)
krig.predict <- matrix(krig$predict, nrow=315, ncol=315)
krig.var <- matrix(krig$krige.var, nrow=315, ncol=315)
image.plot(krig.predict, xlab="x", ylab="y", col=terrain.colors(30), cex.lab=1.4)
contour(krig.predict, xlab="x", ylab="y", cex.lab=1.4)

image.plot(krig.var, xlab="x", ylab="y", cex.lab=1.4)


##d)
control2 <- krige.control(type.krige="OK",cov.model="powered.exponential",cov.pars=c(2500,100),kappa=1.5, trend.d="2nd", trend.l="2nd")
krig2 <- krige.conv(coords=coords, data=data$z, krige=control2, locations=L)

krig2.predict <- matrix(krig2$predict, nrow=315, ncol=315)
krig2.var <- matrix(krig2$krige.var, nrow=315, ncol=315)
image.plot(krig2.predict, xlab="x", ylab="y", col=terrain.colors(30), cex.lab=1.4)
contour(krig2.predict, xlab="x", ylab="y", cex.lab=1.4)

image.plot(krig2.var, xlab="x", ylab="y", cex.lab=1.4)

##e)
x_0 <- 850
mu_0 <- krig.predict[100,100]
sigma_0 <- sqrt(krig.var[100,100])
p850 <- pnorm(x_0, mu_0, sigma_0, lower.tail=FALSE)
p850
upper <- mu_0 + sigma_0 * qnorm(0.1,lower.tail=FALSE)
upper