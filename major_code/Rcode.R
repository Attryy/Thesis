library(elliptic)
library(geneplotter)
library(QRMlib)
library(spatstat)

library(lmf)



n=10^3

#================================================
# grid to define the shape set
u <- owin(c(-1,1),c(-1,1)) # square of side 2
w <- as.mask(u, eps=0.01) # 200 x 200 grid
X <- raster.x(w)
Y <- raster.y(w)
#================================================
#================================================
# radial variable for the density generator f_0(t) ~ exp(-t^2/2) (in analogy with the normal case)
g<-rgamma(n, shape=(3/2))
r<-sqrt(2*g)
#================================================

# normal density
sim.data.norm <- rmnorm(n,mean=c(0,0), varcov=matrix(c(1,0.2,1,0.2),byrow=T,ncol=2))
write.table(sim.data.norm, file="datan")

plot(sim.data.norm, pch=19, cex=0.5, xlab="", ylab="", cex.axis=1.5, xlim=c(-4,4), ylim=c(-4,4))
saveeps("plot_datan", width=8, asp=1)

#================================================
# homothetic density with gauge fn of the skew-normal limit set
a1<--1; a2<-6 #omega=0
E.sn <- owin(c(-1,1), c(-1,1), mask=((X^2 + Y^2 <= 1)&(a1*X+a2*Y>=0)|(X^2 + Y^2 +(a1*X+a2*Y)^2<= 1)&(a1*X+a2*Y<0)))

u <- runifpoint(n, win = E.sn, giveup=1000)
sim.data.Hsn <- r*cbind(u$x,u$y)
write.table(sim.data.Hsn, file="dataHsn")

plot(sim.data.Hsn, pch=19, cex=0.5, xlab="", ylab="", cex.axis=1.5, xlim=c(-4,4), ylim=c(-4,4))
saveeps("plot_dataHsn", width=8, asp=1)

#================================================
# homothetic density with triangular level sets
Tr <- owin(c(-1,1),c(-1,1), poly=list(x=c(1,-1,0),y=c(1,0,-1)))
u <- runifpoint(n, win = Tr, giveup=1000)
sim.data.HT <- r*cbind(u$x,u$y)
write.table(sim.data.HT, file="dataHT")

plot(sim.data.HT, pch=19, cex=0.5, xlab="", ylab="", cex.axis=1.5, xlim=c(-4,4), ylim=c(-4,4))
saveeps("plot_dataHT", width=8, asp=1)

#================================================
# homothetic density with shape set = limit set for the meta t density with normal margins

lam<-1; theta<-2
E.meta <- owin(c(-1,1), c(-1,1), mask=((lam+2)*pmax(abs(X),abs(Y))^theta-(X^2 + Y^2) <= lam))
u <- runifpoint(n, win = E.meta, giveup=1000)
sim.data.Hmeta <- r*cbind(u$x,u$y)
write.table(sim.data.Hmeta, file="dataHmeta")

plot(sim.data.Hmeta, pch=19, cex=0.5, xlab="", ylab="", cex.axis=1.5, xlim=c(-4,4), ylim=c(-4,4))
saveeps("plot_dataHmeta", width=8, asp=1)

