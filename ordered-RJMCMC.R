library(elliptic)
library(geneplotter)
library(QRMlib)
library(spatstat)
library(lmf)
library(mvtnorm)
library(rootSolve)
library(forcats)
library(dplyr)


# gF 
# density generator: weibull distribution
# Should note the different parameterization of the density generators for D-class and homothetic density 
gF = function(x,shape, scale)
{exp(-(x^shape)/(scale^shape))}

gF_t = function(x,p=2,v=4)
{ (1+x/v)^(-0.5*(v+p))}

#change gF according to the estimates from Allen's algorithm estimation.R
#How to determine the constant for the denominator? 
#Note the estimated c is not the actual value since we used the scaled radial variables
#Or does this denominator really influence our final results?
#gF = function(x)
#{ exp(-x^(beta_hat/2)/??)}


# KpriorF
# prior of K (number of knots): Poisson with lambda

KpriorF = function(K, lambda = 5 * m)
{	dpois(K, lambda)}


# BF
# Bernoulli polynomial of degree 4 (i.e., l = 2)

BF = function(x)
{	(x * (1 - x))^2 - 1/30 	}

# RF
# reproducing kernel for l=2

RF = function(theta)
{	-(2 * pi)^3 / 24 * BF(theta/(2*pi))}


# cF
# function to calculate c from K basis 

cF = function(K, omega, v)
{
  omega1 = matrix(-omega, K, K)
  omega2 = matrix(omega, K, K, byrow = T)
  omegaMat = omega1 + omega2
  omegaMat[omegaMat<0] = omegaMat[omegaMat<0] + 2*pi
  Rfull = matrix(1, K+1, K+1)
  Rfull[K+1, 1] = 0
  Rfull[1:K, 2:(K+1)] = RF(omegaMat)
  vfull = c(v, 0)
  c = solve(Rfull, vfull)
  c
}


# polarF
# change data to polar coordinate with respect to center B
# also return the distance of data to B (in square)
# omegaData: polar coordinate of data (angle)
# disData: distance of data of data to B (in square)

polarF = function(data, B)
{
  n = dim(data)[1]
  cosData = data[, 1] - B[1]
  sinData = data[, 2] - B[2]
  tanData = sinData / cosData
  omegaData = atan(tanData)
  index1 = cosData > 0	# cos(theta)
  index2 = sinData > 0	# sin(theta)
  omegaData[!index1] = omegaData[!index1] + pi
  omegaData[index1 & (!index2)] = omegaData[index1 & (!index2)] + 2*pi
  disData = cosData^2 + sinData^2
  list(omegaData = omegaData, disData = disData)
}


# DF
# calculate the dispersion D for data given K basis: omega and v
# input: omegaData is the polar coordinate of data (from polarF)
# output: D (a vector, evalated at data)

DF = function(omegaData, c, K, omega, v)
{
  n = length(omegaData)
  omegaDifMat = matrix(rep(omegaData, K) - rep(omega, each = n), n, K)
  omegaDifMat[omegaDifMat<0] = omegaDifMat[omegaDifMat<0] + 2 * pi
  logD = cbind(1, RF(omegaDifMat)) %*% c
  D = exp(logD)
  D[D > 100] = 100  #how was the choice of the threshold be decided  ?
  D
}


# A fixed-dimension MCMC update for the location and dispersion 
simF1 = function(data, nsim, subset,shape, scale, Bmov = 0.2, vmov = 0.4, K, Bstart = c(0, 0), omegastart, vstart)
{
  n = dim(data)[1]
  Bsim = matrix(NA, 2, nsim)
  Ksim = rep(NA, nsim)
  vsim = list()	
  omegasim = list()
  omegaDatasim = matrix(NA, n, nsim)
  csim = list()
  Dsim = matrix(NA, n, nsim)
  KDsim = matrix(NA, nsim)
  densim = matrix(NA, n, nsim)
  
  # start
  B = Bstart
  omega = omegastart	# omega start
  v = vstart	# vstart
  c = cF(K, omega, v)
  polarData = polarF(data, B)
  omegaData = polarData$omegaData
  disData = polarData$disData
  D = DF(omegaData, c, K, omega, v)
  x = as.vector(disData / D)
  KD = integrate(DF, 0, 2*pi, c = c, K = K, omega = omega, v = v, 
                 stop.on.error = F)$value
  den = gF(x,shape,scale) / KD
  llk = sum(log(den))
  
  # indicators for acceptance rates
  Bacc = 0
  Eacc = Eall = 0
  Facc = Fall = 0
  
  for (i in 1:nsim)
  {	
    for (j in 1:2)
    {
      # sample B1
      Btem = B[j]
      Btemnew = rnorm(1, Btem, Bmov)	# propose new B1: symmetric
      Bnew = B
      Bnew[j] = Btemnew
      priorratio = dnorm(Btemnew, 0, 10) / dnorm(Btem, 0, 10)
      
      polarDatanew = polarF(data, Bnew)
      omegaDatanew = polarDatanew$omegaData
      disDatanew = polarDatanew$disData
      Dnew = DF(omegaDatanew, c, K, omega, v)
      xnew = disDatanew / Dnew
      dennew = sqrt(gF(xnew,shape,scale) / KD)
      llknew = sum(log(dennew)) 
      lkratio = exp(llknew - llk) 
      r = lkratio * priorratio 
      u = runif(1) 
      if (r >= u)
      { B = Bnew; omegaData = omegaDatanew; disData = disDatanew; 
      D = Dnew; x = xnew; den = dennew; llk = llknew; Bacc = Bacc+1 }
    }
    
    # sample D (i.e., omega and v)
    type = runif(1)		# indicator		
    if (type < 0.5) 
    {
      # sample v
      Eall = Eall + 1
      vnew = v
      k = sample(K, 1)
      vk = v[k]
      vknew = rnorm(1, vk, vmov)
      vnew[k] = vknew
      priorratio = dnorm(vknew, 0, 1.4) / dnorm(vk, 0, 1.4)   # vprior
      
      cnew = cF(K, omega, vnew)
      Dnew = DF(omegaData, cnew, K, omega, vnew)
      xnew = disData / Dnew
      KDnew = integrate(DF, 0, 2*pi, c = cnew, K = K, omega = omega, 
                        v = vnew, stop.on.error = FALSE)$value
      dennew = sqrt(gF(xnew,shape,scale) / KDnew)  
      llknew = sum(log(dennew)) 
      lkratio = exp(llknew - llk)   
      r = lkratio * priorratio 
      u = runif(1)
      if (r >= u)
      { v = vnew; c = cnew; D = Dnew; x = xnew; KD = KDnew; 
      den = dennew; llk = llknew; Eacc = Eacc+1 }
    } else
    {
      # sample omega
      Fall = Fall + 1
      omeganew = omega
      k = sample(K, 1)
      omegak = omega[k]
      omegaknew = -1
      while (omegaknew < 0 | omegaknew > 2*pi)
      { omegaknew = rnorm(1, omegak, omegamov) }
      omeganew[k] = omegaknew
      priorratio = 1	# uniform
      cnew = cF(K, omeganew, v)
      Dnew = DF(omegaData, cnew, K, omeganew, v)
      xnew = disData / Dnew
      KDnew = integrate(DF, 0, 2*pi, c = cnew, K = K, omega = omeganew, 
                        v = v, stop.on.error = FALSE)$value
      dennew = sqrt(gF(xnew,shape,scale) / KDnew) 
      llknew = sum(log(dennew)) 
      lkratio = exp(llknew - llk) 
      r = lkratio 
      u = runif(1) 
      if (r >= u)
      { omega = omeganew; c = cnew; D = Dnew; x = xnew; KD = KDnew; 
      den = dennew; llk = llknew; Facc = Facc+1 }
    } 
    print(i)
    if (i %in% subset)
    {
      l = l + 1
      Bsim[, l] = B
      omegasim[[l]] = omega
      vsim[[l]] = v
      csim[[l]] = c
      densim[, l] = den
      KDsim[l] = KD
      Ksim[l] = K
    }
    
  }
  list(B = Bsim, v= vsim, omega = omegasim, omegaData = omegaDatasim, c = csim, D = Dsim, den = densim, K = Ksim, KD = KDsim, rate = c(Bacc/(2*nsim), Eacc/Eall, Facc/Fall))
}





# simF_adap
# simulate MCMC samples using reversive jump MCMC method using a secondary Markov chain to modify proposed moves
# input:
# data: 2*n matrix
# nsim: number of simulation runs (default = 12000)
# nsimk : number of the secondary simulation runs (default = 12000)
# subset: vector, length L, indexes of the MCMC samples saved 
#		  (default = seq(20000, nsim, 10))
# shape: the shape parameter in the weibull density generator
# scale : the scale parameter in the weibull density generator
# sym: indicator: sym=0 (default) not symetric D;
#				  sym=1: symmetric D (not finished)
# Bmov: sd in the proposal of B (default = 0.2)
# vmov: sd in the proposal of v (default = 0.4)
# omegamov: sd in the proposal of omega (default = 0.2)
# vbirthmov: sd in the proposal of new v when K->K+1 (default = 0.2)
# Kstart: starting value of K (default = 5, arbitrary and not important)
# Bstart: starting value of B (default = c(0,0), arbitrary and not important)
# output: a list
# B: 2*L matrix, samples of B (center)
# v: a list, L elements, samples of v (values of knots)
# omega: a list, L elements, samples of omega (locations of knots)
# c: a list, L elements, samples of c (coefficients of smoothing function)
# den: n*L matrix, samples of densities at data
# K: vector of L elements, samples of K (numbers of knots)
# KD: vector of L elements, samples of KD (normalizer)
# rate: c(Bacc, Cacc, Call, Dacc, Dall, Eacc, Eall, Facc, Fall)
# 		acceptance rates for each type of moves

simF = function(data, nsim, subset,shape, scale, sym = 0, Bmov = 0.2, vmov = 0.4, vbirthmov = 0.2, typemov = jkmov, Kstart = 6, Bstart = c(0, 0))
{ 
  
  L = length(subset)
  n = dim(data)[1]
  Bsim = matrix(NA, nrow=2, ncol=L)
  Ksim = rep(NA, L)
  vsim = vector("list", L)	
  omegasim = vector("list", L)	
  csim = vector("list", L)	
  KDsim = rep(NA, L)
  densim = matrix(NA, n, L)
  
  # start
  B = Bstart
  K = Kstart
  omega = 2*(1:K) / (2*K+1) * 2*pi	# omega start:equal-spaced (0, 2*pi)
  v = rep(0, K)	# vstart: all 0
  c = cF(K, omega, v)
  polarData = polarF(data, B)
  omegaData = polarData$omegaData
  disData = polarData$disData
  D = DF(omegaData, c, K, omega, v)#calculate the dispersion D for data given K basis: omega and v
  x = as.vector(disData / D) 
  KD = (1/(2*pi))*integrate(DF, 0, 2*pi, c = c, K = K, omega = omega, v = v, 
                            stop.on.error = F)$value
  
  den = gF(x,shape,scale) / KD 
  llk = sum(log(den)) 
  
  # indicators for acceptance rates
  Bacc = 0 
  Cacc = Call = 0 
  Dacc = Dall = 0 
  Eacc = Eall = 0 
  Facc = Fall = 0
  
  l = 0
  for (i in 1:nsim)
  {
    for (j in 1:2)
    {
      # sample B1
      Btem = B[j]
      Btemnew = rnorm(1, Btem, Bmov)	# propose new B1: symmetric
      Bnew = B
      Bnew[j] = Btemnew
      priorratio = dnorm(Btemnew, 0, 10) / dnorm(Btem, 0, 10)
      
      polarDatanew = polarF(data, Bnew)
      omegaDatanew = polarDatanew$omegaData
      disDatanew = polarDatanew$disData
      Dnew = DF(omegaDatanew, c, K, omega, v)
      xnew = as.vector(disDatanew / Dnew)
      dennew = gF(xnew,shape,scale) / KD
      llknew = sum(log(dennew))
      lkratio = exp(llknew - llk)
      r = lkratio * priorratio
      u = runif(1)
      if (r >= u)
      { B = Bnew; omegaData = omegaDatanew; disData = disDatanew; 
      D = Dnew; x = xnew; den = dennew; llk = llknew; Bacc = Bacc+1 }
    }
    
    
    # sample D (i.e., K, omega and v)
    type = runif(1)		# type determines one of 4 kinds of moves for D
    if (type < typemov[K, 1])
    {
      # death: K -> K-1
      Call = Call + 1
      k = sample(K, 1)
      vk = v[k]
      vnew = v[-k]
      omegak = omega[k]
      omeganew = omega[-k]
      omegakDif = matrix(rep(omegak, K-1) - omeganew, 1, K-1)
      omegakDif[omegakDif<0] = omegakDif[omegakDif<0] + 2*pi
      cnew = cF(K-1, omeganew, vnew)
      logD = cbind(1, RF(omegakDif)) %*% cnew
      Dnew = DF(omegaData, cnew, K-1, omeganew, vnew)
      xnew = as.vector(disData / Dnew)
      KDnew = (1/(2*pi))*integrate(DF, 0, 2*pi, c = cnew, K = K-1, omega = omeganew, 
                                   v = vnew, stop.on.error = FALSE)$value 
      dennew = gF(xnew,shape,scale) / KDnew 
      llknew = sum(log(dennew))  
      lkratio = exp(llknew - llk)
      if(k==K){
        r = (lkratio* dnorm(vk, logD, vbirthmov)*(2*pi-omega[k-1])) / (dnorm(vk, 0, 1.4)*(2*pi-omega[k])*(omega[k]-omega[k-1])*2*(2*K+1))
       } else 
         if (k==1){ 
            r = (lkratio* dnorm(vk, logD, vbirthmov)*(omega[k+1]-0)) / (dnorm(vk, 0, 1.4)*(omega[k+1]-omega[k])*(omega[k]-0)*2*(2*K+1))
      } else
        { r = (lkratio* dnorm(vk, logD, vbirthmov)*(omega[k+1]-omega[k-1])) / (dnorm(vk, 0, 1.4)*(omega[k+1]-omega[k])*(omega[k]-omega[k-1])*2*(2*K+1))}
       u = runif(1)
      if (r >= u)
      { K = K-1; v = vnew; omega = omeganew; c = cnew; D = Dnew; x = xnew;
      KD = KDnew; den = dennew; llk = llknew; Cacc = Cacc+1 }
    } else
      if (type < typemov[K, 2])
      {
        # birth: K -> K+1
        Dall = Dall + 1
        omegak = runif(1, 0, 2*pi)
        k = omegak%/%(2*2*pi/(2*K+1))
        omeganew = c(omega, omegak)
        id  = c( seq_along(omega), k+0.5 )
        omeganew = omeganew[order(id)]
        omegakDif = matrix(rep(omegak, K) - omega, 1, K)
        omegakDif[omegakDif<0] = omegakDif[omegakDif<0] + 2*pi
        logD = cbind(1, RF(omegakDif)) %*% c
        vk = rnorm(1, logD, vbirthmov)   # vbirthmov move
        vnew = c(v, vk)
        id  = c( seq_along(v), k+0.5 )
        vnew = vnew[order(id)]
        cnew = cF(K+1, omeganew, vnew)
        Dnew = DF(omegaData, cnew, K+1, omeganew, vnew)
        xnew = as.vector(disData / Dnew)
        KDnew = (1/(2*pi))* integrate(DF, 0, 2*pi, c = cnew, K = K+1, omega = omeganew, 
                                      v = vnew, stop.on.error = FALSE)$value
        dennew = gF(xnew,shape,scale) / KDnew
        llknew = sum(log(dennew))
        lkratio = exp(llknew - llk)
        if(k == K){
          r = (lkratio * dnorm(vk, 0, 1.4)*2*(2*K+3)*(omegak-omega[k])*(2*pi-omegak)) / ((2*pi-omega[k])*dnorm(vk, logD, vbirthmov))
       } else 
         if(k == 0) {
           r = (lkratio * dnorm(vk, 0, 1.4)*2*(2*K+3)*(omegak-0)*(omega[k+1]-omegak)) / ((omega[k+1]-0)*dnorm(vk, logD, vbirthmov)) 
           } else
           {
             r = (lkratio * dnorm(vk, 0, 1.4)*2*(2*K+3)*(omegak-omega[k])*(omega[k+1]-omegak)) / ((omega[k+1]-omega[k])*dnorm(vk, logD, vbirthmov)) 
           }
         u = runif(1) 
        if (r >= u)
        { K = K+1; v = vnew; omega = omeganew; c = cnew; D = Dnew; x = xnew;
        KD = KDnew; den = dennew; llk = llknew; Dacc = Dacc+1 }
      } else
        if (type < typemov[K, 3]) 
        {
          # sample v (change value)
          Eall = Eall + 1
          vnew = v
          k = sample(K, 1)
          vk = v[k]
          vknew = rnorm(1, vk, vmov)
          vnew[k] = vknew
          priorratio = dnorm(vknew, 0, 1.4) / dnorm(vk, 0, 1.4)   # vprior
          
          cnew = cF(K, omega, vnew)
          Dnew = DF(omegaData, cnew, K, omega, vnew)
          xnew = as.vector(disData / Dnew)
          KDnew = (1/(2*pi))*integrate(DF, 0, 2*pi, c = cnew, K = K, omega = omega, 
                                       v = vnew, stop.on.error = FALSE)$value
          dennew = gF(xnew,shape,scale) / KDnew
          llknew = sum(log(dennew))
          lkratio = exp(llknew - llk)
          r = lkratio * priorratio
          u = runif(1)
          if (r >= u)
          { v = vnew; c = cnew; D = Dnew; x = xnew; KD = KDnew; 
          den = dennew; llk = llknew; Eacc = Eacc+1 }
        } else
        { 
          # sample omega (change location)
          
          
         Fall = Fall + 1
         omeganew = omega
         k = sample(K, 1)
         omegak = omega[k]
         if( k==K) {
           omegaknew = runif(1, omega[k-1], 2*pi)
         } else 
            if (k==1){
                   omegaknew = runif(1,0,omega[k+1])
          } else {
            omegaknew = runif(1, omega[k-1], omega[k+1])
                  }
          omeganew[k] = omegaknew
         if (k==K){
            priorratio = ((2*pi-omegak)*(omegak-omega[k-1]))/((2*pi-omega[k])*(omega[k]-omega[k-1]))
         } else 
            if (k==1){
             priorratio = ((omega[k+1]-omegak)*(omegak-0))/((omega[k+1]-omega[k])*(omega[k]-0))
           } else {
            priorratio = ((omega[k+1]-omegaknew)*(omegaknew-omega[k-1]))/((omega[k+1]-omega[k])*(omega[k]-omega[k-1]))}
          cnew = cF(K, omeganew, v)
          Dnew = DF(omegaData, cnew, K, omeganew, v)
          xnew = as.vector(disData / Dnew)
          KDnew = (1/(2*pi))*integrate(DF, 0, 2*pi, c = cnew, K = K, omega = omeganew, 
                                      v = v, stop.on.error = FALSE)$value
          dennew = gF(xnew,shape,scale) / KDnew
          llknew = sum(log(dennew))
          lkratio = exp(llknew - llk)
          r = lkratio * priorratio
          u = runif(1)
          if (r >= u)
          { omega = omeganew; c = cnew; D = Dnew; x = xnew; KD = KDnew; 
          den = dennew; llk = llknew; Facc = Facc+1 }
        }
    
    print(i)
    if (i %in% subset)
    {
      l = l + 1
      Bsim[, l] = B
      omegasim[[l]] = omega
      vsim[[l]] = v
      csim[[l]] = c
      densim[, l] = den
      KDsim[l] = KD
      Ksim[l] = K
    }
  }
  list(B = Bsim, v = vsim, omega = omegasim, c = csim,
       den = densim, K = Ksim, KD = KDsim, 
       rate = c(Bacc/(2*nsim), Cacc/Call, Dacc/Dall, Eacc/Eall, Facc/Fall))
}  

A<-matrix(c(1,0.5,0.5,1),byrow=T,ncol=2) #covariance matrix
#Check positive definiteness
eigen(A) $values
# normalize the data points
#cov_12<-(svd(A)$v)%*%sqrt(diag(svd(A)$d))
#cov_neg12<- solve(cov_12)

m = 2 #dimension 
l = 2
muNu = rep(0)
sigmaNu = 1.4
delta = 0.1

Kmax = 40
Kprior = KpriorF(1:Kmax)

bk = Kprior[-1] / Kprior[-Kmax]
dk = c(0, 1 / bk)[-Kmax]
bk[bk > 1] = 1
dk[dk > 1] = 1

c = min(0.9 / (dk + bk))
bk = c * bk
dk = c * dk
pik = etak = (1 - bk - dk) / 2

jk = cbind(dk, bk, etak, pik)
jkmov = matrix(c(dk, dk + bk, dk + bk + etak), Kmax - 1, 3)


nsim = 12000
subset = seq(2010, nsim, 10)


sim.data.norm <- rmnorm(n=5000,mean=c(0,0),varcov=A)
#for (i in 1:n){
# sim.data.norm[i,]=cov_neg12%*%(sim.data.norm[i,])
#}
sim.data.norm <- data.frame(x=sim.data.norm[,1],y=sim.data.norm[,2])
sim.data.norm[,1]<-(sim.data.norm[,1]-mean(sim.data.norm[,1]))/sd(sim.data.norm[,1])
sim.data.norm[,2]<-(sim.data.norm[,2]-mean(sim.data.norm[,2]))/sd(sim.data.norm[,2])
sim.data.norm$lsq <- sim.data.norm$x^2+sim.data.norm$y^2

#re-order df by lsq (decreasing)
sim.data.norm <- sim.data.norm [with(sim.data.norm, order(lsq, decreasing = TRUE)), ]
sim.data.norm <- sim.data.norm [,1:2]

#plot(sim.data.norm[,1:2], pch=19, cex=0.5, xlab="", ylab="", cex.axis=1.4, xlim=c(-4,4), ylim=c(-4,4))

#Bmov = 0.05, vmov = 0.05,vbirthmov = 0.08
#Bmov = 0.05, vmov = 0.04,vbirthmov = 0.11
#Bmov = 0.06, vmov = 0.05,vbirthmov = 0.11,
res_normal = simF(sim.data.norm, nsim,subset, shape = 0.9726099, scale =2.150028,sym=0, Bmov = 0.05, vmov = 0.05,vbirthmov = 0.11, typemov = jkmov, Kstart = 6, Bstart = c(0, 0))

save(res_normal, file= "res_normal_test_RJMCMC.Rda")
