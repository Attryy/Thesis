# code
# Modelling directional dispersion through hyperspherical log-splines


# m = 2, l = 2

# gF 
# density generator: weibull distribution
# Should note the different parameterization of the density generators for D-class and homothetic density 
gF = function(x,shape, scale)
{(x/scale)^{shape-1}*exp(-(x^shape)/(scale^shape))}


# gF_t 
# density generator from a bivariate t-distribution with df=4
gF_t = function(x,p=2,v=4)
{ (1+x/v)^(-0.5*(v+p))}


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
  omegaData[which(is.na(tanData))] = 0 
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






# simF
# simulate MCMC samples using reversive jump MCMC method
# input:
# data: 2*n matrix
# nsim: number of simulation runs (default = 120000)
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

simF = function(data, nsim, subset,shape, scale, sym = 0, Bmov = 0.2, vmov = 0.4, omegamov = 0.2, vbirthmov = 1, typemov = jkmov, Kstart = 6, Bstart = c(0, 0))
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
  omega = (1:K - 1/2) / K * 2*pi	# omega start:equal-spaced (0, 2*pi)
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
      if(r == Inf | r == -Inf | is.na(r)) {j =1; next}
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
      r = (lkratio* dnorm(vk, logD, vbirthmov)) / dnorm(vk, 0, 1.4) 
      u = runif(1)
      if(r == Inf | r == -Inf | is.na(r)) { next}
      if (r >= u)
      { K = K-1; v = vnew; omega = omeganew; c = cnew; D = Dnew; x = xnew;
      KD = KDnew; den = dennew; llk = llknew; Cacc = Cacc+1 }
    } else
      if (type < typemov[K, 2])
      {
        # birth: K -> K+1
        Dall = Dall + 1
        omegak = runif(1, 0, 2*pi)
        omegakDif = matrix(rep(omegak, K) - omega, 1, K)
        omegakDif[omegakDif<0] = omegakDif[omegakDif<0] + 2*pi
        logD = cbind(1, RF(omegakDif)) %*% c
        vk = rnorm(1, logD, vbirthmov)   # vbirthmov move
        omeganew = c(omega, omegak)
        vnew = c(v, vk)
        cnew = cF(K+1, omeganew, vnew)
        Dnew = DF(omegaData, cnew, K+1, omeganew, vnew)
        xnew = as.vector(disData / Dnew)
        KDnew = (1/(2*pi))* integrate(DF, 0, 2*pi, c = cnew, K = K+1, omega = omeganew, 
                          v = vnew, stop.on.error = FALSE)$value
        dennew = gF(xnew,shape,scale) / KDnew
        llknew = sum(log(dennew))
        lkratio = exp(llknew - llk)
        r = lkratio * dnorm(vk, 0, 1.4) / (dnorm(vk, logD, vbirthmov)) 
        u = runif(1) 
        if(r == Inf | r == -Inf | is.na(r)) { next}
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
          if(r == Inf | r == -Inf | is.na(r)) { next}
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
          omegaknew = -1
          while (omegaknew < 0 | omegaknew > 2*pi)
          { omegaknew = rnorm(1, omegak, omegamov) }
          omeganew[k] = omegaknew
          priorratio = 1	# uniform
          cnew = cF(K, omeganew, v)
          Dnew = DF(omegaData, cnew, K, omeganew, v)
          xnew = as.vector(disData / Dnew)
          KDnew = (1/(2*pi))*integrate(DF, 0, 2*pi, c = cnew, K = K, omega = omeganew, 
                            v = v, stop.on.error = FALSE)$value
          dennew = gF(xnew,shape,scale) / KDnew
          llknew = sum(log(dennew))
          lkratio = exp(llknew - llk)
          r = lkratio
          u = runif(1)
          if(r == Inf | r == -Inf | is.na(r)) { next}
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


# simF1
# fix K, other parts similar to simF
# Bmov = 0.2, vmov = 0.4, omegamov = 0.2, K = 6, Bstart = c(0, 0

simF1 = function(data, nsim, subset,shape, scale, Bmov , vmov , omegamov , K , Bstart = c(0, 0))
{
  n = dim(data)[1]
  L = length(subset)
  Bsim = matrix(NA, 2, L)
  Ksim = rep(NA, L)
  vsim = list()	
  omegasim = list()
  omegaDatasim = matrix(NA, n, L)
  csim = list()
  Dsim = matrix(NA, n, L)
  KDsim = matrix(NA, L)
  densim = matrix(NA, n, L)
  
  # start
  B = Bstart
  # omega = (1:K - 1/2) / K * 2*pi	# omega start:equal-spaced (0, 2*pi)
  omega = c(res[[3]], runif(1,0,2*pi))
  v = rep(0, K) 	# vstart: all 0
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
      priorratio = dnorm(Btemnew, 0, 1) / dnorm(Btem, 0, 1)
      
      polarDatanew = polarF(data, Bnew)
      omegaDatanew = polarDatanew$omegaData
      disDatanew = polarDatanew$disData
      Dnew = DF(omegaDatanew, c, K, omega, v)
      xnew = disDatanew / Dnew
      dennew = gF(xnew,shape,scale) / KD
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
      dennew = gF(xnew,shape,scale) / KDnew
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


# simF2
# simulate parameters when K, and omega are fixed

simF2 = function(data, nsim, Bmov = 0.2, vmov = 0.4, K = 6, Bstart = c(0, 0))
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
  omega = (1:K - 1/2) / K * 2*pi	# omega start:equal-spaced (0, 2*pi)
  v = rep(0, K)	# vstart: all 0
  c = cF(K, omega, v)
  polarData = polarF(data, B)
  omegaData = polarData$omegaData
  disData = polarData$disData
  D = DF(omegaData, c, K, omega, v)
  x = as.vector(disData / D)
  KD = integrate(DF, 0, 2*pi, c = c, K = K, omega = omega, v = v, stop.on.error = F)$value
  den = gF(x,k,lambda) / KD
  llk = sum(log(den))
  
  # indicators for acceptance rates
  Bacc = 0
  Eacc = 0
  
  for (i in 1:nsim)
  {
    # sample B
    for (j in 1:2)
    {
      Btem = B[j]
      Btemnew = rnorm(1, Btem, Bmov)	# propose new B1: symmetric
      Bnew = B
      Bnew[j] = Btemnew
      priorratio = dnorm(Btemnew, 0, 1) / dnorm(Btem, 0, 1)
      
      polarDatanew = polarF(data, Bnew)
      omegaDatanew = polarDatanew$omegaData
      disDatanew = polarDatanew$disData
      Dnew = DF(omegaDatanew, c, K, omega, v)
      xnew = disDatanew / Dnew
      dennew = gF(xnew,k,lambda) / KD
      llknew = sum(log(dennew))
      lkratio = exp(llknew - llk)
      r = lkratio * priorratio
      u = runif(1)
      if (r >= u)
      { B = Bnew; omegaData = omegaDatanew; disData = disDatanew; 
      D = Dnew; x = xnew; den = dennew; llk = llknew; 
      Bacc = Bacc+1 }
    }
    
    # sample v
    vnew = v
    k = sample(K, 1)
    vk = v[k]
    vknew = rnorm(1, vk, vmov)
    vnew[k] = vknew
    priorratio = dnorm(vknew, 0, 1.4) / dnorm(vk, 0, 1.4)   # vprior
    
    cnew = cF(K, omega, vnew)
    Dnew = DF(omegaData, cnew, K, omega, vnew)
    xnew = disData / Dnew
    KDnew = integrate(DF, 0, 2*pi, c = cnew, K = K, omega = omega, v = vnew, stop.on.error = FALSE)$value
    dennew = gF(xnew,k,lambda) / KDnew
    llknew = sum(log(dennew))
    lkratio = exp(llknew - llk)
    r = lkratio * priorratio
    u = runif(1)
    if (r >= u)
    { v = vnew; c = cnew; D = Dnew; x = xnew; KD = KDnew; 
    den = dennew; llk = llknew; Eacc = Eacc+1 }
    
    Bsim[, i] = B
    omegasim[[i]] = omega
    vsim[[i]] = v
    omegaDatasim[, i] = omegaData
    csim[[i]] = c
    Dsim[, i] = D
    densim[, i] = den
    KDsim[i] = KD
    Ksim[i] = K
  }
  list(B = Bsim, v= vsim, omega = omegasim, omegaData = omegaDatasim, c = csim, D = Dsim, den = densim, K = Ksim, KD = KDsim, rate = c(Bacc, Eacc))
}


# simF_order
# simulate MCMC samples using reversive jump MCMC method using a different structure and proposal for the knot locations
# input:
# data: 2*n matrix
# nsim: number of simulation runs (default = 12000)
# subset: vector, length L, indexes of the MCMC samples saved 
#		  (default = seq(20000, nsim, 10))
# shape: the shape parameter in the weibull density generator
# scale : the scale parameter in the weibull density generator
# Bmov: sd in the proposal of B (default = 0.05)
# vmov: sd in the proposal of v (default = 0.05)
# vbirthmov: sd in the proposal of new v when K->K+1 (default = 0.11)
# Kstart: starting value of K (default = 6)
# Bstart: starting value of B (default = c(0,0))
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

simF_order = function(data, nsim, subset,shape, scale, Bmov , vmov , vbirthmov , typemov = jkmov, Kstart = 6, Bstart = c(0, 0))
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
  for ( i in 1:nsim )
  {
    for (j in 1:2)
    {
      # sample B1
      Btem = B[j]
      Btemnew = rnorm(1, Btem, Bmov)	# propose new B1: symmetric
      Bnew = B
      Bnew[j] = Btemnew
      priorratio = dnorm(Btemnew, 0, 1) / dnorm(Btem, 0, 1)
      
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
      if(r == Inf | r == -Inf| is.na(r)) {j = 1; next}
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
      } else {
          r = (lkratio* dnorm(vk, logD, vbirthmov)*(omega[k+1]-omega[k-1])) / (dnorm(vk, 0, 1.4)*(omega[k+1]-omega[k])*(omega[k]-omega[k-1])*2*(2*K+1))}
      u = runif(1)
      if(r == Inf | r == -Inf| is.na(r)) { next}
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
        } else {
            r = (lkratio * dnorm(vk, 0, 1.4)*2*(2*K+3)*(omegak-omega[k])*(omega[k+1]-omegak)) / ((omega[k+1]-omega[k])*dnorm(vk, logD, vbirthmov)) 
        }
        u = runif(1) 
        if(r == Inf | r == -Inf| is.na(r)) { next}
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
          if(r == Inf | r == -Inf | is.na(r)) { next}
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
            if (k==1) {
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
          if(r == Inf | r == -Inf | is.na(r)) { next}
          if (r >= u)
          { omega = omeganew; c = cnew; D = Dnew; x = xnew; KD = KDnew; 
          den = dennew; llk = llknew; Facc = Facc+1 }
        }
    i = i+1
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



# simF_t
# simulate MCMC samples using reversive jump MCMC method with student t density generator with df=4
# input:
# data: 2*n matrix
# nsim: number of simulation runs (default = 120000)
# subset: vector, length L, indexes of the MCMC samples saved 
#		  (default = seq(20000, nsim, 10))
# k: the shape parameter in the weibull density generator
# lambda : the scale parameter in the weibull density generator
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

simF_t = function(data, nsim, subset, sym = 0, Bmov = 0.2, vmov = 0.4, omegamov = 0.2, vbirthmov = 0.2, typemov = jkmov, Kstart = 6, Bstart = c(0, 0))
{
  L = length(subset)
  n = dim(data)[1]
  Bsim = matrix(NA, 2, L)
  Ksim = rep(NA, L)
  vsim = vector("list", L)	
  omegasim = vector("list", L)	
  csim = vector("list", L)	
  KDsim = rep(NA, L)
  densim = matrix(NA, n, L)
  
  # start
  B = Bstart
  K = Kstart
  omega = (1:K - 1/2) / K * 2*pi	# omega start:equal-spaced (0, 2*pi)
  v = rep(0, K)	# vstart: all 0
  c = cF(K, omega, v)
  polarData = polarF(data, B)
  omegaData = polarData$omegaData
  disData = polarData$disData
  D = DF(omegaData, c, K, omega, v)#calculate the dispersion D for data given K basis: omega and v
  x = as.vector(disData / D) 
  KD = integrate(DF, 0, 2*pi, c = c, K = K, omega = omega, v = v, 
                 stop.on.error = F)$value 
  den = gF_t(x) / KD 
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
      priorratio = dnorm(Btemnew, 0, 1) / dnorm(Btem, 0, 1)
      
      polarDatanew = polarF(data, Bnew)
      omegaDatanew = polarDatanew$omegaData
      disDatanew = polarDatanew$disData
      Dnew = DF(omegaDatanew, c, K, omega, v)
      xnew = as.vector(disDatanew / Dnew)
      dennew = gF_t(xnew) / KD
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
      KDnew = integrate(DF, 0, 2*pi, c = cnew, K = K-1, omega = omeganew, 
                        v = vnew, stop.on.error = FALSE)$value 
      dennew = gF_t(xnew) / KDnew 
      llknew = sum(log(dennew))  
      lkratio = exp(llknew - llk)
      r = lkratio / dnorm(vk, 0, 1.4) * dnorm(vk, logD, vbirthmov)
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
        omegakDif = matrix(rep(omegak, K) - omega, 1, K)
        omegakDif[omegakDif<0] = omegakDif[omegakDif<0] + 2*pi
        logD = cbind(1, RF(omegakDif)) %*% c
        vk = rnorm(1, logD, vbirthmov)   # vbirthmov move
        omeganew = c(omega, omegak)
        vnew = c(v, vk)
        cnew = cF(K+1, omeganew, vnew)
        Dnew = DF(omegaData, cnew, K+1, omeganew, vnew)
        xnew = as.vector(disData / Dnew)
        KDnew = integrate(DF, 0, 2*pi, c = cnew, K = K+1, omega = omeganew, 
                          v = vnew, stop.on.error = FALSE)$value
        dennew = gF_t(xnew) / KDnew
        llknew = sum(log(dennew))
        lkratio = exp(llknew - llk)
        r = lkratio * dnorm(vk, 0, 1.4) / dnorm(vk, logD, vbirthmov)
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
          KDnew = integrate(DF, 0, 2*pi, c = cnew, K = K, omega = omega, 
                            v = vnew, stop.on.error = FALSE)$value
          dennew = gF_t(xnew) / KDnew
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
          omegaknew = -1
          while (omegaknew < 0 | omegaknew > 2*pi)
          { omegaknew = rnorm(1, omegak, omegamov) }
          omeganew[k] = omegaknew
          priorratio = 1	# uniform
          cnew = cF(K, omeganew, v)
          Dnew = DF(omegaData, cnew, K, omeganew, v)
          xnew = as.vector(disData / Dnew)
          KDnew = integrate(DF, 0, 2*pi, c = cnew, K = K, omega = omeganew, 
                            v = v, stop.on.error = FALSE)$value
          dennew = gF_t(xnew) / KDnew
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
  list(B = Bsim, v = vsim, omega = omegasim, c = csim,
       den = densim, K = Ksim, KD = KDsim, 
       rate = c(Bacc/(2*nsim), Cacc/Call, Dacc/Dall, Eacc/Eall, Facc/Fall))
}


# simF_order_t
# simulate MCMC samples using reversive jump MCMC method with student t density generator with df=4 under the different 
#          setup of the knot location
# input:
# data: 2*n matrix
# nsim: number of simulation runs (default = 120000)
# subset: vector, length L, indexes of the MCMC samples saved 
#		  (default = seq(20000, nsim, 10))
# Bmov: sd in the proposal of B (default = 0.2)
# vmov: sd in the proposal of v (default = 0.4)
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

simF_order_t = function(data, nsim, subset, Bmov , vmov , vbirthmov , typemov = jkmov, Kstart = 6, Bstart = c(0, 0))
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
  
  den = gF_t(x) / KD 
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
      priorratio = dnorm(Btemnew, 0, 1) / dnorm(Btem, 0, 1)
      
      polarDatanew = polarF(data, Bnew)
      omegaDatanew = polarDatanew$omegaData
      disDatanew = polarDatanew$disData
      Dnew = DF(omegaDatanew, c, K, omega, v)
      xnew = as.vector(disDatanew / Dnew)
      dennew = gF_t(xnew) / KD
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
      dennew = gF_t(xnew) / KDnew 
      llknew = sum(log(dennew))  
      lkratio = exp(llknew - llk)
      if(k==K){
        r = (lkratio* dnorm(vk, logD, vbirthmov)*(2*pi-omega[k-1])) / (dnorm(vk, 0, 1.4)*(2*pi-omega[k])*(omega[k]-omega[k-1])*2*(2*K+1))
      } else 
        if (k==1) { 
          r = (lkratio* dnorm(vk, logD, vbirthmov)*(omega[k+1]-0)) / (dnorm(vk, 0, 1.4)*(omega[k+1]-omega[k])*(omega[k]-0)*2*(2*K+1))
      } else {
         r = (lkratio* dnorm(vk, logD, vbirthmov)*(omega[k+1]-omega[k-1])) / (dnorm(vk, 0, 1.4)*(omega[k+1]-omega[k])*(omega[k]-omega[k-1])*2*(2*K+1))}
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
        dennew = gF_t(xnew) / KDnew
        llknew = sum(log(dennew))
        lkratio = exp(llknew - llk)
        if(k == K){
          r = (lkratio * dnorm(vk, 0, 1.4)*2*(2*K+3)*(omegak-omega[k])*(2*pi-omegak)) / ((2*pi-omega[k])*dnorm(vk, logD, vbirthmov))
        } else 
          if(k == 0) {
            r = (lkratio * dnorm(vk, 0, 1.4)*2*(2*K+3)*(omegak-0)*(omega[k+1]-omegak)) / ((omega[k+1]-0)*dnorm(vk, logD, vbirthmov)) 
        } else {
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
          dennew = gF_t(xnew) / KDnew
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
          dennew = gF_t(xnew) / KDnew
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




# p4eqn
# the equation used for solving p4 measure of posterior density 
# input: 
# p4: any value of the density at a vector x
# delta: parameter (proportion of prior used)
# m: number of samples of parameters
# den: 1*m vector of densities, each row represent density of one # MCMC sample, getting from the MCMC results


p4eqn = function(p4, delta, m, den)
{
  eta = 1 - delta
  den[den <= 1.e-20] = 1.e-20
  tem = 1 / (delta * p4 + eta * den)
  tem1 = delta * m / eta
  up = tem1 + sum(den * tem)
  down = tem1 * p4 + sum(tem)
  up / down - p4
}

# p4NR
# calculate p4 posterior density using a Newton-Raphason method

p4NR = function(delta, m, den, pstart = 0.1, alpha = 1, tol = 1.e-5, itermax = 100000)
{
  eta = 1 - delta
  p4 = pstart
  if (max(den) <= 1.e-20)
  { return(0)	}
  den[den <= 1.e-20] = 1.e-20		# for numerical computation
  F = 1
  iter = 0
  while (abs(F) >= tol & iter < itermax)
  {
    iter = iter + 1
    tem = 1 / (delta * p4 + eta * den)
    tem1 = delta * m / eta
    tem2 = tem^2
    
    up = tem1 + sum(den * tem) 
    down = tem1 * p4 + sum(tem)
    F = up / down - p4
    term1 = - delta * sum(den * tem2) * down
    term2 = up * (tem1 - delta * sum(tem2))
    deriv = (term1 + term2) / down^2 - 1
    dif = F / deriv * alpha
    new = p4 - dif
    while (new < 0)
    {
      dif = dif / 2
      new = p4 - dif
    }
    p4 = new
  }
  if (iter >= itermax) { cat("Maximum Iterition Exceeded")} else
  { p4 }
}	



# p4postDataF
# calculate the posterior density at given data points
# input:
# delta: parameter
# den: densities at given data points when sampling parameters are given
#      from MCMC chain, already calculated and saved
# range: 1*2 vector, the range of the bisection method for root of p4eqn
#		default (0, 5)
# output:
# post: 1*n vector of the posterior densties at given data points 

p4postDataF = function(delta, den, range = c(0,5))
{
  n = dim(den)[1]
  m = dim(den)[2]
  post = rep(NA, n)
  for (i in 1:n)
  {
    dentem = den[i, ]
    #	post[i] = uniroot(p4eqn, range, delta = delta, m = m, den = dentem)$root
    post[i] = p4NR(delta, m, den = dentem)
  }
  post
}


# p4postF
# general function to calculate posterior density at any given (x, y)
# input:
# z = (x, y)
# res: contains B, c, omega, v, KD, ... from MCMC chain
# output:
# a value f(z): the posterior denstiy at z = (x, y)

p4postF = function(z, res, subset)
{
  B = res$B[, subset]
  v = res$v[subset]
  c = res$c[subset]
  KD = res$KD[subset]
  m = length(KD)
  K = res$K[subset]
  omega = res$omega[subset]
  denz = rep(NA, m)
  for (k in 1:m)
  {
    polar = polarF(z, B[, k])
    omegaz = polar$omegaData
    disz = polar$disData
    Dz = DF(omegaz, c[[k]], K[k], omega[[k]], v[[k]])
    xz = disz / Dz
    denz[k] = gF(xz,shape = 1,scale = 2) / KD[k]
  }	 
  #	post = uniroot(p4eqn, c(0, 1), delta = delta, m = m, den = denz)$root
  post = p4NR(delta = delta, m = m, den = denz)
  post
}

p4postF_limit = function(z, res, subset ,shape, scale )
{
  B = res$B[, subset]
  v = res$v[subset]
  c = res$c[subset]
  KD = res$KD[subset]
  m = length(KD)
  K = res$K[subset]
  omega = res$omega[subset]
  denz = rep(NA, m)
  for (k in 1:m)
  {
    polar = polarF(z, B[, k])
    omegaz = polar$omegaData
    disz = polar$disData
    Dz = DF(omegaz, c[[k]], K[k], omega[[k]], v[[k]])
    xz = disz / Dz
    denz[k] = gF(xz,shape = shape,scale = scale) / KD[k]
  }	 
  #	post = uniroot(p4eqn, c(0, 1), delta = delta, m = m, den = denz)$root
  post = p4NR(delta = delta, m = m, den = denz)
  post
}
# postxyzF
# calculate posterior density at grids (x, y), to be used in contour plot
# input:
# res: results from MCMC run
# subset: subset of the MCMC run used
# xrange, yrange: 1*2 vector, range of xgrids and ygrids
# xngrid, yngrid: number of grids to be evaluated
# output:
# x, y, z (density)

postxyzF = function(res, subset, xrange, yrange, xngrid, yngrid)
{
  x = seq(xrange[1], xrange[2], length.out = xngrid)
  y = seq(yrange[1], yrange[2], length.out = yngrid)
  zres = matrix(NA, xngrid, yngrid)
  for (i in 1:xngrid)
  {
    for (j in 1:yngrid)
    {
      z = matrix(c(x[i], y[j]), 1, 2)
      zres[i, j] = p4postF(z, res, subset)
    }
  }
  list(x = x, y = y, z = zres)
}	

postxyzF_limit = function(res, subset, xrange, yrange, xngrid, yngrid, shape, scale)
{
  x = seq(xrange[1], xrange[2], length.out = xngrid)
  y = seq(yrange[1], yrange[2], length.out = yngrid)
  zres = matrix(NA, xngrid, yngrid)
  for (i in 1:xngrid)
  {
    for (j in 1:yngrid)
    {
      z = matrix(c(x[i], y[j]), 1, 2)
      zres[i, j] = p4postF_limit(z, res, subset, shape, scale)
    }
  }
  list(x = x, y = y, z = zres)
}	

# gridDF
# generate D values at equally spaced omega vector to calculate confidence
# interval for D
# input:
# res from MCMC samples
# subset, the subset of res used
# omeagerange, omegangrid: to generate omega vector
# output: a list
# omega: the grid of omega
# D: n*m matrix, with D values at omega

gridDF = function(res, subset, omegangrid = 100)
{
  omegavec = seq(0, 2*pi, length.out = omegangrid)
  v = res$v[subset]
  c = res$c[subset]
  KD = res$KD[subset]
  omega = res$omega[subset]
  v = res$v[subset]
  m = length(subset)
  K = res$K[subset]
  Dres = matrix(NA, length(omegavec), m)
  for (i in 1:m)
  {
    Dres[, i] = DF(omegavec, c[[i]], K[i], omega[[i]], v[[i]])
  }
  list(omega = omegavec, D = Dres)
}

# q0.025F
# 0.025 quantile

q0.025F = function(x)
{	quantile(x, 0.025)}

# q0.975F
# 0.975 quantile

q0.975F = function(x)
{	quantile(x, 0.975)}

# DciPlotF
# plot confidence interval for D
# input:
# res, subset
# output:
# a plot with mean and 95% CI of D as a function of omega

DciPlotF = function(res, subset, omegangrid = 100, ylim = c(0, 4))
{
  omegares = gridDF(res, subset, omegangrid)
  omega = omegares$omega
  D = omegares$D
  plot(omega, apply(D, 1, mean), type = 'l', ylim = ylim, xlab = "", ylab = "D", xaxt = 'n')
  #lines(omega, apply(D, 1, q0.025F), lty = 2)
  #lines(omega, apply(D, 1, q0.975F), lty = 2)
  axis(1, at = c(pi, 2*pi), c('pi', '2pi'))
}


# contourPlotF
# contour plot of the posterier density
# input:
# res, subset, data
# xngrid, yngrid: default 20 (increase to make finer contours, but slower)
# output:
# a contour plot of the posterior with originial data

contourPlotF = function(res, subset = 1:L, data = data, xngrid = 40, yngrid = 40)
{
  xrange = c(min(data[, 1]) - 0.1, max(data[, 1]) + 0.1)
  yrange = c(min(data[, 2]) - 0.1, max(data[, 2]) + 0.1)
  
  rescontour = postxyzF(res, subset, xrange, yrange, xngrid, yngrid)
  plot(data, pch = 19, cex = 0.5, type = 'p', 
       xlab = '', ylab = '')
  contour(rescontour$x, rescontour$y, rescontour$z, labels = "", 
          xlab = "", ylab = "", add = TRUE)
}

contourPlotF_limit = function(res, subset = 1:L, data = data, xngrid = 40, yngrid = 40, shape, scale)
{
  xrange = c(min(data[, 1]) - 0.1, max(data[, 1]) + 0.1)
  yrange = c(min(data[, 2]) - 0.1, max(data[, 2]) + 0.1)
  
  rescontour = postxyzF_limit(res, subset, xrange, yrange, xngrid, yngrid, shape, scale)
  plot(data, pch = 19, cex = 0.5, type = 'p', 
       xlab = '', ylab = '')
  contour(rescontour$x, rescontour$y, rescontour$z, labels = "", 
          xlab = "", ylab = "", add = TRUE)
}

# EstimatedD
# Estimated dispersion function
# input: omegaX: omega values where the estimated dispersion will be evaluated
#        res: a result list from the Reversible Jump MCMC algorithm
#        subset: the indice of the estimates to be include 
EstimatedD<-function(omegaX,res,subset){
  c=res$c[subset]
  v=res$v[subset]
  KD=res$KD[subset]
  K=res$K[subset]
  omega=res$omega[subset]
  n=length(omegaX)
  m = length(subset)
  D = matrix(rep(0,m*n),ncol=m,nrow=n)
  for(i in 1:m){
    omegadiff=matrix(rep(omegaX,K[i])-rep(omega[[i]],each=n),nrow=n,ncol=K[i],byrow=FALSE)
    omegadiff[omegadiff<0]=omegadiff[omegadiff<0]+2*pi
    logD=cbind(1,RF(omegadiff))%*%c[[i]]
    D[,i]=exp(logD)
  }
  D[D>100]=100
  return(apply(D, 1, mean))
}

# TrueD_normal
# Compute the true directional disperison for normal case
# Input: omega: oemga values where the estimated dispersion will be evaluated
#        cov: covariance matrix for normal distribution
TrueD_normal<-function(omega,cov){
  n=length(omega)
  D<-matrix(rep(0,2*n),ncol=2,nrow=n)
  y=rep(0,n)
  for(i in 1:n){
    D[i,]<-c(cos(omega[i]),sin(omega[i]))
    y[i]<-1/(t(as.matrix(D[i,]))%*%solve(cov)%*%as.matrix(D[i,]))
  }
  return(y)
}


# TrueD_Hsn
# Compute the true directional disperison for a homothetic density with the shape set = limit set of a skew normal density 
# Input: omega: oemga values where the estimated dispersion will be evaluated
#        cov: covariance matrix for skew normal distribution
#        alpha : the skewness parameter for skew normal distribution, "1*d vector"
TrueD_Hsn<-function(omega,cov, alpha){
  n=length(omega)
  D<-matrix(rep(0,2*n),ncol=2,nrow=n)
  y=rep(0,n)
  for(i in 1:n){
    D[i,] <- c(cos(omega[i]),sin(omega[i]))
    y[i] <- 1/(t(as.matrix(D[i,]))%*%solve(cov)%*%as.matrix(D[i,])+(t(as.matrix(alpha))%*%as.matrix(D[i,]))^2*(t(as.matrix(alpha))%*%as.matrix(D[i,])<0))
      }
  return(y)
}

# EllD_normal
# Compute the estimated directional disperison for a homothetic density with an elliptical shape set 
# Input: omega: oemga values where the estimated dispersion will be evaluated
#        estcov: estimated Pearson's correlation coefficient 
EllD_normal<-function(omega,estcov){
  n <- length(omega)
  D <- matrix(rep(0,2*n),ncol=2,nrow=n)
  cov <- matrix (c(1, estcov, estcov, 1), byrow = TRUE, nrow = 2)
  y=rep(0,n)
  for(i in 1:n){
    D[i,]<-c(cos(omega[i]),sin(omega[i]))
    y[i]<-1/(t(as.matrix(D[i,]))%*%solve(cov)%*%as.matrix(D[i,]))
  }
  return(y)
}



#===============================================================
# Functions for computing symmetric difference
#===============================================================
# diff_ normal.limit
# measure the difference between the estimated and the true dispersion function for normal case with the limit set estimation
# Input: omega: omega values where the estimated dispersion will be evaluated
#        cov: covariance matrix for normal distribution
#        res: estimations from the limit set approach 
#        r: the estimated scaling factor
#        rt: the true scaling factor 
diff_normal.limit<-function(omega,cov,res,r,rt){
  #y=(r*sqrt(EstimatedDlimit(omega, res.limit = res))-rt*sqrt(TrueD_normal(omega,cov)))
  y=abs(r*EstimatedDlimit(omega, res.limit = res)-rt*sqrt(TrueD_normal(omega,cov))) 
    return(y)
  } 

# diff_ normal.mc
# measure the difference between the estimated and the true dispersion function for normal case with the MCMC estimation
# Input: omega: omega values where the estimated dispersion will be evaluated
#        cov: covariance matrix for normal distribution
#        res: estimations from the limit set approach 
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        subset: the indices of the MCMC samples
diff_normal.mc<-function(omega,cov,res,r,rt,subset){ 
  #y=abs(r*sqrt(EstimatedD(omega, res = res,subset= subset))-rt*sqrt(TrueD_normal(omega,cov)))
  y=(r*sqrt(EstimatedD(omega, res = res,subset= subset))-rt*sqrt(TrueD_normal(omega,cov)))
  return(y)
 }


# diff_ normal.ell
# measure the difference between the estimated and the true dispersion function with estimator under elliptical assumption 
# Input: omega: omega values where the estimated dispersion will be evaluated
#        cov: covariance matrix for normal distribution
#        estcov: the estimated Pearson's correlation coefficient 
#        r: the estimated scaling factor
#        rt: the true scaling factor 

diff_normal.ell<-function( omega,cov, estcov, r,rt ){
 # y= r*sqrt(EllD_normal(omega, estcov))-rt*sqrt(TrueD_normal(omega,cov))
  y=abs(r*sqrt(EllD_normal(omega, estcov))-rt*sqrt(TrueD_normal(omega,cov)))
  return(y)
} 

#===============================================
# diff_ Hsn.limit
# measure the difference between the estimated and the true dispersion function for Hsn case with the limit set estimation
# Input: omega: omega values where the estimated dispersion will be evaluated
#        cov: covariance matrix for normal distribution
#        res: estimations from the limit set approach 
#        alpha : the skewness parameter for skew normal distribution, "1*d vector"
#        r: the estimated scaling factor
#        rt: the true scaling factor 
diff_Hsn.limit<-function(omega,cov,res, alpha,r,rt){
  y=abs(r*EstimatedDlimit(omega, res.limit = res)-rt*sqrt(TrueD_Hsn(omega,cov, alpha))) 
  return(y)
} 

# diff_Hsn.mc
# measure the difference between the estimated and the true dispersion function for Hsn case with the MCMC estimation
# Input: omega: omega values where the estimated dispersion will be evaluated
#        cov: covariance matrix for normal distribution
#        res: estimations from the limit set approach 
#        alpha : the skewness parameter for skew normal distribution, "d*1 vector"
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        subset: the indices of the MCMC samples
diff_Hsn.mc<-function(omega,cov,res, alpha,r,rt,subset){  
  y=abs(r*sqrt(EstimatedD(omega, res = res,subset= subset))-rt*sqrt(TrueD_Hsn(omega,cov, alpha)))
  return(y)
} 


# diff_Hsn.ell
# measure the difference between the estimated and the true dispersion function with estimator under elliptical assumption 
# Input: omega: omega values where the estimated dispersion will be evaluated
#        cov: covariance matrix for normal distribution
#        estcov: the estimated Pearson's correlation coefficient 
#        alpha : the skewness parameter for skew normal distribution, "d*1 vector"
#        r: the estimated scaling factor
#        rt: the true scaling factor 

diff_Hsn.ell<-function( omega,cov, estcov, alpha, r,rt ){
  y=abs(r*sqrt(EllD_normal(omega, estcov))-rt*sqrt(TrueD_Hsn(omega,cov, alpha)))
  return(y)
} 




# symdiff_normal.limit
# compute the symmetric difference between estimated and true dispersions for normal case with the limit set approach 
# Input: cov: the covariace matrix of the bivariate normal distribution
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        sample_sub: the subset of the sample data from a t distribution with 'distance' and 'angle' that lies in the desired section
#        mean : the location of the proposal distribution
#        res.limit: result list of the limit set approach
#        nsample: the number of sample points generated 
symdiff_normal.limit<-function(sample_sub,cov,r,rt,mean,res.limit,nsample){
  kn=length(res.limit$x)/4
  x=res.limit$x
  t= res.limit$omega
  angle=sample_sub$angle
  n= length(sample_sub$x)
  m= length(t)
  Des <- EstimatedDlimit(angle, res.limit)
  Dest<-TrueD_normal(angle, cov) 
  sample_sub$test1<- (sample_sub$distance>r^2*Des^2)
  sample_sub$test2<- (sample_sub$distance>rt^2*Dest)
  sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
  samplecm <- subset(sample_sub,sym==1)
  if( length(samplecm$x)==0){
    print("No observations sampled in this section!")
    sym =0
  } 
  else  {
    sym<- sum(dmvnorm(samplecm[,1:2],sigma=cov)/dmvnorm(samplecm[,1:2],mean= mean, sigma=cov))/nsample
  }
  return(sym)
}  

# symdiff_normal.mc
# compute the symmetric difference between estimated and true dispersions for normal case with the MCMC approach
# Input: cov: the covariace matrix of the bivariate normal distribution
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        sample_sub: the subset of the sample data from a normal distribution with 'distance' and 'angle' that lies in the desired section
#        mean : the location of the proposal distribution
#        res: the list from RJMCMC algorithm
#        subset: the subset of the MCMC sampls to be used in the estimation of the dispersion
#        nsample: the number of sample points generated 
symdiff_normal.mc<- function(sample_sub, cov, r, rt, mean, res, subset,nsample){
  n = length(sample_sub$x)
  angle <- sample_sub$angle 
  Dest<-TrueD_normal(angle, cov) 
  Des<-EstimatedD(angle,res=res,subset)
  sample_sub$test1<- (sample_sub$distance>r^2*Des)
  sample_sub$test2<- (sample_sub$distance>rt^2*Dest)
  sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
  samplecm <- subset(sample_sub,sym==1)
  if( length(samplecm$x)==0){
    print("No observations sampled in this section!")
    sym = 0 
  }
  else {
    sym<- sum(dmvnorm(samplecm[,1:2],sigma=cov)/dmvnorm(samplecm[,1:2],mean= mean, sigma=cov))/nsample
    
  }
  return(sym)
  
}




# symdiff_normal.ell
# compute the symmetric difference between estimated and true dispersions for normal case ( u= (0,0)) with the estimator under elliptical assumption 
# Input:
#        sample_sub: the subset of the sample data from a normal distribution with 'distance' and 'angle' that lies in the desired section
#        cov : the covariance matrix for the original & the proposal distribution 
#        estcov : the estimated Pearson's correlation coefficient from the estimator under elliptical assumption 
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        mean : the location of the proposal distribution
#        nsample: the number of sample points generated 

symdiff_normal.ell<- function(sample_sub, cov, estcov, r, rt, mean,nsample){
  n = length(sample_sub$x)
  angle <- sample_sub$angle 
  D<-TrueD_normal(angle, cov) 
  D_est<-EllD_normal(angle,estcov)
  sample_sub$test1<- (sample_sub$distance>r^2*D_est)
  sample_sub$test2<- (sample_sub$distance>rt^2*D)
  sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
  samplecm <- subset(sample_sub,sym==1)
  if( length(samplecm$x)==0){
    print("No observations sampled in this section!")
    sym = 0 
    
  }
  else {
    sym<- sum(dmvnorm(samplecm[,1:2],sigma=cov)/dmvnorm(samplecm[,1:2],mean= mean, sigma=cov))/nsample
    
  }
  return(sym)
  
}


#==================================================
# symdiff_Hsn.ell
# compute the symmetric difference between estimated and true dispersions for Hsn case ( u= (0,0)) with the estimator under elliptical assumption 
# Input:
#        sample_sub: the subset of the sample data from a normal distribution with 'distance' and 'angle' that lies in the desired section
#        cov : the covariance matrix for the original & the proposal distribution 
#        estcov : the estimated Pearson's correlation coefficient from the estimator under elliptical assumption
#        alpha : the skewness parameter for skew normal distribution, "d*1 vector"
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        location : the location of the proposal distribution
#        nsample: the number of sample points generated 

symdiff_Hsn.ell<- function(sample_sub, cov, estcov, alpha, r, rt, location,nsample){ 
  n = length(sample_sub$x)
  angle <- sample_sub$angle 
  if(length(angle) == 0){
    
    print("No observations sampled in this section!") 
    sym = 0
  }
  else {
  sample_sub$D<-TrueD_Hsn(angle, cov, alpha) 
  sample_sub$D_est<-EllD_normal(angle,estcov)
  sample_sub$test1<- (sample_sub$distance>r^2*sample_sub$D_est)
  sample_sub$test2<- (sample_sub$distance>rt^2*sample_sub$D)
  sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
  samplecm <- subset(sample_sub,sym==1)
  samplecm$distance_new <- polarF(samplecm,location)$disData
  samplecm$angle_new <- polarF(samplecm, location)$omegaData
  angle <- samplecm$angle_new
  samplecm$D_est_new<-TrueD_Hsn(angle,cov, alpha)
  if( length(samplecm$x)==0){
    print("No observations sampled in this section!") 
    sym = 0
  }
  else { 
   
    sym<- sum(exp(-0.5*(samplecm$distance/samplecm$D))/exp(-0.5*(samplecm$distance_new/samplecm$D_est_new)))/nsample
  }   
  }
  return(sym)
} 


# symdiff_Hsn.limit
# compute the symmetric difference between estimated and true dispersions for Hsn case ( u= (0,0)) with the estimator from limit set approach 
# Input:
#        sample_sub: the subset of the sample data from a normal distribution with 'distance' and 'angle' that lies in the desired section
#        cov : the covariance matrix for the original & the proposal distribution 
#        res : the result list of the limit set approach
#        alpha : the skewness parameter for skew normal distribution, "d*1 vector"
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        location : the location of the proposal distribution
#        nsample: the number of sample points generated 

symdiff_Hsn.limit<- function(sample_sub, cov, res, alpha, r, rt, location,nsample){
  n = length(sample_sub$x)
  angle <- sample_sub$angle 
  if(length(angle) == 0){
    
    print("No observations sampled in this section!") 
    sym = 0
  }
  else {
    sample_sub$D<-TrueD_Hsn(angle, cov, alpha) 
    sample_sub$D_est<-EllD_normal(angle,estcov)
    sample_sub$test1<- (sample_sub$distance>r^2*sample_sub$D_est)
    sample_sub$test2<- (sample_sub$distance>rt^2*sample_sub$D)
    sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
    samplecm <- subset(sample_sub,sym==1)
    samplecm$distance_new <- polarF(samplecm,location)$disData
    samplecm$angle_new <- polarF(samplecm, location)$omegaData
    angle <- samplecm$angle_new
    samplecm$D_est_new<-TrueD_Hsn(angle,cov, alpha)
    if( length(samplecm$x)==0){
      print("No observations sampled in this section!") 
      sym = 0
    }
    else { 
      
      sym<- sum(exp(-0.5*(samplecm$distance/samplecm$D))/exp(-0.5*(samplecm$distance_new/samplecm$D_est_new)))/nsample
    }  
  }
  return(sym)
} 




# symdiff_Hsn.mc
# compute the symmetric difference between estimated and true dispersions for Hsn case with the MCMC approach
# Input: cov: the covariace matrix of the bivariate normal distribution
#        alpha: the skewness parameter for the Hsn case "d*1" vector
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        sample_sub: the subset of the sample data from a normal distribution with 'distance' and 'angle' that lies in the desired section
#        location : the location of the proposal distribution
#        res: the list from RJMCMC algorithm
#        subset: the subset of the MCMC sampls to be used in the estimation of the dispersion
#        nsample: the number of sample points generated 

symdiff_Hsn.mc<- function(sample_sub, cov, alpha, r, rt, location, res, subset,nsample){
  n = length(sample_sub$x)
  angle1 <- sample_sub$angle1
  angle <- sample_sub$angle
  if (length(angle1) == 0){
    print("No observations sampled in this section!")
    sym = 0
  }
  else {
    sample_sub$D<-TrueD_Hsn(angle, cov, alpha) 
    sample_sub$D_est<-EstimatedD(angle, res = res,subset= subset)
    sample_sub$test1<- (sample_sub$distance>r^2*sample_sub$D_est)
    sample_sub$test2<- (sample_sub$distance>rt^2*sample_sub$D)
    sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
    samplecm <- subset(sample_sub,sym==1)
    samplecm$distance_new <- polarF(samplecm,location)$disData
    samplecm$angle_new <- polarF(samplecm, location)$omegaData
    angle <- samplecm$angle_new
    samplecm$D_est_new<-TrueD_Hsn(angle,cov, alpha)
    if( length(samplecm$x)==0){
      print("No observations sampled in this section!") 
      sym = 0
    }
    else { 
      
      sym<- sum(exp(-0.5*(samplecm$distance/samplecm$D))/exp(-0.5*(samplecm$distance_new/samplecm$D_est_new)))/nsample
    }  
  }
  return(sym)
}


#=================================================================

#==================================================
# symdiff_Ell.ell
# compute the symmetric difference between estimated and true dispersions for Elliptical case ( u= (0,0)) with the estimator under elliptical assumption 
# Input:
#        sample_sub: the subset of the sample data from a normal distribution with 'distance' and 'angle' that lies in the desired section
#        cov : the covariance matrix for the original & the proposal distribution 
#        estcov : the estimated Pearson's correlation coefficient from the estimator under elliptical assumption
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        location : the location of the proposal distribution
#        nsample: the number of sample points generated 

symdiff_Ell.ell<- function(sample_sub, cov, estcov, r, rt, location,nsample){ 
  n = length(sample_sub$x)
  angle <- sample_sub$angle 
  sample_sub$D<-TrueD_normal(angle, cov) 
  sample_sub$D_est<-EllD_normal(angle,estcov)
  sample_sub$test1<- (sample_sub$distance>r^2*sample_sub$D_est)
  sample_sub$test2<- (sample_sub$distance>rt^2*sample_sub$D)
  sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
  samplecm <- subset(sample_sub,sym==1)
  samplecm$distance_new <- polarF(samplecm,location)$disData
  samplecm$angle_new <- polarF(samplecm, location)$omegaData
  angle <- samplecm$angle_new
  samplecm$D_est_new<-TrueD_normal(angle,cov)
  if( length(samplecm$x)==0){
    print("No observations sampled in this section!") 
    sym = 0
  }
  else {
    sym<- sum(((sqrt(samplecm$distance/samplecm$D))^(-1)*exp(-sqrt(samplecm$distance/samplecm$D)))/((sqrt(samplecm$distance_new/samplecm$D_est_new))^(-1)*exp(-sqrt(samplecm$distance_new/samplecm$D_est_new))))/nsample

  }   
  return(sym)
} 



# symdiff_Ell.limit
# compute the symmetric difference between estimated and true dispersions for elliptical case ( u= (0,0)) with the estimator from limit set approach 
# Input:
#        sample_sub: the subset of the sample data from a normal distribution with 'distance' and 'angle' that lies in the desired section
#        cov : the covariance matrix for the original & the proposal distribution 
#        res : the result list of the limit set approach
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        location : the location of the proposal distribution
#        nsample: the number of sample points generated 

symdiff_Ell.limit<- function(sample_sub, cov, res, r, rt, location,nsample){
  n = length(sample_sub$x)
  angle <- sample_sub$angle 
  sample_sub$D<-TrueD_normal(angle, cov) 
  sample_sub$D_est<-EstimatedDlimit(angle, res.limit = res)
  sample_sub$test1<- (sample_sub$distance>r^2*sample_sub$D_est^2)
  sample_sub$test2<- (sample_sub$distance>rt^2*sample_sub$D)
  sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
  samplecm <- subset(sample_sub,sym==1)
  samplecm$distance_new <- polarF(samplecm,location)$disData
  samplecm$angle_new <- polarF(samplecm, location)$omegaData
  angle <- samplecm$angle_new
  samplecm$D_est_new<-TrueD_normal(angle,cov)
  if( length(samplecm$x)==0){
    print("No observations sampled in this section!") 
    sym = 0
  }
  else { 
    sym<- sum(((sqrt(samplecm$distance/samplecm$D))^(-1)*exp(-sqrt(samplecm$distance/samplecm$D)))/((sqrt(samplecm$distance_new/samplecm$D_est_new))^(-1)*exp(-sqrt(samplecm$distance_new/samplecm$D_est_new))))/nsample
  } 
  return(sym)
} 



# symdiff_Ell.mc
# compute the symmetric difference between estimated and true dispersions for elliptical case with the MCMC approach
# Input: cov: the covariace matrix of the bivariate normal distribution
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        sample_sub: the subset of the sample data from a normal distribution with 'distance' and 'angle' that lies in the desired section
#        location : the location of the proposal distribution
#        res: the list from RJMCMC algorithm
#        subthe subset of the MCMC sampls to be used in the estimation of the dispersion
#        nsample: the number of sample points generated 
symdiff_Ell.mc<- function(sample_sub, cov, r, rt, location, res, subset,nsample){
  n = length(sample_sub$x)
  angle <- sample_sub$angle 
  if( length(angle) ==0){
    print("No observations sampled in this section!")
    sym = 0
  }
  else{
  sample_sub$D<-TrueD_normal(angle, cov) 
  sample_sub$Dest<-EstimatedD(angle,res=res,subset)
  sample_sub$test1<- (sample_sub$distance>r^2*sample_sub$Dest)
  sample_sub$test2<- (sample_sub$distance>rt^2*sample_sub$D)
  sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
  samplecm <- subset(sample_sub,sym==1)
  samplecm$distance_new <- polarF(samplecm,location)$disData
  samplecm$angle_new <- polarF(samplecm, location)$omegaData
  angle <- samplecm$angle_new
  samplecm$D_est_new<-TrueD_normal(angle,cov)
  
  if( length(samplecm$x)==0){
    print("No observations sampled in this section!")
    sym = 0 
  }
  else {
    
    sym<- sum(((sqrt(samplecm$distance/samplecm$D))^(-1)*exp(-sqrt(samplecm$distance/samplecm$D)))/((sqrt(samplecm$distance_new/samplecm$D_est_new))^(-1)*exp(-sqrt(samplecm$distance_new/samplecm$D_est_new))))/nsample
  }
  }
  return(sym)
  
}




#==================================================
# symdiff_Hsn1.ell
# compute the symmetric difference between estimated and true dispersions for Hsn1 case ( u= (0,0)) with the estimator under elliptical assumption 
# Input:
#        sample_sub: the subset of the sample data from a normal distribution with 'distance' and 'angle' that lies in the desired section
#        cov : the covariance matrix for the original & the proposal distribution 
#        estcov : the estimated Pearson's correlation coefficient from the estimator under elliptical assumption
#        alpha : the skewness parameter for skew normal distribution, "d*1 vector"
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        location : the location of the proposal distribution
#        nsample: the number of sample points generated 

symdiff_Hsn1.ell<- function(sample_sub, cov, estcov, alpha, r, rt, location,nsample){  
  n = length(sample_sub$x)
  angle <- sample_sub$angle 
  if(length(angle) == 0){
    print("No obs!")
    sym =0
    
  }
  else{
  sample_sub$D<-TrueD_Hsn(angle, cov, alpha) 
  sample_sub$D_est<-EllD_normal(angle,estcov)
  sample_sub$test1<- (sample_sub$distance>r^2*sample_sub$D_est)
  sample_sub$test2<- (sample_sub$distance>rt^2*sample_sub$D)
  sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
  samplecm <- subset(sample_sub,sym==1)
  samplecm$distance_new <- polarF(samplecm,location)$disData
  samplecm$angle_new <- polarF(samplecm, location)$omegaData
  angle <- samplecm$angle_new
  samplecm$D_est_new<-TrueD_Hsn(angle,cov, alpha) 
  if( length(samplecm$x)==0){
    print("No observations sampled in this section!") 
    sym = 0
  }
  else { 
   
    sym<- sum((sqrt(samplecm$distance/samplecm$D))^(-1)*exp(-sqrt(samplecm$distance/samplecm$D))/(((sqrt(samplecm$distance_new/samplecm$D_new))^(-1)*exp(-sqrt(samplecm$distance_new/samplecm$D_new)))))/nsample
  }  
  }
  return(sym)
} 



# symdiff_Hsn1.limit
# compute the symmetric difference between estimated and true dispersions for Hsn1 case ( u= (0,0)) with the estimator from limit set approach 
# Input:
#        sample_sub: the subset of the sample data from a normal distribution with 'distance' and 'angle' that lies in the desired section
#        cov : the covariance matrix for the original & the proposal distribution 
#        res : the result list of the limit set approach
#        alpha : the skewness parameter for skew normal distribution, "d*1 vector"
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        location : the location of the proposal distribution
#        nsample: the number of sample points generated 

symdiff_Hsn1.limit<- function(sample_sub, cov, res, alpha, r, rt, location,nsample){
  n = length(sample_sub$x)
  angle <- sample_sub$angle 
  sample_sub$D<-TrueD_Hsn(angle, cov, alpha) 
  sample_sub$D_est<-EstimatedDlimit(angle, res.limit = res)
  sample_sub$test1<- (sample_sub$distance>r^2*sample_sub$D_est^2)
  sample_sub$test2<- (sample_sub$distance>rt^2*sample_sub$D)
  sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
  samplecm <- subset(sample_sub,sym==1)
  samplecm$distance_new <- polarF(samplecm,location)$disData
  samplecm$angle_new <- polarF(samplecm, location)$omegaData
  angle <- samplecm$angle_new
  samplecm$D_est_new<-TrueD_Hsn(angle,cov, alpha) 
  if( length(samplecm$x)==0){
    print("No observations sampled in this section!") 
    sym = 0
  }
  else { 
    sym<- sum(((sqrt(samplecm$distance/samplecm$D))^(-1)*exp(-sqrt(samplecm$distance/samplecm$D)))/((sqrt(samplecm$distance_new/samplecm$D_est_new))^(-1)*exp(-sqrt(samplecm$distance_new/samplecm$D_est_new))))/nsample
    
  } 
  return(sym)
} 


# symdiff_Hsn1.mc
# compute the symmetric difference between estimated and true dispersions for Hsn1 case with the MCMC approach
# Input: cov: the covariace matrix of the bivariate normal distribution
#        alpha: the skewness parameter for the Hsn case "d*1" vector
#        r: the estimated scaling factor
#        rt: the true scaling factor 
#        sample_sub: the subset of the sample data from a normal distribution with 'distance' and 'angle' that lies in the desired section
#        location : the location of the proposal distribution
#        res: the list from RJMCMC algorithm
#        subset: the subset of the MCMC sampls to be used in the estimation of the dispersion
#        nsample: the number of sample points generated 

symdiff_Hsn1.mc<- function(sample_sub, cov, alpha, r, rt, location, res, subset,nsample){ 
  n = length(sample_sub$x)
  angle = sample_sub$angle
  if( n ==0){
    print("No observations sampled in this section!")
    sym = 0
  }
  else{
  sample_sub$D<-TrueD_Hsn(angle, cov, alpha) 
  sample_sub$Dest<-EstimatedD(angle,res=res,subset)
  sample_sub$test1<- (sample_sub$distance>r^2*sample_sub$Dest)
  sample_sub$test2<- (sample_sub$distance>rt^2*sample_sub$D)
  sample_sub$sym <- (sample_sub$test1+sample_sub$test2)
  samplecm <- subset(sample_sub,sym==1)
  samplecm$distance_new <- polarF(samplecm,location)$disData
  samplecm$angle_new <- polarF(samplecm, location)$omegaData
  angle <- samplecm$angle_new
  samplecm$D_est_new<-TrueD_Hsn(angle,cov, alpha) 
  if( length(samplecm$x)==0){
    print("No observations sampled in this section!")
    sym = 0
  }
  else {
    sym<- sum(((sqrt(samplecm$distance/samplecm$D))^(-1)*exp(-sqrt(samplecm$distance/samplecm$D)))/((sqrt(samplecm$distance_new/samplecm$D_est_new))^(-1)*exp(-sqrt(samplecm$distance_new/samplecm$D_est_new))))/nsample
  }
  }
  return(sym)
}




#=================================================================
#Functions for estimating extreme quantiles
#=================================================================
#theta_hat as a function of k and X & n
theta_hat<-function(data,k,n){
  num=(1/k)*(sum(log(data[(n-k+1):n]/data[n-k])))
  den=(1/k)*(sum(log(log(n/(seq(1,k))))))-log(log((n+1)/(k+1)))
  theta<-num/den
  return(theta)
}

#quantile estimator as a function of k, X, p & n
est_quantile<-function(data,k,p,n){
  x<-data[n-k+1]*((log(1/p))/(log(n/k)))^(theta_hat(data,k,n))
  return(x)
}


# shiftedmean_normal : function to compute the shifted means for proposal 
# Input: root: the intersection points
#       cov: covariance matrix
#       rt: the true scaling factor for the pre-specified prob level  
#Output: u: a length(root) by 2 matrix with each row representing a shifted mean 
shiftedmean_normal <- function(root,cov,rt){        
  n <- length(root)
  u <- matrix (0, ncol=2, nrow = n)
  w1<-mean(c(root[1]+2*pi,root[n]))
  r1<- sqrt(TrueD_normal(w1,cov))
  u[1,] = c(rt*r1*cos(w1),rt*r1*sin(w1))
  for (i in 2:n){
    w<- mean (c(root[i-1], root[i]))
    r<- sqrt(TrueD_normal(w, cov))
    u[i,] <- c(rt*r*cos(w),rt*r*sin(w))
  } 
  return(u)
}

# shiftedmean_Hsn : function to compute the shifted means for proposal for Hsn case
# Input: root: the intersection points
#       cov: covariance matrix
#       alpha: 
#       rt: the true scaling factor for the pre-specified prob level  
#Output: u: a length(root) by 2 matrix with each row representing a shifted mean 
shiftedmean_Hsn <- function(root,cov,alpha,rt){        
  n <- length(root)
  u <- matrix (0, ncol=2, nrow = n)
  w1<-mean(c(root[1]+2*pi,root[n]))
  r1<- sqrt(TrueD_Hsn(w1,cov, alpha))
  u[1,] = c(rt*r1*cos(w1),rt*r1*sin(w1))
  for (i in 2:n){
    w<- mean (c(root[i-1], root[i]))
    r<- sqrt(TrueD_Hsn(w, cov, alpha))
    u[i,] <- c(rt*r*cos(w),rt*r*sin(w))
  } 
  return(u)
}

