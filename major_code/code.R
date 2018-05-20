# code
# Modelling directional dispersion through hyperspherical log-splines

# m = 2, l = 2

# gF
# density generator: normal

gF = function(x)
{ exp(- x/2) } 

#gF = function(x,k, lamda)
#{exp(-(x/lambda)^k)}

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
	D[D > 100] = 100
	D
}






# simF
# simulate MCMC samples using reversive jump MCMC method
# input:
# data: 2*n matrix
# nsim: number of simulation runs (default = 120000)
# subset: vector, length L, indexes of the MCMC samples saved 
#		  (default = seq(20000, nsim, 10))
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

simF = function(data, nsim, subset, sym = 0, Bmov = 0.2, vmov = 0.4, omegamov = 0.2, vbirthmov = 1, typemov = jkmov, Kstart = 6, Bstart = c(0, 0))
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
	x = disData / D
	KD = integrate(DF, 0, 2*pi, c = c, K = K, omega = omega, v = v, 
			stop.on.error = F)$value
	den = gF(x) / KD
	llk = sum(log(den))#??
	
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
			xnew = disDatanew / Dnew
			dennew = gF(xnew) / KD
			llknew = sum(log(dennew))
			lkratio = exp(llknew - llk)
			r = lkratio * priorratio
			u = runif(1)
			if (r >= u)
			{ B = Bnew; omegaData = omegaDatanew; disData = disDatanew; 
			  D = Dnew; x = xnew; den = dennew; llk = llknew; Bacc = Bacc+1 }
		}
	  
	  #sample the parameter gamma in the density generator
	  
	  
				
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
			xnew = disData / Dnew
			KDnew = integrate(DF, 0, 2*pi, c = cnew, K = K-1, omega = omeganew, 
						v = vnew, stop.on.error = FALSE)$value
			dennew = gF(xnew) / KDnew
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
			xnew = disData / Dnew
			KDnew = integrate(DF, 0, 2*pi, c = cnew, K = K+1, omega = omeganew, 
						v = vnew, stop.on.error = FALSE)$value
			dennew = gF(xnew) / KDnew
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
			xnew = disData / Dnew
			KDnew = integrate(DF, 0, 2*pi, c = cnew, K = K, omega = omega, 
						v = vnew, stop.on.error = FALSE)$value
			dennew = gF(xnew) / KDnew
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
			xnew = disData / Dnew
			KDnew = integrate(DF, 0, 2*pi, c = cnew, K = K, omega = omeganew, 
						v = v, stop.on.error = FALSE)$value
			dennew = gF(xnew) / KDnew
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


# simF1
# fix K, other parts similar to simF

simF1 = function(data, nsim, Bmov = 0.2, vmov = 0.4, omegamov = 0.2, K = 6, Bstart = c(0, 0))
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
	x = disData / D
	KD = integrate(DF, 0, 2*pi, c = c, K = K, omega = omega, v = v, 
			stop.on.error = F)$value
	den = gF(x) / KD
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
			dennew = gF(xnew) / KD
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
			dennew = gF(xnew) / KDnew
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
			dennew = gF(xnew) / KDnew
			llknew = sum(log(dennew))
			lkratio = exp(llknew - llk)
			r = lkratio
			u = runif(1)
			if (r >= u)
			{ omega = omeganew; c = cnew; D = Dnew; x = xnew; KD = KDnew; 
			  den = dennew; llk = llknew; Facc = Facc+1 }
		}
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
	list(B = Bsim, v= vsim, omega = omegasim, omegaData = omegaDatasim, c = csim, D = Dsim, den = densim, K = Ksim, KD = KDsim, rate = c(Bacc, Eacc, Eall, Facc, Fall))
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
	x = disData / D
	KD = integrate(DF, 0, 2*pi, c = c, K = K, omega = omega, v = v, stop.on.error = F)$value
	den = gF(x) / KD
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
			priorratio = dnorm(Btemnew, 0, 10) / dnorm(Btem, 0, 10)
		
			polarDatanew = polarF(data, Bnew)
			omegaDatanew = polarDatanew$omegaData
			disDatanew = polarDatanew$disData
			Dnew = DF(omegaDatanew, c, K, omega, v)
			xnew = disDatanew / Dnew
			dennew = gF(xnew) / KD
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
		dennew = gF(xnew) / KDnew
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

p4NR = function(delta, m, den, pstart = 0.1, alpha = 1, tol = 1.e-5, itermax = 1000)
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
# den: densities at give data points when sampling parameters are given
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
#		post[i] = uniroot(p4eqn, range, delta = delta, m = m, den = dentem)$root
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
		denz[k] = gF(xz) / KD[k]
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

contourPlotF = function(res, subset = 1:L, data = data, xngrid = 20, yngrid = 20)
{
	xrange = c(min(data[, 1]) - 0.1, max(data[, 1]) + 0.1)
	yrange = c(min(data[, 2]) - 0.1, max(data[, 2]) + 0.1)
	
	rescontour = postxyzF(res, subset, xrange, yrange, xngrid, yngrid)
	plot(data, pch = 19, cex = 0.5, type = 'p', 
		xlab = '', ylab = '')
	contour(rescontour$x, rescontour$y, rescontour$z, labels = "", 
		xlab = "", ylab = "", add = TRUE)
}


