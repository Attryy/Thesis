#returns dataframe of simulated random variables with angle and radius
getMaxPoints = function(df, kn = 10){
  
  df$lsq = df$x^2 + df$y^2 #length squared
  df$t = atan(df$y/df$x) #theta
  df[nrow(df)+1,] = 0 #placeholder row  ???
  df$max = FALSE #TRUE if point has the greatest radius in that segments, FALSE otherwise
  
  #adjust for negative r. The cycle of the tan(x) function is pi.
  df = within(df, t[(x < 0) & (y < 0)] <- (t[(x < 0) & (y < 0)] + pi)) #4 quadrant
  df = within(df, t[(x < 0) & (y > 0)] <- (t[(x < 0) & (y > 0)] + pi)) #3 quadrant
  
  #normalize t to [0, 2pi]
  df = within(df, t[t < 0] <- t[t < 0] + 2*pi)
  df = within(df, t[t > 2*pi] <- t[t > 2*pi] - 2*pi)
  
  #re-order df by lsq (decreasing)
  df = df[with(df, order(lsq, decreasing = TRUE)), ]
  
  k = 2*pi/kn#angle interval per segment
  
  #vector for each segment's max value
  t = numeric(kn)
  lsq = numeric(kn)
  x = numeric(kn)
  y = numeric(kn)
  l = numeric(kn)
  
  df$max = F
  
  for (i in 1:kn){
    max = 0
    pos = nrow(df)
    
    for (j in 1:(nrow(df)-1)){
      if (df$t[j] >= k*(i-1) & df$t[j] <= k*i){
        if (df$lsq[j] > max){
          max = df$lsq[j]
          pos = j
          break #df is ordered by lsq, so first point in the segment has the greatest length
        }
      }
    }
    
    df$max[pos] = TRUE
  }
  
  
  df = df[1:(nrow(df)-1),] #delete placeholder row ???
  
  df$col = "black"
  df$col[df$max] = "red"
  
  return(df)
}

#returns data frame of max points in each segment
interpolateDataFrame = function(df){
  df2 = df[df$max,] #extract maximum values
  df2 = df2[with(df2, order(t)), ] #reorder by theta, t
  df2[nrow(df2)+1,] = df2[1,]
  df2$t[nrow(df2)] = df2$t[1] + 2*pi
  df2$l = sqrt(df2$lsq)
  
  return(df2)
}

#returns vector polynomials from interpolating through max points
interpolatePolynomials = function(df2, kn = 10, DIFF = 20){
  #spline interpolation
  n = nrow(df2)
  A = matrix(0, 4*(n-1), 4*(n-1))
  b = matrix(0, 4*(n-1), 1)
  
  #linear equations for spline
  #f continuous
  for (i in 1:(n-1)){ #f_i(x_i) = y_i, i = 1,...,n-1
    cl = (4*(i-1)+1)
    xi = df2$t[i]
    
    A[i,cl] = 1
    A[i,cl+1] = xi
    A[i,cl+2] = xi^2
    A[i,cl+3] = xi^3
    
    b[i,1] = df2$l[i]
  }
  
  for (i in 1:(n-1)){ #f_i(x_{i+1}) = y_{i+1}, i = 1,...,n-1
    cl = (4*(i-1)+1)
    rw = i + (n-1)
    xi = df2$t[i+1]
    
    A[rw,cl] = 1
    A[rw,cl+1] = xi
    A[rw,cl+2] = xi^2
    A[rw,cl+3] = xi^3
    
    b[rw,1] = df2$l[i+1]
  }
  
  #f' continuous
  for (i in 1:(n-2)){ #f_{i}'(x_{i+1}) - f_{i+1}'(x_{i+1}) = 0
    cl = (4*(i-1)+1)
    rw = i + (2*(n-1))
    xi = df2$t[i+1]
    
    A[rw,cl+1] = 1
    A[rw,cl+2] = 2*xi
    A[rw,cl+3] = 3*xi^2
    
    A[rw,cl+1+4] = -1
    A[rw,cl+2+4] = -(2*xi)
    A[rw,cl+3+4] = -(3*xi^2)
    
    b[rw,1] = 0
  }
  
  #f'' continuous
  for (i in 1:(n-2)){ #f_{i}''(x_{i+1}) - f_{i+1}(x_{i+1}) = 0
    cl = (4*(i-1)+1)
    rw = i + (2*(n-1) + (n-2))
    xi = df2$t[i+1]
    
    A[rw,cl+2] = 2
    A[rw,cl+3] = 6*xi
    
    A[rw,cl+2+4] = -2
    A[rw,cl+3+4] = -(6*xi)
    
    b[rw,1] = 0
  }
  
  #f''(0) = f''(x_n)
  A[4*(n-1)-1,3] = 2
  A[4*(n-1)-1,4] = 6*df2$t[1]
  A[4*(n-1)-1,4*(n-1)-1] = -2
  A[4*(n-1)-1,4*(n-1)] = -6*df2$t[n]
  
  #f'(0) = f'(x_n)
  A[4*(n-1),2] = 1
  A[4*(n-1),3] = 2*df2$t[1]
  A[4*(n-1),4] = 3*df2$t[1]^2
  A[4*(n-1),4*(n-1)-2] = -1
  A[4*(n-1),4*(n-1)-1] = -2*df2$t[n]
  A[4*(n-1),4*(n-1)] = -3*df2$t[n]^2
  
  x = solve(A,b) #vector of coefficients
  
  df2$kink = FALSE
  df2$diff = 0
  
  for (i in 2:(nrow(df2)-2)){
    left = (df2$l[i]-df2$l[i-1])/(df2$t[i]-df2$t[i-1])
    right = (df2$l[i+1]-df2$l[i])/(df2$t[i+1]-df2$t[i])
    
    rw = 4*(i-1)+1
    p0 = x[rw]
    p1 = x[rw+1]
    p2 = x[rw+2]
    p3 = x[rw+3]
    xi = df2$t[i]
    
    diff = 2*p2+6*p3*xi
    
    df2$diff[i] = diff
    
    if (sign(left) != sign(right) & abs(diff) >= DIFF & !df2$kink[i-1]){
      df2$kink[i] = TRUE
    }
  }
  
  for (i in 1:nrow(df2)){
    if (df2$kink[i]){
      rw = (i-1)+2*(n-1)
      cl = (4*(i-2)+1)
      A[rw,] = 0
      b[rw,1] = 0
      
      A[rw,cl] = 0
      A[rw,cl+1] = 1
      A[rw,cl+2] = 2*df2$t[i-1]
      A[rw,cl+3] = 3*df2$t[i-1]^2
      
      b[rw,1] = (df2$l[i]-df2$l[i-1])/(df2$t[i]-df2$t[i-1])
      
      A[rw,cl] = 0
      A[rw,cl+1] = 0
      A[rw,cl+2] = 2
      A[rw,cl+3] = 6*df2$t[i-1]
      
      b[rw,1] = 0
    }		
  }
  
  for (i in 1:nrow(df2)){
    if (df2$kink[i]){
      cl = (4*(i-1)+1)
      rw = (i-1) + (2*(n-1) + (n-2))
      A[rw,] = 0
      b[rw,1] = 0
      
      A[rw,cl] = 0
      A[rw,cl+1] = 1
      A[rw,cl+2] = 2*(df2$t[i])
      A[rw,cl+3] = 3*(df2$t[i])^2
      
      b[rw,1] = (df2$l[i+1]-df2$l[i])/(df2$t[i+1]-df2$t[i])
      
      A[rw,cl] = 0
      A[rw,cl+1] = 0
      A[rw,cl+2] = 2
      A[rw,cl+3] = 6*df2$t[i-1]
      
      b[rw,1] = 0
    }
  }
  
  x = solve(A,b) #vector of coefficients
  return(x)
}

getScale = function(df2, x, delta = 250){
  #Here the parameter delta can be chosen differently. Large delta means finer gird.
  maxY = 0
  
  maxYs = sapply(seq(1, (nrow(df2)-1)),
                 getMaxYSegments,
                 x = x, 
                 delta = delta,df2 = df2)
  
  return(1/max(maxYs))
}

getMaxYSegments = function(j, x, delta,df2){
  maxYSpline = sapply(seq(0, delta),
                      getMaxYSplines, 
                      delta = delta,
                      diff = (df2$t[j+1]-df2$t[j]),
                      theta0 = df2$t[j],
                      rw = 4*(j-1)+1,
                      x = x)
  
  return (max(maxYSpline))
}

getMaxYSplines = function(k, delta, diff, theta0, rw, x){
  p0 = x[rw]
  p1 = x[rw+1]
  p2 = x[rw+2]
  p3 = x[rw+3]
  theta = theta0 + k*(diff/delta)
  
  r = p0 + p1*theta + p2*theta^2 + p3*theta^3
  
  return (max(r*sin(theta), r*cos(theta)))
}



#LimitSetRadial
#to estimate the limit set without reference to the density generator.
#Input: data--a data frame
#Input: iter--number of iterations for ??
#Input: kn--number of segments to consider when estimating the shape set, namely the number of basis for cubic spline interpolation
#Output: estimates for beta and c

LimitSetRadial<-function(data,kn){
  n=dim(data)[1]
  data = getMaxPoints(df= data, kn = kn)
  data2 = interpolateDataFrame(df = data)  #obtain dataframe for max points in each segment
  x = interpolatePolynomials(df2 = data2) #obtain vector of polynomials for cubic spline
  scale = getScale(df2 = data2, x, delta = 250) #obtain the scaling factor 
  x = scale*x
  data2$x = scale*data2$x #scale the limit set
  data2$y = scale*data2$y
  #transform the data into polar coordinate with respect to the center?? Or data should be normalized with respect to the sample mean??
  angle<-polarF(data=cbind(data$x,data$y),B=c(0,0))$omegaData
  dist<-polarF(data=cbind(data$x,data$y),B=c(0,0))$disData
  #compute the radial components
  radial<-NULL
  for(j in 1: n){ 
    for (i in 1:kn){
      rw = 4*(i-1)+1
      p0 = x[rw]
      p1 = x[rw+1]
      p2 = x[rw+2]
      p3 = x[rw+3]
      
      if(angle[j]>=data2$t[i]& angle[j]<=data2$t[i+1]) {
        radial[j]<-sqrt(dist[j])/(p0+p1*angle[j]+p2*angle[j]^2 +p3*angle[j]^3)
      } 
      
      else if (angle[j]<data2$t[1]) {
        radial[j]<-sqrt(dist[j])/(x[(4*kn)-7]+x[(4*kn)-6]*(angle[j]+2*pi)+x[(4*kn)-5]*(angle[j]+2*pi)^2 + x[(4*kn-4)]*(angle[j]+2*pi)^3) 
        
      } 
    }
  }
  
  return(list(radial=radial,x=x,omega=data2$t,data= data))
} 

# EstimatedDlimit
# Estimated  function for the boundary of the shape set with the limit set estimation 
# r(omega)^2 = D(d_x)
# input: omegaX: omega values where the estimated dispersion will be evaluated
#        res.limit: a result list from the Limit set estimation 
EstimatedDlimit<-function(omegaX,res.limit){
  x=res.limit$x
  omega=res.limit$omega
  n=length(omegaX)
  m = length(omega)
  D = rep(0,n)
  for(j in 1: n){ 
    for (i in 1:(m-1)){
      rw = 4*(i-1)+1
      p0 = x[rw]
      p1 = x[rw+1]
      p2 = x[rw+2]
      p3 = x[rw+3]
      
      if(omegaX[j]>=omega[i]& omegaX[j]<=omega[i+1]) {
        D[j]<-(p0+p1*omegaX[j]+p2*omegaX[j]^2 +p3*omegaX[j]^3)
      } 
    }
      
       if (omegaX[j]<omega[1]) {
        D[j]<-(x[(4*m)-7]+x[(4*m)-6]*(omegaX[j]+2*pi)+x[(4*m)-5]*(omegaX[j]+2*pi)^2 + x[(4*m-4)]*(omegaX[j]+2*pi)^3) 
        
      } 
    }
 return(D) 
} 


#==========================================================================================
 int<-function(t,d=2,cov){
   y<-(t^(d-1))*(1/(sqrt(4*(pi^2)*det(cov))))*exp(-0.5*t^2)
   return(y)
 }
 
 
 prob<-function(r,d=2,cov,p){
   y<-integrate(int,lower=0,upper=r,d,cov)$value
   D_volume=2*pi/(sqrt(4*(4/3)*(4/3)-(4/3)^2)) #depends on the cov and n_D
   pr=1-y*d*D_volume-p
   return(pr)
 }

#==========================================================================================
#approach 1.1: fitting a Weibull density using MLE
#==========================================================================================
#log.llk
#log likelihood function of the Weibull distribution with two parameters: shape and scale
log.llk<-function(x,parameter){
  -1*(log(parameter[1])-parameter[1]*log(parameter[2])+(parameter[1]-1)*mean(log(x))-mean((x/parameter[2])^parameter[1]))
}



#==========================================================================================
#approach 1.2: fitting a Gamma density using MLE
#==========================================================================================
#log.llk
#log likelihood function of the Weibull distribution with two parameters: shape-k and scale-theta
log.llk.gamma<-function(x,shape, scale){
  n <- length(x)
 y <- (shape-1)*sum(log(x))-(1/scale)*sum(x)-n*log(gamma(shape))-n*shape*log(theta)
 return (y)}
#Use a function called 'optim'to maximize the log likelihood function
#optim(par = rep(1,2), log.llk, x = sqrt(r_true)) 

source("EstimationDispersionFunctions.R")
 

#==========================================================================================
#approach 2: construct a density generator by asymptotic assumptions
#==========================================================================================


#Estimatebc
#compute the estimated beta & c with raidal component and chosen threshold t_n
Estimatebc<-function(radial,tn){
  rn<-sort(radial,decreasing = TRUE)
  n<-length(rn)
  num<-(1/tn)*sum(log(log(n/(seq(1,tn)))))-log(log(n/tn))
  den<-(1/tn)*sum(log(rn[1:tn]))-log(rn[tn+1])
  beta_hat<-num/den 
  c_hat<-(1/tn)*sum(log(n/(seq(1,tn)))/((rn[1:tn])^beta_hat))  
  return(list(beta=beta_hat,c=c_hat))
}















