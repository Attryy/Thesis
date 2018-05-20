simulateNormal = function(n, rho = 0){
     #covariance matrix
     S = diag(2)
     S[1,2] = S[2,1] = rho
     
     #simulated points
     z = rmnorm(n,mean=rep(0,2), S)
     
     #create dataframe for the points
     x = z[,1]
     y = z[,2]
     return (data.frame(x,y))
}

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
                    delta = delta)
     
     return(1/max(maxYs))
}

getMaxYSegments = function(j, x, delta){
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
     
     return (r*sin(theta))
}
