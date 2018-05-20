library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

tn = 2000 #the threshold for estimations of the von mises functions 
n = 5*10^3 # sample size
A = matrix(c(1,0.5,0.5,1),byrow=T,ncol=2) #covariance matrix

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

result<-foreach(i = 1:20, .combine=rbind) %dopar% {
  sim.data.norm <- rmnorm(n,mean = c(0,0),varcov = A)
 
  sim.data.norm <- data.frame(x = sim.data.norm[,1],y = sim.data.norm[,2])
  sim.data.norm[,1]<- (sim.data.norm[,1]-mean(sim.data.norm[,1]))/sd(sim.data.norm[,1])
  sim.data.norm[,2]<- (sim.data.norm[,2]-mean(sim.data.norm[,2]))/sd(sim.data.norm[,2])
  sim.data.norm$lsq <- sim.data.norm$x^2+sim.data.norm$y^2
  
  #re-order df by lsq (decreasing)
  sim.data.norm <- sim.data.norm [with(sim.data.norm, order(lsq, decreasing = TRUE)), ]
  sim.data.norm <- sim.data.norm [,1:2]
  
  res_normal = simF(sim.data.norm, nsim,subset, shape = 1 , scale = 2,sym = 0, Bmov = 0.05, vmov = 0.1,omegamov = 0.1,vbirthmov = 0.08, typemov = jkmov, Kstart = 6, Bstart = c(0, 0))
  res_t = simF_t(sim.data.norm,nsim,subset,sym = 0, Bmov = 0.05, vmov = 0.1,omegamov = 0.2,vbirthmov = 1,typemov = jkmov, Kstart = 6, Bstart = c(0, 0))
  
  
  omegadat = polarF(as.matrix(sim.data.norm[,1:2]),c(0,0))$omegaData
  distancedat = polarF(as.matrix(sim.data.norm[,1:2]),c(0,0))$disData
  r_hat = distancedat/EstimatedD(omegadat,res_normal,1:L)# squared radial component 
  r_hatt = distancedat/EstimatedD(omegadat,res_t,1:L)
  radial_normal = sqrt(r_hat)
  radial_t = sqrt(r_hatt)
  
  # Estimation of the Von Mises function with asymptotic results
  estimates = Estimatebc(radial_normal, tn)
  estimates_t = Estimatebc(radial_t, tn)
  beta = estimates$beta
  c = estimates$c
  
  beta_t = estimates_t$beta
  c_t = estimates_t$c
  
  root.mc = uniroot.all(diff_normal.mc,interval=c(0,2*pi),tol=.Machine$double.eps^0.25,cov=A,res=res_normal,r=1,rt=1,subset=1:L)
  root.mc_t = uniroot.all(diff_normal.mc,interval=c(0,2*pi),tol=.Machine$double.eps^0.25,cov=A,res=res_t,r=1,rt=1,subset=1:L)
  
  if(length(root.mc)==0){
    root.mc=seq(0,2*pi,length.out = 6)[2:5]
  }
  u<- shiftedmean(root.mc, cov=A, rt=1)
  if(length(root.mc_t)==0){
    root.mc_t=seq(0,2*pi,length.out = 6)[2:5]
  }
  ut<-shiftedmean(root.mc_t, cov = A, rt = 1)
  
  nsample = 10000
  sample = list()
  sample_t = list()
  subset = list()
  subset_t = list()
  nlist = nrow(u)
  nlist_t = nrow(ut)
  sym = 0
  sym_t = 0
  for (j in 1: nlist){ 
    print(j)
    sample[[j]] = rmnorm(nsample,mean = u[j,], varcov= A)
    sample[[j]] = data.frame(x = sample[[j]][,1],y = sample[[j]][,2])
    sample[[j]]$distance = polarF(sample[[j]],c(0,0))$disData
    sample[[j]]$angle = polarF(sample[[j]],B=c(0,0))$omegaData
    sample[[j]]$id = cut(sample[[j]]$angle, c(0,root.mc,2*pi))
    level = levels(sample[[j]]$id)
    sample[[j]]$id = sample[[j]]$id %>% fct_collapse( "first"= c(level[1],level[nlist+1]))
    level = levels(sample[[j]]$id) # the first level is the shifted mean between the first and the last intersection points
    subset[[j]] = subset(sample[[j]], id ==level[j])
    sym = sum(sym[j],symdiff_normal.mc(subset[[j]], cov = A, r = 1, rt = 1, mean = u[j,], res = res_normal, subset = 1:L,nsample = nsample))
    
  } 
  
  for (j in 1: nlist_t){ 
    print(j)
    sample_t[[j]]= rmnorm(nsample,mean= ut[j,], varcov= A)
    sample_t[[j]]= data.frame(x=sample_t[[j]][,1],y= sample_t[[j]][,2])
    sample_t[[j]]$distance = polarF(sample_t[[j]],c(0,0))$disData
    sample_t[[j]]$angle = polarF(sample_t[[j]],B=c(0,0))$omegaData
    sample_t[[j]]$id = cut(sample_t[[j]]$angle, c(0,root.mc_t,2*pi))
    level = levels(sample_t[[j]]$id)
    sample_t[[j]]$id = sample_t[[j]]$id %>% fct_collapse( "first"= c(level[1],level[nlist+1]))
    level = levels(sample_t[[j]]$id) # the first level is the shifted mean between the first and the last intersection points
    subset_t[[j]]= subset(sample_t[[j]], id==level[j])
    sym_t = sum(sym_t[j],symdiff_normal.mc(subset_t[[j]], cov = A, r = 1, rt = 1, mean = ut[j,], res = res_t, subset = 1:L,nsample = nsample))
    
  } 
  print(sprintf("Interation: %d", i)) 
  flush.console() 
  Sys.sleep(1)
  c(beta,c,beta_t,c_t, sym, sym_t)
}

save(result, file= "density_generator.Rda")