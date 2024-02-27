#Stuff i miiiiiight need later

#banned suckas banned



#Using H. Joe (2006) notation
# Incomplete, gives slightly different results than Joe's own library "clusterGeneration" with func rcorrmatrix
UnifCovMatSamp<-function(nsamp, covvec){
  d<-length(covvec)
  mylist<- list()
  i<-1
  param<-1+((d-2)/2) #parameter in beta distribution
  while(i<=nsamp){
    cormat<-matrix(NA, nrow=d, ncol=d)
    cormat[row(cormat)==col(cormat)]<-rep(1, d) #diagonal equal to 1
    initsamp<-2*rbeta(d-1, param, param) - 1 #changing support from (0,1) to (-1, 1)
    cormat[row(cormat)==1+col(cormat)]<-initsamp
    #cormat[row(cormat)==1+col(cormat)]<-runif((d-1), min=-1, max=1) #sample (i, i-1) uniformly
    
    cormat[row(cormat)+1==col(cormat)]<-cormat[row(cormat)==1+col(cormat)]
    
    for(k in 2:(d-1)){
      for(j in 1:(d-k)){
        r_1<-cormat[j,(j+1):(j+k-1)]
        r_3<-cormat[(j+k),(j+1):(j+k-1)]
        R_j_jk<-cormat[(j:(j+k)),(j:(j+k))]
        R_2 <- R_j_jk[2:(k), 2:(k)]
        R_2_inv <- solve(R_2)
        D_jk_2 <- (1 - t(r_1) %*% R_2_inv %*% r_1) %*% (1 - t(r_3) %*% R_2_inv %*% r_3)
        D_jk <- sqrt(abs(D_jk_2))
        #Sometimes D_jk squared is negative (not possible), but will be rejected in space of positive definite matrices i think
        cormat[j, (j+k)] <- t(r_1) %*% R_2_inv %*% r_3 + runif(1, min=-1, max=1) * D_jk
        cormat[(j+k), j] <- cormat[j, (j+k)]
      }
    }
    cond<-ifelse(eigen(cormat)$val>0, 1, 0) #Counting number of positive eigenvalues (=d <=> p. def.)
    if(sum(cond)==d){
      mylist[[i]] <- cormat
      i <- i + 1
    }
  }
  return(mylist)
}

#check of marginal distribution (self-built slightly wrong)
a<-UnifCovMatSamp(10000, c(1,1,1,1,1))
rowind<-1 #indexes
colind<-5
b<-rep(NA, 10000)
d<-rep(NA, 10000)
for(i in 1:10000){
  b[i]<-a[[i]][colind,rowind]
  d[i]<-rcorrmatrix(5, 1)[colind,rowind]
}
hist(d)


# Test of "the shape of partial correlation matrices" corollary 2 (artner, wellingerhof)
corgen<-function(dim, nsamp){
  mylist<-list()
  i<-1
  while(i<=nsamp){
    samped<-runif((dim*(dim-1)/2), min=-1, max=1)
    mat<- - diag(1, nrow=dim, ncol=dim)
    mat[lower.tri(mat)==TRUE]<- samped
    mat[upper.tri(mat)==TRUE] <- mat[lower.tri(mat)==TRUE] #partial correlation matrix constructed
    
    cormat<- - mat
    cond<-ifelse(eigen(cormat)$val>0, 1, 0)
    if(sum(cond)==dim){
      mylist[[i]] <- cormat
      i <- i + 1
    }
  }
  return(mylist)
}
a<-corgen(3, 10000)
b<-rep(NA, 10000)
for(i in 1:10000){
  b[i]<-a[[i]][2,3]
}
hist(b)

#Jarle's code
set.seed(1)
n <- 5
Sigma <- diag(n)
partialcorr <- replicate(10000,{
  repeat{
    u <- runif(n*(n-1)/2, -1, 1)
    Sigma[lower.tri(Sigma)] <- u
    Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
    if (sum(eigen(Sigma)$val>0)==n)  
      break()
  }
  P <- solve(Sigma)
  D <- diag(1/sqrt(diag(P)))
  (D %*% P %*% D)[lower.tri(P)]
})

toeplitz(c(-1,partialcorr[,1]))


n <- 3
Sigma <- diag(n)
partialcorr <- replicate(10000,{
  repeat{
    u <- runif(n*(n-1)/2, -1, 1)
    Sigma[lower.tri(Sigma)] <- u
    Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
    if (sum(eigen(Sigma)$val>0)==n)  
      break()
  }
  P <- solve(Sigma)
  D <- diag(1/sqrt(diag(P)))
  (D %*% P %*% D)[lower.tri(P)]
})



#for(i in 1:n){
#block <- Data_mat[((i-1)*m+1):(i*m),]
#print(i)
#Skipping iteration if error - happens rarely (~ 1% chance per iteration)
#Assuming this should not change the samples drastically
#tryCatch({
#paramat[(1 + (i-1)*mcmcsamps):(i*mcmcsamps),] <- MCMCmetrop1R(target_dens, theta.init=init, x=block, mcmc=mcmcsamps)
#}, error=function(e){})
#}


#for(i in 1:n){
#block <- Data_mat[((i-1)*m+1):(i*m),]
#print(i)
#Skipping iteration if error - happens rarely (~ 1% chance per iteration)
#Assuming this should not change the samples drastically
#tryCatch({
#paramat[(1 + (i-1)*mcmcsamps):(i*mcmcsamps),] <- MCMCmetrop1R(target_dens, theta.init=init, x=block, mcmc=mcmcsamps)
#}, error=function(e){})
#}


target_dens<-function(theta, x){
  d <- (1/2) * (sqrt(8 * length(theta) + 9) - 3) #integer solution to equation len(theta) = d/2 * (3+d)
  mu<-theta[1:d]
  margvar<-theta[(d+1):(d*2)]
  if(any(margvar<=0)){
    -Inf
  }
  else{
    parcorrs<-theta[((d*2)+1):((d*2) + d*(d-1)/2)]
    parcorrmat <- diag(-1, nrow=d)
    parcorrmat[lower.tri(parcorrmat)==TRUE] <- parcorrs
    parcorrmat <- parcorrmat + t(parcorrmat) - diag(diag(parcorrmat))
    
    
    if ((sum(eigen(parcorrmat)$val<0)==d)){ # if partial corr mat is negative definite
      S <- - parcorrmat
      S_inv <- solve(S)
      D_S_inv <- solve(diag(sqrt(diag(S_inv))))
      corrmat <-(D_S_inv %*% S_inv %*% D_S_inv)
      Sigma <- diag(margvar) %*% corrmat %*% diag(margvar)
      
      logdens <- dmvnorm(x, mean=mu, sigma=Sigma, log=TRUE)
      logprior <- -3/2 * sum(log(margvar))
      sum(logdens) + logprior
    }
    
    else{-Inf}
  }
}

for(i in 1:n){
  theta <- paramat[i,]
  margvar<-theta[(d+1):(d*2)]
  parcorrs<-theta[((d*2)+1):((d*2) + d*(d-1)/2)]
  parcorrmat <- diag(-1, nrow=d)
  parcorrmat[lower.tri(parcorrmat)==TRUE] <- parcorrs
  parcorrmat <- parcorrmat + t(parcorrmat) - diag(diag(parcorrmat))
  
  S <- - parcorrmat
  S_inv <- solve(S)
  D_S_inv <- solve(diag(sqrt(diag(S_inv))))
  corrmat <- (D_S_inv %*% S_inv %*% D_S_inv)
  paramat[i, ((d*2)+1):((d*2) + d*(d-1)/2)] <- corrmat[lower.tri(corrmat)==TRUE]
}

resid_vec_ourmod<-rep(NA, mcmcsamps*n*npredsamp)
for(i in 1:n){
  for(j in (1 + (1-i)*mcmcsamps):(i*mcmcsamps)){
    print(i)
    theta <- paramat[j,] # 1 sample from posterior
    
    corrmat <- diag(1, nrow=d)
    corrmat[lower.tri(corrmat)==TRUE] <- theta[(2*d + 1):para_len]
    corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
    Sigma <- diag(theta[(d+1):(2*d)]) %*% corrmat %*% diag(theta[(d+1):(2*d)])
    
    Sigma_11 <- Sigma[1,1]
    Sigma_12 <- Sigma[1, 2:d]
    Sigma_21 <- Sigma[2:d, 1]
    Sigma_22_inv <- solve(Sigma[2:d, 2:d])
    
    mu_1 <- theta[1]
    mu_2 <- theta[2:d]
    
    condmean<-as.numeric(mu_1 + Sigma_12 %*% Sigma_22_inv %*% (mth_obs[i,2:d] - mu_2))
    condvar<-as.numeric(Sigma_11 - Sigma_12 %*% Sigma_22_inv %*%  Sigma_21)
    
    # 1 sample of predicted y given x_1, ... on observation not included in model
    predsim <- rnorm(1, mean = condmean, sd <- sqrt(condvar)) 
    resid_vec_ourmod[j] <- (mth_obs[i,1] - predsim)^2
  }
}

#1.1.4 predictive simulations
npredsamp <- 1
predsamps <- matrix(NA, nrow=(mcmcsamps*n*npredsamp), ncol=d)

for(i in 1:(mcmcsamps*n)){
  print(i)
  theta <- paramat[i,]
  
  corrmat <- diag(1, nrow=d)
  corrmat[lower.tri(corrmat)==TRUE] <- theta[(2*d + 1):para_len]
  corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
  
  Sigma <- diag(theta[(d+1):(2*d)]) %*% corrmat %*% diag(theta[(d+1):(2*d)])
  
  a <- rmvnorm(1, mean=theta[1:d], sigma=Sigma)
  predsamps[i,] <- a
}

for(i in 1:n){
  y <- Data_mat[(((i-1)*m)+1):(i*m - holdout), 1]
  X <- Data_mat[(((i-1)*m)+1):(i*m - holdout), 2:d]
  mod_obj <- blasso(X, y, thin=thinning, T=(burninit+samps))
  betas<-mod_obj$beta[(burninit:(burninit + samps)),] 
  mu_blasso<-mod_obj$mu[(burninit:(burninit + samps))]
  sigmasq_blasso <- mod_obj$s2[(burninit:(burninit + samps))]
  for(j in 1:samps){
    for(k in 1:holdout){
      pred <- rnorm(1, mean = mu_blasso[j] + mth_obs[(k+((i-1)*holdout)),2:d] %*% betas[j,], 
                    sd = sqrt(sigmasq_blasso[j]))
      resid_vec_blasso[k + (holdout*((i*j)-1))] <- (mth_obs[(k+((i-1)*holdout)),1] - pred)
      print(k + (holdout*((i*j)-1)))
    }
  }
}
