library(MCMCpack)
library(mvtnorm)
library(corpcor)
library(clusterGeneration)
library(monomvn)

# To observarsjoner fra regresjonsmodellen
x <- rbind(c(0,0),c(1,1))
# posterioritetthet, logistic prior på glogit(rho),
# uniform improper på mu1, mu2 og log(sigma_1), log(sigma_2)
target <- function(theta, x, scale) {
  logPrior <- dlogis(theta[5], scale=scale, log=TRUE)
  mu <- theta[1:2]
  sigma <- exp(theta[3:4])
  rho <- 2*plogis(theta[5]) - 1
  Sigma <- diag(sigma) %*% toeplitz(c(1,rho)) %*% diag(sigma)
  logLik <- dmvnorm(x, mu, Sigma, log=TRUE)
  sum(logLik) + logPrior
}

# Generisk Metropolis-Hastings
s <- .4 # (scale til logistic prior til transformert rho)
chain <- MCMCmetrop1R(target, theta.init=rep(0,5), x=x, mcmc=1e+6, thin=100, scale=s)
# summary stats, merk at mange parametere ikke har endelig forventning
summary(chain)
# Tilbaketransformasjon av posteriori samples av rho
rho <- 2*plogis(chain[,5])-1
hist(rho, breaks=100, prob=TRUE)
# prior til rho
curve(2*((1-rho)/(rho+1))^(1/s)/(s*(1+((1-rho)/(rho+1))^(1/s))^2*(1-rho^2)), xname="rho", add=TRUE)
# posteriorifordeling til abs(rho)
hist(abs(rho),breaks=100, prob=TRUE)
# prior til abs(rho)
curve(4*((1-rho)/(rho+1))^(1/s)/(s*(1+((1-rho)/(rho+1))^(1/s))^2*(1-rho^2)), xname="rho", add=TRUE)
# trace-plot
plot(chain)
predictivesample <- function(theta) {
  mu <- theta[1:2]
  sigma <- exp(theta[3:4])
  rho <- 2*plogis(theta[5]) - 1
  Sigma <- diag(sigma) %*% toeplitz(c(1,rho)) %*% diag(sigma)
  c(x=rmvnorm(1, mu, Sigma), beta=c(mu[2] - Sigma[2,1]/Sigma[1,1]*mu[1],Sigma[2,1]/Sigma[1,1]))
}
xpredict <- t(apply(chain, 1, predictivesample))
# Predictiv tetthet til x_3
plot(xpredict[,1:2], pch=".", xlim=c(-14,15), ylim=c(-14,15))
# Observasjonene
points(x[,1],x[,2],col="red",pch=16)
# Sjekk av konsistens
mean(xpredict[,1]<0)
mean(xpredict[,1]<1)
mean(xpredict[,2]<0)
mean(xpredict[,2]<1)
hist.trunc <- function(x, trunc=20, by=1) {
  hist(x, breaks=seq(0, max(x)+by, by=by), xlim=c(0,trunc), prob=TRUE)
}
# posteriori til sigma_1^2
hist.trunc(exp(chain[,3])^2, 10, .1)
# lik invers gamma som forventet
curve(dinvgamma(x, shape=1/2, scale=1/4), add=TRUE)


improved_target_dens<-function(theta,x){
  Samp_cov <- cov(x)
  d <- -1/2 + sqrt(1/4 + 2 * length(theta))
  sigma <- theta[1:d]
  parcorrs <- theta[(d+1):length(theta)]
  if(any(sigma<=0)){-Inf}
  else{
    parcorrmat <- diag(-1, nrow=d)
    parcorrmat[lower.tri(parcorrmat)==TRUE] <- parcorrs
    parcorrmat <- parcorrmat + t(parcorrmat) - diag(diag(parcorrmat))
    
    if ((sum(eigen(parcorrmat)$val<0)==d)){ # if partial corr mat is negative definite
      S <- - parcorrmat
      S_inv <- solve(S)
      #calculating precision matrix
      Sigma_inv <- diag(1/sigma) %*% diag(sqrt(diag(S_inv))) %*% S %*% diag(sqrt(diag(S_inv))) %*% diag(1/sigma)
      -((n-1)/2) * (log(1/det(Sigma_inv)) + sum(diag(Sigma_inv %*% Samp_cov))) - sum(log(sigma))
    }
    else{-Inf}
  }
}


corr_from_pcor<-function(n, d){
  Sigma <- - diag(d)
  replicate(n,{
    repeat{
      u <- runif(d*(d-1)/2, -1, 1)
      Sigma[lower.tri(Sigma)] <- u #partial corrs
      Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
      if ((sum(eigen(Sigma)$val<0)==d)) 
        break()
    }
    # Lemma 2 from Artner 2022 space of partial correlation matrices to convert to correlation matrix
    S <- - Sigma
    S_inv <- solve(S)
    D_S_inv <- solve(diag(sqrt(diag(S_inv))))
    (D_S_inv %*% S_inv %*% D_S_inv)[lower.tri(S_inv)==TRUE]
  })
}

metrop_samp <- function(n, m, para_len, init, Data_mat, mcmcsamps, target_dens, burn=500, holdout=FALSE){
  paramat <- matrix(NA, nrow=n*mcmcsamps, ncol=para_len)
  if(holdout==FALSE){t <- m} #use all data to sample from posterior
  else{t <- m-1} #do not use mth observation
  for(i in 1:n){
    block <- Data_mat[((i-1)*m+1):(i*t),]
    tryCatch({
      sigma_pcor <- MCMCmetrop1R(improved_target_dens, theta.init=init, burnin = burn, x=block, mcmc=mcmcsamps)
      mu <- matrix(NA, nrow=mcmcsamps, ncol=d)
      for(j in 1:mcmcsamps){
        parcorrs <- sigma_pcor[j,(d+1):ncol(sigma_pcor)]
        sigmas<-sigma_pcor[j,1:d]
        parcorrmat <- diag(-1, nrow=d)
        parcorrmat[lower.tri(parcorrmat)==TRUE] <- parcorrs
        parcorrmat <- parcorrmat + t(parcorrmat) - diag(diag(parcorrmat))
        S <- - parcorrmat
        S_inv <- solve(S)
        D_S_inv <- solve(diag(sqrt(diag(S_inv))))
        corrmat <- (D_S_inv %*% S_inv %*% D_S_inv)
        Sigma <- diag(sigmas) %*% corrmat %*% diag(sigmas)
        mu[j,] <- rmvnorm(1, mean=colMeans(block), sigma=Sigma)
      }
      paramat[(1 + (i-1)*mcmcsamps):(i*mcmcsamps),] <- cbind(mu, sigma_pcor)
    }, error=function(e){})
  }
  return(paramat)
}


#Berger & Sun 2008, Accept-Reject bivariate normal
alg<-function(datamat, nsamples){
  x <- t(datamat)
  n <- ncol(x)
  s11 <- sum((x[1,] - mean(x[1,]))^2)
  s22 <- sum((x[2,] - mean(x[2,]))^2)
  s12 <- sum( t(x[1,] - mean(x[1,])) %*% (x[2,] - mean(x[2,])) )
  mean_1<-mean(x[1,])
  mean_2<-mean(x[2,])
  S_mat <- matrix(data=c(s11,s12,s12,s22), nrow=2)
  
  samp_rho<-c()
  samp_sigma1<-c()
  samp_sigma2<-c()
  samp_mu1<-c()
  samp_mu2<-c()
  
  z<-0
  while(z<nsamples){
    prop <- riwish(n-1, S_mat) #Uses different parametrization than in Berger & Sun (input S instead of S^-1)
    rho <- prop[1,2]/(sqrt(prop[1,1]) * sqrt(prop[2,2]))
    rej_bound <- (1 - (rho^2))^(3/2)
    m<-runif(1)
    if(m <= rej_bound){
      samp_rho<-append(samp_rho, rho)
      samp_sigma1<-append(samp_sigma1, sqrt(prop[1,1]))
      samp_sigma2<-append(samp_sigma2, sqrt(prop[2,2]))
      z <- z+1
    }
  }
  
  rhoest<-mean(samp_rho) #why am i taking the mean??
  sigma1est<-mean(samp_sigma1) #update: just set nsamp to 1
  sigma2est<-mean(samp_sigma2)
  Sigma<-matrix(data=c(sigma1est^2, rhoest*sigma1est*sigma2est, rhoest*sigma1est*sigma2est, sigma2est^2), nrow=2)
  mu<-rmvnorm(nsamples, c(mean_1, mean_2), sigma=Sigma/nsamples)
  
  return(cbind(samp_rho, samp_sigma1, samp_sigma2, mu[,1], mu[,2]))
}

predsamp<-function(theta){
  var_11<-theta[2]^2  
  var_21<-theta[1]*theta[2]*theta[3]
  var_22<-theta[3]^2
  Sigma<-matrix(data=c(var_11,var_21,var_21,var_22), nrow=2)
  rmvnorm(predsims, c(theta[4], theta[5]), sigma=Sigma)
}

n<-100000
#data<-rmvnorm(3, mean=rep(0, 2), sigma=matrix(data=c(1,0.5,0.5,1), nrow=2))
Sigma<-matrix(data=c(10, 5, 5, 100), nrow=2)
a<-sqrt(solve(matrix(data=c(Sigma[1,1], 0, 0, Sigma[2,2]), nrow=2)))
cor_mat<-a %*% Sigma %*% a
cor_mat
data<-rmvnorm(n, mean=c(5,-5), sigma=Sigma)
a<-alg(data, n) # simulates rho, sigma1, sigma2, mu1, mu2
plot(data, col="red")


hist(a[,1]) #rho
hist(a[,2], breaks=1000) #sigma1
hist(a[,3], breaks=1000) #sigma2
hist(a[,4], breaks=1000) #mu1
hist(a[,5], breaks=1000) #mu2

#Check means
mean(a[,1])
mean(a[,2])
mean(a[,3])
mean(a[,4])
mean(a[,5])

#Simulating sample from predictive distribution (not yet reparametrized)
predsims<-10000
predsamp<-function(theta){
  var_11<-theta[2]^2  
  var_21<-theta[1]*theta[2]*theta[3]
  var_22<-theta[3]^2
  Sigma<-matrix(data=c(var_11,var_21,var_21,var_22), nrow=2)
  rmvnorm(predsims, c(theta[4], theta[5]), sigma=Sigma)
}

#Prediktiv tetthet
b<-t(apply(a, 1, FUN=predsamp))
plot(b[,1], b[,2], xlab="x", ylab="y", ylim=c(-20, 20), xlim=c(-20, 20))
points(data, col="red")

sorted1<-sort(data[,1], decreasing=TRUE)
sorted2<-sort(data[,2], decreasing=TRUE)

mean(ifelse(b[,1]>sorted1[1], 1, 0))
mean(ifelse(b[,1]>sorted1[2], 1, 0))
mean(ifelse(b[,1]>sorted1[3], 1, 0))
mean(ifelse(b[,1]>sorted1[4], 1, 0))

mean(ifelse(b[,2]>sorted2[1], 1, 0))
mean(ifelse(b[,2]>sorted2[2], 1, 0))
mean(ifelse(b[,2]>sorted2[3], 1, 0))
mean(ifelse(b[,2]>sorted2[4], 1, 0))

#Sample rho and other parameters
n<-10000

#holder <- rlogis(n,0,1) #prior from berger sun
#rho <- 2*plogis(holder)-1
rho<-runif(n, -0.99, 0.99)

#rho<-seq(from=-0.99, to=0.99, length.out=n)
#rho<-runif(n, min=-0.0000001, max=0.0000001)

#rho<-seq(from=-0.99, to=0.99, length.out=n)
hist(rho)

mu_1 <- 15
mu_2 <- 15
sigma_1 <- 10
sigma_2 <- 1

#Sampling data

meanvec<-c(mu_1, mu_2)
m<-4
x_1<-c()
x_2<-c()

for(i in 1:n){
  print(i)
  Sigma <- matrix(data=c(sigma_1^2,sigma_1*sigma_2*rho[i],sigma_1*sigma_2*rho[i],sigma_2^2), nrow=2)
  a <- rmvnorm(m, mean=meanvec, sigma=Sigma)
  x_1 <- append(a[,1], x_1)
  x_2 <- append(a[,2], x_2)
}

x <- cbind(x_1, x_2)

a_1<-0
a_2<-0
a_3<-0
a_4<-0
a_5<-0
b_1<-0
b_2<-0
b_3<-0
b_4<-0
b_5<-0

parasims<-10
predsims<-1
it<-floor(n/m)

for (i in 1:(it)){
  newdata<-x[((i-1)*m+1):(i*m),] #Block of m of the data for each iteration
  wadup<-alg(newdata,parasims) #simulating parameters for this block
  sims<-apply(as.matrix(wadup), 1, FUN=predsamp) #simulating predictive samples for simulated parameters
  
  sorted_1<-sort(newdata[,1], decreasing=TRUE)
  sorted_2<-sort(newdata[,2], decreasing=TRUE)
  
  a_1 <- a_1 + sum(ifelse(sims[1:predsims,]>sorted_1[1], 1, 0))
  a_2 <- a_2 + sum(ifelse(sims[1:predsims,]>sorted_1[2], 1, 0))
  a_3 <- a_3 + sum(ifelse(sims[1:predsims,]>sorted_1[3], 1, 0))
  a_4 <- a_4 + sum(ifelse(sims[1:predsims,]>sorted_1[4], 1, 0))
  a_5 <- a_5 + sum(ifelse(sims[1:predsims,]>sorted_1[5], 1, 0))
  
  b_1 <- b_1 + sum(ifelse(sims[(predsims+1):(predsims*2),]>sorted_2[1], 1, 0))
  b_2 <- b_2 + sum(ifelse(sims[(predsims+1):(predsims*2),]>sorted_2[2], 1, 0))
  b_3 <- b_3 + sum(ifelse(sims[(predsims+1):(predsims*2),]>sorted_2[3], 1, 0))
  b_4 <- b_4 + sum(ifelse(sims[(predsims+1):(predsims*2),]>sorted_2[4], 1, 0))
  b_5 <- b_5 + sum(ifelse(sims[(predsims+1):(predsims*2),]>sorted_2[5], 1, 0))
  
  print(i)
}

tot_comb<-it*parasims*predsims
c(a_1, a_2, a_3, a_4, a_5)/tot_comb #amount of x_pred_n+1 above each observed x_i up to n
c(b_1, b_2, b_3, b_4, b_5)/tot_comb #amount of y_pred_n+1 above each observed y_i up to n
seq(from=1/(m+1), to=m/(m+1), by=1/(m+1)) #Expected under t dist

#General check of predictive distribution for >= 3 dimensions

set.seed(1)

#1: Sample priors
n <- 100 #samples
d <- 3 #dimension

muvec<-rep(0,d)
sigmavec<-rep(1,d)

#Sampling n realizations of correlations with uniform marginals
corrs <- corr_from_pcor(n,d)


#2 Sample Data given priors

m <- 5 #how many datapoints we use to estimate the posterior sample
Data_mat<-matrix(0, nrow=(n*m), ncol=d)


for(i in 1:n){
  corrmat <- diag(1, nrow=d)
  corrmat[lower.tri(corrmat)==TRUE] <- corrs[,i]
  corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
  Sigma <- diag(sigmavec) %*% corrmat %*% diag(sigmavec)
  #Output
  a <- rmvnorm(m, mean=muvec, sigma=Sigma)
  Data_mat[(((i-1)*m)+1):(i*m),] <- a
}

#3 Estimate parameters

init <- c(rep(1,d), rep(0, d*(d-1)/2))
para_len <- length(init) + d
mcmcsamps <- 10

#3.1 Diagnostics (compute ESS, trace plots, etc... for different blocks)

diagnostic_obj_3 <- MCMCmetrop1R(improved_target_dens, theta.init=init, 
                               burnin = 1000, x=Data_mat[7:9,], mcmc=100000)
plot(diagnostic_obj_3) #trace plots not convincing for correlations with m=d

diagnostic_obj_4 <- MCMCmetrop1R(improved_target_dens, theta.init=init, 
                                 burnin = 2000, x=Data_mat[1:4,], mcmc=100000)

plot(diagnostic_obj_4) #trace plots very convincing for m >= d + 1

diagnostic_obj_50 <- MCMCmetrop1R(improved_target_dens, theta.init=init, burnin = 2000,
                                  x=rmvnorm(50, mean=muvec, sigma=diag(sigmavec)), mcmc=100000)

plot(diagnostic_obj_50)

#3.2 

paramat <- metrop_samp(n, m, para_len, init, Data_mat, mcmcsamps, improved_target_dens)

#4: simulate predictive samples given samples from posterior

mu_1<-paramat[,1]
n_missing_rows<-sum(ifelse(is.na(mu_1[is.na(mu_1)]), 1, 0))
n <- n - n_missing_rows/mcmcsamps
paramat<-paramat[complete.cases(paramat),]

#4.1 Transforming partial corrs to corrs

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

#5: Counting amount of times above observations

countmat<-matrix(0, nrow=m, ncol=d)

for(i in 1:n){
  obs<-Data_mat[(((i-1)*m)+1):(i*m),] #block of data from step 2
  for(j in 1:d){
    orderstat<-sort(obs[,j], decreasing=TRUE) #empirical order statistics for jth dimension
    for(k in 1:m){
      datablock<-predsamps[((i-1)*mcmcsamps + 1): (i*mcmcsamps),j]
      countmat[k,j]<- countmat[k,j] + sum(ifelse(datablock<orderstat[k], 1, 0)) #add count of simulations above observed value
    }
  }
}

probmat<-countmat/(n*mcmcsamps)
probmat

rowMeans(probmat)

#1: Comparison Bayesian Lasso and reparametrized model


# From our prior

set.seed(1)

#1.1.1: Sample n priors
n<-100
d<-3 #dimension

muvec<-rep(0,d)
sigmavec<-rep(1,d)

corrs <- corr_from_pcor(n,d)

#1.1.2 Sample Data given priors

m <- 4 #how many datapoints per iteration.
#We use m-1 observations to fit the model and then compare
#sample from predicted (from model) with remaining obs

Data_mat<-matrix(0, nrow=(n*m), ncol=d)
mth_obs<-matrix(NA, nrow = n, ncol=d) #holdout observations

for(i in 1:n){
  corrmat <- diag(1, nrow=d)
  corrmat[lower.tri(corrmat)==TRUE] <- corrs[,i]
  corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
  Sigma <- diag(sigmavec) %*% corrmat %*% diag(sigmavec)
  
  Data_mat[(((i-1)*m)+1):(i*m),] <- rmvnorm(m, mean=muvec, sigma=Sigma)
  mth_obs[i,]<-Data_mat[(i*m),]
}


#1.1.3 posterior sampling

init <- c(rep(1,d), rep(0, d*(d-1)/2)) #initialization for c(sigma, parcorr) in MCMC
para_len <- length(init) + d
mcmcsamps <- 10
burnin_mcmc <- 2000

paramat <- metrop_samp(n, m, para_len, init, Data_mat, mcmcsamps, 
                       improved_target_dens, burn=burnin_mcmc, holdout=TRUE)


mu_1<-paramat[,1]
#removing crashed samples
cond<-is.na(mu_1[seq(1, n*mcmcsamps, by=mcmcsamps)])
mth_obs<-mth_obs[!cond,]

n_missing_rows <- sum(ifelse(cond, 1, 0))
n <- n - n_missing_rows
paramat <- paramat[complete.cases(paramat),]


# Transforming partial corrs to corrs

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

#1.1.4 predictive simulations
npredsamp <- 1
predsamps <- matrix(NA, nrow=(mcmcsamps*n*npredsamp), ncol=d)

for(i in 1:(mcmcsamps*n)){
  theta <- paramat[i,]
  
  corrmat <- diag(1, nrow=d)
  corrmat[lower.tri(corrmat)==TRUE] <- theta[(2*d + 1):para_len]
  corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
  
  Sigma <- diag(theta[(d+1):(2*d)]) %*% corrmat %*% diag(theta[(d+1):(2*d)])
  
  a <- rmvnorm(1, mean=theta[1:d], sigma=Sigma)
  predsamps[i,] <- a
}

#1.1.5 residuals in 1 dimension

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
    predsim <- rnorm(1, mean = condmean, sd <- condvar) 
    resid_vec_ourmod[j] <- (mth_obs[i,1] - predsim)^2
  }
}

mean(resid_vec_ourmod[complete.cases(resid_vec_ourmod)])

#hist
hist(paramat[,7])
hist(paramat[,8])
hist(paramat[,9])

mean(paramat[,7])
mean(paramat[,8])
mean(paramat[,9])

#1.1.6 Bayesian lasso

resid_vec_blasso<-rep(NA, n*npredsamp)
burninit<-10000
for(i in 1:n){
  y <- Data_mat[(((i-1)*m)+1):(i*m - 1), 1]
  X <- Data_mat[(((i-1)*m)+1):(i*m - 1), 2:d]
  mod_obj <- blasso(X, y, T=burninit)
  betas<-mod_obj$beta[burninit,] 
  mu_blasso<-mod_obj$mu[burninit]
  sigmasq_blasso <- mod_obj$s2[burninit]
  pred <- rnorm(1, mean = mu_blasso + mth_obs[i,2:d] %*% betas, 
                sd = sqrt(sigmasq_blasso))
  resid_vec_blasso[i] <- (mth_obs[i,1] - pred)
}

mean(resid_vec_blasso^2)

# Generating data from bayesian lasso and then comparing

#1.2.1 Generating priors

set.seed(2)

n <- 500
d <- 3
p <- (d-1) #for notational purposes, denote vector (y, x_1 , ..., x_p)

r <- 1
delta <- 1 
lambda_sq <- rgamma(n, shape=delta, rate=r)

sigma_sq <- 1/rgamma(n, shape=3, scale=3) 
#inverse gamma to maintain conjugacy as stated in article
tau_sq <- matrix(NA, nrow=n, ncol=p)
beta_samp <- matrix(NA, nrow = n, ncol = p)


for(i in 1:n){
  for(j in 1:(d-1)){
    tau_sq[i,j] <- rexp(1, rate = lambda_sq[i]/2)
  }
}

for(i in 1:n){
  beta_samp[i,] <- rmvnorm(1, rep(0, p), sigma_sq[i] * diag(tau_sq[i,]))
}

# Generating X
m <- 10

muvec<-rep(0,p)
sigmavec<-rep(1,p)

corrs <- corr_from_pcor(n, d)

fulldata <- matrix(NA, nrow=(m*n), ncol = p)
mth_obs <- matrix(NA, nrow=n, ncol = d)

for(i in 1:n){
  
  corrmat <- diag(1, nrow=p)
  corrmat[lower.tri(corrmat)==TRUE] <- if(p==2){corrs[i]}else{corrs[i,]}
  corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
  Sigma <- diag(sigmavec) %*% corrmat %*% diag(sigmavec)
  fulldata[(1 + ((i-1)*m)):(i*m),] <- rmvnorm(m, muvec, Sigma)
}

#Standardizing X, as assumed for B. Lasso
for(i in 1:p){
  fulldata[,i] <- (fulldata[,i] - mean(fulldata[,i]))/sd(fulldata[,i])
}

#Generate y conditional on X, \Beta, \sigma^2 (\mu not includedn not sure how to and if i should)
y <- rep(NA, n*m)
for(i in 1:(n)){
  y[(1 + (i-1)*m): (i*m)] <- rmvnorm(1, mean=fulldata[(1 + ((i-1)*m)):((i*m)),] %*% beta_samp[i,], sigma_sq[i] * diag(m))
}

Data_mat <- cbind(y,fulldata)

for(i in 1:n){
  mth_obs[i,] <- Data_mat[(i*m),]
}

#1.2.2 Bayesian lasso

resid_vec_blasso<-rep(NA, n)
burnin_blasso<-5000
for(i in 1:n){
  y <- Data_mat[(((i-1)*m)+1):(i*m - 1), 1]
  X <- Data_mat[(((i-1)*m)+1):(i*m - 1), 2:d]
  mod_obj <- blasso(X, y, T=burnin_blasso)
  betas<-mod_obj$beta[burnin_blasso,] #accepting only last sample of posterior
  mu_blasso<-mod_obj$mu[burnin_blasso]
  sigmasq_blasso <- mod_obj$s2[burnin_blasso]
  pred <- rnorm(1, mean = mu_blasso + mth_obs[i,2:d] %*% betas, 
                sd = sqrt(sigmasq_blasso))
  resid_vec_blasso[i] <- (mth_obs[i,1] - pred)
}

mean(resid_vec_blasso^2)



#1.2.3 Our model

init <- c(rep(1,d), rep(0, d*(d-1)/2))
para_len <- length(init) + d
mcmcsamps <- 10
burnin_mcmc <- 2000

paramat <- metrop_samp(n, m, para_len, init, Data_mat, mcmcsamps, 
                       improved_target_dens, burn=burnin_mcmc, holdout=TRUE)


mu_1<-paramat[,1]
#removing crashed samples
cond<-is.na(mu_1[seq(1, n*mcmcsamps, by=mcmcsamps)])
mth_obs<-mth_obs[!cond,]

n_missing_rows <- sum(ifelse(cond, 1, 0))
n <- n - n_missing_rows
paramat <- paramat[complete.cases(paramat),]

# Transforming partial corrs to corrs

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


npredsamp <- 1
predsamps <- matrix(NA, nrow=(mcmcsamps*n*npredsamp), ncol=d)

for(i in 1:(mcmcsamps*n)){
  theta <- paramat[i,]
  
  corrmat <- diag(1, nrow=d)
  corrmat[lower.tri(corrmat)==TRUE] <- theta[(2*d + 1):para_len]
  corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
  
  Sigma <- diag(theta[(d+1):(2*d)]) %*% corrmat %*% diag(theta[(d+1):(2*d)])
  
  predsamps[i,] <- rmvnorm(1, mean=theta[1:d], sigma=Sigma)
}

#1.2.4 residuals in 1 dimension

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
    resid_vec_ourmod[j] <- (mth_obs[i,1] - predsim)
  }
}

resid_vec_ourmod_cleaned <- resid_vec_ourmod[complete.cases(resid_vec_ourmod)]
mean(resid_vec_ourmod_cleaned^2)



par(mfrow=c(3,3))
hist(paramat[,1], breaks=100, main="mu_1")
hist(paramat[,2], breaks=100, main="mu_2")
hist(paramat[,3], breaks=100, main="mu_3")
hist(paramat[,4], breaks=100, main="sigma_1")
hist(paramat[,5], breaks=100, main="sigma_2")
hist(paramat[,6], breaks=100, main="sigma_3")
hist(paramat[,7], breaks=100, main="rho_12")
hist(paramat[,8], breaks=100, main="rho_13")
hist(paramat[,9], breaks=100, main="rho_23")

colMeans(paramat)

