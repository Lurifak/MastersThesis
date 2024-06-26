#Loading packages
library(MCMCpack)
library(mvtnorm)
library(corpcor)
library(clusterGeneration)
library(monomvn)
library(coda)
library(ggplot2)
library(ggalt)
library(latex2exp)
library(gridExtra)
library(datasets)



#Function for adjusting tuning parameter in MCMCpack as a function of acceptance probability
tune_adjust<-function(a_prob){
  if(a_prob>0.25){
    4*(a_prob-0.25) + 1
  } else{
      3 * (a_prob) +  1/4
  }
}


# Expression proportional to posterior of covariance matrix in eq. 3.3
improved_target_dens <- function(theta,x){
  L <- length(theta)
  Samp_cov <- cov(x)
  n <- nrow(x)
  d <- ncol(x)
  sigma <- theta[1:d]
  parcorrs <- theta[(d+1):L]
  if(any(sigma<=0)){
    -Inf
  } else {
    parcorrmat <- diag(-1, nrow=d)
    parcorrmat[lower.tri(parcorrmat)==TRUE] <- parcorrs
    parcorrmat <- parcorrmat + t(parcorrmat) - diag(diag(parcorrmat))
    
    if((sum(eigen(parcorrmat)$val<0)==d)){ # if partial corr mat is negative definite => Corr. mat positive def
      
      Sigma_inv <- solve(diag(sigma) %*% pcor2cor(parcorrmat+diag(2, nrow=d)) %*% diag(sigma)) #pcor2cor defines p. corr mat having unit diagonal instead of negative unit diag
      -((n-1)/2) * (log(1/det(Sigma_inv)) + sum(diag(Sigma_inv %*% Samp_cov))) - sum(log(sigma))
    } else {
      -Inf
    }
  }
}

#Function for transforming MVN variables to Regression coefficients

#PMAt = parameter matrix, d = dim of mvn vector
#Assumed rows in PMat: mu_s, sigmas and corrs
betafrommvn <- function(PMat, d){
  n <- nrow(PMat)
  len <- ncol(PMat)
  betas <- matrix(NA, nrow=n, ncol=d)
  for(i in 1:n){
    theta <- PMat[i,]
  
    corrmat <- diag(1, nrow=d)
    corrmat[lower.tri(corrmat)==TRUE] <- theta[(2*d + 1):len]
    corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
    Sigma <- diag(theta[(d+1):(2*d)]) %*% corrmat %*% diag(theta[(d+1):(2*d)])
    
    Sigma_12 <- Sigma[1, 2:d]
    Sigma_22_inv <- solve(Sigma[2:d, 2:d])
    betas[i,2:d] <- Sigma_12 %*% Sigma_22_inv #slopes
    betas[i,1] <- theta[1] - betas[i,2:d] %*% theta[2:d] #intercept
  }
  betas
}


#Generating correlations from uniform marginals on partial correlation matrix
corr_from_pcor<-function(n, d){
  Sigma <- - diag(d)
  replicate(n,{
    repeat{
      u <- runif(d*(d-1)/2, -1, 1)
      Sigma[lower.tri(Sigma)] <- u #partial corrs
      Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
      if ((sum(eigen(Sigma)$val<0)==d)){break()}
    }
    # Lemma 2 from Artner 2022 space of partial correlation matrices to convert to correlation matrix
    S <- - Sigma
    S_inv <- solve(S)
    D_S_inv <- diag(1/sqrt(diag(S_inv))) # D^-1 = diag (1 / diag(D)) when D diagonal matrix
    (D_S_inv %*% S_inv %*% D_S_inv)[lower.tri(S_inv)==TRUE]
  })
}


#Sampling from full posterior of posterior in Eq. 3.4 given data
#Variables: n is amount of chains
# m is amount of observed values per chain
# Data_mat is the matrix of the response and the design matrix
# holdout is amount of samples to not use to fit the model 
metrop_samp <- function(n, m, para_len, Data_mat, mcmcsamps, target_dens, burn=500, holdout=0){
  paramat <- matrix(NA, nrow=n*mcmcsamps, ncol=para_len)
  if(holdout==0){
    t <- m
  } else {
    t <- m - holdout
  }
  
  pattern <- "[0-9]+\\.?[0-9]*" #for extracting acceptance prob
  
  tot_acc_post <- 0
  count_post <- 0
  
  ESS_mat_holder <- matrix(0, nrow=n, ncol=((d) + (d)*(d-1)/2))
  
  for(i in 1:n){
    block <- Data_mat[((i-1)*m+1):(((i-1)*m)+t),]
    
    #Determining initialization cheaply
    covmat <- cov(block)
    pcorrmat <- cor2pcor(cov2cor(covmat))
    init <- c(sqrt(diag(covmat)), pcorrmat[lower.tri(pcorrmat)==TRUE])
    
    #For estimated correlations almost on edge of parameter space it can crash - thus skip these
    
    tuning <- 1
    acc_prob_pre <- 0
    ess_test <- rep(1, (d) + (d)*(d-1)/2)
    
    #test-chain
    #sink captures acceptance rate string printed
    
    tryCatch({
      
      sink(file="test.txt")
      
      testchain <- MCMCmetrop1R(improved_target_dens, 
                                theta.init=init, burnin = 500, x=block, mcmc=250)
      
      out <- readLines("test.txt")
      
      acc_prob_pre <- regmatches(out, regexpr(pattern, out))
      
      acc_prob_pre <- as.numeric(acc_prob_pre)
      
      tuning <- tune_adjust(acc_prob_pre)
      
    }, error=function(e){})
    
    sink()
    
    #full chain
    
    tryCatch({
      
      sink(file="test.txt")
      
      sigma_pcor <- MCMCmetrop1R(improved_target_dens, theta.init=init, 
                                 burnin = (burn), x=block, mcmc=mcmcsamps, tune=tuning)
      
      out <- readLines("test.txt")
      
      acc_prob_post <- regmatches(out, regexpr(pattern, out))
      
      acc_prob_post <- as.numeric(acc_prob_post)
      
      tot_acc_post <- tot_acc_post + acc_prob_post
      
      count_post <- count_post + 1
      
      ESS_mat_holder[i, ] <- effectiveSize(sigma_pcor)
      
      mu <- matrix(NA, nrow=mcmcsamps, ncol=d)
      
      for(j in 1:mcmcsamps){
        parcorrs <- sigma_pcor[j,(d+1):ncol(sigma_pcor)]
        sigmas <- sigma_pcor[j,1:d]
        parcorrmat <- diag(1, nrow=d) #pcor2cor assumes positive unit diagonal in partial correlation matrix
        parcorrmat[lower.tri(parcorrmat)==TRUE] <- parcorrs
        parcorrmat <- parcorrmat + t(parcorrmat) - diag(diag(parcorrmat))
        corrmat <- pcor2cor(parcorrmat)
        Sigma <- (1/m) * diag(sigmas) %*% corrmat %*% diag(sigmas)
        mu[j,] <- rmvnorm(1, mean=(colMeans(block)), sigma = (Sigma))
      }
      paramat[(1 + (i-1)*mcmcsamps):(i*mcmcsamps),] <- cbind(mu, sigma_pcor)
    }, error=function(e){})
    sink()
    cat(count_post, "of", i, " runs accepted \n")
  }
  
  ESS_mat <<- ESS_mat_holder
  accrate <<- tot_acc_post/count_post
  cat("Average acceptance post tuning probability was ", tot_acc_post/count_post, "\n")
  cat("Count Post was ", count_post, "\n")
  
  return(paramat)
}

#transforms the partial correlation columns to correlations in a matrix
pcors_to_corrs<-function(paramat, d){
  n <- nrow(paramat)
  for(i in 1:n){
    theta <- paramat[i,]
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
  return(paramat)
}


#Berger & Sun 2008, Accept-Reject bivariate normal described in their article
Sun_reject<-function(datamat, nsamples){
  x <- t(datamat)
  n <- ncol(x)
  s11 <- sum((x[1,] - mean(x[1,]))^2)
  s22 <- sum((x[2,] - mean(x[2,]))^2)
  s12 <- sum( t(x[1,] - mean(x[1,])) %*% (x[2,] - mean(x[2,])) )
  mean_1<-mean(x[1,])
  mean_2<-mean(x[2,])
  S_mat <- matrix(data=c(s11,s12,s12,s22), nrow=2) #Sample covariance
  
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
  
  rhoest<-samp_rho
  sigma1est<-samp_sigma1
  sigma2est<-samp_sigma2
  Sigma<-matrix(data=c(sigma1est^2, rhoest*sigma1est*sigma2est, rhoest*sigma1est*sigma2est, sigma2est^2), nrow=2)
  mu<-rmvnorm(nsamples, c(mean_1, mean_2), sigma=Sigma/nsamples)
  
  return(cbind(samp_rho, samp_sigma1, samp_sigma2, mu[,1], mu[,2]))
}


#Function to sample realizations of x given parameters for the bivariate normal
xsim <- function(thetamat, nsims){
  n <- nrow(thetamat)
  x <- rep(NA, (n*nsims))
  for(i in 1:n){
    x[((i-1)*nsims + 1):(i*nsims)] <- rnorm(nsims, mean=para_sims[i,5], sd=para_sims[i,3])
  }
  x
}


#Function to sample realizations of y given parameters, x for the bivariate normal
ysim<-function(theta, x){
  n <- nrow(theta)
  L <- length(x)
  y <- rep(NA, L)
  for(i in 1:n){
    rho <- theta[i,1]
    sigmas <- theta[i,2:3]
    mu <- theta[i,4:5]
    condvar <- (1-rho^2)*sigmas[1]^2
    for(j in 1:(L/n)){
      condmean <- mu[1] + rho*(sigmas[1]/sigmas[2]) * (x[((i-1)*(L/n) + j)] - mu[2])
      y[((i-1)*(L/n) + j)] <- rnorm(1, mean=condmean, sd=sqrt(condvar))
    }
  }
  y
}



#############################################################################################
### Sampling from bivariate normal
#############################################################################################


rm(list = setdiff(ls(), lsf.str())) #Clear variables

#Sample rho and other parameters
set.seed(1)

n<-100

rho<-runif(n, -0.99, 0.99)

mu_1 <- 0
mu_2 <- 0
sigma_1 <- 1
sigma_2 <- 1

#Sampling data

meanvec<-c(mu_1, mu_2)
m <- 5
holdout <- 2
t <- m + holdout #total samples
y<-rep(NA, (n*m))
x<-rep(NA, (n*m))
holdout_x <- rep(NA, (n*holdout))
holdout_y <- rep(NA, (n*holdout))

for(i in 1:n){
  Sigma <- matrix(data=c(sigma_1^2, sigma_1*sigma_2*rho[i], sigma_1*sigma_2*rho[i], sigma_2^2), nrow=2)
  a <- rmvnorm((m+holdout), mean=meanvec, sigma=Sigma)
  y[((i-1)*m + 1): ((i-1)*m + m)] <- a[1:m, 1]
  x[((i-1)*m + 1): ((i-1)*m + m)] <- a[1:m, 2]
  holdout_y[((i-1)*holdout + 1):(i*holdout)] <- a[(m+1):(m+holdout), 1]
  holdout_x[((i-1)*holdout + 1):(i*holdout)] <- a[(m+1):(m+holdout), 2]
}

z <- cbind(y, x)

count <- rep(0, m)

parasims <- 1
predsims <- 3

residvec <- rep(NA, n*holdout)

paramat <- matrix(NA, nrow=parasims*n, ncol=5)

for (i in 1:(n)){
  newdata <- z[((i-1)*m + 1): (i*m),] #Block of m of the data for each iteration
  para_sims <- Sun_reject(newdata, parasims) #simulating parameters for this block
  x_sim <- xsim(para_sims, predsims) #simulating x
  y_x <- ysim(para_sims, x_sim)  #simulating y conditional on x
  
  sorted <- sort(newdata[,1], decreasing=TRUE)
  
  for(j in 1:m){
    count[j] <- count[j] + sum(ifelse(y_x>sorted[j], 1, 0)) #counting times ysim above observed vals
  }
  
  hold <- cbind(holdout_y[((i-1)*holdout + 1):(i*holdout)], holdout_x[((i-1)*holdout + 1):(i*holdout)])
  y_sim <- ysim(para_sims, hold[,2])
  
  for(k in 1:holdout){
    residvec[(i-1)*(holdout) + k] <- y_sim[k] - hold[k,1]
  }
  
  paramat[((i-1)*parasims + 1): (i*parasims),1:2] <- para_sims[,4:5]
  paramat[((i-1)*parasims + 1): (i*parasims),3:4] <- para_sims[,2:3]
  paramat[((i-1)*parasims + 1): (i*parasims),5] <- para_sims[,1]
  
}

mean(residvec^2)

tot_comb<-n*parasims*predsims
count/tot_comb #amount of y_pred_n+1 above each observed y_i up to n
seq(from=1/(m+1), to=m/(m+1), by=1/(m+1)) #Expected under t dist

par(mfrow=c(2,3))
hist(paramat[,1], main="mu_1")
hist(paramat[,2], main="mu_2")
hist(paramat[,3], main="sigma_1", breaks=30, xlim=c(0,10))
hist(paramat[,4], main="sigma_2", breaks=30, xlim=c(0,10))
hist(paramat[,5], main="rho")

#Visualizing betas for our prior

#B_1 and B_2

n <- 3000
d <- 3
muvec<-rep(0,d)
corrs <- corr_from_pcor(n,d)

par_mat_111 <- cbind(matrix(muvec, nrow=n, ncol=d), matrix(c(1,1,1), nrow=n, ncol=d), t(corrs))
par_mat_122 <- cbind(matrix(muvec, nrow=n, ncol=d), cbind(rep(1, n), rep(2, n), rep(2,n)), t(corrs))
par_mat_123 <- cbind(matrix(muvec, nrow=n, ncol=d), cbind(rep(1, n), rep(2, n), rep(3,n)), t(corrs))
par_mat_211 <- cbind(matrix(muvec, nrow=n, ncol=d), cbind(rep(2, n), rep(1, n), rep(1,n)), t(corrs))


prior_betas_111 <- betafrommvn(par_mat_111, d)
prior_betas_122 <- betafrommvn(par_mat_122, d)
prior_betas_123 <- betafrommvn(par_mat_123, d)
prior_betas_211 <- betafrommvn(par_mat_211, d)

colnames(prior_betas_111) <- c("Beta_0", "Beta_1", "Beta_2")
colnames(prior_betas_122) <- c("Beta_0", "Beta_1", "Beta_2")
colnames(prior_betas_123) <- c("Beta_0", "Beta_1", "Beta_2")
colnames(prior_betas_211) <- c("Beta_0", "Beta_1", "Beta_2")

par(mfrow=c(2,2))
plot(prior_betas_111[,2], prior_betas_111[,3], cex=0.1, ylim=c(-3,3), xlim=c(-3,3), xlab=TeX(r"($\Beta_1)"), ylab=TeX(r"($\Beta_2)"), main = expression(sigma == "[1, 1, 1]"))
plot(prior_betas_122[,2], prior_betas_122[,3], cex=0.1, ylim=c(-3,3), xlim=c(-3,3), xlab=TeX(r"($\Beta_1)"), ylab=TeX(r"($\Beta_2)"), main = expression(sigma == "[1, 2, 2]"))
plot(prior_betas_123[,2], prior_betas_123[,3], cex=0.1, ylim=c(-3,3), xlim=c(-3,3), xlab=TeX(r"($\Beta_1)"), ylab=TeX(r"($\Beta_2)"), main = expression(sigma == "[1, 2, 3]"))
plot(prior_betas_211[,2], prior_betas_211[,3], cex=0.1, ylim=c(-3,3), xlim=c(-3,3), xlab=TeX(r"($\Beta_1)"), ylab=TeX(r"($\Beta_2)"), main = expression(sigma == "[2, 1, 1]"))


#############################################################################################
### General check of consistency for predictive distribution for 3 or more dimensions
#############################################################################################


rm(list = setdiff(ls(), lsf.str()))

set.seed(1)

#1: Sample priors
n <- 1000 #samples
d <- 4 #dimension

muvec<-rep(0,d)
sigmavec<-rep(1,d)

#Sampling n realizations of correlations with uniform marginals
corrs <- corr_from_pcor(n,d)

#2 Sample Data given priors

m <- 20 #how many datapoints we use to estimate the posterior sample
Data_mat<-matrix(NA, nrow=(n*m), ncol=d)


for(i in 1:n){
  corrmat <- diag(1, nrow=d)
  corrmat[lower.tri(corrmat)==TRUE] <- corrs[,i]
  corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
  Sigma <- diag(sigmavec) %*% corrmat %*% diag(sigmavec)
  Data_mat[(((i-1)*m)+1):(i*m),] <- rmvnorm(m, mean=muvec, sigma=Sigma)
}

#3 Estimate parameters

para_len <- 2*d + (d*(d-1))/2
mcmcsamps <- 1000
paramat_pcor <- metrop_samp(n, m, para_len, Data_mat, mcmcsamps, improved_target_dens, 
                            burn=1000)


ESS_mat_red <- matrix(ESS_mat[ESS_mat!=0], ncol=(para_len-d)) #Removing NA runs
ESS_mat_red
colMeans(ESS_mat_red)
a <- sum(ESS_mat_red)/(ncol(ESS_mat_red)*nrow(ESS_mat_red)) #average ESS total
a

#4: simulate predictive samples given samples from posterior


#4 Remove NA rows and redefine parameters
mu_1<-paramat_pcor[,1]
n_missing_rows<-sum(ifelse(is.na(mu_1[is.na(mu_1)]), 1, 0))
n <- n - n_missing_rows/mcmcsamps
paramat_pcor<-paramat_pcor[complete.cases(paramat_pcor),]

#4.1 Transforming partial corrs to corrs

paramat <- pcors_to_corrs(paramat_pcor, d)

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

countmat <- matrix(0, (n*m), ncol=d)
for(i in 1:n){
  print(i)
  obs <- Data_mat[(((i-1)*m)+1):(i*m),] #block of data from step 2
  for(j in 1:d){
    orderstat <- sort(obs[,j], decreasing=TRUE) #empirical order statistics for jth dimension
    for(k in 1:m){
      datablock <- predsamps[((i-1)*mcmcsamps + 1): (i*mcmcsamps),j]
      countmat[(i-1)*m + k,j] <- sum(ifelse(datablock<orderstat[k], 1, 0)) #add count of simulations above observed value
    }
  }
}

side_by_side <- matrix(NA, nrow=m, ncol=(d*n))
for(i in 1:m){
  for(j in 1:n){
    side_by_side[i,((j-1)*d+1):(j*d)] <- countmat[i + (j-1)*m,]
  }
}

sd_est <- apply(side_by_side/n, 1, FUN=sd)
sd_est

seq(from=m, to=1)/(m+1) #expected
rowMeans(side_by_side)/mcmcsamps #realization


#############################################################################################
### Comparing MSPE with data from MVN priors
#############################################################################################

rm(list = setdiff(ls(), lsf.str())) #removes all variables except functions

set.seed(1)

#1.1.1: Sample n priors
n <- 10000 #amount of samples to run
d <- 4 #dimension

muvec<-rep(0,d)
sigmavec<-rep(1,d)

corrs <- corr_from_pcor(n,d)

#1.1.2 Sample Data given priors

m <- 60 #how many datapoints per iteration.
holdout <- 20
#We use m-holdout observations to fit the model and then compare
#sample from predicted (from model) with remaining observations

Data_mat<-matrix(NA, nrow=(n*m), ncol=d)
mth_obs<-matrix(NA, nrow = (n*holdout), ncol=d) #holdout observations

#Generate data
for(i in 1:n){
  corrmat <- diag(1, nrow=d)
  corrmat[lower.tri(corrmat)==TRUE] <- corrs[,i]
  corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
  Sigma <- diag(sigmavec) %*% corrmat %*% diag(sigmavec)
  Data_mat[(((i-1)*m)+1):(i*m),] <- rmvnorm(m, mean=muvec, sigma=Sigma)
}

#Standardizing, as assumed for B. Lasso
for(i in 1:n){
  for(j in 1:d){
    Data_mat[((i-1)*m+1):(i*m),j] <- (Data_mat[((i-1)*m+1):(i*m),j] - mean(Data_mat[((i-1)*m+1):(i*m),j]))/sd(Data_mat[((i-1)*m+1):(i*m),j])
  }
}

for(i in 1:n){
  mth_obs[((i-1)*holdout+1):(i*holdout),]<-Data_mat[((i*m) - (holdout-1)):(i*m),]
}

#1.1.3 posterior sampling

para_len <- 2*d + (d*(d-1)/2)
mcmcsamps <- 1000
burnin_mcmc <- 1000

paramat_pcor <- metrop_samp(n, m, para_len, Data_mat, mcmcsamps, 
                       improved_target_dens, burn=burnin_mcmc, holdout=holdout)


ESS_mat_red <- matrix(ESS_mat[ESS_mat!=0], ncol=(para_len-d))
colMeans(ESS_mat_red)
sum(ESS_mat_red)/(ncol(ESS_mat_red)*nrow(ESS_mat_red)) #average ESS per parameter



#removing crashed samples
mu_1<-paramat_pcor[,1]
cond <- is.na(mu_1[seq(1, ((n-1)*mcmcsamps+1), by=mcmcsamps)])
removals<-rep(TRUE, holdout*n)
for(i in 1:n){
  if(cond[i]==TRUE){
    removals[(((i-1)*holdout)+1):(i*holdout)] <- rep(FALSE, holdout)
  }
}
mth_obs_1 <- mth_obs[removals,]
n_missing_rows <- sum(ifelse(cond, 1, 0))
n_1 <- n - n_missing_rows
paramat_pcor <- paramat_pcor[complete.cases(paramat_pcor),]


# Transforming partial corrs to corrs
paramat <- pcors_to_corrs(paramat_pcor, d)

betamat <- betafrommvn(paramat, d)

#1.1.5 residuals in 1 dimension

batch_size <- round(sqrt(mcmcsamps*holdout))
infomat <- matrix(NA, nrow = holdout*mcmcsamps*n_1, ncol=5)
colnames(infomat) <- c("Predictions", "Actual", "Residual", "Estimated SE", "Corrected MSE")

for(i in 1:n_1){
  for(j in 1:mcmcsamps){
    theta <- paramat[((i-1)*mcmcsamps + j),] # 1 sample from posterior
    
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
    condvar<-as.numeric(Sigma_11 - Sigma_12 %*% Sigma_22_inv %*%  Sigma_21)
    
    for(k in 1:holdout){
      condmean <- as.numeric(mu_1 + Sigma_12 %*% Sigma_22_inv %*% (mth_obs_1[(k+((i-1)*holdout)),2:d] - mu_2))
      
      # 1 sample of predicted y given x_1, ... on observation not included in model
      predsim <- rnorm(1, mean = condmean, sd <- sqrt(condvar))
      
      infomat[((i-1)*(mcmcsamps*holdout) + (j-1)*holdout + k),1] <- predsim
      infomat[((i-1)*(mcmcsamps*holdout) + (j-1)*holdout + k),2] <- mth_obs_1[(k+((i-1)*holdout)),1]
      infomat[(i-1)*(mcmcsamps*holdout) + (j-1)*holdout + k, 3] <- (mth_obs_1[(k+((i-1)*holdout)),1] - predsim)
    }
  }
}

placeholder <- matrix(data=NA, nrow=mcmcsamps*holdout, ncol=n)

for(i in 1:n_1){
  placeholder[1:(mcmcsamps*holdout), i] <- infomat[((i-1)*(mcmcsamps*holdout)+1):(i*mcmcsamps*holdout),1]
}

temp <- batchSE(as.mcmc(placeholder), batchSize = batch_size)

for(i in 1:n_1){
  infomat[((i-1)*(mcmcsamps*holdout)+1):(i*mcmcsamps*holdout),4] <- rep(temp[i], (mcmcsamps*holdout))
}

infomat[,5] <- infomat[,3]^2 - infomat[,4]^2

#1.1.6 Bayesian lasso

burninit<-1000
samps<-1000
batch_size <- round(sqrt(samps*holdout))
thinning <- NULL
betamat_b <- matrix(NA, nrow=n*samps, ncol=d)
infomat_b <- matrix(NA, nrow = holdout*samps*n, ncol=5)
colnames(infomat_b) <- c("Predictions", "Actual", "Residual", "Estimated SE", "Corrected MSE")


for(i in 1:n){
  y <- Data_mat[(((i-1)*m)+1):(i*m - holdout), 1]
  X <- Data_mat[(((i-1)*m)+1):(i*m - holdout), 2:d]
  mod_obj <- blasso(X, y, thin=thinning, T=(burninit+samps), normalize=FALSE)
  betas<-mod_obj$beta[(burninit:(burninit + samps - 1)),] 
  mu_blasso<-mod_obj$mu[(burninit:(burninit + samps - 1))]
  betamat_b[((i-1)*samps+1):(i*samps),] <- cbind(mu_blasso, betas)
  sigmasq_blasso <- mod_obj$s2[(burninit:(burninit + samps))]
  for(j in 1:samps){
    for(k in 1:holdout){
      pred <- rnorm(1, mean = mu_blasso[j] + mth_obs[(k+((i-1)*holdout)),2:d] %*% betas[j,], 
                    sd = sqrt(sigmasq_blasso[j]))
      actual <- mth_obs[k + (i-1)*holdout,1]
      infomat_b[((i-1)*(samps*holdout) + (j-1)*holdout + k), 1] <- pred
      infomat_b[((i-1)*(samps*holdout) + (j-1)*holdout + k), 2] <- mth_obs[(k+((i-1)*holdout)),1]
      infomat_b[((i-1)*(samps*holdout) + (j-1)*holdout + k), 3] <- pred - actual
    }
  }
}

placeholder <- matrix(data=NA, nrow=samps*holdout, ncol=n)

for(i in 1:n){
  placeholder[1:(samps*holdout), i] <- infomat_b[((i-1)*(samps*holdout)+1):(i*samps*holdout),1]
}

temp <- batchSE(as.mcmc(placeholder), batchSize = batch_size)

for(i in 1:n){
  infomat_b[((i-1)*(samps*holdout)+1):(i*samps*holdout),4] <- rep(temp[i], (samps*holdout))
}

infomat_b[,5] <- infomat_b[,3]^2 - infomat_b[,4]^2

mean(infomat[,3]^2) #Estimated MSPE MVN
mean(infomat_b[,3]^2) #Estimated MSPE B. Lasso

accrate #acceptance rate
quantile(ESS_mat_red, 0.01) #quantile of ESS
n_1 #amount of total accepted runs

par(mfrow=c(2,1))
hist(infomat[,3], breaks=100)
hist(infomat_b[,3], breaks=100)

#############################################################################################
### Comparing MSPE with data from B. Lasso priors
#############################################################################################

#1.2.1 Generating priors

rm(list = setdiff(ls(), lsf.str())) #removes all variables except functions

set.seed(1)

n <- 10000 #Amount of mcmc chains
d <- 4
p <- (d-1)
m <- 60 #amount of total observations (including nonobserved values)
holdout <- 20

r <- 1
delta <- 1 
lambda_sq <- rgamma(n, shape=r, rate=delta)

sigma_sq <- 1 #has prior 1/sigma_sq^2 (improper) so we fix it, see notation for B. Lasso
tau_sq <- matrix(NA, nrow=n, ncol=p)
beta_samp <- matrix(NA, nrow = n, ncol = p)


for(i in 1:n){
  for(j in 1:(d-1)){
    tau_sq[i,j] <- rexp(1, rate = lambda_sq[i]/2)
  }
}

for(i in 1:n){
  beta_samp[i,] <- rmvnorm(1, rep(0, p), sigma_sq * diag(tau_sq[i,]))
}

# Generating X

muvec<-rep(0,p)
sigmavec<-rep(1,p)

corrs <- corr_from_pcor(n, p)

fulldata <- matrix(NA, nrow=(m*n), ncol = p)
mth_obs<-matrix(NA, nrow = (n*holdout), ncol=d) #holdout observations

for(i in 1:n){
  corrmat <- diag(1, nrow=p)
  corrmat[lower.tri(corrmat)==TRUE] <- if(p==2){corrs[i]}else{corrs[,i]}
  corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
  Sigma <- diag(sigmavec) %*% corrmat %*% diag(sigmavec)
  fulldata[(1 + ((i-1)*m)):(i*m),] <- rmvnorm(m, muvec, Sigma)
}

#Generate y conditional on X, \Beta, \sigma^2 (mu is improper and set to 0)
y <- rep(NA, n*m)
for(i in 1:(n)){
  print(i)
  y[(1 + (i-1)*m): (i*m)] <- rmvnorm(1, mean=fulldata[(1 + ((i-1)*m)):((i*m)),] %*% beta_samp[i,], sigma_sq * diag(m))
}

Data_mat <- cbind(y,fulldata) #Response and design matrix in one matrix

#Standardizing, as assumed for B. Lasso
for(i in 1:n){
  for(j in 1:d){
    Data_mat[((i-1)*m+1):(i*m),j] <- (Data_mat[((i-1)*m+1):(i*m),j] - mean(Data_mat[((i-1)*m+1):(i*m),j]))/sd(Data_mat[((i-1)*m+1):(i*m),j])
  }
}

for(i in 1:n){
  mth_obs[((i-1)*holdout+1):(i*holdout),]<-Data_mat[((i*m) - (holdout-1)):(i*m),]
}

#1.2.2 Bayesian lasso
burninit<-1000
samps<-1000
batch_size <- round(sqrt(samps*holdout))
thinning <- NULL
betamat_b <- matrix(NA, nrow=n*samps, ncol=d)
infomat_b <- matrix(NA, nrow = holdout*samps*n, ncol=5)
colnames(infomat_b) <- c("Predictions", "Actual", "Residual", "Estimated SE", "Corrected MSE")


for(i in 1:n){
  y <- Data_mat[(((i-1)*m)+1):(i*m - holdout), 1]
  X <- Data_mat[(((i-1)*m)+1):(i*m - holdout), 2:d]
  mod_obj <- blasso(X, y, thin=thinning, T=(burninit+samps), normalize=FALSE)
  betas<-mod_obj$beta[(burninit:(burninit + samps - 1)),] 
  mu_blasso<-mod_obj$mu[(burninit:(burninit + samps - 1))]
  betamat_b[((i-1)*samps+1):(i*samps),] <- cbind(mu_blasso, betas)
  sigmasq_blasso <- mod_obj$s2[(burninit:(burninit + samps))]
  for(j in 1:samps){
    for(k in 1:holdout){
      pred <- rnorm(1, mean = mu_blasso[j] + mth_obs[(k+((i-1)*holdout)),2:d] %*% betas[j,], 
                    sd = sqrt(sigmasq_blasso[j]))
      actual <- mth_obs[k + (i-1)*holdout,1]
      
      infomat_b[((i-1)*(samps*holdout) + (j-1)*holdout + k), 1] <- pred
      infomat_b[((i-1)*(samps*holdout) + (j-1)*holdout + k), 2] <- mth_obs[(k+((i-1)*holdout)),1]
      infomat_b[((i-1)*(samps*holdout) + (j-1)*holdout + k), 3] <- pred - actual
    }
  }
}

placeholder <- matrix(data=NA, nrow=samps*holdout, ncol=n)

for(i in 1:n){
  placeholder[1:(samps*holdout), i] <- infomat_b[((i-1)*(samps*holdout)+1):(i*samps*holdout),1]
}

temp <- batchSE(as.mcmc(placeholder), batchSize = batch_size)

for(i in 1:n){
  infomat_b[((i-1)*(samps*holdout)+1):(i*samps*holdout),4] <- rep(temp[i], (samps*holdout))
}

infomat_b[,5] <- infomat_b[,3]^2 - infomat_b[,4]^2

mean(infomat_b[,5])


#1.2.3 Our model

para_len <- 2*d + (d*(d-1)/2)
mcmcsamps <- 1000
burnin_mcmc <- 1000

paramat_pcor <- metrop_samp(n, m, para_len, Data_mat, mcmcsamps, 
                            improved_target_dens, burn=burnin_mcmc, holdout=holdout)


ESS_mat_red <- matrix(ESS_mat[ESS_mat!=0], ncol=(para_len-d))
colMeans(ESS_mat_red)
a <- sum(ESS_mat_red)/(ncol(ESS_mat_red)*nrow(ESS_mat_red)) #average ESS per parameter
a

mu_1 <- paramat_pcor[,1]
#removing crashed samples
cond <- is.na(mu_1[seq(1, ((n-1)*mcmcsamps+1), by=mcmcsamps)])

removals<-rep(TRUE, holdout*n)

for(i in 1:n){
  if(cond[i]==TRUE){
    removals[(((i-1)*holdout)+1):(i*holdout)] <- rep(FALSE, holdout)
  }
}

mth_obs_1 <- mth_obs[removals,]

n_missing_rows <- sum(ifelse(cond, 1, 0))
n_1 <- n - n_missing_rows
paramat_pcor <- paramat_pcor[complete.cases(paramat_pcor),]


# Transforming partial corrs to corrs
paramat <- pcors_to_corrs(paramat_pcor, d)

betamat <- betafrommvn(paramat, d)


#1.2.4 residuals in 1 dimension
batch_size <- round(sqrt(mcmcsamps*holdout))
infomat <- matrix(NA, nrow = holdout*mcmcsamps*n_1, ncol=5)
colnames(infomat) <- c("Predictions", "Actual", "Residual", "Estimated SE", "Corrected MSE")

for(i in 1:n_1){
  print(i)
  for(j in 1:mcmcsamps){
    theta <- paramat[((i-1)*mcmcsamps + j),] # 1 sample from posterior
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
    condvar<-as.numeric(Sigma_11 - Sigma_12 %*% Sigma_22_inv %*%  Sigma_21)
    
    for(k in 1:holdout){
      condmean <- as.numeric(mu_1 + Sigma_12 %*% Sigma_22_inv %*% (mth_obs_1[(k+((i-1)*holdout)),2:d] - mu_2))
      
      # 1 sample of predicted y given x_1, ... on observation not included in model
      predsim <- rnorm(1, mean = condmean, sd <- sqrt(condvar))
      
      infomat[((i-1)*(mcmcsamps*holdout) + (j-1)*holdout + k),1] <- predsim
      infomat[((i-1)*(mcmcsamps*holdout) + (j-1)*holdout + k),2] <- mth_obs_1[(k+((i-1)*holdout)),1]
      infomat[(i-1)*(mcmcsamps*holdout) + (j-1)*holdout + k, 3] <- (mth_obs_1[(k+((i-1)*holdout)),1] - predsim)
    }
  }
}

placeholder <- matrix(data=NA, nrow=mcmcsamps*holdout, ncol=n_1)

for(i in 1:n_1){
  placeholder[1:(mcmcsamps*holdout), i] <- infomat[((i-1)*(mcmcsamps*holdout)+1):(i*mcmcsamps*holdout),1]
}

temp <- batchSE(as.mcmc(placeholder), batchSize = batch_size)

for(i in 1:n_1){
  infomat[((i-1)*(mcmcsamps*holdout)+1):(i*mcmcsamps*holdout),4] <- rep(temp[i], (mcmcsamps*holdout))
}

infomat[,5] <- infomat[,3]^2 - infomat[,4]^2

mean(infomat[,3]^2)
mean(infomat[,5])
mean(infomat_b[,5])


accrate
n_1
quantile(ESS_mat_red, 0.01)


#Plots, not included in thesis

#betamat d=3
par(mfrow=c(2,3))
hist(betamat[,1], breaks=50, main="Beta_0")
hist(betamat[,2], breaks=50, main="Beta_1")
hist(betamat[,3], breaks=50, main="Beta_2")
hist(betamat_b[,1], breaks=50, main="Beta_0")
hist(betamat_b[,2], breaks=50, main="Beta_1")
hist(betamat_b[,3], breaks=50, main="Beta_2")

#betamat d=4
par(mfrow=c(3,3))
hist(betamat[,1], breaks=50, main="Beta_0")
hist(betamat[,2], breaks=50, main="Beta_1")
hist(betamat[,3], breaks=50, main="Beta_2")
hist(betamat[,4], breaks=50, main="Beta_3")

hist(betamat_b[,1], breaks=50, main="Beta_0")
hist(betamat_b[,2], breaks=50, main="Beta_1")
hist(betamat_b[,3], breaks=50, main="Beta_2")
hist(betamat_b[,4], breaks=50, main="Beta_3")

#ourmod d=3
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

#ourmod d=4
par(mfrow=c(4,4))
hist(paramat[,1], breaks=100, main="mu_1")
hist(paramat[,2], breaks=100, main="mu_2")
hist(paramat[,3], breaks=100, main="mu_3")
hist(paramat[,4], breaks=100, main="mu_4")
hist(paramat[,5], breaks=100, main="sigma_1")
hist(paramat[,6], breaks=100, main="sigma_2")
hist(paramat[,7], breaks=100, main="sigma_3")
hist(paramat[,8], breaks=100, main="sigma_4")
hist(paramat[,9], breaks=100, main="rho_12")
hist(paramat[,10], breaks=100, main="rho_13")
hist(paramat[,11], breaks=100, main="rho_14")
hist(paramat[,12], breaks=100, main="rho_23")
hist(paramat[,13], breaks=100, main="rho_24")
hist(paramat[,14], breaks=100, main="rho_34")

#############################################################################################
### Example: Longley dataset (Orange trees circumference (cm) vs age (days))
#############################################################################################

# Code not meant to be general, only change n_obs

rm(list = setdiff(ls(), lsf.str())) #Clear variables

set.seed(4)

samps <- 100000 #number of posterior samples
n_obs <- 20

data <- Orange
data <- data[,-1] #Remove type of tree from dataset
data <- data[sample(35,n_obs),] #select n_obs observations for estimating parameters

summary(data)
dim(data)

lmcoef <- lm(circumference ~ age, data=data)$coefficients
lmcoef
thinning <- 1000

bcoef <- blasso(data$age, data$circumference, T=samps/thinning, thin=thinning)
bcoefs <- cbind(bcoef$mu, bcoef$beta)
colMeans(bcoefs) #bayesian lasso posterior mean estimates

covmat <- cov(data)
pcorrmat <- cor2pcor(cov2cor(covmat))
init <- c(sqrt(diag(covmat)), pcorrmat[lower.tri(pcorrmat)==TRUE])

chain <- MCMCmetrop1R(improved_target_dens, theta.init=init, thin=thinning,
                      burnin = 10000, x=data, mcmc=samps)
colMeans(chain)
effectiveSize(chain)

#sample mu
mu <- matrix(NA, nrow=samps/thinning, ncol=2)
beta_0 <- rep(NA, nrow=samps/thinning)
for(j in 1:(samps/thinning)){
  sigmas <- c(chain[j,1], chain[j,2])
  corrmat <- toeplitz(c(1,chain[j,3]))
  Sigma <- (1/n_obs) * diag(sigmas) %*% corrmat %*% diag(sigmas)
  mu[j,] <- rmvnorm(1, mean=(colMeans(data)), sigma = (Sigma))
}

beta_1 <- (chain[,2]/chain[,1]) * chain[,3]
for(j in 1:(samps/thinning)){
  beta_0[j] <- mu[j,2] - mu[j,1] * beta_1[j]
}

lmcoef
blas <- colMeans(bcoefs)
blas
mvnmod <- colMeans(cbind(beta_0, beta_1))
mvnmod

plot(data, main="Orange tree estimates", cex=0.8)
abline(a=lmcoef[1], b=lmcoef[2])
abline(a=blas[1], b=blas[2], col="blue")
abline(a=mvnmod[1], b=mvnmod[2], col="red")
legend("topleft", legend=c("MLE", "B.Lasso", "MVN"), col=c("black", "blue", "red"), lty=1, cex=0.8)


par(mfrow=c(2,2))
plot(data, main="Regression lines, n=20", cex=1.4, xlim=c(0,2000), ylim=c(0,300), type="n", pch=19)
for(i in 91:100){
  abline(a=bcoefs[i,1], b=bcoefs[i,2], col="forestgreen", lwd=1.5)
  abline(a=beta_0[i], b=beta_1[i], col="goldenrod3", lwd=1.5)
}
points(data, pch=19)