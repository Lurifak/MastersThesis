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

tune_adjust<-function(a_prob){
  if(a_prob>0.25){4*(a_prob-0.25) + 1}
  else{3 * (a_prob) +  1/4}
}

#test 2 using parametrization p(sigma) = 1/sigma

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
    
    if((sum(eigen(parcorrmat)$val<0)==d)){ # if partial corr mat is negative definite <=> corrmat positive def
      
      Sigma_inv <- solve(diag(sigma) %*% pcor2cor(parcorrmat+diag(2, nrow=d)) %*% diag(sigma)) #pcor2cor defines p. corr mat having unit diagonal instead of negative unit diag
      -((n-1)/2) * (log(1/det(Sigma_inv)) + sum(diag(Sigma_inv %*% Samp_cov))) - sum(log(sigma))
    } else {
      -Inf
    }
  }
}

betafrommvn <- function(PMat, d){ #PMAt = parameter matrix, d = dim of mvn vector
  #Assumed rows in PMat: mu_s, sigmas and corrs
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
    block <- Data_mat[((i-1)*m+1):(i*t),]
    
    #Determining initialization cheaply
    covmat <- cov(block)
    pcorrmat <- cor2pcor(cov2cor(covmat))
    init <- c(sqrt(diag(covmat)), pcorrmat[lower.tri(pcorrmat)==TRUE])
    
    #For estimated correlations almost on edge of parameter space it can crash - thus skip these
    
    tuning <- 1
    acc_prob_pre <- 0
    ess_test <- rep(1, (d) + (d)*(d-1)/2)
    
    #test-chain
    
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
    
    
    tryCatch({
      
      sink(file="test.txt")
      
      sigma_pcor <- MCMCmetrop1R(improved_target_dens, theta.init=init, 
                                 burnin = burn, x=block, mcmc=mcmcsamps, tune=tuning)
      
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

xsim <- function(thetamat, nsims){
  n <- nrow(thetamat)
  x <- rep(NA, (n*nsims))
  for(i in 1:n){
    x[((i-1)*nsims + 1):(i*nsims)] <- rnorm(nsims, mean=wadup[i,5], sd=wadup[i,3])
  }
  x
}

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
n<-100000

rho<-runif(n, -0.99, 0.99)

hist(rho)

mu_1 <- 0
mu_2 <- 0
sigma_1 <- 1
sigma_2 <- 1

#Sampling data

meanvec<-c(mu_1, mu_2)
m <- 5
holdout <- 2
y<-c()
x<-c()
holdout_x <- c()
holdout_y <- c()

for(i in 1:n){
  print(i)
  Sigma <- matrix(data=c(sigma_1^2, sigma_1*sigma_2*rho[i], sigma_1*sigma_2*rho[i], sigma_2^2), nrow=2)
  a <- rmvnorm((m+holdout), mean=meanvec, sigma=Sigma)
  y <- append(a[1:(m), 1], y)
  x <- append(a[1:(m), 2], x)
  holdout_y <- append(a[(m+1):(m+holdout), 1], holdout_y)
  holdout_x <- append(a[(m+1):(m+holdout), 2], holdout_x)
}

z <- cbind(y, x)

count <- rep(0, m)

parasims <- 60
predsims <- 3

residvec <- rep(NA, n*holdout)

paramat <- matrix(NA, nrow=parasims*n, ncol=5)

for (i in 1:(n)){
  print(i)
  newdata <- z[((i-1)*m + 1): (i*m),] #Block of m of the data for each iteration
  wadup <- alg(newdata, parasims) #simulating parameters for this block
  x_sim <- xsim(wadup, predsims) #simulating x
  y_x <- ysim(wadup, x_sim)  #simulating y conditional on x
  
  sorted <- sort(newdata[,1], decreasing=TRUE)
  
  for(j in 1:m){
    count[j] <- count[j] + sum(ifelse(y_x>sorted[j], 1, 0)) #counting times ysim above observed vals
  }
  
  hold <- cbind(holdout_y[((i-1)*holdout + 1):(i*holdout)], holdout_x[((i-1)*holdout + 1):(i*holdout)])
  y_sim <- ysim(wadup, hold[,2])
  
  for(k in 1:holdout){
    residvec[(i-1)*(holdout) + k] <- y_sim[k] - hold[k,1]
  }
  
  paramat[((i-1)*parasims + 1): (i*parasims),1:2] <- wadup[,4:5]
  paramat[((i-1)*parasims + 1): (i*parasims),3:4] <- wadup[,2:3]
  paramat[((i-1)*parasims + 1): (i*parasims),5] <- wadup[,1]
  
}

mean(residvec^2)

tot_comb<-n*parasims*predsims
count/tot_comb #amount of y_pred_n+1 above each observed y_i up to n
seq(from=1/(m+1), to=m/(m+1), by=1/(m+1)) #Expected under t dist

par(mfrow=c(1,2))
hist(betas[,1], main="Beta_0")
hist(betas[,2], main="Beta_1")

par(mfrow=c(2,3))
hist(paramat[,1], main="mu_1")
hist(paramat[,2], main="mu_2")
hist(paramat[,3], main="sigma_1", breaks=30, xlim=c(0,10))
hist(paramat[,4], main="sigma_2", breaks=30, xlim=c(0,10))
hist(paramat[,5], main="rho")


#############################################################################################
### Diagnostics (compute ESS, trace plots, etc... to determine tuning, burnin)
#############################################################################################


# First we determine tuning parameter 'tune' in mcmcmetrop1r

# We have 2 steps: first determine norm of tuning vector, 
# then determine weight between marginal variances and p. correlations
# measured by how well our ESS is

# 1: Pregenerating data from correlations near edge of parameter space
# i.e. partial corrs will also lie close to parameter space



rm(list = setdiff(ls(), lsf.str()))

set.seed(1)


samp_corrs <- corr_from_pcor(5000, 3)

# 1.1 pick out m observations close to parameter edge

d <- 3
j <- 0
i <- 0
m <- 20
par_list <- list()
while((j<m) & (i<5000)){
  i <- i + 1
  corrmat <- diag(1, nrow=d)
  corrmat[lower.tri(corrmat)==TRUE] <- samp_corrs[,i]
  corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
  if((any(abs(eigen(corrmat)$values) < 10^(-1))) & (any(abs(eigen(corrmat)$values) > 10^(-2)))){
    par_list[[j+1]] <- corrmat
    j <- j+1
  }
}

par_list #all correlations near edge of parameter space

n_obs <- 4

data_list <- list()
for(i in 1:m){
  data_list[[i]] <- rmvnorm(n_obs,mean=rep(0,d), sigma=par_list[[i]])
}

# 1.2 checking best absolute tune param size

n_tune_vec <- 4
tot_ESS <- matrix(0, nrow=n_tune_vec, ncol=2*d)

for(i in 1:n_tune_vec){
  tune_vec <- i/2  #1/2, 1, 3/2, ... for each tuning param
  for(j in 1:m){
    #Determining initialization cheaply
    covmat <- cov(data_list[[j]])
    pcorrmat <- cor2pcor(cov2cor(covmat))
    init <- c(sqrt(diag(covmat)), pcorrmat[lower.tri(pcorrmat)==TRUE])
    
    tryCatch({
      chain <- MCMCmetrop1R(improved_target_dens, theta.init=init, 
                          burnin = 500, x=data_list[[j]],tune=tune_vec, mcmc=3000)
    
      tot_ESS[i,] <- tot_ESS[i,] + effectiveSize(chain)
    }, error=function(e){})
  }
}

tot_ESS
rowMeans(tot_ESS)

#Checking increasing tuning parameter for marginal variances

tot_ESS <- matrix(0, nrow=n_tune_vec, ncol=2*d)

for(i in 1:n_tune_vec){
  tune_vec <- c(rep(1, d)*i/2, rep(1,d))
  for(j in 1:m){
    #Determining initialization cheaply
    covmat <- cov(data_list[[j]])
    pcorrmat <- cor2pcor(cov2cor(covmat))
    init <- c(sqrt(diag(covmat)), pcorrmat[lower.tri(pcorrmat)==TRUE])
    
    tryCatch({
      chain <- MCMCmetrop1R(improved_target_dens, theta.init=init, 
                            burnin = 500, x=data_list[[j]], tune=tune_vec, mcmc=2000)
      
      tot_ESS[i,] <- tot_ESS[i,] + effectiveSize(chain)
    }, error=function(e){})
  }
}

tot_ESS
rowMeans(tot_ESS) #tuning params for sigmas = 1 best

#Checking tuning parameter for partial correlations

tot_ESS <- matrix(0, nrow=n_tune_vec, ncol=2*d)

for(i in 1:n_tune_vec){
  tune_vec <- c(rep(1, d), rep(1,d)*(i))
  for(j in 1:m){
    #Determining initialization cheaply
    covmat <- cov(data_list[[j]])
    pcorrmat <- cor2pcor(cov2cor(covmat))
    init <- c(sqrt(diag(covmat)), pcorrmat[lower.tri(pcorrmat)==TRUE])
    
    tryCatch({
      chain <- MCMCmetrop1R(improved_target_dens, theta.init=init, 
                            burnin = 1000, x=data_list[[j]],tune=tune_vec, mcmc=1000)
      
      tot_ESS[i,] <- tot_ESS[i,] + effectiveSize(chain)
    }, error=function(e){})
  }
}

tot_ESS
rowMeans(tot_ESS)


#Testing adapting tuning by ESS
tot_ESS_tune <- matrix(0, nrow=m, ncol=2*d)
tot_ESS <- matrix(0, nrow=m, ncol=2*d)


for(j in 1:m){
  #Determining initialization cheaply
  covmat <- cov(data_list[[j]])
  pcorrmat <- cor2pcor(cov2cor(covmat))
  init <- c(sqrt(diag(covmat)), pcorrmat[lower.tri(pcorrmat)==TRUE])
    
  tryCatch({
    testchain <- MCMCmetrop1R(improved_target_dens, theta.init=init, 
                              burnin = 500, x=data_list[[j]], mcmc=200)
    
    ess <- effectiveSize(testchain)
    tot <- sum(ess)
    tune_corrected <- c(rep(1,d) * sum(ess[1:d])/tot, rep(1,d) * sum(ess[(d+1):(2*d)])/tot)
    
    chain_tuned <- MCMCmetrop1R(improved_target_dens, theta.init=init,
                                burnin = 500, x=data_list[[j]],tune=tune_corrected, mcmc=2000)
    
    chain <- MCMCmetrop1R(improved_target_dens, theta.init=init, 
                            burnin = 500, x=data_list[[j]], mcmc=2000)
      
    tot_ESS_tune[j,] <- effectiveSize(chain_tuned)
    tot_ESS[j,] <- effectiveSize(chain)
    
  }, error=function(e){})
}

tot_ESS_tune
tot_ESS
colMeans(tot_ESS_tune)
colMeans(tot_ESS)



#calculating intialization
set.seed(3)
examp_data <- rmvnorm(4, mean=rep(0,d), sigma = diag(1, nrow=3))

covmat <- cov(examp_data)
pcorrmat <- cor2pcor(cov2cor(covmat))
init_3 <- c(sqrt(diag(covmat)), pcorrmat[lower.tri(pcorrmat)==TRUE])
init_3



diagnostic_obj_3 <- MCMCmetrop1R(improved_target_dens, theta.init=init_3, 
                                 burnin = 1, x=examp_data,#tune=c(1,1,1,12,12,12), 
                                 mcmc=1700)

sink(file="test.txt")

sink()

out <- readLines("test.txt")

# Extract the number using regex
numbers <- regmatches(out, regexpr(pattern, out))

# Convert the extracted number to numeric
numbers <- as.numeric(numbers)

# Print the extracted number
print(numbers)

effectiveSize(diagnostic_obj_3)

set.seed(7)

examp_data <- rmvnorm(4, mean=rep(0,d), sigma = par_list[[14]])

covmat <- cov(examp_data)
pcorrmat <- cor2pcor(cov2cor(covmat))
cov2cor(covmat)
init_3 <- c(sqrt(diag(covmat)), pcorrmat[lower.tri(pcorrmat)==TRUE])
init_3


for(i in 1:100){
  diagnostic_obj_3 <- MCMCmetrop1R(improved_target_dens, theta.init=init_3, 
                                 burnin = 500, seed=i,x=examp_data, tune=1, 
                                 mcmc=200)
}

ø <- MCMCmetrop1R(improved_target_dens, theta.init=init_3, 
             burnin = 500, seed=6,x=examp_data, tune=1, 
             mcmc=10000)

a <- effectiveSize(diagnostic_obj_3)
a

diagnostic_obj <- MCMCmetrop1R(improved_target_dens, theta.init=init_3, 
                                 burnin = 500, x=examp_data, tune=0.32/0.25, 
                                 mcmc=3000)

b <- effectiveSize(diagnostic_obj)
b
sum(a)
sum(b)


plot(diagnostic_obj_3) #trace plots convincing for m >= d + 1

batchSE(as.mcmc(cbind(diagnostic_obj_3, rep(1, 100000))))

diagnostic_obj_4 <- MCMCmetrop1R(improved_target_dens, theta.init=c(rep(1,4), rep(0,6)), 
                                 burnin = 1,
                                 x=rmvnorm(10, mean=rep(0,4), sigma=diag(rep(1,4))), 
                                 mcmc=10000)

effectiveSize(diagnostic_obj_4)

plot(diagnostic_obj_4)

diagnostic_obj_5 <- MCMCmetrop1R(improved_target_dens, theta.init=c(rep(1,5), rep(0,10)), 
                                 burnin = 1,
                                 x=rmvnorm(6, mean=rep(0,5), sigma=diag(rep(1,5))),
                                   mcmc=100000)

effectiveSize(diagnostic_obj_5)

plot(diagnostic_obj_5)

#Visualizing betas for our prior

#B_0 and B_1

n <- 10000

rho <- runif(n, -1, 1)
params_1 <- cbind(0, 1, 1, 1, rho) #mu_y, mu_x, s_y, s_x, rho
params_2 <- cbind(1, 2, 1, 1, rho)
params_3 <- cbind(1, 1, 2, 1, rho)
params_4 <- cbind(1, 1, 2, 2, rho)

prior_beta_1 <- betafrommvn(params_1,2)
prior_beta_2 <- betafrommvn(params_2,2)
prior_beta_3 <- betafrommvn(params_3,2)
prior_beta_4 <- betafrommvn(params_4,2)

colnames(prior_beta_1) <- c("Beta_0", "Beta_1")
colnames(prior_beta_2) <- c("Beta_0", "Beta_1")
colnames(prior_beta_3) <- c("Beta_0", "Beta_1")
colnames(prior_beta_4) <- c("Beta_0", "Beta_1")

plt_1 <-  ggplot(data = prior_beta_1, mapping = aes(x = Beta_0, y = Beta_1)) +
  xlim(-6, 6) + ylim(-6,6) + xlab(TeX(r"($\Beta_0)")) + ylab(TeX(r"($\Beta_1)")) +
  geom_hex(bins=150) + scale_fill_continuous(type = "viridis")

plt_2 <-  ggplot(data = prior_beta_2, mapping = aes(x = Beta_0, y = Beta_1)) +
  xlim(-6, 6) + ylim(-6,6) + xlab(TeX(r"($\Beta_0)")) + ylab(TeX(r"($\Beta_1)")) +
  geom_hex(bins=150) + scale_fill_continuous(type = "viridis")

plt_3 <-  ggplot(data = prior_betas_111, mapping = aes(x = Beta_1, y = Beta_2)) +
  xlim(-6, 6) + ylim(-6,6) + xlab(TeX(r"($\Beta_1)")) + ylab(TeX(r"($\Beta_2)")) +
  geom_hex(bins=150) + scale_fill_continuous(type = "viridis")




#B_1 and B_2

n <- 1000000
d <- 3
muvec<-rep(0,d)
corrs <- corr_from_pcor(n,d)

par_mat_111 <- cbind(matrix(muvec, nrow=n, ncol=d), matrix(c(1,1,1), nrow=n, ncol=d), t(corrs))
par_mat_122 <- cbind(matrix(muvec, nrow=n, ncol=d), cbind(rep(1, n), rep(2, n), rep(2,n)), t(corrs))
par_mat_123 <- cbind(matrix(muvec, nrow=n, ncol=d), cbind(rep(1, n), rep(2, n), rep(3,n)), t(corrs))
par_mat_311 <- cbind(matrix(muvec, nrow=n, ncol=d), cbind(rep(3, n), rep(1, n), rep(1,n)), t(corrs))


prior_betas_111 <- betafrommvn(par_mat_111, d)
prior_betas_122 <- betafrommvn(par_mat_122, d)
prior_betas_123 <- betafrommvn(par_mat_123, d)
prior_betas_311 <- betafrommvn(par_mat_311, d)

colnames(prior_betas_111) <- c("Beta_0", "Beta_1", "Beta_2")
colnames(prior_betas_122) <- c("Beta_0", "Beta_1", "Beta_2")
colnames(prior_betas_123) <- c("Beta_0", "Beta_1", "Beta_2")
colnames(prior_betas_311) <- c("Beta_0", "Beta_1", "Beta_2")



plt_111 <-  ggplot(data = prior_betas_111, mapping = aes(x = Beta_1, y = Beta_2)) +
  xlim(-6, 6) + ylim(-6,6) + xlab(TeX(r"($\Beta_1)")) + ylab(TeX(r"($\Beta_2)")) +
  geom_hex(bins=150) + scale_fill_continuous(type = "viridis")

plt_122 <-  ggplot(data = prior_betas_122, mapping = aes(x = Beta_1, y = Beta_2)) +
  xlim(-6, 6) + ylim(-6,6) + xlab(TeX(r"($\Beta_1)")) + ylab(TeX(r"($\Beta_2)")) +
  geom_hex(bins=150) + scale_fill_continuous(type = "viridis")

plt_123 <- ggplot(data = prior_betas_123, mapping = aes(x = Beta_1, y = Beta_2)) +
  xlim(-6, 6) + ylim(-6,6) + xlab(TeX(r"($\Beta_1)")) + ylab(TeX(r"($\Beta_2)")) +
  geom_hex(bins=150) + scale_fill_continuous(type = "viridis")

plt_311 <- ggplot(data = prior_betas_311, mapping = aes(x = Beta_1, y = Beta_2)) +
  xlim(-6, 6) + ylim(-6,6) + xlab(TeX(r"($\Beta_1)")) + ylab(TeX(r"($\Beta_2)")) +
  geom_hex(bins=150) + scale_fill_continuous(type = "viridis")


grid.arrange(plt_111, plt_122, plt_123, plt_311, ncol=2, nrow=2)

plt_111_dens <-  ggplot(data = prior_betas_111, mapping = aes(x = Beta_1, y = Beta_2)) +
  xlim(-4, 4) + ylim(-4,4) + xlab(TeX(r"($\Beta_1)")) + ylab(TeX(r"($\Beta_2)")) +
  stat_density_2d(geom = "raster",aes(fill = after_stat(density)),contour = FALSE) + 
  scale_fill_viridis_c(option="turbo")

plt_122_dens <-  ggplot(data = prior_betas_122, mapping = aes(x = Beta_1, y = Beta_2)) +
  xlim(-4, 4) + ylim(-4,4) + xlab(TeX(r"($\Beta_1)")) + ylab(TeX(r"($\Beta_2)")) +
  stat_density_2d(geom = "raster",aes(fill = after_stat(density)),contour = FALSE) + 
  scale_fill_viridis_c(option="turbo")

plt_123_dens <-  ggplot(data = prior_betas_123, mapping = aes(x = Beta_1, y = Beta_2)) +
  xlim(-4, 4) + ylim(-4,4) + xlab(TeX(r"($\Beta_1)")) + ylab(TeX(r"($\Beta_2)")) +
  stat_density_2d(geom = "raster",aes(fill = after_stat(density)),contour = FALSE) + 
  scale_fill_viridis_c(option="turbo")

plt_311_dens <-  ggplot(data = prior_betas_311, mapping = aes(x = Beta_1, y = Beta_2)) +
  xlim(-4, 4) + ylim(-4,4) + xlab(TeX(r"($\Beta_1)")) + ylab(TeX(r"($\Beta_2)")) +
  stat_density_2d(geom = "raster",aes(fill = after_stat(density)),contour = FALSE) + 
  scale_fill_viridis_c(option="turbo")

grid.arrange(plt_111_dens, plt_122_dens, plt_123_dens, plt_311_dens, ncol=2, nrow=2)


#General check of predictive distribution for >= 3 dimensions

rm(list = setdiff(ls(), lsf.str()))

set.seed(1)

#1: Sample priors
n <- 5 #samples
d <- 3 #dimension

muvec<-rep(0,d)
sigmavec<-rep(1,d)

#Sampling n realizations of correlations with uniform marginals
corrs <- corr_from_pcor(n,d)

#2 Sample Data given priors

m <- 4 #how many datapoints we use to estimate the posterior sample
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


#Standardizing (y, x_1, ..., x_p)
for(i in 1:d){
  Data_mat[,i] <- (Data_mat[,i] - mean(Data_mat[,i]))/sd(Data_mat[,i])
}

#3 Estimate parameters

para_len <- 2*d + (d*(d-1))/2
mcmcsamps <- 1000
paramat_pcor <- metrop_samp(n, m, para_len, Data_mat, mcmcsamps, improved_target_dens, 
                            burn=1000)


ESS_mat_red <- matrix(ESS_mat[ESS_mat!=0], ncol=(para_len-d))
ESS_mat_red
colMeans(ESS_mat_red)
a <- sum(ESS_mat_red)/(ncol(ESS_mat_red)*nrow(ESS_mat_red)) #average ESS per parameter
a

#4: simulate predictive samples given samples from posterior

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


rowMeans(side_by_side)/mcmcsamps

#1: Comparison Bayesian Lasso and reparametrized model


# From our prior

rm(list = setdiff(ls(), lsf.str())) #removes all variables except functions

set.seed(1)

#1.1.1: Sample n priors
n<-50
d<-3 #dimension

muvec<-rep(0,d)
sigmavec<-rep(1,d)

corrs <- corr_from_pcor(n,d)

#1.1.2 Sample Data given priors

m <- 2020 #how many datapoints per iteration.
holdout <- 20
#We use m-holdout observations to fit the model and then compare
#sample from predicted (from model) with remaining observations

Data_mat<-matrix(0, nrow=(n*m), ncol=d)
mth_obs<-matrix(NA, nrow = (n*holdout), ncol=d) #holdout observations

for(i in 1:n){
  corrmat <- diag(1, nrow=d)
  corrmat[lower.tri(corrmat)==TRUE] <- corrs[,i]
  corrmat <- corrmat + t(corrmat) - diag(diag(corrmat))
  Sigma <- diag(sigmavec) %*% corrmat %*% diag(sigmavec)
  Data_mat[(((i-1)*m)+1):(i*m),] <- rmvnorm(m, mean=muvec, sigma=Sigma)
}

#Standardizing (y, x_1, ..., x_p)
for(i in 1:d){
  Data_mat[,i] <- (Data_mat[,i] - mean(Data_mat[,i]))/sd(Data_mat[,i])
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
a <- sum(ESS_mat_red)/(ncol(ESS_mat_red)*nrow(ESS_mat_red)) #average ESS per parameter
a

mu_1<-paramat_pcor[,1]
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

#1.1.5 residuals in 1 dimension

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

mean(infomat[,3]^2)
mean(infomat_b[,3]^2)

mean(infomat[,5])
mean(infomat_b[,5])

par(mfrow=c(2,1))
hist(infomat[,3], breaks=100)
hist(infomat_b[,3], breaks=100)

#Looking at first realization
realization1<-Data_mat[1:m,]
ourmod_mvn <- paramat[1:mcmcsamps,]
ourmod_beta <- betamat[1:mcmcsamps,]
blasso_beta <- betamat_b[1:samps,]

colnames(realization1) <- c("y","x_1","x_2")
lmcoef <- lm(y ~., data=as.data.frame(realization1))$coef #lm coefficients

lmcoef
colMeans(ourmod_beta) #Our model coefficients
colMeans(blasso_beta) #Bayesian Lasso coefficients

summary(betamat[1:mcmcsamps,]) 
summary(betamat_b[1:samps,])

par(mfrow=c(2,3)) #top row our model
hist(ourmod_beta[,1], breaks=20, main="Beta_0")
hist(ourmod_beta[,2], breaks=20, main="Beta_1")
hist(ourmod_beta[,3], breaks=20, main="Beta_2")
hist(blasso_beta[,1], breaks=20, main="Beta_0")
hist(blasso_beta[,2], breaks=20, main="Beta_1")
hist(blasso_beta[,3], breaks=20, main="Beta_2")

summary(ourmod_mvn)

par(mfrow=c(3,3))
hist(ourmod_mvn[,1])
hist(ourmod_mvn[,2])
hist(ourmod_mvn[,3])
hist(ourmod_mvn[,4])
hist(ourmod_mvn[,5])
hist(ourmod_mvn[,6])
hist(ourmod_mvn[,7])
hist(ourmod_mvn[,8])
hist(ourmod_mvn[,9])

# Generating data from bayesian lasso and then comparing

#1.2.1 Generating priors

rm(list = setdiff(ls(), lsf.str())) #removes all variables except functions

set.seed(2)

n <- 100
d <- 3
p <- (d-1) #for notational purposes, denote vector (y, x_1 , ..., x_p)
m <- 70
holdout <- 20

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

#Generate y conditional on X, \Beta, \sigma^2 (mu not included not sure how to and if i should)
y <- rep(NA, n*m)
for(i in 1:(n)){
  print(i)
  y[(1 + (i-1)*m): (i*m)] <- rmvnorm(1, mean=fulldata[(1 + ((i-1)*m)):((i*m)),] %*% beta_samp[i,], sigma_sq[i] * diag(m))
}

Data_mat <- cbind(y,fulldata)

#Standardizing, as assumed for B. Lasso
for(i in 1:d){
  Data_mat[,i] <- (Data_mat[,i] - mean(Data_mat[,i]))/sd(Data_mat[,i])
}

for(i in 1:n){
  mth_obs[((i-1)*holdout+1):(i*holdout),]<-Data_mat[((i*m) - (holdout-1)):(i*m),]
}

#1.2.2 Bayesian lasso
burninit<-1000
samps<-200
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



#1.2.3 Our model

para_len <- 2*d + (d*(d-1)/2)
mcmcsamps <- 3000
burnin_mcmc <- 1000

paramat_pcor <- metrop_samp(n, m, para_len, Data_mat, mcmcsamps, 
                            improved_target_dens, burn=burnin_mcmc, holdout=holdout)


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
infomat <- matrix(NA, nrow = holdout*mcmcsamps*n, ncol=5)
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

placeholder <- matrix(data=NA, nrow=mcmcsamps*holdout, ncol=n)

for(i in 1:n_1){
  placeholder[1:(mcmcsamps*holdout), i] <- infomat[((i-1)*(mcmcsamps*holdout)+1):(i*mcmcsamps*holdout),1]
}

temp <- batchSE(as.mcmc(placeholder), batchSize = batch_size)

for(i in 1:n_1){
  print(i)
  infomat[((i-1)*(mcmcsamps*holdout)+1):(i*mcmcsamps*holdout),4] <- rep(temp[i], (mcmcsamps*holdout))
}

infomat[,5] <- infomat[,3]^2 - infomat[,4]^2


mean(infomat[,5])
mean(infomat_b[,5])

#plots

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

#betamat d=5
par(mfrow=c(3,4))
hist(betamat[,1], breaks=50, main="Beta_0")
hist(betamat[,2], breaks=50, main="Beta_1")
hist(betamat[,3], breaks=50, main="Beta_2")
hist(betamat[,4], breaks=50, main="Beta_3")
hist(betamat[,5], breaks=50, main="Beta_4")

hist(betamat_b[,1], breaks=50, main="Beta_0")
hist(betamat_b[,2], breaks=50, main="Beta_1")
hist(betamat_b[,3], breaks=50, main="Beta_2")
hist(betamat_b[,4], breaks=50, main="Beta_3")
hist(betamat_b[,5], breaks=50, main="Beta_4")

# ourmod d=3
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

#for d=4
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

colMeans(paramat)


