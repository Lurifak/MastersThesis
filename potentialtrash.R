#Stuff i miiiiiight need later

#banned suckas banned


# To observarsjoner fra regresjonsmodellen
x <- rbind(c(0,0),c(1,1))
# posterioritetthet, logistic prior på glogit(rho),
# uniform improper på mu1, mu2 og log(sigma_1), log(sigma_2)
target <- function(theta, x, scale){
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







preds <- rep(NA, holdout*mcmcsamps*n)
resid_vec_ourmod<-rep(NA, holdout*mcmcsamps*n)

infomat <- matrix(NA, nrow = holdout*mcmcsamps*n, ncol=3)
colnames(infomat) <- c("Predictions", "Squared Residual", "Estimated SE")

for(i in 1:n){
  for(j in 1:mcmcsamps){
    print(j)
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
      condmean <- as.numeric(mu_1 + Sigma_12 %*% Sigma_22_inv %*% (mth_obs[(k+((i-1)*holdout)),2:d] - mu_2))
      
      # 1 sample of predicted y given x_1, ... on observation not included in model
      predsim <- rnorm(1, mean = condmean, sd <- sqrt(condvar))
      
      infomat[((i-1)*(mcmcsamps*holdout) + (j-1)*holdout + k),1] <- predsim
      infomat[(i-1)*(mcmcsamps*holdout) + (j-1)*holdout + k, 2] <- (mth_obs[(k+((i-1)*holdout)),1] - predsim)
    }
  }
}

burninit<-2000
samps<-20
thinning <- 10
infomat_b <- matrix(NA, nrow = holdout*samps*n, ncol=5)
colnames(infomat) <- c("Predictions", "Actual", "Residual", "Estimated SE", "Corrected MSE")


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
      resid <- mth_obs[k + (i-1)*holdout,1] - pred
      
      infomat_b[((i-1)*(samps*holdout) + (j-1)*holdout + k), 1] <- pred
      infomat_b[((i-1)*(samps*holdout) + (j-1)*holdout + k), 2] <- mth_obs[(k+((i-1)*holdout)),1]
      infomat_b[((i-1)*(samps*holdout) + (j-1)*holdout + k), 2] <- pred - resid
    }
  }
}

placeholder <- matrix(data=NA, nrow=mcmcsamps*holdout, ncol=n)

for(i in 1:n){
  placeholder[1:(mcmcsamps*holdout), i] <- infomat_b[((i-1)*(samps*holdout)+1):(i*samps*holdout),1]
}

temp <- batchSE(as.mcmc(placeholder), batchSize = batch_size)

for(i in 1:n){
  infomat_b[((i-1)*(mcmcsamps*holdout)+1):(i*mcmcsamps*holdout),4] <- rep(temp[i], (samps*holdout))
}

infomat_b[,5] <- infomat_b[,3]^2 - infomat_b[,4]^2

#Test 1

set.seed(2)
n <- 2
d <- 3

muvec<-rep(0,d)
sigmavec<-rep(1,d)

corrs <- corr_from_pcor(n,d)

m <- 100000 #how many datapoints per iteration.
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

#Standardizing (x_1, ..., x_p)
for(i in 2:d){
  Data_mat[,i] <- (Data_mat[,i] - mean(Data_mat[,i]))/sd(Data_mat[,i])
}

#1.1.3 posterior sampling

para_len <- 2*d + (d*(d-1)/2)
mcmcsamps <- 1000
burnin_mcmc <- 5000

paramat_pcor <- metrop_samp(n,m,para_len, Data_mat, mcmcsamps,improved_target_dens, burn=burnin_mcmc,holdout=0)

# Transforming partial corrs to corrs
paramat <- pcors_to_corrs(paramat_pcor, d)

betamat <- betafrommvn(paramat, d)

#blasso
burninit<-5000
samps<-1000
thinning <- NULL
betamat_b <- matrix(NA, nrow=(n-1)*samps, ncol=d)
infomat_b <- matrix(NA, nrow = holdout*samps*n, ncol=5)
colnames(infomat) <- c("Predictions", "Actual", "Residual", "Estimated SE", "Corrected MSE")


for(i in 1:(n-1)){
  y <- Data_mat[(((i-1)*m)+1):(i*m - holdout), 1]
  X <- Data_mat[(((i-1)*m)+1):(i*m - holdout), 2:d]
  mod_obj <- blasso(X, y, thin=thinning, T=(burninit+samps))
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

colnames(Data_mat) <- c("y", "x_1", "x_2")

#betamat
par(mfrow=c(2,3))
hist(betamat[,1], breaks=50, main="Beta_0")
hist(betamat[,2], breaks=50, main="Beta_1")
hist(betamat[,3], breaks=50, main="Beta_2")
hist(betamat_b[,1], breaks=50, main="Beta_0")
hist(betamat_b[,2], breaks=50, main="Beta_1")
hist(betamat_b[,3], breaks=50, main="Beta_2")
colMeans(betamat)
colMeans(betamat_b)
lm(y ~ x_1 + x_2, as.data.frame(Data_mat[1:100000,]))$coefficients


improved_target_dens<-function(theta,x){
  L <- length(theta)
  Samp_cov <- cov(x)
  d <- -1/2 + sqrt(1/4 + 2 * L)
  sigma <- theta[1:d]
  parcorrs <- theta[(d+1):L]
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

improved_target_dens<-function(theta,x){
  L <- length(theta)
  Samp_cov <- cov(x)
  n <- nrow(x)
  d <- ncol(x)
  sigma <- theta[1:d]
  parcorrs <- theta[(d+1):L]
  if(any(sigma<=0)){-Inf}
  else{
    parcorrmat <- diag(-1, nrow=d)
    parcorrmat[lower.tri(parcorrmat)==TRUE] <- parcorrs
    parcorrmat <- parcorrmat + t(parcorrmat) - diag(diag(parcorrmat))
    
    if ((sum(eigen(parcorrmat)$val<0)==d)){ # if partial corr mat is negative definite
      #S <- - parcorrmat
      #S_inv <- solve(S)
      #calculating precision matrix
      #Sigma_inv <- diag(1/sigma) %*% diag(sqrt(diag(S_inv))) %*% S %*% diag(sqrt(diag(S_inv))) %*% diag(1/sigma)
      Sigma_inv <- solve(diag(sigma) %*% pcor2cor(parcorrmat+diag(2, nrow=d)) %*% diag(sigma))
      -((n-1)/2) * (log(1/det(Sigma_inv)) + sum(diag(Sigma_inv %*% Samp_cov))) - sum(log(sigma))
    }
    else{-Inf}
  }
}

for(i in 1:n){
  print(i)
  Sigma <- matrix(data=c(sigma_1^2,sigma_1*sigma_2*rho[i],sigma_1*sigma_2*rho[i],sigma_2^2), nrow=2)
  a <- rmvnorm(m, mean=meanvec, sigma=Sigma)
  y <- append(a[,1], x_1)
  x <- append(a[,2], x_2)
}

z <- cbind(y, x)

k <- 2
holdouts <- matrix(NA, nrow=n*k, ncol=2)
holdouts[1:(n*k),2] <- rnorm(n*k, mean=mu_2, sd=sigma_2)

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

parasims<-1
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
