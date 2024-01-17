library(MCMCpack)
library(mvtnorm)
library(corpcor)
library(clusterGeneration)

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
  
  rhoest<-mean(samp_rho)
  sigma1est<-mean(samp_sigma1)
  sigma2est<-mean(samp_sigma2)
  Sigma<-matrix(data=c(sigma1est^2, rhoest*sigma1est*sigma2est, rhoest*sigma1est*sigma2est, sigma2est^2), nrow=2)
  mu<-rmvnorm(nsamples, c(mean_1, mean_2), sigma=Sigma/nsamples)
  
  return(cbind(samp_rho, samp_sigma1, samp_sigma2, mu[,1], mu[,2]))
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
n<-1000

#holder <- rlogis(n,0,1) #prior from berger sun
#rho <- 2*plogis(holder)-1
#rho[rho>0.99]

#rho<-seq(from=-0.99, to=0.99, length.out=n)
rho<-runif(n, min=-0.0000001, max=0.0000001)

#rho<-seq(from=-0.99, to=0.99, length.out=n)
hist(rho)

mu_1 <- 0
mu_2 <- 0
sigma_1 <- 10
sigma_2 <- 10

#Sampling data

meanvec<-c(mu_1, mu_2)
m<-5
x_1<-c()
x_2<-c()

for(i in 1:n){
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

parasims<-40
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

#General check of predictive distribution for >= 3 dimensions

#1: Sample n priors
n<-10000
d<-3 #dimension

muvec<-rep(0,d)
sigmavec<-rep(1,d)

Sigma <- diag(d)
corr_from_pcor <- replicate(n,{
  repeat{
    u <- runif(d*(d-1)/2, -1, 1)
    Sigma[lower.tri(Sigma)] <- u
    Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
    if (sum(eigen(Sigma)$val>0)==d)  
      break()
  }
  P <- solve(Sigma)
  D <- diag(1/sqrt(diag(P)))
  corr <- - (D %*% P %*% D)[lower.tri(P)]
  # valid correlations from Artner (2022) "The space of partial correlation matrices" corollary 2
  corr
})

x_1<-c()
x_2<-c()
x_3<-c()
m<-4


for(i in 1:n){
  corrmat<- diag(1, nrow=d)
  corrmat[lower.tri(corrmat)==TRUE]<-corr_from_pcor[,i]
  corrmat[upper.tri(corrmat)==TRUE]<-corrmat[lower.tri(corrmat)==TRUE]
  Sigma <- diag(sigmavec) %*% corrmat %*% diag(sigmavec)
  a <- rmvnorm(m, mean=muvec, sigma=Sigma)
  x_1 <- append(a[,1], x_1)
  x_2 <- append(a[,2], x_2)
  x_3 <- append(a[,3], x_3)
}
