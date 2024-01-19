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
