library(faraway)

set.seed(2)

mypredict <- function(mod, new_data, alpha=0.05){
  
  #Gives the first nonzero index of a vector
  first_nonz_ind<-function(example_vec){
    i <- 1
    while (i<length(example_vec)){
      if(example_vec[i]!=0)
      {break}
      else
      {i <- i+1}
    }
    return(i)
  } 
  
  X <- mod$model
  new_form <- update(formula(mod), NULL ~.)
  x_0 <- model.matrix(new_form, new_data, xlev=mod$xlevels)
  ind <- first_nonz_ind(x_0)
  pred<-predict.glm(mod, newdata=new_data, type="link", se.fit=T)
  
  mle_eta_0 <- as.numeric(pred$fit)
  se<-qnorm(1-alpha/2)*pred$se.fit
  
  #Defining function giving the value of eta for the desired quantile by the inverted LRT 
  cutoff<-function(x)
  {return(proflik(x) - proflik(mle_eta_0) + qchisq(1-alpha, 1)/2)}
  
  #Extracting response and design matrix
  mat <- model.matrix(formula(mod), X)
  
  #Name of response variable
  resp_name <- all.vars(formula(mod))[1]
  
  #Response column
  response <- X[,resp_name]
  
  #In case of the response variable being binary it has to be converted from (1 or 2) to (0 or 1):
  if(is.factor(response)){response <- as.numeric(as.factor(response)) - 1}
  
  #Defining X_star from stackexchange
  X_star<-subset(mat, select = -(ind+1)) - (1/x_0[ind+1])*mat[,ind+1]%*%t(x_0[-(ind+1)])
  
  #Defining the function for the profile likelihood given eta_0
  proflik<-function(eta_0){
    offset_star <- as.vector((1/x_0[ind+1]) * eta_0 * mat[,ind+1])
    new_mod<-glm.fit(X_star, response, offset=offset_star, family = family(mod))
    log_lik_new<-length(new_mod$coef) - new_mod$aic/2 #AIC := 2p - 2 ln L <=> ln L = p - AIC/2
    return(log_lik_new)
  }
  
  # Todo - select interval in a smarter way.
  left<-uniroot(cutoff, lower=mle_eta_0 - 2 * se, upper=mle_eta_0, extendInt="upX")$root
  right<-uniroot(cutoff, lower=mle_eta_0, upper=mle_eta_0 + 2 * se, extendInt="downX")$root
  return(c(left, right))
}

#Test when n=10000
library(ISLR)
dataset<-ISLR::Default
mod_test <- glm(default~., family=binomial(link = "logit"), data=as.data.frame(dataset))
summary(mod_test)
pred <- predict.glm(mod_test, newdata=data.frame(student="Yes", balance=2000, income=40000), type="link", se.fit=T)
conf_int <- c(pred$fit-1.96*pred$se.fit, pred$fit+1.96*pred$se.fit)

#95 % Wald for eta_0
conf_int

#Supposed 95% profile conf. int for eta_0
mypredict(mod_test, new_data=data.frame(student="Yes", balance=2000, income=40000), alpha=0.05)

#Test when n=1000
dataset1<-dataset[sample(10000, size=1000, replace=F),]
mod_test_1 <- glm(default~., family=binomial(link = "logit"), data=as.data.frame(dataset1))
pred_1 <- predict.glm(mod_test_1, newdata=data.frame(student="Yes", balance=2000, income=40000), type="link", se.fit=T)
conf_int_1 <- c(pred_1$fit-1.96*pred_1$se.fit, pred_1$fit+1.96*pred_1$se.fit)
conf_int_1
mypredict(mod_test_1, new_data=data.frame(student="Yes", balance=2000, income=40000), alpha=0.05)


#Test when n=200
dataset2<-dataset[sample(10000, size=200, replace=F),]
mod_test_2 <- glm(default~., family=binomial(link = "logit"), data=as.data.frame(dataset2))
pred_2 <- predict.glm(mod_test_2, newdata=data.frame(student="Yes", balance=2000, income=40000), type="link", se.fit=T)
conf_int_2 <- c(pred_2$fit-1.96*pred_2$se.fit, pred_2$fit+1.96*pred_2$se.fit)
conf_int_2
mypredict(mod_test_2, new_data=data.frame(student="Yes", balance=2000, income=40000), alpha=0.05)

#Test when n=50
dataset3<-dataset[sample(10000, size=50, replace=F),]
mod_test_3 <- glm(default~., family=binomial(link = "logit"), data=as.data.frame(dataset3))
pred_3 <- predict.glm(mod_test_3, newdata=data.frame(student="Yes", balance=2000, income=40000), type="link", se.fit=T)
conf_int_3 <- c(pred_3$fit-1.96*pred_3$se.fit, pred_3$fit+1.96*pred_3$se.fit)
conf_int_3
mypredict(mod_test_3, new_data=data.frame(student="Yes", balance=2000, income=40000), alpha=0.05)





#Coverage for n=1000

alpha<-0.05
se<-qnorm(1-alpha/2)

m<-1000
design<-subset(dataset1, select=-1)
y<-simulate(mod_test_1, nsim=m)
colnames(y)<-rep("default", m)
test_proflik<-rep(NA, m)
test_wald<-rep(NA,m)
oldpred<-predict(mod_test_1, newdata=data.frame(student="Yes", balance=2000, income=40000), type="link", se.fit=T)
for (i in 1:m){
  mod_ci<-glm(formula(mod_test), family=family(mod_test), data=cbind(y[i],design))
  pred<-predict(mod_ci, newdata=data.frame(student="Yes", balance=2000, income=40000), type="link", se.fit=T)
  conf_int_prof<-mypredict(mod_ci, new_data=data.frame(student="Yes", balance=2000, income=40000), alpha=0.1)
  conf_int_wald<-c(pred$fit-se*pred$se.fit, pred$fit+se*pred$se.fit)
  test_proflik[i] <-ifelse(conf_int_prof[1]<=oldpred & oldpred<=conf_int_prof[2], 1, 0)
  test_wald[i] <-ifelse(conf_int_wald[1]<=oldpred$fit & oldpred$fit<=conf_int_wald[2], 1, 0)
}
sum(test_proflik)
sum(test_wald)


#Coverage for n=200
m<-100
design<-subset(dataset2, select=-1)
y<-simulate(mod_test_2, nsim=m)
colnames(y)<-rep("default", m)
test_proflik<-rep(NA, m)
test_wald<-rep(NA,m)
oldpred<-predict(mod_test_2, newdata=data.frame(student="Yes", balance=2000, income=40000), type="link", se.fit=T)
for (i in 1:m){
  mod_ci<-glm(formula(mod_test), family=family(mod_test), data=cbind(y[i],design))
  pred<-predict(mod_ci, newdata=data.frame(student="Yes", balance=2000, income=40000), type="link", se.fit=T)
  conf_int_prof<-mypredict(mod_ci, new_data=data.frame(student="Yes", balance=2000, income=40000), alpha=0.01)
  conf_int_wald<-c(pred$fit-se*pred$se.fit, pred$fit+se*pred$se.fit)
  test_proflik[i] <-ifelse(conf_int_prof[1]<=oldpred & oldpred<=conf_int_prof[2], 1, 0)
  test_wald[i] <-ifelse(conf_int_wald[1]<=oldpred$fit & oldpred$fit<=conf_int_wald[2], 1, 0)
}
sum(test_proflik)
sum(test_wald)

#Coverage for n=50
m<-100
design<-subset(dataset3, select=-1)
y<-simulate(mod_test_3, nsim=m)
colnames(y)<-rep("default", m)
test_proflik<-rep(NA, m)
test_wald<-rep(NA,m)
oldpred<-predict(mod_test_3, newdata=data.frame(student="Yes", balance=2000, income=40000), type="link", se.fit=T)
for (i in 1:m){
  mod_ci<-glm(formula(mod_test), family=family(mod_test), data=cbind(y[i],design))
  pred<-predict(mod_ci, newdata=data.frame(student="Yes", balance=2000, income=40000), type="link", se.fit=T)
  conf_int_prof<-mypredict(mod_ci, new_data=data.frame(student="Yes", balance=2000, income=40000), alpha=0.1)
  conf_int_wald<-c(pred$fit-se*pred$se.fit, pred$fit+se*pred$se.fit)
  test_proflik[i] <-ifelse(conf_int_prof[1]<=oldpred & oldpred<=conf_int_prof[2], 1, 0)
  test_wald[i] <-ifelse(conf_int_wald[1]<=oldpred$fit & oldpred$fit<=conf_int_wald[2], 1, 0)
}
sum(test_proflik)
sum(test_wald)




#Single covariate testing
library(ISLR)
dataset_s<-subset(ISLR::Default, select=c(1,3))

#n=1000
dataset_s_1000<-dataset_s[sample(10000, size=1000, replace=F),]
mod_s_1000 <- glm(default~., family=binomial(link = "logit"), data=as.data.frame(dataset_s_1000))
alpha<-0.1
se<-qnorm(1-alpha/2)

m<-100
design<-subset(dataset_s_1000, select=-1)
y<-simulate(mod_s_1000, nsim=m)
colnames(y)<-rep("default", m)
test_proflik<-rep(NA, m)
test_wald<-rep(NA,m)
oldpred<-predict(mod_s_1000, newdata=data.frame(balance=3000), type="link", se.fit=T) #p_0 roughly 0.995
for (i in 1:m){
  mod_ci<-glm(formula(mod_s_1000), family=family(mod_s_1000), data=cbind(y[i],design))
  pred<-predict(mod_ci, newdata=data.frame(balance=3000), type="link", se.fit=T)
  conf_int_prof<-mypredict(mod_ci, new_data=data.frame(balance=3000), alpha=0.1)
  conf_int_wald<-c(pred$fit-se*pred$se.fit, pred$fit+se*pred$se.fit)
  test_proflik[i] <-ifelse(conf_int_prof[1]<=oldpred$fit & oldpred$fit<=conf_int_prof[2], 1, 0)
  test_wald[i] <-ifelse(conf_int_wald[1]<=oldpred$fit & oldpred$fit<=conf_int_wald[2], 1, 0)
}
sum(test_proflik)
sum(test_wald)

#n=500
dataset_s_500<-dataset_s[sample(10000, size=500, replace=F),]
mod_s_500 <- glm(default~., family=binomial(link = "logit"), data=as.data.frame(dataset_s_500))

m<-100
design<-subset(dataset_s_500, select=-1)
y<-simulate(mod_s_500, nsim=m)
colnames(y)<-rep("default", m)
test_proflik<-rep(NA, m)
test_wald<-rep(NA,m)
oldpred<-predict(mod_s_500, newdata=data.frame(balance=3000), type="link", se.fit=T)
for (i in 1:m){
  mod_ci<-glm(formula(mod_s_500), family=family(mod_s_500), data=cbind(y[i],design))
  pred<-predict(mod_ci, newdata=data.frame(balance=3000), type="link", se.fit=T)
  conf_int_prof<-mypredict(mod_ci, new_data=data.frame(balance=3000), alpha=0.05)
  conf_int_wald<-c(pred$fit-se*pred$se.fit, pred$fit+se*pred$se.fit)
  test_proflik[i] <-ifelse(conf_int_prof[1]<=oldpred$fit & oldpred$fit<=conf_int_prof[2], 1, 0)
  test_wald[i] <-ifelse(conf_int_wald[1]<=oldpred$fit & oldpred$fit<=conf_int_wald[2], 1, 0)
}
sum(test_proflik)/m
sum(test_wald)/m

# n=1000

dataset_s_1000<-dataset_s[sample(10000, size=1000, replace=F),]
mod_s_1000 <- glm(default~., family=binomial(link = "logit"), data=as.data.frame(dataset_s_1000))
alpha<-0.05
se<-qnorm(1-alpha/2)

m<-100
design<-subset(dataset_s_1000, select=-1)
y<-simulate(mod_s_1000, nsim=m)
colnames(y)<-rep("default", m)
test_proflik<-rep(NA, m)
test_wald<-rep(NA,m)
oldpred<-predict(mod_s_1000, newdata=data.frame(balance=3000), type="link", se.fit=T)
for (i in 1:m){
  mod_ci<-glm(formula(mod_s_1000), family=family(mod_s_1000), data=cbind(y[i],design))
  pred<-predict(mod_ci, newdata=data.frame(balance=3000), type="link", se.fit=T)
  conf_int_prof<-mypredict(mod_ci, new_data=data.frame(balance=3000), alpha=0.05)
  conf_int_wald<-c(pred$fit-se*pred$se.fit, pred$fit+se*pred$se.fit)
  test_proflik[i] <-ifelse(conf_int_prof[1]<=oldpred$fit & oldpred$fit<=conf_int_prof[2], 1, 0)
  test_wald[i] <-ifelse(conf_int_wald[1]<=oldpred$fit & oldpred$fit<=conf_int_wald[2], 1, 0)
}
sum(test_proflik)/m
sum(test_wald)/m

#Poisson test
n<-1000
x_1 <- rnorm(n)
lam <- exp(1 + x_1)
y<-rep(NA, n) 
for(i in 1:n)
  {y[i]<-rpois(n=1, lambda=lam[i])}
name<-c("y", "x_1")
dataset_pois<-as.data.frame(cbind(y,x_1), col.names=name)
mod_pois <- glm(y~., data=dataset_pois, family=poisson)

m<-100
alpha<-0.05
se<-qnorm(1-alpha/2)
design<-subset(dataset_pois, select=-1)
y<-simulate(mod_pois, nsim=m)
colnames(y)<-rep("y", m)
test_proflik<-rep(NA, m)
test_wald<-rep(NA,m)
oldpred<-predict(mod_pois, newdata=data.frame(x_1=4), type="link", se.fit=T)
for (i in 1:m){
  mod_ci<-glm(formula(mod_pois), family=family(mod_pois), data=cbind(y[i],design))
  pred<-predict(mod_ci, newdata=data.frame(x_1=4), type="link", se.fit=T)
  conf_int_prof<-mypredict(mod_ci, new_data=data.frame(x_1=4), alpha=0.05)
  conf_int_wald<-c(pred$fit-se*pred$se.fit, pred$fit+se*pred$se.fit)
  test_proflik[i] <-ifelse(conf_int_prof[1]<=oldpred$fit & oldpred$fit<=conf_int_prof[2], 1, 0)
  test_wald[i] <-ifelse(conf_int_wald[1]<=oldpred$fit & oldpred$fit<=conf_int_wald[2], 1, 0)
}
sum(test_proflik)/m
sum(test_wald)/m

 
#Poisson test 2
n<-10
x_1 <- rnorm(n)
lam <- exp(5 * x_1)
y<-rep(NA, n) 
for(i in 1:n)
{y[i]<-rpois(n=1, lambda=lam[i])}
name<-c("y", "x_1")
dataset_pois<-as.data.frame(cbind(y,x_1), col.names=name)
mod_pois <- glm(y~., data=dataset_pois, family=poisson)

m<-100
alpha<-0.05
se<-qnorm(1-alpha/2)
design<-subset(dataset_pois, select=-1)
y<-simulate(mod_pois, nsim=m)
colnames(y)<-rep("y", m)
test_proflik<-rep(NA, m)
test_wald<-rep(NA,m)
oldpred<-predict(mod_pois, newdata=data.frame(x_1=4), type="link", se.fit=T)
for (i in 1:m){
  mod_ci<-glm(formula(mod_pois), family=family(mod_pois), data=cbind(y[i],design))
  pred<-predict(mod_ci, newdata=data.frame(x_1=4), type="link", se.fit=T)
  conf_int_prof<-mypredict(mod_ci, new_data=data.frame(x_1=4), alpha=0.05)
  conf_int_wald<-c(pred$fit-se*pred$se.fit, pred$fit+se*pred$se.fit)
  test_proflik[i] <-ifelse(conf_int_prof[1]<=oldpred$fit & oldpred$fit<=conf_int_prof[2], 1, 0)
  test_wald[i] <-ifelse(conf_int_wald[1]<=oldpred$fit & oldpred$fit<=conf_int_wald[2], 1, 0)
}
sum(test_proflik)/m
sum(test_wald)/m

#Gamma test
library(faraway)
data(wafer)
m3 <- glm(formula = resist ~ x1 + x2 + x3 + x4,
                         family  = Gamma(link = "log"),
                         data    = wafer)
m<-100
alpha<-0.05
se<-qnorm(1-alpha/2)
design<-subset(m3$model, select=-1)
y<-simulate(m3, nsim=m)
colnames(y)<-rep("resist", m)
test_proflik<-rep(NA, m)
test_wald<-rep(NA,m)
oldpred<-predict(m3, newdata=data.frame("x1"= '+',"x2"= '+',"x3"='+',"x4"='+'), type="link", se.fit=T)
for (i in 1:m){
  mod_ci<-glm(formula(m3), family=family(m3), data=cbind(y[i], design))
  pred<-predict(mod_ci, newdata=data.frame("x1"= '+',"x2"= '+',"x3"='+',"x4"='+'), type="link", se.fit=T)
  conf_int_prof<-mypredict(mod_ci, new_data=data.frame("x1"= '+',"x2"= '+',"x3"='+',"x4"='+'), alpha=0.05)
  conf_int_wald<-c(pred$fit-se*pred$se.fit, pred$fit+se*pred$se.fit)
  test_proflik[i] <-ifelse(conf_int_prof[1]<=oldpred$fit & oldpred$fit<=conf_int_prof[2], 1, 0)
  test_wald[i] <-ifelse(conf_int_wald[1]<=oldpred$fit & oldpred$fit<=conf_int_wald[2], 1, 0)
}
sum(test_proflik)/m
sum(test_wald)/m


#Poisson test 3 cauchy
n<-10
x_1 <- rcauchy(n)
lam <- exp(x_1)
y<-rep(NA, n) 
for(i in 1:n)
{y[i]<-rpois(n=1, lambda=lam[i])}
name<-c("y", "x_1")
dataset_pois<-as.data.frame(cbind(y,x_1), col.names=name)
mod_pois <- glm(y~., data=dataset_pois, family=poisson)

m<-100
alpha<-0.05
se<-qnorm(1-alpha/2)
design<-subset(dataset_pois, select=-1)
y<-simulate(mod_pois, nsim=m)
colnames(y)<-rep("y", m)
test_proflik<-rep(NA, m)
test_wald<-rep(NA,m)
oldpred<-predict(mod_pois, newdata=data.frame(x_1=4), type="link", se.fit=T)
for (i in 1:m){
  mod_ci<-glm(formula(mod_pois), family=family(mod_pois), data=cbind(y[i],design))
  pred<-predict(mod_ci, newdata=data.frame(x_1=4), type="link", se.fit=T)
  conf_int_prof<-mypredict(mod_ci, new_data=data.frame(x_1=4), alpha=0.05)
  conf_int_wald<-c(pred$fit-se*pred$se.fit, pred$fit+se*pred$se.fit)
  test_proflik[i] <-ifelse(conf_int_prof[1]<=oldpred$fit & oldpred$fit<=conf_int_prof[2], 1, 0)
  test_wald[i] <-ifelse(conf_int_wald[1]<=oldpred$fit & oldpred$fit<=conf_int_wald[2], 1, 0)
}
sum(test_proflik)/m
sum(test_wald)/m




# Manuell test intercept-modell

alpha<-0.05

expit<-function(x){return(1/(1+exp(-x)))}

logit<-function(x){return(log(x/(1-x)))}

loglik<-function(beta){return(k * beta - n * log(1 + exp(beta)))}

prof <- function(x){
  result <- loglik(x) - loglik(log(k/(n-k))) + qchisq(1-alpha, 1)/2
  return(result)
}

wald <- function(n,k){
  p_hat<-k/n
  width <- (1/n) * qnorm(alpha/2) * sqrt(k * ((n - k) / n))
  return(c(p_hat-width, p_hat+width))
}


n<-20
k_obs<-rep(1:(n-1))
p_lower<-rep(NA, n-1)
p_upper<-rep(NA, n-1)
k<-1
for(i in 1:(n-1)){
  p_lower[i]<-expit(uniroot(prof, lower=-10, upper=logit(k/n), extendInt="upX")$root)
  p_upper[i]<-expit(uniroot(prof, lower=logit(k/n), upper=10, extendInt="downX")$root)
  k<-k+1
}

x_wald<-rep(0:n)
p_wald_lower<-rep(NA, (n-1))
p_wald_upper<-rep(NA, (n-1))
for(i in 1:(n-1)){
  p_wald_upper[i]<-wald(n,i)[1]
  p_wald_lower[i]<-wald(n,i)[2]
}

plot(NULL, xlab="k", ylab="p_hat", main="n=20, red=prof", ylim=c(-0.1,1.1),xlim=c(0,n))
points(k_obs, p_wald_upper)
points(k_obs, p_wald_lower)
points(k_obs, p_lower, col="red") #Prof
points(k_obs, p_upper, col="red")

#Coverage test
p<-5/13
abline(h=p) #prof: observing 4 to 12 gives p in conf.int. Wald: 5 to 11 gives p in conf.int
m<-10000000
x<-rbinom(m, n, p)
sum(ifelse((x>=4)&(x<=12), 1, 0))/m #Nominal coverage for prof
sum(ifelse((x>=5)&(x<=11), 1, 0))/m #Nominal coverage for wald

#Coverage test 2
p<-12/13
abline(h=p) #prof: observing 16 to 19 gives p in conf.int. Wald: 15 to 19 gives p in conf.int
m<-10000000
x<-rbinom(m, n, p)
sum(ifelse((x>=16)&(x<=19), 1, 0))/m #Nominal coverage for prof
sum(ifelse((x>=15)&(x<=19), 1, 0))/m #Nominal coverage for wald

#Coverage test 3
p<-1/13
abline(h=p) #prof: observing 1 to 4 gives p in conf.int. Wald: 1 to 5 gives p in conf.int
m<-10000000
x<-rbinom(m, n, p)
sum(ifelse((x>=1)&(x<=4), 1, 0))/m #Nominal coverage for prof
sum(ifelse((x>=1)&(x<=5), 1, 0))/m #Nominal coverage for wald

# plot of coverage for different values of p
p_vals<-seq(0.01,1, by=0.01)
wald_cov<-c()
prof_cov<-c()
a<-c()
b<-c()
for(j in 1:length(p_vals)){
  a<-ifelse((p_lower<=p_vals[j]) & (p_upper>=p_vals[j]), 1, 0)
  b<-ifelse((p_wald_lower<=p_vals[j]) & (p_wald_upper>=p_vals[j]), 1, 0)
  dens<-dbinom(rep(1:(n-1)),n,p_vals[j])
  prof_cov[j]<-sum(dens[a==1])
  wald_cov[j]<-sum(dens[b==1])
}

plot(NULL, xlim=c(0,1), ylim=c(0,1), main="exact coverage for each value of p, n=20, prof=red", ylab="coverage", xlab="p")
lines(p_vals, wald_cov)
lines(p_vals, prof_cov, col="red")
abline(h=0.95)
