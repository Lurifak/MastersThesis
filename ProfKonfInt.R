mypredict <- function(mod, newdata, alpha=0.05){
  
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
  ind <- first_nonz_ind(x_0)
  a <- as.numeric(predict.glm(mod, newdata=data.frame(student="No", balance=2000, income=40000), type="link"))
  #Defining function giving the value of eta for the desired quantile by the inverted LRT 
  cutoff<-function(x)
  #{return(proflik(x) - as.numeric(logLik(mod)) + qchisq(alpha, 1)/2)}
  {return(proflik(x) - proflik(a) + qchisq(alpha, 1)/2)}
  
  #Extracting response and design matrix
  mat <- model.matrix(formula(mod), X)
  x_0 <- model.matrix(formula(mod), newdata)
  
  #Name of response variable
  resp_name <- all.vars(formula(mod))[1]
  
  #Response column
  response <- X[,resp_name]
  
  x_0 <- as.numeric(newdata)
  
  #In case of the response variable being binary it has to be converted from (1 or 2) to (0 or 1):
  if(is.factor(response)){response <- as.numeric(as.factor(response)) - 1}
  
  #Defining X_star from stackexchange
  X_star<-subset(mat, select = -(ind+1)) - (1/x_0[ind])*mat[,ind+1]%*%t(c(1,x_0[-ind]))
  
  #Defining the function for the profile likelihood given eta_0
  proflik<-function(eta_0){
    offset_star <- as.vector((1/x_0[ind]) * eta_0 * mat[,ind+1])
    new_mod<-glm.fit(X_star, response, offset=offset_star, family = family(mod))
    log_lik_new<-length(new_mod$coef) - new_mod$aic/2 #AIC := 2p - 2 ln L <=> ln L = p - AIC/2
    return(log_lik_new)
  }
  
  # Todo - select interval in a smarter way. Probably select as some function of the wald interval endpoints
  # Currently just extending the interval as a test
  left<-uniroot(cutoff, lower=-5, upper=0)$root
  right<-uniroot(cutoff, lower=0, upper=10)$root
  return(c(left, right))
}

#Test
library(ISLR)
dataset<-ISLR::Default
mod_test <- glm(default~., family=binomial(link = "logit"), data=as.data.frame(dataset))
summary(mod_test)
pred <- predict.glm(mod_test, newdata=data.frame(student="Yes", balance=2000, income=40000), type="link", se.fit=T)
conf_int <- c(pred$fit-1.96*pred$se.fit, pred$fit+1.96*pred$se.fit)

#95 % Wald for eta_0
conf_int

#Supposed 95% profile conf. int for eta_0 - not working
mypredict(mod_test, newdata=c(0,2000,40000), alpha=0.05)

# Test plot for the root - not achieving anything close to 0 in a much larger interval than the wald 95% quantiles
# Have to run line 62 to 80 and then everything inside the mypredict function to be able to run the lines below
mod<-mod_test; alpha=0.05; newdata<-data.frame(student="Yes", balance=2000, income=40000)

å<-lapply(seq(-10, 10, by=0.1), cutoff)
plot(seq(-10, 10, by=0.1), å, xlab="eta_0", ylab="f(eta_0)")
abline(h=0)
  