rnormgamma <- function(n, mu, lambda, alpha, beta) {
  if (length(n) > 1) 
    n <- length(n)
  tau <- rgamma(n, alpha, beta)
  x <- rnorm(n, mu, sqrt(1/(lambda*tau)))
  data.frame(tau = tau, x = x)
}

n <- 1000
m <- 5
mu_0 <- 0
kappa_0 <- 1
alpha_0 <- 1
beta_0 <- 1

count <- rep(0,m)

for(i in 1:n){
  lambda <- rgamma(1, shape=alpha_0, rate=beta_0)
  mu <- rnorm(mu_0, sd = sqrt((kappa_0*lambda)^-1))
  data <- rnorm(m, mu, sd=(sqrt(lambda^-1)))
  sorted <- sort(data)
  
  mu_n <- ((kappa_0 * mu_0) + m * mean(data))/(kappa_0 + m)
  kappa_n <- kappa_0 + m
  alpha_n <- alpha_0 + m/2
  beta_n <- beta_0 + (1/2) * sum((data - mean(data))^2) + ((kappa_0 * n) * (mean(data) - mu_0)^2)/(2*(kappa_0+m))
  
  predsims <- rt(m, 2*alpha_n)
}