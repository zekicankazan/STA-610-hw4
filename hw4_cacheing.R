## Code to produce the CSV files of posterior samples for every combination 
## of rho and eta

library(mvtnorm)
library(brms)
library(readr)

rho_vec <- c(seq(-0.95, -0.15, 0.05), seq(-0.1, 0.95, 0.05))
eta_vec <- c(0.1, seq(0.2, 0.8, 0.2), 1:5, 10)
rho_eta_combs <- expand.grid(eta = eta_vec, rho = rho_vec)

for(i in 1:nrow(rho_eta_combs)){
  J <- 8 #Number of groups
  n <- rep(25,J) #Observations per group
  beta0 <- 5 #Fixed effect intercept
  beta1 <- 2 #Fixed effect slope
  tau0 <- 2 #Random intercept standard deviation
  tau1 <- 1 #Random slope standard deviation
  sigma <- 1 #Residual standard deviation
  rho <- rho_eta_combs$rho[i]
  eta <- rho_eta_combs$eta[i]
  
  tau_mat <- diag(c(tau0, tau1))
  Omega <- matrix(c(1,rho,rho,1),ncol=2)
  Sigma <- tau_mat %*% Omega %*% tau_mat
  
  set.seed(1)
  beta <- rmvnorm(J, mean = c(beta0, beta1),  sigma = Sigma) #Sample coefficients
  
  x <- c(); y <- c()
  for(j in 1:J){
    xj <- rnorm(n[j])
    yj <- rnorm(n[j], beta[j,1] + beta[j,2]*xj,sigma)
    x <- c(x,xj); y <- c(y,yj)
  }
  
  data <- data.frame(x, y, group = rep(1:J, n))
  
  Int_prior <- paste0("normal(",as.character(mean(y)),", ",
                      as.character(10*sd(y)),")")
  b_prior <- paste0("normal(0, ", as.character(2.5*sd(y)/sd(x)),")")
  L_prior <- paste0("lkj_corr_cholesky(", as.character(eta), ")")
  
  prior <- c(set_prior(Int_prior, class = "Intercept"),
             set_prior(b_prior, class = "b"),
             set_prior(L_prior, class = "L"))
  
  mod <- brm(y ~ x + (x | group), data = data, cores = 4, chains = 4, 
             warmup = 1000, iter = 5000, prior = prior, seed = 1)
  
  draws <- as_draws_df(mod)
  post <- round(draws[,"cor_group__Intercept__x"],3) #Round to make CSV smaller
  
  write_csv(post, paste0("STA-610-hw4/hw4_cache/eta_",
                         as.character(eta), "_rho_", as.character(rho), ".csv"))
}