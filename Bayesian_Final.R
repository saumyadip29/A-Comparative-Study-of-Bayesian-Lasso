rm(list=ls())
library(mvtnorm)
library(truncdist)
library(tmvmixnorm)
library(monomvn)
library(lubridate)
library(lars)
library(lasso2)





sigma2 <- 1  #chosen \sigma^2
beta <- c(3, 1.5, 2, 0, 1, 0, 0, 0)  #Chosen beta
n <- 400  #sample size
p <- length(beta)
mu <- rep(0,p)
sig <- matrix(rep(0.95, p*p), nrow = p) + 0.05*diag(p) #variance covariance matrix
X <- rmvnorm(n, mu, sig)
eps <- rmvnorm(1, rep(0, n), sigma2*diag(n))
y <- X%*%beta + t(eps) #generated y



lambda_in <- 2 #initial choice of \lambda
sigma2_in <- 4 #initial choice of \sigma^2
a <- 1; b <- 10; #parameter for the \lambda prior
R <- 11000 #MCMC iterations


#Function to draw samples from truncated normal distribution
t_norm <- function(X, y, sigma2,uj){
  beta_hat <- as.vector(solve(t(X)%*%X)%*%t(X)%*%y)
  s_hat <- sigma2*solve(t(X)%*%X)
  samp <-  rtmvn(1, Mean = beta_hat, Sigma = s_hat,
                 D=diag(1,length(beta_hat)),lower = -sqrt(sigma2)*uj, 
                 upper = sqrt(sigma2)*uj,int=rep(0,length(beta_hat)))
  #samp <- rmvnorm(1, beta_hat, s_hat)
  # while (sum(abs(samp) < sigma*uj) != length(samp)){
  #   samp <- t(rmvnorm(1, beta_hat,  s_hat))
  #}
  return(samp)
}


#Function to draw samples from truncated gamma distribution
t_gam <- function(X, y, beta, uj){
  sh <- (n-1+p)/2
  ra <- (t(y-X%*%beta)%*%(y-X%*%beta))/2
  samp <- rtrunc(1, "gamma", b = 1/max((beta^2)/(uj^2)), shape = sh, rate = ra)
  
  # while(samp > 1/max((beta^2)/(uj^2))){
  #   samp <- rgamma(1, sh, ra)  
  # }
  return(1/samp)
}

burnin <- 1000
#Main MCMC iteration function
iter <- function(X, y, lambda_in, beta_in, sigma2_in, a, b, R, p){
  l <- list(sigma2_1 = vector(length = R-burnin), beta_1 = matrix(0, nrow = p, ncol = R-burnin), lambda_1 = vector( length = R-burnin))
  lambda <- lambda_in
  beta <- beta_in
  sigma2 <- sigma2_in 
  
  for(i in 1:R){
    if(i<=burnin){
      ujp <- rexp(p, lambda)
      uj <- ujp + abs(beta)/sqrt(sigma2)
      beta <- as.vector(t_norm(X, y, sigma2, uj))
      sigma2 <- t_gam(X, y, beta, uj)
      lambda <- rgamma(1, a+2*p, b+sum(abs(beta)))
    }
    else{
      ujp <- rexp(p, lambda)
      uj <- ujp + abs(beta)/sqrt(sigma2)
      
      beta <- as.vector(t_norm(X, y, sigma2, uj))
      l[[2]][,i-burnin] <- beta
      
      #print(1/max(beta^2/uj^2))
      
      sigma2 <- t_gam(X, y, beta, uj)
      l[[1]][i-burnin] <- sigma2 
      
      
      lambda <- rgamma(1, a+2*p, b+sum(abs(beta)))
      l[[3]][i-burnin] <- lambda
    }
  }
  return(l)
}


#New Bayesian Lasso Step
t1 <- 100
n_train <- 0.75*nrow(X)
lambda_m <- array(0)
sigma2_m <- array(0)
beta_m <- matrix(0, nrow = p, ncol = t1)
mse <- array(0)
shell('cls')
print("Running New Bayesian Lasso:")
pb = txtProgressBar(min = 0, max = t1, style = 3, width = t1, char = '=')
init <- numeric(t1)
end <- numeric(t1)
for(i in 1:t1){
  init[i] <- Sys.time()
  t <- sample(1:n, n_train)
  X_train <- X[t,]
  X_test <- X[-t,]
  y_train <- y[t]
  y_train <- y_train - mean(y_train)
  y_test <- y[-t]
  
  #Scaling the training and test design matrix
  X_sc_train <- scale(X_train)
  X_mu <- apply(X_train, 2, mean)
  X_sd <- apply(X_train, 2, sd)
  X_sc_test <- matrix(0, nrow = nrow(X_test), ncol = ncol(X_test))
  for(j in 1:ncol(X_test)){
    X_sc_test[,j] <- (X_test[,j]-X_mu[j])/X_sd[j]
  }
  
  #initial beta estimated from the OLS estimate 
  beta_in <- as.vector(solve(t(X_sc_train)%*%X_sc_train)%*%t(X_sc_train)%*%y_train)
  l <- 0
  l <- iter(X_sc_train, y_train, lambda_in, beta_in, sigma2_in, a, b, R, p)
  lambda_m[i] <- mean(l$lambda_1)
  sigma2_m[i] <- mean(l$sigma2_1)
  beta_hat <- apply(l$beta_1, 1, mean)
  beta_m[,i] <- beta_hat
  yhat <- as.vector(X_sc_test%*%beta_hat)
  mse[i] <- mean((y_test-yhat)^2)
  end[i] <- Sys.time()
  setTxtProgressBar(pb, i)
  time <- round(seconds_to_period(sum(end - init)), 0)
  est <- t1 * (mean(end[end != 0] - init[init != 0])) - time
  remainining <- round(seconds_to_period(est), 0)
  
  cat(paste(" // Execution time:", time,
            " // Estimated time remaining:", remainining), "")
}
#close(pb)
par(mfrow = c(3,3))
for(i in 1:8){
  plot(l$beta_1[i,], type = 'l', main = paste('Trace plot for beta_', i))
}

shell('cls')



t1 <- 100
# Old Bayesian Lasso
mse2 <- array(0)
print("Running Old Bayesian Lasso: ")
pb = txtProgressBar(min = 0, max = t1, style = 3, width = t1, char = '=')
init <- numeric(t1)
end <- numeric(t1)
for (i in 1:t1){
  init[i] <- Sys.time()
  t <- sample(1:n, n_train)
  X_train <- X[t,]
  X_test <- X[-t,]
  y_train <- y[t]
  y_train <- y_train - mean(y_train)
  y_test <- y[-t]
  
  #Scaling the training and test design matrix
  X_sc_train <- scale(X_train)
  X_mu <- apply(X_train, 2, mean)
  X_sd <- apply(X_train, 2, sd)
  X_sc_test <- matrix(0, nrow = nrow(X_test), ncol = ncol(X_test))
  for(j in 1:ncol(X_test)){
    X_sc_test[,j] <- (X_test[,j]-X_mu[j])/X_sd[j]
  }
  beta_in <- as.vector(solve(t(X_sc_train)%*%X_sc_train)%*%t(X_sc_train)%*%y_train)
  
  model <- blasso(X, y, T = 100, beta = beta_in, lambda2 = 0.01, RJ = F, rd = c(1, 4))
  beta_hat2 <- apply(model$beta, 2, mean)
  yhat2 <- X_sc_test%*%beta_hat2
  mse2[i] <- mean((y_test-yhat2)^2)
  end[i] <- Sys.time()
  setTxtProgressBar(pb, i)
  time <- round(seconds_to_period(sum(end - init)), 0)
  est <- t1 * (mean(end[end != 0] - init[init != 0])) - time
  remainining <- round(seconds_to_period(est), 0)
  
  cat(paste(" // Execution time:", time,
            " // Estimated time remaining:", remainining), "")
}
shell('cls')


#LARS
t1 <- 10
mse3 <- array(0)
print("Running Lasso: ")
pb = txtProgressBar(min = 0, max = t1, style = 3, width = t1, char = '=')
init <- numeric(t1)
end <- numeric(t1)
for(i in 1:t1){
  init[i] <- Sys.time()
  t <- sample(1:n, n_train)
  X_train <- X[t,]
  X_test <- X[-t,]
  y_train <- y[t]
  y_train <- y_train - mean(y_train)
  y_test <- y[-t]
  
  #Scaling the training and test design matrix
  X_sc_train <- scale(X_train)
  X_mu <- apply(X_train, 2, mean)
  X_sd <- apply(X_train, 2, sd)
  X_sc_test <- matrix(0, nrow = nrow(X_test), ncol = ncol(X_test))
  for(j in 1:ncol(X_test)){
    X_sc_test[,j] <- (X_test[,j]-X_mu[j])/X_sd[j]
  }
  model2 = lars(
    X_sc_train, y_train, type = "lasso",
    normalize = F, intercept =  F
  )
  lars_reg <- cv.lars(X_sc_train, y_train, K=10, type = 'lasso', intercept = F, mode = 'step')
  beta_hat_3 = predict.lars(model2, s = 8, type = "coef", mode = "step")$coefficients
  mse3[i] = mean((X_test%*%beta_hat_3 - y_test)^2)
  end[i] <- Sys.time()
  setTxtProgressBar(pb, i)
  time <- round(seconds_to_period(sum(end - init)), 0)
  est <- t1 * (mean(end[end != 0] - init[init != 0])) - time
  remainining <- round(seconds_to_period(est), 0)
  
  cat(paste(" // Execution time:", time,
            " // Estimated time remaining:", remainining), "")
}

shell('cls')



#Printing the MMSE values
print(paste("MMSE value (NBLASSO): ",round(median(mse),6)))

print(paste("MMSE value (OBLASSO): ",round(median(mse2),6)))

print(paste("MMSE value (LASSO): ",round(median(mse3),6)))

##########################################################################

#Diabetes Data Analysis

data("diabetes")
attach(diabetes)
y<-diabetes$y
X<-diabetes$x
data.all = data.frame(cbind(X, y ))
n <- dim(data.all)[1] 
set.seed(1306)
test <- sample(n, round(n/4))
data.train <- data.all[-test,]
data.test <- data.all[test,]
X_train <- X[-test,]
X_test <- X[test,]
y <- data.all$y-mean(data.all$y)
y_train <- y[-test]
y_test <- y[test]
n_train <- dim(data.train)[1] 
n_test <- dim(data.test)[1]
p<-dim(X_train)[2]

X_sc_train<-scale(X_train)
X_sc_test<-scale(X_test)
#Setting Initial Values and prior parameters
lambda_in <- 0.01
beta_in <- as.vector(solve(t(X_sc_train)%*%X_sc_train)%*%t(X_sc_train)%*%y_train)
sigma_in <- 1
a <- 1; b <- 10;
R <- 1100
burnin <- 100

t <- 1
lambda_m <- array(0)
sigma_m <- array(0)
beta_m <- matrix(0, nrow = p, ncol = t)
mse <- array(0)
pb <- txtProgressBar(min = 0,
                     max = t,
                     style = 3,
                     width = t, # Needed to avoid multiple printings
                     char = "#") 

for(i in 1:t){
  l <- 0
  l <- iter(X_sc_train, y_train, lambda_in, beta_in, sigma_in, a, b, R, p)
  lambda_m[i] <- mean(l$lambda_1)
  sigma_m[i] <- mean(l$sigma_1)
  beta_hat <- apply(l$beta_1, 1, mean)
  beta_m[,i] <- beta_hat
  yhat <- as.vector(X_sc_test%*%beta_hat)
  mse[i] <- mean((y_test-yhat)^2)
  
  setTxtProgressBar(pb, i)
}

#Finding the median values of the obtained parameters
median(lambda_m)
median(sigma_m)
apply(beta_m, 1, median)
median(mse)

par(mfrow = c(5,2))
for(i in 1:10){
  hist(l$beta_1[i,], freq = F, main = paste('Histogram plot for', colnames(X)[i]), xlab = paste(colnames(X)[i]))
  lines(density(l$beta_1[i,]))
}

par(mfrow = c(5,2))
for(i in 1:10){
  plot(l$beta_1[i,], type = 'l', main = paste('Trace plot for', colnames(X)[i]),xlab = paste(colnames(X)[i]))
}
###########################################################################

#Prostate Data Analysis

data("Prostate")
attach(Prostate)
y<-Prostate[,'lpsa']
X<-Prostate[,-9]
data.all = data.frame(cbind(X, y ))
n <- dim(data.all)[1] 
set.seed(1306)
test <- sample(n, round(n/4))
data.train <- data.all[-test,]
data.test <- data.all[test,]
X_train <- X[-test,]
X_test <- X[test,]
y <- data.all$y-mean(data.all$y)
y_train <- y[-test]
y_test <- y[test]
n_train <- dim(data.train)[1] 
n_test <- dim(data.test)[1]
p<-dim(X_train)[2]
X_sc_train<-scale(X_train)
X_sc_test<-scale(X_test)
#Setting Initial Values and prior parameters
lambda_in <- 0.01
beta_in <- as.vector(solve(t(X_sc_train)%*%X_sc_train)%*%t(X_sc_train)%*%y_train)
sigma_in <- 1
a <- 1; b <- 10;
R <-11000
burnin <- 1000

#MSE of OLS
mean(sum(y_test-X_sc_test%*%beta_in)^2)

#Iterating the Gibb's Sampler t times
t <- 1
lambda_m <- array(0)
sigma_m <- array(0)
beta_m <- matrix(0, nrow = p, ncol = t)
mse <- array(0)
pb <- txtProgressBar(min = 0,
                     max = t,
                     style = 3,
                     width = t, # Needed to avoid multiple printings
                     char = "#") 

for(i in 1:t){
  l <- 0
  l <- iter(X_sc_train, y_train, lambda_in, beta_in, sigma_in, a, b, R, p)
  lambda_m[i] <- mean(l$lambda_1)
  sigma_m[i] <- mean(l$sigma_1)
  beta_hat <- apply(l$beta_1, 1, mean)
  beta_m[,i] <- beta_hat
  yhat <- as.vector(X_sc_test%*%beta_hat)
  mse[i] <- mean((y_test-yhat)^2)
  
  setTxtProgressBar(pb, i)
}

#Finding the median values of the obtained parameters
median(lambda_m)
median(sigma_m)
apply(beta_m, 1, median)
median(mse)


# Old Bayesian Lasso
mse2 <- array(0)
for (i in 1:t){
  model <- blasso(X, y, T = 10000, beta = beta_in, lambda2 = 9, RJ = F, rd = c(1, 10))
  beta_hat2 <- apply(model$beta, 2, mean)
  yhat2 <- X_sc_test%*%beta_hat2
  mse2[i] <- mean((y_test-yhat2)^2)
}
median(mse2)

#LASSO
lars_reg <- cv.lars(X_sc_train, y_train, K=10, type = 'lasso', intercept = F, mode = 'step')
model2 <- lars(X_sc_train, y_train, type = 'lasso', normalize = F, intercept = F)
beta_hat4 <- predict.lars(model2, s = 6, type = 'coef', mode = 'step')$coefficients
beta_hat4
mean((X_sc_test%*%beta_hat4-y_test)^2)

par(mfrow = c(4,2))
for(i in 1:8){
  hist(l$beta_1[i,], freq = F, main = paste('Histogram plot for', colnames(X)[i]), xlab = paste(colnames(X)[i]))
  lines(density(l$beta_1[i,]))
}

par(mfrow = c(4,2))
for(i in 1:8){
  plot(l$beta_1[i,], type = 'l', main = paste('Trace plot for', colnames(X)[i]),xlab = paste(colnames(X)[i]))
}
