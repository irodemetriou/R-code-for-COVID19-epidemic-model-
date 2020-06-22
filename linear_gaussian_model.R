# Xn = rho*Xn-1 + tau*Vn 
# Yn = Xn + sigma*Wn 
install.packages("ggplot")
install.packages("invgamma")
library(ggplot2)
library(invgamma)

rho <- 0.95
tau <- 1
sigma <- 1 
T <- 100 
N <- 1000

Y <- vector("numeric")
X <- vector("numeric")

X[1] <- rnorm(1)
Y[1] <- X[1] + sigma*rnorm(1)

for (t in 1:T){
  X[t+1] <- rho*X[t] + tau*rnorm(1)
  Y[t+1] <- X[t+1] + sigma*rnorm(1) 
}

#####################
#KALMAN FILTER

muk <- vector("numeric") 
sigmak <- vector("numeric")
mk <- vector("numeric")
sk <- vector("numeric")
klm <- vector("numeric")
kf_mean <- vector("numeric")
kf_var <- vector("numeric") 
marginal_kf <- vector("numeric")

kf_mean[1] <- 0
kf_var[1] <- 1
marginal_kf[1] <- dnorm(Y[1], mean = 0, sd = sqrt(sigma^2+1))

set.seed(2020)
for (j in 1:T){
  muk[j] <- rho*kf_mean[j]
  sigmak[j] <- (rho^2)*kf_var[j] + tau^2
  mk[j] <- muk[j]
  sk[j] <- sigmak[j] + sigma^2
  klm[j] <- sigmak[j]/sk[j]
  
  kf_mean[j+1] <- muk[j] + klm[j]*(Y[j+1] - mk[j])
  kf_var[j+1] <- sigmak[j] - klm[j]*sigmak[j]
  marginal_kf[j+1] <- dnorm(Y[j+1], mean = mk[j], sd = sqrt(sk[j]))*marginal_kf[j]
}

x <- 1:(T+1)
df <- data.frame(x, X, kf_mean)
df2 <- data.frame(x,Y)

# TRAJECTORY AND KF MEAN 
ggplot(df, aes(x)) +
  geom_line(aes(y=X, colour="True trajectory")) + 
  geom_line(aes( y=kf_mean, colour="KF mean")) +
  ggtitle("True state trajectory and Kalman Filter mean") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("values") + theme(legend.position = "bottom")

ggplot(df2, aes(x)) +
  geom_line(aes( y=Y), colour="black") +
  ggtitle("Observations") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("Y") 



  
######################################
# SIS implementation using bootstrap proposal 
w <- matrix(0, T+1 , N)
W <- matrix(0, T+1, N)
XS <- matrix(0, T+1, N)
ESS <-vector("numeric")
sis_mean <- vector("numeric")
sis_var <- vector("numeric")
marginal_sisb <- matrix(0, T+1)

set.seed(2020)
for (i in 1:(T+1)){
  if (i == 1){
    XS[i,] <- rnorm(N)
    #logw[i,] <- -0.5*log(2*pi*sigma^2)*matrix(1,N) -0.5*((Y[i]*matrix(1,N) - XS[i,])^2)/(sigma^2)
    w[i,] <- dnorm(Y[i], mean = XS[i,], sd = sigma)
    W[i,] <- w[i,]/sum(w[i,])
    marginal_sisb[i] <- sum(w[i,])/N
  } else{
    XS[i,] <- rnorm(N, mean = rho*XS[i-1,], sd = tau)
    #logw[i,] <- logw[i-1,] -0.5*log(2*pi*sigma^2)*matrix(1,N) -0.5*((Y[i]*matrix(1,N) - XS[i,])^2)/(sigma^2)
    
    wx <- dnorm(Y[i], mean = XS[i,], sd = sigma) 
    w[i,] <- wx*W[i-1,]
    W[i,] <- w[i,]/sum(w[i,])
    marginal_sisb[i] <- marginal_sisb[i-1]*sum(W[i,])/N
  }
  # convert the weights and normalise  
  ESS[i] <- 1/(sum(W[i,]^2))
  sis_mean[i] <- sum(W[i,]*XS[i,])
  sis_var[i] <- sum(W[i,]*XS[i,]^2) - sis_mean[i]^2
}



######################################
# SIS implementation using optimal proposal 
w2 <- matrix(0, T+1 , N)
W2 <- matrix(0, T+1, N)
XS2 <- matrix(0, T+1, N)
marginal <- matrix(0, T+1)
sislikelihood2 <- vector("numeric")
ESS2 <-vector("numeric")
sis_mean2 <- vector("numeric")
sis_var2 <- vector("numeric")

ss <- (sigma^2)/(1+sigma^2)
set.seed(2020)
for (i in 1:(T+1)){
  if (i == 1){
    XS2[i,] <- rnorm(N, mean = Y[i]/(sigma^2+1), sd = sqrt(sigma^2/(sigma^2+1)))
    #logw2[i,] <- -0.5*log(2*pi*(tau^2 + sigma^2))*matrix(1,N) - 0.5*((Y[i]*matrix(1,N))^2)/(tau^2 + sigma^2)
    w2[i,] <- (dnorm(XS2[i,], mean = 0, sd = 1)*dnorm(Y[i], mean = XS2[i,], sd = sigma))/dnorm(XS2[i,], mean = Y[i]/(sigma^2+1), sd = sqrt(sigma^2/(sigma^2+1)))
    W2[i,] <- w2[i,]/sum(w2[i,])
    marginal[i] <- sum(w2[i,])/N
  } else{ 
    XS2[i,] <- rnorm(N, mean = (rho*sigma^2*XS2[i-1,] + tau^2*Y[i])/(sigma^2+tau^2), sd = sqrt((sigma^2*tau^2)/(sigma^2+tau^2)))
    #XS2[i,] <- ss*(rho*XS[i-1,] + Y[i]*matrix(1,N)/(sigma^2) ) + sqrt(ss)*rnorm(N)
    #logw2[i,] <- logw2[i-1,] -0.5*log(2*pi*(tau^2 + sigma^2))*matrix(1,N) -0.5*((Y[i]*matrix(1,N) - rho*XS2[i-1,])^2)/(tau^2 + sigma^2)
    w2x <- dnorm(Y[i], rho*XS2[i-1,], sd = sqrt(sigma^2+tau^2))
    w2[i,] <- W2[i-1,]*w2x
    W2[i,] <- w2[i,]/sum(w2[i,])
    marginal[i] <- marginal[i-1]*sum(W2[i,])/N
  }
  ESS2[i] <- 1/(sum(W2[i,]^2))
  sis_mean2[i] <- sum(W2[i,]*XS2[i,])
  sis_var2[i] <- sum(W2[i,]*XS2[i,]^2) - sis_mean2[i]^2
}

x <- 1:(T+1)

# mean plots for SIS prior and optimal proposal 
df_means <- data.frame(x, kf_mean, sis_mean, sis_mean2)
ggplot(df_means, aes(x)) +
  geom_line(aes( y=kf_mean, colour="KF mean"), linetype = "dashed") + 
  geom_line(aes( y=sis_mean, colour="bootstrap proposal")) +
  geom_line(aes( y=sis_mean2, colour="optimal proposal"), linetype = "dashed") + 
  ggtitle("Estimates of filter mean") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") + theme(legend.position = "bottom")


# variance plots for SIS prior and optimal proposal 
df_var <- data.frame(x, kf_var, sis_var, sis_var2)
ggplot(df_var, aes(x)) +
  geom_line(aes(y=kf_var, colour="KF variance")) + 
  geom_line(aes(y=sis_var, colour="bootstrap proposal")) +
  geom_line(aes(y=sis_var2, colour="optimal proposal"), linetype = "dashed") + 
  ggtitle("Estimates of filter variance") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("V(X)") + theme(legend.position = "bottom")


# ESS plot for SIS prior and optimal proposal 
dfess <- data.frame(x, SIS_ESS = ESS, SIS2_ESS = ESS2)
ggplot(dfess, aes(x)) +
  geom_line(aes( y=ESS, colour="bootstrap proposal")) + 
  geom_line(aes( y=ESS2, colour="optimal proposal"), linetype = "dashed") +
  ggtitle("Estimates of Effective Sample Size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("values") + theme(legend.position = "bottom")


# Marginals of KF and SIS for bootstrap and optimal proposal 
df_margs <- data.frame(x, log(marginal_kf), log(marginal_sisb), log(marginal))
ggplot(df_margs, aes(x)) +
  geom_line(aes( y=log(marginal_sisb), colour="Bootstrap proposal")) + 
  geom_line(aes( y=log(marginal_kf), colour="Kalman Filter"), linetype = "dashed") + 
  geom_line(aes( y=log(marginal), colour="Optimal Proposal"), linetype = "dashed") +
  ggtitle("Estimates of Marginal") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("values") + theme(legend.position = "bottom")


  


######################################
# SIR Implementation using bootstrap proposal 
wr <- matrix(0, T+1 , N)
wnew <- matrix(0, T+1, N)
XR <- matrix(0, T+1, N)
xbar <- matrix(0, T+1, N)
logmarginal_lik <- vector("numeric")
RESS <-vector("numeric")
sir_mean <- vector("numeric")
sir_var <- vector("numeric")

for (i in 1:(T+1)){
  if (i == 1){
    XR[i,] <- rnorm(N)
    # unormalised weights
    wr[i,] <- dnorm(Y[i], mean = XR[i,], sd = sigma)
    w <- wr[i,]/(sum(wr[i,]))
    wnew[i,] <- w
    
    logmarginal_lik[i] <- log(sum(wnew)/N)
    
    # resampling step 
    resampling <- sample(XR[i,], size = N, replace = TRUE, prob = w)
    xbar[i,] <- resampling
    
  } else{
    XR[i,] <- rnorm(N, mean = rho*XR[i-1,], sd= tau)
    # unormalised weights
    wr[i,] <- dnorm(Y[i], mean = XR[i,], sd = sigma)
    w <- wr[i,]/(sum(wr[i,]))
    wnew[i,] <- w
    logmarginal_lik[i] <- logmarginal_lik[i-1] + log(sum(wnew)/N)
    
    for (j in 1:i){
      # resampling step 
      resampling <- sample(XR[j,], size = N, replace = TRUE, prob = w)
      XR[j,] <- resampling 
    }
    xbar[i,] <- XR[i,]
  }
  RESS[i] = 1/sum(wnew[i,]^2)
  sir_mean[i] <-sum(wnew[i,]*XR[i,]) 
  sir_var[i] <- sum(wnew[i,]*(XR[i,])^2) - sir_mean[i]^2 
}

# plot of simulated trajectories for bootstrap 
plot(XR[,1], type = "l", xlab = "Time", ylab = "X", ylim = c(-11,11), 
     col = "lightgrey", main = "Bootstrap Proposal for SIR")
for (i in 2:ncol(XR)){
  lines(XR[,i], lwd = 0.3, col = "lightgrey")
}
lines(kf_mean, lwd = 1, col="purple")
plot(X, type = "l", main = "True trajectories")





######################################
# SIR Implementation using optimal proposal
wr2 <- matrix(0, T+1 , N)
wnew2 <- matrix(0, T+1, N)
XR2 <- matrix(0, T+1, N)
xbar2 <- matrix(0, T+1, N)
logmarginal_lik2 <- vector("numeric")
RESS2 <-vector("numeric")
sir_mean2 <- vector("numeric")
sir_var2 <- vector("numeric")

for (i in 1:(T+1)){
  if (i == 1){
    XR2[i,] <- rnorm(N, mean = ss*Y[i]/sigma^2, sd = sqrt(ss))
    # unormalised weights
    wr2[i,] <- dnorm(XR2[i,], mean = 0, sd = 1)*dnorm(Y[i], mean= XR2[i,], sd = sigma)/dnorm(XR2[i,], mean = ss*Y[i]/sigma^2, sd = sqrt(ss))
    w <- wr2[i,]/(sum(wr2[i,]))
    wnew2[i,] <- w
    logmarginal_lik2[i] <- log(sum(wnew2)/N)
    
    # resampling step 
    resampling <- sample(XR2[i,], size = N, replace = TRUE, prob = w)
    xbar2[i,] <- resampling
    
  } else{
    XR2[i,] <- rnorm(N, mean = ss*(rho*XR2[i-1,] + Y[i]/sigma^2), sd = sqrt(ss))
    # unormalised weights
    wr2[i,] <- dnorm(Y[i], mean = rho*XR2[i-1,], sd = sqrt(sigma^2 + tau^2))
    w <- wr2[i,]/(sum(wr2[i,]))
    wnew2[i,] <- w
    logmarginal_lik2[i] <- logmarginal_lik2[i-1] + log(sum(wnew2)/N)
    
    for (j in 1:i){
      # resampling step 
      resampling <- sample(XR2[j,], size = N, replace = TRUE, prob = w)
      XR2[j,] <- resampling 
    }
    xbar2[i,] <- XR2[i,]
  }
  RESS2[i] = 1/sum(wnew2[i,]^2)
  sir_mean2[i] <-sum(wnew2[i,]*XR2[i,]) 
  sir_var2[i] <- sum(wnew2[i,]*(XR2[i,])^2) - sir_mean2[i]^2 
}

# plot of simulated trajectories for optimal 
plot(XR2[,1], type = "l", xlab = "Time", ylab = "X", ylim = c(-11,11), 
     col = "lightgrey", main = "Optimal Proposal for SIS")
for (i in 2:ncol(XR2)){
  lines(XR2[,i], lwd = 0.3, col = "lightgrey")
}
lines(kf_mean, lwd = 1, col="purple")

# FILTER MEAN, VARIANCE AND ESS PLOT FOR SIR OPTIMAL AND BOOTSTRAP 
data <- 1:(T+1)
meandframe <- data.frame(data, kf_mean, sir_mean, sir_mean2)
ggplot(meandframe, aes(data)) +
  geom_line(aes( y=kf_mean, colour="KF mean"), linetype = "dashed") + 
  geom_line(aes( y=sir_mean, colour="bootstrap proposal")) +
  geom_line(aes( y=sir_mean2, colour="optimal proposal"), linetype = "dashed") + 
  ggtitle("Estimates of filter mean (SIR)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") + theme(legend.position = "bottom")


sirdframe <- data.frame(data, kf_var, sir_var, sir_var2)
ggplot(sirdframe, aes(data)) +
  geom_line(aes( y=kf_var, colour="KF mean"), linetype = "dashed") + 
  geom_line(aes( y=sir_var, colour="bootstrap proposal")) +
  geom_line(aes( y=sir_var2, colour="optimal proposal"), linetype = "dashed") + 
  ggtitle("Estimates of filter variance (SIR)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("V(X)") + theme(legend.position = "bottom")


lik <- data.frame(data, marginal_kf ,logmarginal_lik, logmarginal_lik2)
ggplot(lik, aes(data)) + 
  geom_line(aes(y = logmarginal_lik, colour = "bootstrap"), linetype = "dashed") +
  geom_line(aes(y = logmarginal_lik2, colour = "optimal")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Estimates for Marginal Likelihoods") +
  xlab("time") + ylab("Log Marginal") + theme(legend.position = "bottom")
  
  
essdframe <- data.frame(data, RESS, RESS2)
ggplot(essdframe, aes(data)) +
  geom_line(aes( y=RESS, colour="bootstrap proposal")) + 
  geom_line(aes( y=RESS2, colour="optimal proposal"), linetype = "dashed") +
  ggtitle("Estimates of Effective Sample Size (SIR)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("ESS") + theme(legend.position = "bottom")
  
  
  


################################################

# SIR and Particle Marginal Metropolis Hasting (PMMH)

SIR <- function(Y, N, T, rho, tau, sigma){
  
  wr2 <- matrix(0, T+1 , N)
  wnew2 <- matrix(0, T+1, N)
  XR2 <- matrix(0, T+1, N)
  xbar2 <- matrix(0, T+1, N)
  marginal_lik2 <- vector("numeric")
  logmarginal_lik2 <- vector("numeric")
  RESS2 <-vector("numeric")
  sir_mean2 <- vector("numeric")
  sir_var2 <- vector("numeric")
  ss <- (sigma^2)/(1+sigma^2)

  for (i in 1:(T+1)){
    if (i == 1){
      XR2[i,] <- rnorm(N, mean = ss*Y[i]/sigma^2, sd = sqrt(ss))
      # unormalised weights
      wr2[i,] <- dnorm(XR2[i,], mean = 0, sd = 1)*dnorm(Y[i], mean= XR2[i,], sd = sigma)/dnorm(XR2[i,], mean = ss*Y[i]/sigma^2, sd = sqrt(ss))
      w <- wr2[i,]/(sum(wr2[i,]))
      wnew2[i,] <- w
      marginal_lik2[i] <- sum(wr2[i,])/N
      logmarginal_lik2[i] <- log(sum(wr2[i,])/N)
      
      # resampling step 
      resampling <- sample(XR2[i,], size = N, replace = TRUE, prob = w)
      xbar2[i,] <- resampling
      
    } else{
      XR2[i,] <- rnorm(N, mean = ss*(rho*XR2[i-1,] + Y[i]/sigma^2), sd = sqrt(ss))
      # unormalised weights
      wr2[i,] <- dnorm(Y[i], mean = rho*XR2[i-1,], sd = sqrt(sigma^2 + tau^2))
      w <- wr2[i,]/(sum(wr2[i,]))
      wnew2[i,] <- w
      marginal_lik2[i] <- marginal_lik2[i-1]*sum(wr2[i,])/N
      logmarginal_lik2[i] <- logmarginal_lik2[i-1] + log(sum(wr2[i,])/N)
      
      for (j in 1:i){
        # resampling step 
        resampling <- sample(XR2[j,], size = N, replace = TRUE, prob = w)
        XR2[j,] <- resampling 
      }
      xbar2[i,] <- XR2[i,]
    }
    RESS2[i] = 1/sum(wnew2[i,]^2)
    sir_mean2[i] <-sum(wnew2[i,]*XR2[i,]) 
    sir_var2[i] <- sum(wnew2[i,]*(XR2[i,])^2) - sir_mean2[i]^2 
  }
list("Mean" = sir_mean2, "Variance" = sir_var2,"X" = XR2, "X.bar" = xbar2, "ESS" = RESS2, "Weights" = wnew2,
    "Marginal_Likelihood "= marginal_lik2  ,"Log_Marginal_Likelihood" = logmarginal_lik2)
}

sir_ex <- SIR(Y, N, T, rho, tau, sigma)



#### implement PMMH 
# n: iterations 
# b: burn in times for convergence 
# N: number of samples 
# thetas: vector of parameters
# lamda: step size 
# ap: acceptance ratio for Metropolis Hastings 
n = 200
b = 100
lamda = 0.2
theta0 =c(runif(1,-1,1), sqrt(rinvgamma(1,1,1)))

pmmh <- function(N, T, n, b, Y, theta0, lamda){
  ap <- vector("numeric")
  thetas <- matrix(0, n+b, 2)
  thetas[1,] <- theta0
  X <- matrix(0, n+b, T+1)
  Xbar <- matrix(0, n+b, T+1)
  Likelihood <- vector("numeric")
  Loglikelihood <- vector("numeric")
  
  for (i in 1:(n+b)){
    if (i == 1){
      sir <- SIR(Y, N, T, rho = theta0[1], tau = theta0[2], sigma = 0.5)
      res <- sample(N, 1, prob = sir$Weights[T+1,])
      X[i,] <- sir$X[,res]
      Xbar[i,] <- sir$X.bar[,res]
      Likelihood[i] <- sir$Marginal_Likelihood[i]
      Loglikelihood[i] <- sir$Log_Marginal_Likelihood[i]
    } else {
      # sample proposal s.t. theta' = theta(n-1) + lamda*N(0,1)
      tdash <- thetas[i-1,] + rnorm(2, mean = 0, sd = lamda)
      sir <- SIR(Y, N, T, rho = tdash[1], tau = tdash[2], sigma = 0.5)
      res <- sample(N, 1, prob = sir$Weights[T+1,])
      
      Likelihood[i] <- sir$Marginal_Likelihood[T+1]
      Loglikelihood[i] <- sir$Log_Marginal_Likelihood[T+1]
      
      # MH acceptance probability 
      num <- Likelihood[i]*dunif(tdash[1], -1, 1)*dinvgamma(tdash[2]^2, 1, 1)*
        dnorm(thetas[i-1,1], mean =  tdash[1], sd = lamda)*dnorm(thetas[i-1,2], mean =  tdash[2], sd = lamda)
      den <- Likelihood[i-1]*dunif(thetas[i-1,1], -1, 1)*dinvgamma(thetas[i-1,2]^2, 1, 1)*
        dnorm(tdash[1], mean =  thetas[i-1,1], sd = lamda)*dnorm(tdash[2], mean =  thetas[i-1,2], sd = lamda)
     
      p <- ifelse(is.na(num/den), 1, num/den)
      a <- min(1, p)
      ap <- c(ap, a)
      
      U <- runif(1)
      if (U <= a){
        thetas[i,] <- tdash
        res <- sample(N, 1, prob = sir$Weights[T+1,])
        X[i,] <- sir$X[,res]
        Xbar[i,] <- sir$X.bar[,res]
      } else {
        thetas[i,] <- thetas[i-1,] 
        X[i,] <- X[i-1,]
        Xbar[i,] <- Xbar[i-1,]
      }
    }
  }
  # discard first b iterations 
  X <- X[c((b+1):(n+b)),]
  Xbar <- Xbar[c((b+1):(n+b)),]
  thetas <- thetas[c((b+1):(n+b)),]^2
  ap <- ap[c((b+1):(n+b))]
  Likelihood <- Likelihood[c((b+1):(n+b))]
  Loglikelihood <- Loglikelihood[c((b+1):(n+b))]
  
  list("X"= X, "Xbar" = Xbar, "thetas" = thetas, "ap" =  ap, "Likelihood"=Likelihood, "Loglikelihood" = Loglikelihood)
}

pex <- pmmh(N, T, n, b, Y, theta0, lamda)

hist(pex$thetas[,1],col=rgb(0,0,1,1/4),breaks = 40,freq = FALSE, xlab = expression(paste(rho)),
     main = expression(paste("Histogram of ", rho)))

hist(pex$thetas[,2],col=rgb(0,0,1,1/4),breaks = 40,freq = FALSE, xlab = expression(paste(tau^2)),
     main = expression(paste("Histogram of ", tau^2)), xlim = c(0,3))



### Homework page 37
# fix rho = 1 and let tau^2, sigma^2 have prior distirbution Inverse Gamma 

n = 1000
b = 100
lamda = 0.2
rho = 1 
theta0 =c(sqrt(rinvgamma(1,1,1)), sqrt(rinvgamma(1,1,1))) # tau^2, sigma^2 

pmmh2 <- function(N, T, n, b, Y, theta0, lamda){
  ap <- vector("numeric")
  thetas <- matrix(0, n+b, 2)
  thetas[1,] <- theta0
  X <- matrix(0, n+b, T+1)
  Xbar <- matrix(0, n+b, T+1)
  Likelihood <- vector("numeric")
  Loglikelihood <- vector("numeric")
  
  for (i in 1:(n+b)){
    if (i == 1){
      sir <- SIR(Y, N, T, rho = 1, tau = theta0[1], sigma = theta0[2])
      res <- sample(N, 1, prob = sir$Weights[T+1,])
      X[i,] <- sir$X[,res]
      Xbar[i,] <- sir$X.bar[,res]
      Likelihood[i] <- sir$Marginal_Likelihood[i]
      Loglikelihood[i] <- sir$Log_Marginal_Likelihood[i]
    } else {
      # sample proposal s.t. theta' = theta(n-1) + lamda*N(0,1)
      tdash <- thetas[i-1,] + rnorm(2, mean = 0, sd = lamda)
      sir <- SIR(Y, N, T, rho = 1, tau = tdash[1], sigma = abs(tdash[2]))
      res <- sample(N, 1, prob = sir$Weights[T+1,])
      
      Likelihood[i] <- sir$Marginal_Likelihood[T+1]
      Loglikelihood[i] <- sir$Log_Marginal_Likelihood[T+1]
      
      # MH acceptance probability 
      num <- Likelihood[i]*dinvgamma(tdash[1]^2, 1, 1)*dinvgamma(tdash[2]^2, 1, 1)*
        dnorm(thetas[i-1,1], mean =  tdash[1], sd = lamda)*dnorm(thetas[i-1,2], mean =  tdash[2], sd = lamda)
      den <- Likelihood[i-1]*dinvgamma(thetas[i-1,1]^2, 1, 1)*dinvgamma(thetas[i-1,2]^2, 1, 1)*
        dnorm(tdash[1], mean =  thetas[i-1,1], sd = lamda)*dnorm(tdash[2], mean =  thetas[i-1,2], sd = lamda)
      
      p <- ifelse(is.na(num/den), 1, num/den)
      a <- min(1, p)
      ap <- c(ap, a)
      
      U <- runif(1)
      if (U <= a){
        thetas[i,] <- tdash
        res <- sample(N, 1, prob = sir$Weights[T+1,])
        X[i,] <- sir$X[,res]
        Xbar[i,] <- sir$X.bar[,res]
      } else {
        thetas[i,] <- thetas[i-1,] 
        X[i,] <- X[i-1,]
        Xbar[i,] <- Xbar[i-1,]
      }
    }
  }
  # discard first b iterations 
  X <- X[c((b+1):(n+b)),]
  Xbar <- Xbar[c((b+1):(n+b)),]
  thetas <- thetas[c((b+1):(n+b)),]^2
  ap <- ap[c((b+1):(n+b))]
  Likelihood <- Likelihood[c((b+1):(n+b))]
  Loglikelihood <- Loglikelihood[c((b+1):(n+b))]
  
  list("X"= X, "Xbar" = Xbar, "thetas" = thetas, "ap" =  ap, "Likelihood"=Likelihood, "Loglikelihood" = Loglikelihood)
}


pex2 <- pmmh2(N, T, n, b, Y, theta0, lamda)

hist(pex2$thetas[,1],col=rgb(0,0,1,1/4),breaks = 40,freq = FALSE, xlab = expression(paste(tau^2)),
     main = expression(paste("Histogram of ", tau)))

hist(pex2$thetas[,2],col=rgb(0,0,1,1/4),breaks = 40,freq = FALSE, xlab = expression(paste(sigma^2)),
     main = expression(paste("Histogram of ", sigma^2)), xlim = c(0,5))


thetasdframe <- data.frame(ind = 1:1000, taunew = pex2$thetas[,1], sigmanew = pex2$thetas[,2])
ggplot(thetasdframe, aes(x = taunew)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", tau^2))) +
  xlab(expression(paste(tau^2))) + ylab("Density") + theme(legend.position = "bottom")

ggplot(thetasdframe, aes(x = sigmanew)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", sigma^2))) +
  xlab(expression(paste(sigma^2))) + ylab("Density") + theme(legend.position = "bottom")


acf(pex2$thetas[,1], main = expression(paste("ACF for ", tau^2)))
acf(pex2$thetas[,2], main = expression(paste("ACF for ", sigma^2)))



