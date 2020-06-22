install.packages("rjags")
library(rjags)
install.packages("devtools")
library(devtools)
install_github("lilywang1988/eSIR")
library(eSIR) 

setRepositories()
install.packages("latex2exp")
install.packages("ggpubr")
library(invgamma)
library(ggplot2)
library(latex2exp)
library(ggpubr)

# load data 
data = as.data.frame(read.csv("/Users/irodemetriou/Desktop/Year 4/SIR /covid/WHO-COVID-19-global-data.csv"))
data = data[c(-2,-4)]
rec = as.data.frame(read.csv("/Users/irodemetriou/Desktop/Year 4/SIR /covid/recovered.csv"))

Italy = data[data$Country.Name == "Italy",]
UK = data[data$Country.Name == "The United Kingdom",]
France = data[data$Country.Name == "France",]
Germany = data[data$Country.Name == "Germany",]
Spain = data[data$Country.Name == "Spain",]

Italy = Italy[4:64,]
UK = UK[2:62,]
France = France[9:69,]
Germany = Germany[5:65,]
Spain = Spain[2:62,]

Italy = data.frame(Italy,Day = 1:61 ,Cumulative.Recovered = rec$Italy.Cumulative.Recovered[1:61], Recovered = rec$Italy.Recovered[1:61])
UK = data.frame(UK,Day = 1:61 ,Cumulative.Recovered = rec$UK.Cumulative.Recovered[1:61], Recovered = rec$UK.Recovered[1:61])
France = data.frame(France, Day = 1:61 ,Cumulative.Recovered = rec$France.Cumulative.Recovered[1:61], Recovered = rec$France.Recovered[1:61])
Germany = data.frame(Germany, Day = 1:61 ,Cumulative.Recovered = rec$Germany.Cumulative.Recovered[1:61], Recovered = rec$Germany.Recovered[1:61])
Spain = data.frame(Germany,Day = 1:61 ,Cumulative.Recovered = rec$Spain.Cumulative.Recovered[1:61], Recovered = rec$Spain.Recovered[1:61])

# cure rate
lambda <- function(t, lam0, lam1) {lam0*(1- exp(-lam1*t))}

# mortality rate 
kappa <- function(t, kap0, kap1) {kap0*exp(-kap1*t)}


# transmission model
seir.extended <- function(t,dt,alpha,beta,gamma,delta,X,J,lam0,lam1,kap0,kap1){
  
  if (J == 1){
    S = X[1]
    E = X[2]
    I = X[3]
    Q = X[4]
    R = X[5]
    D = X[6]
    P = X[7]
  } else {
    S = X[,1]
    E = X[,2]
    I = X[,3]
    Q = X[,4]
    R = X[,5]
    D = X[,6]
    P = X[,7]
  }
  
  N = S+E+I+Q+R+D+P
  v.se = I*(beta*dt/N)
  v.sp = alpha*dt
  v.ei = gamma*dt
  v.iq = delta*dt
  v.qr = lambda(t, lam0, lam1)*dt
  v.qd = kappa(t, kap0, kap1)*dt
  

  probs = c(v.iq, v.qr, v.qd)
  
  #rbinom(N, size  = XP[i-1,,1], prob = 1-(1-p)^XP[i-1,,2])
  dN.se = rpois(J, lambda = v.se*S*dt)
  dN.sp = rpois(J, lambda = v.sp*S*dt)
  dN.ei = rpois(J, lambda = v.ei*E*dt)
  dN.iq = rpois(J, lambda = v.iq*I*dt)
  dN.qr = rpois(J, lambda = v.qr*Q*dt)
  dN.qd = rpois(J, lambda = v.qd*Q*dt)
  
  
  if (J == 1){

    X[1] = ifelse(S - dN.se - dN.sp <= 0, 0, S - dN.se - dN.sp)
    X[2] = ifelse(E + dN.se - dN.ei <= 0, 0, E + dN.se - dN.ei)
    X[3] = ifelse(I + dN.ei - dN.iq <= 0, 0, I + dN.ei - dN.iq)
    X[4] = ifelse(Q + dN.iq - dN.qr - dN.qd <= 0, 0, Q + dN.iq - dN.qr - dN.qd) 
    X[5] = R + dN.qr 
    X[6] = D + dN.qd 
    X[7] = P + dN.sp
  } else {
    
    dN <- matrix(0, nrow = J, ncol = 7)
    
    X[,1] = ifelse(S - dN.se - dN.sp <= 0, 0, S - dN.se - dN.sp)
    X[,2] = ifelse(E + dN.se - dN.ei <= 0, 0, E + dN.se - dN.ei)
    X[,3] = ifelse(I + dN.ei - dN.iq <= 0, 0, I + dN.ei - dN.iq)
    X[,4] = ifelse(Q + dN.iq - dN.qr - dN.qd <= 0, 0, Q + dN.iq - dN.qr - dN.qd) 
    X[,5] = R + dN.qr 
    X[,6] = D + dN.qd 
    X[,7] = P + dN.sp 
    
    dN[,1] = dN.se
    dN[,2] = dN.ei
    dN[,3] = dN.iq
    dN[,4] = dN.qr
    dN[,5] = dN.qd
    dN[,6] = dN.sp
    
  }
  return(list(X = X, dN = dN, probs = probs))
}



# simulate y trajectory
true.y <- function(tstar,J,T,S0,E0,I0,Q0,R0,D0,P0,alpha,beta,gamma,delta,lam0,lam1,kap0,kap1){
  x.true <- array(rep(0,7*J*(T+1)),dim=c(T+1,J,7))
  proposed_y = array(0,dim = c(T,1,3))
  
  x.true[1,,1] = S0
  x.true[1,,2] = E0
  x.true[1,,3] = I0
  x.true[1,,4] = Q0
  x.true[1,,5] = R0
  x.true[1,,6] = D0
  x.true[1,,7] = P0
  
  set.seed(2020)
  for (i in 1:T) {
    if (length(alpha) > 1){
      alphax <- ifelse(i <= tstar[1], alpha[1], 
                       ifelse(i <= tstar[2], alpha[2], 
                              ifelse(i <= tstar[3], alpha[3], 
                                     ifelse(i <= tstar[4], alpha[4], alpha[5]))))
    } else {
      alphax <- alpha 
    }
    
      myseir <- seir.extended(i,dt,alphax,beta,gamma,delta,x.true[i,1,],J=1,lam0,lam1,kap0,kap1)
    x.true[i+1,1,] = myseir$X
    proposed_y[i,1,1] = rnbinom(1, size = 7, mu = myseir$probs[1]*x.true[i+1,1,4])
    proposed_y[i,1,2] = rnbinom(1, size = 7, mu = myseir$probs[2]*x.true[i+1,1,5])
    proposed_y[i,1,3] = rnbinom(1, size = 7, mu = myseir$probs[3]*x.true[i+1,1,6])

  }
  return(list(proposed_y = proposed_y, xtrue = x.true))
}






# Sequential Monte Carlo
smc <- function(tstar,T,J,dt,alpha,beta,gamma,delta,S0,E0,I0,Q0,R0,D0,P0,lam0,lam1,kap0,kap1,proposed_y, ytrue){
  w <- matrix(0, T, J)
  wnorm <- matrix(0, T, J)
  X <- array(rep(0,7*J*(T+1)),dim=c(T+1,J,7))
  Xbar <- array(rep(0,7*J*(T+1)),dim=c(T+1,J,7))
  marginal_lik <- vector("numeric")
  logmarginal_lik <- vector("numeric")
  ESS <-vector("numeric")
  DN <- array(rep(0,7*J*(T)),dim=c(T,J,7))
  DNbar <- array(rep(0,7*J*(T)),dim=c(T,J,7))
  rep.no = matrix(0, T, 1)
  
  
  # initialise states 
  X[1,,1] = S0
  X[1,,2] = E0
  X[1,,3] = I0
  X[1,,4] = Q0
  X[1,,5] = R0
  X[1,,6] = D0
  X[1,,7] = P0
  
  Xbar = X
  set.seed(2020)
  for (i in 1:T){
    if (length(alpha) > 1){
      alphax <- ifelse(i <= tstar[1], alpha[1], 
                       ifelse(i <= tstar[2], alpha[2], 
                              ifelse(i <= tstar[3], alpha[3], 
                                     ifelse(i <= tstar[4], alpha[4], alpha[5]))))
    } else {
      alphax <- alpha 
    }
    
    
    model <- seir.extended(i,dt,alphax,beta,gamma,delta,X[i,,],J,lam0,lam1,kap0,kap1)
    X[i+1,,] <- model$X
    DN[i,,] <- model$dN
    
    # Measurement process
    y <- dnbinom(proposed_y[i,1,3], mu = model$probs[3]*X[i+1,,6], size = 7)
    w[i,] <- y
    
    
    for (j in 1:J){
      wnorm[i,j] = ifelse(w[i,j] == 0, 1/J , w[i,j]/sum(w[i,]))
      
      if (wnorm[i,j]  > 1/J){
        wnorm[i,j] = wnorm[i,j]
      } else {
        wnorm[i,j] = 1/J}
         }
    
    ESS[i] =  1/sum(wnorm[i,]^2)
    
    if (i==1){
      marginal_lik[i] <- sum(w[i,])/J
      logmarginal_lik[i] <- log(sum(w[i,])/J)
    } else {
      marginal_lik[i] <- marginal_lik[i-1]*sum(w[i,])/J
      logmarginal_lik[i] <- logmarginal_lik[i-1] + log(sum(w[i,])/J)
    }
    
    #Resampling step using Multinomial Resampling 
    number = as.numeric(rmultinom(1,J,wnorm[i,]))
    number[0] = 0
    for(ii in 1:J){
      if((number[ii]==0)==FALSE){
        for(iii in 1:number[ii]){
          k = sum(number[0:(ii-1)])
          X[i+1,iii+k,]=X[i+1,ii,]
          DN[i, iii+k,]= DN[i, ii,]
        }}}
    Xbar[i+1,,] <- X[i+1,,]
    DNbar[i,,] <- DN[i,,]
    
   
    #ind <- sample(1:J, size = J, replace = TRUE, prob = wnorm[i,])
    #DN[1:(i),,] <- DN[1:(i), ind,]
    #X[1:(i+1),,] <-  X[1:(i+1), ind,]
    #Xbar[i+1,,] <- X[(i+1),,]
    #DNbar[i,,] <- DN[(i),,]

    # Reproduction Number 
    rep.no[i] <- beta*(1/delta)*(1-alphax)^i
  }
  return(list(X = Xbar, dN = DNbar, Weights = wnorm, Likelihood = marginal_lik, Loglikelihood = logmarginal_lik, ESS = ESS, Rep.no = rep.no))
}

#36 until 8 march

# define fixed parameters 
J1 = 1
J2 = 250
T = 210

#alpha <- c(0.0001, 0.025)
theta0 <- c(1.3, 0.2, 0.6, 0.05, 0.9, 0.06, 0.05)
alpha_it <- c(0.0001, 0.0005, 0.001, 0.025, 0.01)
tstar_it <- c(33, 37, 39, 120)

theta_uk <- c(1.3, 0.2, 0.6, 0.05, 0.9, 0.05, 0.055)
alpha_uk <- c(0.0001, 0.0005, 0.001, 0.025, 0.01)
tstar_uk <- c(44, 49, 52, 120)

theta_gm <- c(1.3, 0.2, 0.6, 0.05, 0.9, 0.04, 0.07)
alpha_gm <-c(0.0001, 0.0005, 0.001, 0.025, 0.01)
tstar_gm <- c(40, 42, 50, 120)

theta_fr <- c(1.3, 0.2, 0.6, 0.05, 0.9, 0.04, 0.05)
alpha_fr <- c(0.0001, 0.0005, 0.001, 0.025, 0.01)
tstar_fr <- c(41, 42, 44, 120)

alpha_sp <- c(0.0001, 0.0005, 0.001, 0.025, 0.01)
tstar_sp <- c(37, 41, 42, 120)


set.seed(2020)
# Proposed trajectories 
it_truestates <- true.y(tstar=tstar_it,J1,T,S0=60000000-170,E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,
                        alpha=alpha_it,beta=theta0[1],gamma=theta0[2],delta=theta0[3],
                        lam0=theta0[4],lam1=theta0[5],kap0=theta0[6],kap1=theta0[7]) # Italy
prop_it <- it_truestates$proposed_y


uk_truestates <- true.y(tstar=tstar_uk,J1,T,S0=60660000-170,E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,
                    alpha=alpha_uk,beta=theta_uk[1],gamma=theta_uk[2],delta=theta_uk[3],
                    lam0=theta_uk[4],lam1=theta_uk[5],kap0=theta_uk[6],kap1=theta_uk[7]) # UK
prop_uk <- uk_truestates$proposed_y


gm_truestates <- true.y(tstar=tstar_gm,J1,T,S0=83000000-170,E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,
                        alpha=alpha_gm,beta=theta_gm[1],gamma=theta_gm[2],delta=theta_gm[3],
                        lam0=theta_gm[4],lam1=theta_gm[5],kap0=theta_gm[6],kap1=theta_gm[7])
prop_germ <- gm_truestates$proposed_y


fr_truestates <- true.y(tstar=tstar_fr,J1,T,S0=66990000-170,E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,
                        alpha=alpha_fr,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],
                        lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7])
prop_fr <- fr_truestates$proposed_y


sp_truestates <- true.y(tstar=tstar_sp,J1,T,S0=46940000-170,E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,
                    alpha=alpha_sp,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],
                    lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7])
prop_spain <- sp_truestates$proposed_y




# Italy coronavirus simulation for simulated and true data 
smc_it <- smc(tstar=tstar_it,T,J2,dt,alpha=alpha_it,beta=theta0[1],gamma=theta0[2],delta=theta0[3],S0=60000000-170,
         E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta0[4],lam1=theta0[5],kap0=theta0[6],kap1=theta0[7],prop_it)
smc_it$X[,2,]


truevals <- array(0,dim = c(T=61,1,3))
truevals[,1,1] <- Italy$Cumulative.Confirmed
truevals[,1,2] <- Italy$Cumulative.Recovered
truevals[,1,3] <- Italy$Cumulative.Deaths
ittruevals <- smc(tstar=tstar_it,T=61,J2,dt,alpha=alpha_it,beta=theta0[1],gamma=theta0[2],delta=theta0[3],S0=60000000-170,
                E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta0[4],lam1=theta0[5],kap0=theta0[6],kap1=theta0[7],truevals)
ittruevals$X[,2,]




# UK coronavirus simulation for simulated and true data 
smc_uk <- smc(tstar=tstar_uk,T,J2,dt,alpha=alpha_uk,beta=theta_uk[1],gamma=theta_uk[2],delta=theta_uk[3],S0=60660000-170,
         E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_uk[4],lam1=theta_uk[5],kap0=theta_uk[6],kap1=theta_uk[7],prop_uk)
smc_uk$X[,2,]


truevals2 <- array(0,dim = c(T=61,1,3))
truevals2[,1,1] <- UK$Cumulative.Confirmed
truevals2[,1,2] <- UK$Cumulative.Recovered
truevals2[,1,3] <- UK$Cumulative.Deaths
uktruevals <- smc(tstar=tstar_uk,T=61,J2,dt,alpha=alpha_uk,beta=theta_uk[1],gamma=theta_uk[2],delta=theta_uk[3],S0=60000000-170,
                  E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_uk[4],lam1=theta_uk[5],kap0=theta_uk[6],kap1=theta_uk[7],truevals2)
uktruevals$X[,2,]
                  
                  
# Germany coronavirus simulation for simulated and true data 
smc_germ <- smc(tstar=tstar_gm,T,J2,dt,alpha=alpha_gm,beta=theta_gm[1],gamma=theta_gm[2],delta=theta_gm[3],S0=83000000-170,
         E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_gm[4],lam1=theta_gm[5],kap0=theta_gm[6],kap1=theta_gm[7],prop_germ)
smc_germ$X[,2,]

truevals3 <- array(0,dim = c(T=61,1,3))
truevals3[,1,1] <- Germany$Cumulative.Confirmed
truevals3[,1,2] <- Germany$Cumulative.Recovered
truevals3[,1,3] <- Germany$Cumulative.Deaths
gemtruevals <- smc(tstar=tstar_gm,T=61,J2,dt,alpha=alpha_gm,beta=theta_gm[1],gamma=theta_gm[2],delta=theta_gm[3],S0=60000000-170,
                  E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_gm[4],lam1=theta_gm[5],kap0=theta_gm[6],kap1=theta_gm[7],truevals3)
gemtruevals$X[,2,]


# France coronavirus simulation for simulated and true data 
smc_fr <- smc(tstar=tstar_fr,T,J2,dt,alpha=alpha_fr,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],S0=66990000-170,
         E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7],prop_fr)
smc_fr$X[,2,]

truevals4 <- array(0,dim = c(T=61,1,3))
truevals4[,1,1] <- France$Cumulative.Confirmed
truevals4[,1,2] <- France$Cumulative.Recovered
truevals4[,1,3] <- France$Cumulative.Deaths
frtruevals <- smc(tstar=tstar_fr,T=61,J2,dt,alpha=alpha_fr,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],S0=60000000-170,
                   E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7],truevals4)
frtruevals$X[,2,]



# Spain coronavirus simulation for simulated and true data 
smc_sp <- smc(tstar=tstar_sp,T,J2,dt,alpha=alpha_sp,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],S0=46940000-170,
         E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7],prop_spain)
smc_sp$X[,2,]

truevals5 <- array(0,dim = c(T=61,1,3))
truevals5[,1,1] <- Spain$Cumulative.Confirmed
truevals5[,1,2] <- Spain$Cumulative.Recovered
truevals5[,1,3] <- Spain$Cumulative.Deaths
sptruevals <- smc(tstar=tstar_sp,T=61,J2,dt,alpha=alpha_sp,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],S0=60000000-170,
                  E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7],truevals5)
sptruevals$X[,2,]









###### PLOTS OF SMC 
### ITALY
italymeans <- data.frame(time = 1:(T+1), mean1 = rowMeans(smc_it$X[,,1]), mean2 = rowMeans(smc_it$X[,,2]),
                         mean3 = rowMeans(smc_it$X[,,3]), mean4 = rowMeans(smc_it$X[,,4]),
                         mean5 = rowMeans(smc_it$X[,,5]), mean6 = rowMeans(smc_it$X[,,6]),
                         mean7 = rowMeans(smc_it$X[,,7]))

italymeans2 <- data.frame(time = 1:61, m1 = rowMeans(ittruevals$X[1:61,,1]), m2 = rowMeans(ittruevals$X[1:61,,2]),
                          m3 = rowMeans(ittruevals$X[1:61,,3]), m4 = rowMeans(ittruevals$X[1:61,,4]), m5 = rowMeans(ittruevals$X[1:61,,5]), 
                          m6 = rowMeans(ittruevals$X[1:61,,6]), m7 = rowMeans(ittruevals$X[1:61,,7]),
                          truecases = Italy$Cumulative.Confirmed, trueds = Italy$Cumulative.Deaths,
                          truerec = Italy$Cumulative.Recovered)


it1 <- ggplot() +
  geom_line(data =italymeans,  aes(x = time, y=mean1), col = 1)  + 
  geom_point(data =italymeans2,  aes(x = time, y=m1), col = 5)  + 
  ggtitle("Italy Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

it2 <- ggplot() +
  geom_line(data = italymeans, aes(x = time, y=mean2), col = 1)  +
  geom_point(data =italymeans2,  aes(x = time, y=m2), col = 5) +
  ggtitle("Italy Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

it3 <- ggplot() +
  geom_line(data = italymeans, aes(x = time, y=mean3), col = 1)  +
  geom_point(data =italymeans2,  aes(x = time, y=m3), col = 5) +
  ggtitle("Italy Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

it4 <- ggplot() +
  geom_line(data = italymeans, aes(x = time, y=mean4), col = 1)  +
  geom_point(data =italymeans2,  aes(x = time, y=m4), col = 5) +
  geom_point(data =italymeans2,  aes(x = time, y= truecases), col = 2)  +
  ggtitle("Italy Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

it5 <- ggplot() +
  geom_line(data = italymeans, aes(x = time, y=mean5), col = 1)  +
  geom_point(data =italymeans2,  aes(x = time, y=m5), col = 5) +
  geom_point(data =italymeans2,  aes(x = time, y= truerec), col = 2)  +
  ggtitle("Italy Recovered Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

it6 <- ggplot() +
  geom_line(data = italymeans, aes(x = time, y=mean6), col = 1)  +
  geom_point(data =italymeans2,  aes(x = time, y=m6), col = 5) +
  geom_point(data =italymeans2,  aes(x = time, y= trueds), col = 2)  +
  ggtitle("Italy Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


it7 <- ggplot() +
  geom_line(data = italymeans, aes(x = time, y=mean7), col = 1)  +
  geom_point(data =italymeans2,  aes(x = time, y=m7), col = 5) +
  ggtitle("Italy Insusceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

figure <- ggarrange(it1, it2, it3, it4, it5, it6, it7,
                    ncol = 2, nrow = 4)





##################### UK
ukmeans <- data.frame(time = 1:(T+1), mean1 = rowMeans(smc_uk$X[,,1]), mean2 = rowMeans(smc_uk$X[,,2]),
                         mean3 = rowMeans(smc_uk$X[,,3]), mean4 = rowMeans(smc_uk$X[,,4]),
                         mean5 = rowMeans(smc_uk$X[,,5]), mean6 = rowMeans(smc_uk$X[,,6]),
                         mean7 = rowMeans(smc_uk$X[,,7]))

ukmeans2 <- data.frame(time = 1:61, m1 = rowMeans(uktruevals$X[1:61,,1]), m2 = rowMeans(uktruevals$X[1:61,,2]),
                          m3 = rowMeans(uktruevals$X[1:61,,3]), m4 = rowMeans(uktruevals$X[1:61,,4]), m5 = rowMeans(uktruevals$X[1:61,,5]), 
                          m6 = rowMeans(uktruevals$X[1:61,,6]), m7 = rowMeans(uktruevals$X[1:61,,7]),
                          truecases = UK$Cumulative.Confirmed, trueds = UK$Cumulative.Deaths,
                          truerec = UK$Cumulative.Recovered)


uk1 <- ggplot() +
  geom_line(data =ukmeans,  aes(x = time, y=mean1), col = 1)  + 
  geom_point(data =ukmeans2,  aes(x = time, y=m1), col = 5)  + 
  ggtitle("UK Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

uk2 <- ggplot() +
  geom_line(data =ukmeans,  aes(x = time, y=mean2), col = 1)  + 
  geom_point(data =ukmeans2,  aes(x = time, y=m2), col = 5)  + 
  ggtitle("UK Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

uk3 <- ggplot() +
  geom_line(data =ukmeans,  aes(x = time, y=mean3), col = 1)  + 
  geom_point(data =ukmeans2,  aes(x = time, y=m3), col = 5)  + 
  ggtitle("UK Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

uk4 <- ggplot() +
  geom_line(data =ukmeans,  aes(x = time, y=mean4), col = 1)  + 
  geom_point(data =ukmeans2,  aes(x = time, y=m4), col = 5)  + 
  ggtitle("UK Quarantined Population") +
  geom_point(data =ukmeans2,  aes(x = time, y= truecases), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

uk5 <- ggplot() +
  geom_line(data =ukmeans,  aes(x = time, y=mean5), col = 1)  + 
  geom_point(data =ukmeans2,  aes(x = time, y=m5), col = 5)  + 
  geom_point(data =ukmeans2,  aes(x = time, y= truerec), col = 2)  + 
  ggtitle("UK Recovered Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

uk6 <- ggplot() +
  geom_line(data =ukmeans,  aes(x = time, y=mean6), col = 1)  + 
  geom_point(data =ukmeans2,  aes(x = time, y=m6), col = 5)  + 
  geom_point(data =ukmeans2,  aes(x = time, y= trueds), col = 2)  + 
  ggtitle("UK Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

uk7 <- ggplot() +
  geom_line(data =ukmeans,  aes(x = time, y=mean7), col = 1)  + 
  geom_point(data =ukmeans2,  aes(x = time, y=m7), col = 5)  + 
  ggtitle("UK Insusceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


figure2 <- ggarrange(uk1, uk2, uk3, uk4, uk5, uk6, uk7,
                    ncol = 2, nrow = 4)




############# GERMANY
grmeans <- data.frame(time = 1:(T+1), mean1 = rowMeans(smc_germ$X[,,1]), mean2 = rowMeans(smc_germ$X[,,2]),
                      mean3 = rowMeans(smc_germ$X[,,3]), mean4 = rowMeans(smc_germ$X[,,4]),
                      mean5 = rowMeans(smc_germ$X[,,5]), mean6 = rowMeans(smc_germ$X[,,6]),
                      mean7 = rowMeans(smc_germ$X[,,7]))

grmeans2 <- data.frame(time = 1:61, m1 = rowMeans(gemtruevals$X[1:61,,1]), m2 = rowMeans(gemtruevals$X[1:61,,2]),
                       m3 = rowMeans(gemtruevals$X[1:61,,3]), m4 = rowMeans(gemtruevals$X[1:61,,4]), m5 = rowMeans(gemtruevals$X[1:61,,5]), 
                       m6 = rowMeans(gemtruevals$X[1:61,,6]), m7 = rowMeans(gemtruevals$X[1:61,,7]),
                       truecases = Germany$Cumulative.Confirmed, trueds = Germany$Cumulative.Deaths,
                       truerec = Germany$Cumulative.Recovered)


gr1 <- ggplot() +
  geom_line(data =grmeans,  aes(x = time, y=mean1, colour="M1"), col = 1)  + 
  geom_point(data =grmeans2,  aes(x = time, y=m1, colour="M2"), col = 5)  + 
  ggtitle("Germany Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") +  theme(legend.position = "bottom")

gr2 <-ggplot() +
  geom_line(data =grmeans,  aes(x = time, y=mean2), col = 1)  + 
  geom_point(data =grmeans2,  aes(x = time, y=m2), col = 5)  + 
  ggtitle("Germany Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

gr3 <-ggplot() +
  geom_line(data =grmeans,  aes(x = time, y=mean3), col = 1)  + 
  geom_point(data =grmeans2,  aes(x = time, y=m3), col = 5)  + 
  ggtitle("Germany Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

gr4 <-ggplot() +
  geom_line(data =grmeans,  aes(x = time, y=mean4), col = 1)  + 
  geom_point(data =grmeans2,  aes(x = time, y=m4), col = 5)  + 
  geom_point(data =grmeans2,  aes(x = time, y= truecases), col = 2)  + 
  ggtitle("Germany Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

gr5 <-ggplot() +
  geom_line(data =grmeans,  aes(x = time, y=mean5), col = 1)  + 
  geom_point(data =grmeans2,  aes(x = time, y=m5), col = 5)  + 
  ggtitle("Germany Recovered Population") +
  geom_point(data =grmeans2,  aes(x = time, y= truerec), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

gr6 <-ggplot() +
  geom_line(data =grmeans,  aes(x = time, y=mean6), col = 1)  + 
  geom_point(data =grmeans2,  aes(x = time, y=m6), col = 5)  + 
  geom_point(data =grmeans2,  aes(x = time, y= trueds), col = 2)  + 
  ggtitle("Germany Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


gr7 <-ggplot() +
  geom_line(data =grmeans,  aes(x = time, y=mean7), col = 1)  + 
  geom_point(data =grmeans2,  aes(x = time, y=m7), col = 5)  + 
  ggtitle("Germany Insusceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

figure3 <- ggarrange(gr1, gr2, gr3, gr4, gr5, gr6, gr7,
                     ncol = 2, nrow = 4)





############# FRANCE
frmeans <- data.frame(time = 1:(T+1), mean1 = rowMeans(smc_fr$X[,,1]), mean2 = rowMeans(smc_fr$X[,,2]),
                      mean3 = rowMeans(smc_fr$X[,,3]), mean4 = rowMeans(smc_fr$X[,,4]),
                      mean5 = rowMeans(smc_fr$X[,,5]), mean6 = rowMeans(smc_fr$X[,,6]),
                      mean7 = rowMeans(smc_fr$X[,,7]))
frmeans2 <- data.frame(time = 1:61, m1 = rowMeans(frtruevals$X[1:61,,1]), m2 = rowMeans(frtruevals$X[1:61,,2]),
                       m3 = rowMeans(frtruevals$X[1:61,,3]), m4 = rowMeans(frtruevals$X[1:61,,4]), m5 = rowMeans(frtruevals$X[1:61,,5]), 
                       m6 = rowMeans(frtruevals$X[1:61,,6]), m7 = rowMeans(frtruevals$X[1:61,,7]),
                       truecases = France$Cumulative.Confirmed, trueds = France$Cumulative.Deaths,
                       truerec = France$Cumulative.Recovered)

fr1 <- ggplot() +
  geom_line(data =frmeans,  aes(x = time, y=mean1), col = 1)  + 
  geom_point(data =frmeans2,  aes(x = time, y=m1), col = 5)  + 
  ggtitle("France Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

fr2 <- ggplot() +
  geom_line(data =frmeans,  aes(x = time, y=mean2), col = 1)  + 
  geom_point(data =frmeans2,  aes(x = time, y=m2), col = 5)  + 
  ggtitle("France Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

fr3 <- ggplot() +
  geom_line(data =frmeans,  aes(x = time, y=mean3), col = 1)  + 
  geom_point(data =frmeans2,  aes(x = time, y=m3), col = 5)  + 
  ggtitle("France Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

fr4 <- ggplot() +
  geom_line(data =frmeans,  aes(x = time, y=mean4), col = 1)  + 
  geom_point(data =frmeans2,  aes(x = time, y=m4), col = 5)  + 
  geom_point(data =frmeans2,  aes(x = time, y=truecases), col = 2)  + 
  ggtitle("France Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

fr5 <- ggplot() +
  geom_line(data =frmeans,  aes(x = time, y=mean5), col = 1)  + 
  geom_point(data =frmeans2,  aes(x = time, y=m5), col = 5)  + 
  ggtitle("France Recovered Population") +
  geom_point(data =frmeans2,  aes(x = time, y=truerec), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

fr6 <- ggplot() +
  geom_line(data =frmeans,  aes(x = time, y=mean6), col = 1)  + 
  geom_point(data =frmeans2,  aes(x = time, y=m6), col = 5)  + 
  ggtitle("France Deaths") +
  geom_point(data =frmeans2,  aes(x = time, y=trueds), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

fr7 <- ggplot() +
  geom_line(data =frmeans,  aes(x = time, y=mean7), col = 1)  + 
  geom_point(data =frmeans2,  aes(x = time, y=m7), col = 5)  + 
  ggtitle("France Insusceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

figure4 <- ggarrange(fr1, fr2, fr3, fr4, fr5, fr6, fr7,
                     ncol = 2, nrow = 4)





############# FRANCE
spmeans <- data.frame(time = 1:(T+1), mean1 = rowMeans(smc_sp$X[,,1]), mean2 = rowMeans(smc_sp$X[,,2]),
                      mean3 = rowMeans(smc_sp$X[,,3]), mean4 = rowMeans(smc_sp$X[,,4]),
                      mean5 = rowMeans(smc_sp$X[,,5]), mean6 = rowMeans(smc_sp$X[,,6]),
                      mean7 = rowMeans(smc_sp$X[,,7]))

spmeans2 <- data.frame(time = 1:61, m1 = rowMeans(sptruevals$X[1:61,,1]), m2 = rowMeans(sptruevals$X[1:61,,2]),
                       m3 = rowMeans(sptruevals$X[1:61,,3]), m4 = rowMeans(sptruevals$X[1:61,,4]), m5 = rowMeans(sptruevals$X[1:61,,5]), 
                       m6 = rowMeans(sptruevals$X[1:61,,6]), m7 = rowMeans(sptruevals$X[1:61,,7]),
                       truecases = Spain$Cumulative.Confirmed, trueds = Spain$Cumulative.Deaths,
                       truerec = Spain$Cumulative.Recovered)

sp1 <- ggplot() +
  geom_line(data =spmeans,  aes(x = time, y=mean1), col = 1)  + 
  geom_point(data =spmeans2,  aes(x = time, y=m1), col = 5)  + 
  ggtitle("Spain Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

sp2 <- ggplot() +
  geom_line(data =spmeans,  aes(x = time, y=mean2), col = 1)  + 
  geom_point(data =spmeans2,  aes(x = time, y=m2), col = 5)  + 
  ggtitle("Spain Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

sp3 <- ggplot() +
  geom_line(data =spmeans,  aes(x = time, y=mean3), col = 1)  + 
  geom_point(data =spmeans2,  aes(x = time, y=m3), col = 5)  +
  ggtitle("Spain Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

sp4 <- ggplot() +
  geom_line(data =spmeans,  aes(x = time, y=mean4), col = 1)  + 
  geom_point(data =spmeans2,  aes(x = time, y=m4), col = 5)  + 
  geom_point(data =spmeans2,  aes(x = time, y=truecases), col = 2)  + 
  ggtitle("Spain Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

sp5 <- ggplot() +
  geom_line(data =spmeans,  aes(x = time, y=mean5), col = 1)  + 
  geom_point(data =spmeans2,  aes(x = time, y=m5), col = 5)  + 
  geom_point(data =spmeans2,  aes(x = time, y=truerec), col = 2)  + 
  ggtitle("Spain Recovered Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


sp6 <- ggplot() +
  geom_line(data =spmeans,  aes(x = time, y=mean6), col = 1)  + 
  geom_point(data =spmeans2,  aes(x = time, y=m6), col = 5)  + 
  geom_point(data =spmeans2,  aes(x = time, y=trueds), col = 2)  + 
  ggtitle("Spain Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

sp7 <- ggplot() +
  geom_line(data =spmeans,  aes(x = time, y=mean7), col = 1)  + 
  geom_point(data =spmeans2,  aes(x = time, y=m7), col = 5)  + 
  ggtitle("Spain Insusceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


figure5 <- ggarrange(sp1, sp2, sp3, sp4, sp5, sp6, sp7,
                     ncol = 2, nrow = 4)

#### Reproduction Number 
rpnumber <- data.frame(time = 1:(T), repit = smc_it$Rep.no, repuk = smc_uk$Rep.no,
                       repgm = smc_germ$Rep.no, repfr = smc_fr$Rep.no, repsp = smc_sp$Rep.no)
# Italy
itrepfig <- ggplot() +
  geom_line(data = rpnumber,  aes(x = time, y=repit), col = 1)  + 
  ggtitle("Italy Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 

# UK
ukrepfig <- ggplot() +
  geom_line(data = rpnumber,  aes(x = time, y=repuk), col = 1)  + 
  ggtitle("UK Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 


# Germany
germanyrepfig <- ggplot() +
  geom_line(data = rpnumber,  aes(x = time, y=repgm), col = 1)  + 
  ggtitle("Germany Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 


# France
francerepfig <- ggplot() +
  geom_line(data = rpnumber,  aes(x = time, y=repfr), col = 1)  + 
  ggtitle("France Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 



# Spain
spainrepfig <- ggplot() +
  geom_line(data = rpnumber,  aes(x = time, y=repsp), col = 1)  + 
  ggtitle("Spain Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN")


reps <- ggarrange(itrepfig, ukrepfig, germanyrepfig, 
                  francerepfig, spainrepfig, 
                  ncol = 2, nrow = 4)


### WITHOUT ITNERVENTIONS 
# Proposed trajectories 
set.seed(2020)
woutit_truestates <- true.y(tstar=tstar_it,J1,T,S0=60000000-170,E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,
                        alpha=0,beta=theta0[1],gamma=theta0[2],delta=theta0[3],
                        lam0=theta0[4],lam1=theta0[5],kap0=theta0[6],kap1=theta0[7]) # Italy
woutprop_it <- it_truestates$proposed_y


woutuk_truestates <- true.y(tstar=tstar_uk,J1,T,S0=60660000-170,E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,
                        alpha=0,beta=theta_uk[1],gamma=theta_uk[2],delta=theta_uk[3],
                        lam0=theta_uk[4],lam1=theta_uk[5],kap0=theta_uk[6],kap1=theta_uk[7]) # UK
woutprop_uk <- uk_truestates$proposed_y


woutgm_truestates <- true.y(tstar=tstar_gm,J1,T,S0=83000000-170,E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,
                        alpha=0,beta=theta_gm[1],gamma=theta_gm[2],delta=theta_gm[3],
                        lam0=theta_gm[4],lam1=theta_gm[5],kap0=theta_gm[6],kap1=theta_gm[7])
woutprop_germ <- gm_truestates$proposed_y


woutfr_truestates <- true.y(tstar=tstar_fr,J1,T,S0=66990000-170,E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,
                        alpha=0,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],
                        lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7])
woutprop_fr <- fr_truestates$proposed_y


woutsp_truestates <- true.y(tstar=tstar_sp,J1,T,S0=46940000-170,E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,
                        alpha=0,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],
                        lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7])
woutprop_spain <- sp_truestates$proposed_y




# Italy coronavirus simulation for simulated and true data 
woutsmc_it <- smc(tstar=tstar_it,T,J2,dt,alpha=0,beta=theta0[1],gamma=theta0[2],delta=theta0[3],S0=60000000-170,
              E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta0[4],lam1=theta0[5],kap0=theta0[6],kap1=theta0[7],prop_it)
woutsmc_it$X[,2,]


wtruevals <- array(0,dim = c(T=61,1,3))
wtruevals[,1,1] <- Italy$Cumulative.Confirmed
wtruevals[,1,2] <- Italy$Cumulative.Recovered
wtruevals[,1,3] <- Italy$Cumulative.Deaths
woutittruevals <- smc(tstar=tstar_it,T=61,J2,dt,alpha=0,beta=theta0[1],gamma=theta0[2],delta=theta0[3],S0=60000000-170,
                  E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta0[4],lam1=theta0[5],kap0=theta0[6],kap1=theta0[7],truevals)
woutittruevals$X[,2,]




# UK coronavirus simulation for simulated and true data 
woutsmc_uk <- smc(tstar=tstar_uk,T,J2,dt,alpha=0,beta=theta_uk[1],gamma=theta_uk[2],delta=theta_uk[3],S0=60660000-170,
              E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_uk[4],lam1=theta_uk[5],kap0=theta_uk[6],kap1=theta_uk[7],prop_uk)
woutsmc_uk$X[,2,]


wtruevals2 <- array(0,dim = c(T=61,1,3))
wtruevals2[,1,1] <- UK$Cumulative.Confirmed
wtruevals2[,1,2] <- UK$Cumulative.Recovered
wtruevals2[,1,3] <- UK$Cumulative.Deaths
woutuktruevals <- smc(tstar=tstar_uk,T=61,J2,dt,alpha=0,beta=theta_uk[1],gamma=theta_uk[2],delta=theta_uk[3],S0=60000000-170,
                  E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_uk[4],lam1=theta_uk[5],kap0=theta_uk[6],kap1=theta_uk[7],truevals2)
woutuktruevals$X[,2,]


# Germany coronavirus simulation for simulated and true data 
woutsmc_germ <- smc(tstar=tstar_gm,T,J2,dt,alpha=0,beta=theta_gm[1],gamma=theta_gm[2],delta=theta_gm[3],S0=83000000-170,
                E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_gm[4],lam1=theta_gm[5],kap0=theta_gm[6],kap1=theta_gm[7],prop_germ)
woutsmc_germ$X[,2,]

wtruevals3 <- array(0,dim = c(T=61,1,3))
wtruevals3[,1,1] <- Germany$Cumulative.Confirmed
wtruevals3[,1,2] <- Germany$Cumulative.Recovered
wtruevals3[,1,3] <- Germany$Cumulative.Deaths
woutgemtruevals <- smc(tstar=tstar_gm,T=61,J2,dt,alpha=0,beta=theta_gm[1],gamma=theta_gm[2],delta=theta_gm[3],S0=60000000-170,
                   E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_gm[4],lam1=theta_gm[5],kap0=theta_gm[6],kap1=theta_gm[7],truevals3)
woutgemtruevals$X[,2,]


# France coronavirus simulation for simulated and true data 
woutsmc_fr <- smc(tstar=tstar_fr,T,J2,dt,alpha=0,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],S0=66990000-170,
              E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7],prop_fr)
woutsmc_fr$X[,2,]

wtruevals4 <- array(0,dim = c(T=61,1,3))
wtruevals4[,1,1] <- France$Cumulative.Confirmed
wtruevals4[,1,2] <- France$Cumulative.Recovered
wtruevals4[,1,3] <- France$Cumulative.Deaths
woutfrtruevals <- smc(tstar=tstar_fr,T=61,J2,dt,alpha=0,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],S0=60000000-170,
                  E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7],truevals4)
woutfrtruevals$X[,2,]



# Spain coronavirus simulation for simulated and true data 
woutsmc_sp <- smc(tstar=tstar_sp,T,J2,dt,alpha=0,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],S0=46940000-170,
              E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7],prop_spain)
woutsmc_sp$X[,2,]

wtruevals5 <- array(0,dim = c(T=61,1,3))
wtruevals5[,1,1] <- Spain$Cumulative.Confirmed
wtruevals5[,1,2] <- Spain$Cumulative.Recovered
wtruevals5[,1,3] <- Spain$Cumulative.Deaths
woutsptruevals <- smc(tstar=tstar_sp,T=61,J2,dt,alpha=0,beta=theta_fr[1],gamma=theta_fr[2],delta=theta_fr[3],S0=60000000-170,
                  E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,lam0=theta_fr[4],lam1=theta_fr[5],kap0=theta_fr[6],kap1=theta_fr[7],truevals5)
woutsptruevals$X[,2,]



###### PLOTS WITHOUT INTERVENTIONS
###### PLOTS OF SMC 
### ITALY
woutitalymeans <- data.frame(time = 1:(T+1), mean1 = rowMeans(woutsmc_it$X[,,1]), mean2 = rowMeans(woutsmc_it$X[,,2]),
                         mean3 = rowMeans(woutsmc_it$X[,,3]), mean4 = rowMeans(woutsmc_it$X[,,4]),
                         mean5 = rowMeans(woutsmc_it$X[,,5]), mean6 = rowMeans(woutsmc_it$X[,,6]),
                         mean7 = rowMeans(woutsmc_it$X[,,7]))

woutitalymeans2 <- data.frame(time = 1:61, m1 = rowMeans(woutittruevals$X[1:61,,1]), m2 = rowMeans(woutittruevals$X[1:61,,2]),
                          m3 = rowMeans(woutittruevals$X[1:61,,3]), m4 = rowMeans(woutittruevals$X[1:61,,4]), m5 = rowMeans(woutittruevals$X[1:61,,5]), 
                          m6 = rowMeans(woutittruevals$X[1:61,,6]), m7 = rowMeans(woutittruevals$X[1:61,,7]),
                          truecases = Italy$Cumulative.Confirmed, trueds = Italy$Cumulative.Deaths,
                          truerec = Italy$Cumulative.Recovered)


wit1 <- ggplot(data =woutitalymeans, aes(x = time)) +
  geom_line(data = woutitalymeans, aes(x = time, y=mean1), col = 1)  +
  geom_point(data =woutitalymeans2,  aes(x = time, y=m1), col = 5) +
  ggtitle("Italy Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wit2 <- ggplot() +
  geom_line(data = woutitalymeans, aes(x = time, y=mean2), col = 1)  +
  geom_point(data =woutitalymeans2,  aes(x = time, y=m2), col = 5) +
  ggtitle("Italy Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wit3 <- ggplot() +
  geom_line(data = woutitalymeans, aes(x = time, y=mean3), col = 1)  +
  geom_point(data =woutitalymeans2,  aes(x = time, y=m3), col = 5) +
  ggtitle("Italy Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wit4 <- ggplot() +
  geom_line(data = woutitalymeans, aes(x = time, y=mean4), col = 1)  +
  geom_point(data =woutitalymeans2,  aes(x = time, y=m4), col = 5) +
  geom_point(data =woutitalymeans2,  aes(x = time, y= truecases), col = 2)  +
  ggtitle("Italy Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wit5 <- ggplot() +
  geom_line(data = woutitalymeans, aes(x = time, y=mean5), col = 1)  +
  geom_point(data =woutitalymeans2,  aes(x = time, y=m5), col = 5) +
  geom_point(data =woutitalymeans2,  aes(x = time, y= truerec), col = 2)  +
  ggtitle("Italy Recovered Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wit6 <- ggplot() +
  geom_line(data = woutitalymeans, aes(x = time, y=mean6), col = 1)  +
  geom_point(data =woutitalymeans2,  aes(x = time, y=m6), col = 5) +
  geom_point(data =woutitalymeans2,  aes(x = time, y= trueds), col = 2)  +
  ggtitle("Italy Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 



wfigure <- ggarrange(wit1, wit2, wit3, wit4, wit5, wit6,
                    ncol = 2, nrow = 4)





##################### UK
woutukmeans <- data.frame(time = 1:(T+1), mean1 = rowMeans(woutsmc_uk$X[,,1]), mean2 = rowMeans(woutsmc_uk$X[,,2]),
                      mean3 = rowMeans(woutsmc_uk$X[,,3]), mean4 = rowMeans(woutsmc_uk$X[,,4]),
                      mean5 = rowMeans(woutsmc_uk$X[,,5]), mean6 = rowMeans(woutsmc_uk$X[,,6]),
                      mean7 = rowMeans(woutsmc_uk$X[,,7]))

woutukmeans2 <- data.frame(time = 1:61, m1 = rowMeans(woutuktruevals$X[1:61,,1]), m2 = rowMeans(woutuktruevals$X[1:61,,2]),
                       m3 = rowMeans(woutuktruevals$X[1:61,,3]), m4 = rowMeans(woutuktruevals$X[1:61,,4]), m5 = rowMeans(woutuktruevals$X[1:61,,5]), 
                       m6 = rowMeans(woutuktruevals$X[1:61,,6]), m7 = rowMeans(woutuktruevals$X[1:61,,7]),
                       truecases = UK$Cumulative.Confirmed, trueds = UK$Cumulative.Deaths,
                       truerec = UK$Cumulative.Recovered)


wuk1 <- ggplot() +
  geom_line(data =woutukmeans,  aes(x = time, y=mean1), col = 1)  + 
  geom_point(data =woutukmeans2,  aes(x = time, y=m1), col = 5)  + 
  ggtitle("UK Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wuk2 <- ggplot() +
  geom_line(data =woutukmeans,  aes(x = time, y=mean2), col = 1)  + 
  geom_point(data =woutukmeans2,  aes(x = time, y=m2), col = 5)  + 
  ggtitle("UK Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wuk3 <- ggplot() +
  geom_line(data =woutukmeans,  aes(x = time, y=mean3), col = 1)  + 
  geom_point(data =woutukmeans2,  aes(x = time, y=m3), col = 5)  + 
  ggtitle("UK Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wuk4 <- ggplot() +
  geom_line(data =woutukmeans,  aes(x = time, y=mean4), col = 1)  + 
  geom_point(data =woutukmeans2,  aes(x = time, y=m4), col = 5)  + 
  ggtitle("UK Quarantined Population") +
  geom_point(data =woutukmeans2,  aes(x = time, y= truecases), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wuk5 <- ggplot() +
  geom_line(data =woutukmeans,  aes(x = time, y=mean5), col = 1)  + 
  geom_point(data =woutukmeans2,  aes(x = time, y=m5), col = 5)  + 
  geom_point(data =woutukmeans2,  aes(x = time, y= truerec), col = 2)  + 
  ggtitle("UK Recovered Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wuk6 <- ggplot() +
  geom_line(data =woutukmeans,  aes(x = time, y=mean6), col = 1)  + 
  geom_point(data =woutukmeans2,  aes(x = time, y=m6), col = 5)  + 
  geom_point(data =woutukmeans2,  aes(x = time, y= trueds), col = 2)  + 
  ggtitle("UK Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


wfigure2 <- ggarrange(wuk1, wuk2, wuk3, wuk4, wuk5, wuk6, 
                     ncol = 2, nrow = 4)




############# GERMANY
woutgrmeans <- data.frame(time = 1:(T+1), mean1 = rowMeans(woutsmc_germ$X[,,1]), mean2 = rowMeans(woutsmc_germ$X[,,2]),
                      mean3 = rowMeans(woutsmc_germ$X[,,3]), mean4 = rowMeans(woutsmc_germ$X[,,4]),
                      mean5 = rowMeans(woutsmc_germ$X[,,5]), mean6 = rowMeans(woutsmc_germ$X[,,6]),
                      mean7 = rowMeans(woutsmc_germ$X[,,7]))

woutgrmeans2 <- data.frame(time = 1:61, m1 = rowMeans(woutgemtruevals$X[1:61,,1]), m2 = rowMeans(woutgemtruevals$X[1:61,,2]),
                       m3 = rowMeans(woutgemtruevals$X[1:61,,3]), m4 = rowMeans(woutgemtruevals$X[1:61,,4]), m5 = rowMeans(woutgemtruevals$X[1:61,,5]), 
                       m6 = rowMeans(woutgemtruevals$X[1:61,,6]), m7 = rowMeans(woutgemtruevals$X[1:61,,7]),
                       truecases = Germany$Cumulative.Confirmed, trueds = Germany$Cumulative.Deaths,
                       truerec = Germany$Cumulative.Recovered)


wgr1 <- ggplot() +
  geom_line(data =woutgrmeans,  aes(x = time, y=mean1, colour="M1"), col = 1)  + 
  geom_point(data =woutgrmeans2,  aes(x = time, y=m1, colour="M2"), col = 5)  + 
  ggtitle("Germany Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") +  theme(legend.position = "bottom")

wgr2 <-ggplot() +
  geom_line(data =woutgrmeans,  aes(x = time, y=mean2), col = 1)  + 
  geom_point(data =woutgrmeans2,  aes(x = time, y=m2), col = 5)  + 
  ggtitle("Germany Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wgr3 <-ggplot() +
  geom_line(data =woutgrmeans,  aes(x = time, y=mean3), col = 1)  + 
  geom_point(data =woutgrmeans2,  aes(x = time, y=m3), col = 5)  + 
  ggtitle("Germany Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wgr4 <-ggplot() +
  geom_line(data =woutgrmeans,  aes(x = time, y=mean4), col = 1)  + 
  geom_point(data =woutgrmeans2,  aes(x = time, y=m4), col = 5)  + 
  geom_point(data =woutgrmeans2,  aes(x = time, y= truecases), col = 2)  + 
  ggtitle("Germany Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wgr5 <-ggplot() +
  geom_line(data =woutgrmeans,  aes(x = time, y=mean5), col = 1)  + 
  geom_point(data =woutgrmeans2,  aes(x = time, y=m5), col = 5)  + 
  ggtitle("Germany Recovered Population") +
  geom_point(data =woutgrmeans2,  aes(x = time, y= truerec), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wgr6 <-ggplot() +
  geom_line(data =woutgrmeans,  aes(x = time, y=mean6), col = 1)  + 
  geom_point(data =woutgrmeans2,  aes(x = time, y=m6), col = 5)  + 
  geom_point(data =woutgrmeans2,  aes(x = time, y= trueds), col = 2)  + 
  ggtitle("Germany Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 



wfigure3 <- ggarrange(wgr1, wgr2, wgr3, wgr4, wgr5, wgr6,
                     ncol = 2, nrow = 4)





############# FRANCE
woutfrmeans <- data.frame(time = 1:(T+1), mean1 = rowMeans(woutsmc_fr$X[,,1]), mean2 = rowMeans(woutsmc_fr$X[,,2]),
                      mean3 = rowMeans(woutsmc_fr$X[,,3]), mean4 = rowMeans(woutsmc_fr$X[,,4]),
                      mean5 = rowMeans(woutsmc_fr$X[,,5]), mean6 = rowMeans(woutsmc_fr$X[,,6]),
                      mean7 = rowMeans(woutsmc_fr$X[,,7]))
woutfrmeans2 <- data.frame(time = 1:61, m1 = rowMeans(woutfrtruevals$X[1:61,,1]), m2 = rowMeans(woutfrtruevals$X[1:61,,2]),
                       m3 = rowMeans(woutfrtruevals$X[1:61,,3]), m4 = rowMeans(woutfrtruevals$X[1:61,,4]), m5 = rowMeans(woutfrtruevals$X[1:61,,5]), 
                       m6 = rowMeans(woutfrtruevals$X[1:61,,6]), m7 = rowMeans(woutfrtruevals$X[1:61,,7]),
                       truecases = France$Cumulative.Confirmed, trueds = France$Cumulative.Deaths,
                       truerec = France$Cumulative.Recovered)

wfr1 <- ggplot() +
  geom_line(data =woutfrmeans,  aes(x = time, y=mean1), col = 1)  + 
  geom_point(data =woutfrmeans2,  aes(x = time, y=m1), col = 5)  + 
  ggtitle("France Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wfr2 <- ggplot() +
  geom_line(data =woutfrmeans,  aes(x = time, y=mean2), col = 1)  + 
  geom_point(data =woutfrmeans2,  aes(x = time, y=m2), col = 5)  + 
  ggtitle("France Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

wfr3 <- ggplot() +
  geom_line(data =woutfrmeans,  aes(x = time, y=mean3), col = 1)  + 
  geom_point(data =woutfrmeans2,  aes(x = time, y=m3), col = 5)  + 
  ggtitle("France Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

wfr4 <- ggplot() +
  geom_line(data =woutfrmeans,  aes(x = time, y=mean4), col = 1)  + 
  geom_point(data =woutfrmeans2,  aes(x = time, y=m4), col = 5)  + 
  geom_point(data =woutfrmeans2,  aes(x = time, y=truecases), col = 2)  + 
  ggtitle("France Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

wfr5 <- ggplot() +
  geom_line(data =woutfrmeans,  aes(x = time, y=mean5), col = 1)  + 
  geom_point(data =woutfrmeans2,  aes(x = time, y=m5), col = 5)  + 
  ggtitle("France Recovered Population") +
  geom_point(data =woutfrmeans2,  aes(x = time, y=truerec), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

wfr6 <- ggplot() +
  geom_line(data =woutfrmeans,  aes(x = time, y=mean6), col = 1)  + 
  geom_point(data =woutfrmeans2,  aes(x = time, y=m6), col = 5)  + 
  ggtitle("France Deaths") +
  geom_point(data =woutfrmeans2,  aes(x = time, y=trueds), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")


wfigure4 <- ggarrange(wfr1, wfr2, wfr3, wfr4, wfr5, wfr6,
                     ncol = 2, nrow = 4)





############# SPAIN
woutspmeans <- data.frame(time = 1:(T+1), mean1 = rowMeans(woutsmc_sp$X[,,1]), mean2 = rowMeans(woutsmc_sp$X[,,2]),
                      mean3 = rowMeans(woutsmc_sp$X[,,3]), mean4 = rowMeans(woutsmc_sp$X[,,4]),
                      mean5 = rowMeans(woutsmc_sp$X[,,5]), mean6 = rowMeans(woutsmc_sp$X[,,6]),
                      mean7 = rowMeans(woutsmc_sp$X[,,7]))

woutspmeans2 <- data.frame(time = 1:61, m1 = rowMeans(woutsptruevals$X[1:61,,1]), m2 = rowMeans(woutsptruevals$X[1:61,,2]),
                       m3 = rowMeans(woutsptruevals$X[1:61,,3]), m4 = rowMeans(woutsptruevals$X[1:61,,4]), m5 = rowMeans(woutsptruevals$X[1:61,,5]), 
                       m6 = rowMeans(woutsptruevals$X[1:61,,6]), m7 = rowMeans(woutsptruevals$X[1:61,,7]),
                       truecases = Spain$Cumulative.Confirmed, trueds = Spain$Cumulative.Deaths,
                       truerec = Spain$Cumulative.Recovered)

wsp1 <- ggplot() +
  geom_line(data =woutspmeans,  aes(x = time, y=mean1), col = 1)  + 
  geom_point(data =woutspmeans2,  aes(x = time, y=m1), col = 5)  + 
  ggtitle("Spain Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wsp2 <- ggplot() +
  geom_line(data =woutspmeans,  aes(x = time, y=mean2), col = 1)  + 
  geom_point(data =woutspmeans2,  aes(x = time, y=m2), col = 5)  + 
  ggtitle("Spain Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wsp3 <- ggplot() +
  geom_line(data =woutspmeans,  aes(x = time, y=mean3), col = 1)  + 
  geom_point(data =woutspmeans2,  aes(x = time, y=m3), col = 5)  +
  ggtitle("Spain Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wsp4 <- ggplot() +
  geom_line(data =woutspmeans,  aes(x = time, y=mean4), col = 1)  + 
  geom_point(data =woutspmeans2,  aes(x = time, y=m4), col = 5)  + 
  geom_point(data =woutspmeans2,  aes(x = time, y=truecases), col = 2)  + 
  ggtitle("Spain Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wsp5 <- ggplot() +
  geom_line(data =woutspmeans,  aes(x = time, y=mean5), col = 1)  + 
  geom_point(data =woutspmeans2,  aes(x = time, y=m5), col = 5)  + 
  geom_point(data =woutspmeans2,  aes(x = time, y=truerec), col = 2)  + 
  ggtitle("Spain Recovered Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


wsp6 <- ggplot() +
  geom_line(data =woutspmeans,  aes(x = time, y=mean6), col = 1)  + 
  geom_point(data =woutspmeans2,  aes(x = time, y=m6), col = 5)  + 
  geom_point(data =woutspmeans2,  aes(x = time, y=trueds), col = 2)  + 
  ggtitle("Spain Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

wfigure5 <- ggarrange(wsp1, wsp2, wsp3, wsp4, wsp5, wsp6, 
                     ncol = 2, nrow = 4)

#### Reproduction Number 
woutrpnumber <- data.frame(time = 1:(T), repit = woutsmc_it$Rep.no, repuk = woutsmc_uk$Rep.no,
                       repgm = woutsmc_germ$Rep.no, repfr = woutsmc_fr$Rep.no, repsp = woutsmc_sp$Rep.no)
# Italy
woutitrepfig <- ggplot() +
  geom_line(data = woutrpnumber,  aes(x = time, y=repit), col = 1)  + 
  ggtitle("Italy Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 

# UK
woutukrepfig <- ggplot() +
  geom_line(data = woutrpnumber,  aes(x = time, y=repuk), col = 1)  + 
  ggtitle("UK Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 


# Germany
woutgermanyrepfig <- ggplot() +
  geom_line(data = woutrpnumber,  aes(x = time, y=repgm), col = 1)  + 
  ggtitle("Germany Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 


# France
woutfrancerepfig <- ggplot() +
  geom_line(data = woutrpnumber,  aes(x = time, y=repfr), col = 1)  + 
  ggtitle("France Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 



# Spain
woutspainrepfig <- ggplot() +
  geom_line(data = woutrpnumber,  aes(x = time, y=repsp), col = 1)  + 
  ggtitle("Spain Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN")


woutreps <- ggarrange(woutitrepfig, woutukrepfig, woutgermanyrepfig, 
                      woutfrancerepfig, woutspainrepfig, 
                  ncol = 2, nrow = 4)




























################### implement PMCMC 
# n: iterations 
# nburn: burn in times for convergence 
# J: number of particles
# thetas: vector of parameters
# lamda: step size 
# ap: acceptance ratio for Metropolis Hastings 
n = 1000
nburn = 100
lamda = 0.1/n 



PMCMC <- function(tstar,T,J,n,nburn,alpha,theta0,lamda,S0,E0,I0,Q0,R0,D0,P0,proposed_y){
  ap <- vector("numeric")
  thetas <- matrix(0, n+nburn, 7)
  thetas[1,] <- theta0
  Loglikelihood <- vector("numeric")
  X <- array(rep(0,7*(n+nburn)*(T+1)),dim=c(T+1,n+nburn,7))
  Xbar <- array(rep(0,7*(n+nburn)*(T+1)),dim=c(T+1,n+nburn,7))
  dN <- array(rep(0,7*(n+nburn)*(T)),dim=c(T,n+nburn,7))
  Rep.no <- array(rep(0,(n+nburn)*(T)),dim=c(T,n+nburn))
  
  SMC <- smc(tstar,T,J,dt,alpha = alpha ,beta=theta0[1],gamma=theta0[2],delta=theta0[3],
             S0,E0,I0,Q0,R0,D0,P0,lam0=theta0[4],lam1=theta0[5],kap0=theta0[6],kap1=theta0[7],proposed_y)
  res <- sample(J, 1, prob = SMC$Weights[T,]) # resampling 
  X[,1,] <- SMC$X[,res,]
  Xbar[,1,] <- X[,1,]
  
  Loglikelihood[1] <- SMC$Loglikelihood[T] # Marginal Likelihood 
  Rep.no[,1] <- SMC$Rep.no # Reproduction number 
  
  
  for (i in 2:(n+nburn)){
      # sample proposal s.t. theta' = theta(n-1) + lamda*N(0,1)
      thetanew <- thetas[i-1,] + rnorm(7, mean = 0, sd = lamda)
      
      for (jj in 1:7){
        thetanew[jj] <- ifelse(thetanew[jj] < 0, 0.0001, thetanew[jj])
      }
      
      tdash = thetanew
      
      SMC <- smc(tstar,T,J,dt,alpha = alpha ,beta=tdash[1],gamma=tdash[2],delta=tdash[3],S0,E0,I0,Q0,R0,D0,P0,
                 lam0=tdash[4],lam1=tdash[5],kap0=tdash[6],kap1=tdash[7],proposed_y)
      
      Loglikelihood[i] <- SMC$Loglikelihood[T]
      
      xmarginal = Loglikelihood[i]
      ymarginal = Loglikelihood[i-1]
      
      num <- xmarginal*dunif(tdash[1], 1, 1.8)*dunif(tdash[2], 0.1, 0.7)*
        dunif(tdash[3], 0.1,0.7)*dunif(tdash[4], 0.02,0.05)*dunif(tdash[5], 0.02,0.9)*
        dunif(tdash[6], 0.017,0.06)*dunif(tdash[7], 0.017,0.06)*dnorm(thetas[i-1,1], mean =  tdash[1], sd = lamda)*
        dnorm(thetas[i-1,2], mean =  tdash[2], sd = lamda)*dnorm(thetas[i-1,3], mean =  tdash[3], sd = lamda)*
        dnorm(thetas[i-1,4], mean =  tdash[4], sd = lamda)*dnorm(thetas[i-1,5], mean =  tdash[5], sd = lamda)*
        dnorm(thetas[i-1,6], mean =  tdash[6], sd = lamda)*dnorm(thetas[i-1,7], mean =  tdash[7], sd = lamda)
      
      den <- ymarginal*dunif(thetas[i-1,1], 1, 1.8)*dunif(thetas[i-1,2], 0.1,0.7)*
        dunif(thetas[i-1,3], 0.1,0.7)*dunif(thetas[i-1,4], 0.02,0.05)*dunif(thetas[i-1,5], 0.02,0.9)*
        dunif(thetas[i-1,6], 0.017,0.06)*dunif(thetas[i-1,7], 0.017,0.06)*dnorm(tdash[1], mean =  thetas[i-1,1], sd = lamda)*
        dnorm(tdash[2], mean =  thetas[i-1,2], sd = lamda)*dnorm(tdash[3], mean =  thetas[i-1,3], sd = lamda)*
        dnorm(tdash[4], mean =  thetas[i-1,4], sd = lamda)*dnorm(tdash[5], mean =  thetas[i-1,5], sd = lamda)*
        dnorm(tdash[6], mean =  thetas[i-1,6], sd = lamda)*dnorm(tdash[7], mean =  thetas[i-1,7], sd = lamda)
      
      
      prob <- ifelse(is.na(num/den), 1, num/den)
      a <- min(1, prob)
      ap <- c(ap, a)
      
      U <- runif(1)
      if (U <= a){
        #thetas[i,] <- c(tdash[1], tdash[c(-1)]^2)
        thetas[i,] <- thetanew
        res <- sample(J, 1, prob = SMC$Weights[T,])
        X[,i,] <- SMC$X[,res,]
        dN[,i,] <- SMC$dN[,res,]
        Xbar[,i,] <- X[,i,]
        Rep.no[,i] <- SMC$Rep.no # Reproduction number 
        
      } else { 
        thetas[i,] <- thetas[i-1,]
        X[,i,] <- X[,i-1,]
        Xbar[,i,] <- Xbar[,i-1,]
        dN[,i,] <- matrix(0, T, 7)
        Rep.no[,i] <- Rep.no[,i-1]
      }
  }
  
  # discard first nburn iterations 
  Rep.no <- Rep.no[,c((nburn+1):(n+nburn))]
  X <- X[,c((nburn+1):(n+nburn)),]
  Xbar <- Xbar[,c((nburn+1):(n+nburn)),]
  thetas <- thetas[c((nburn+1):(n+nburn)),]
  ap <- ap[c((nburn+1):(n+nburn))]
  Loglikelihood <- Loglikelihood[c((nburn+1):(n+nburn))]
  list("X"= X, "dN" = dN, "thetas" = thetas, "ap" =  ap, "Loglikelihood" = Loglikelihood, "Rep.no" = Rep.no)
}



# Initialise parameters 
T = 210
theta0 <- c(1.3, 0.2, 0.6, 0.05, 0.9, 0.06, 0.05)


# Italy states and parameter inference 
italysim <-PMCMC(tstar=tstar_it,T,J=250,n,nburn,alpha_it,theta0,lamda,S0=60000000-170,E0=166,
                 I0=4,Q0=0,R0=0,D0=0,P0=0,prop_it)
italysim$thetas
italysim$X[,100,]


# UK states and parameter inference 
uksim <- PMCMC(tstar=tstar_uk,T,J=250,n,nburn,alpha = alpha_uk,theta_uk,lamda,S0=60660000-170,E0=166,
              I0=4,Q0=0,R0=0,D0=0,P0=0,prop_uk)
uksim$thetas
uksim$X[,100,]


# Germany states and parameter inference
germanysim <- PMCMC(tstar=tstar_gm,T,J=250,n,nburn,alpha_gm,theta_gm,lamda,S0=83000000-170,
                   E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,proposed_germ)
germanysim$thetas
germanysim$X[,100,]


# France states and parameter inference
francesim <- PMCMC(tstar=tstar_fr,T,J=250,n,nburn,alpha=alpha_fr,theta_fr,lamda,S0=66990000-170,
                   E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,prop_fr)
               francesim$thetas
francesim$X[,100,]


# Spain states and parameter inference
spainsim <- PMCMC(tstar=tstar_sp,T,J=250,n,nburn,alpha=alpha_sp,theta_fr,lamda,S0=46940000-170,
                  E0=166,I0=4,Q0=0,R0=0,D0=0,P0=0,prop_spain)
spainsim$thetas
spainsim$X[,100,]


############# PLOTS FOR PMCMC
########### ITALY PARAMETERS 
itthetas <- data.frame(ind = 1:1000, beta = italysim$thetas[,1], gamma = italysim$thetas[,2],
                       delta = italysim$thetas[,3], lamda0 = italysim$thetas[,4], lamda1 = italysim$thetas[,5],
                       kappa0 = italysim$thetas[,6], kappa1 = italysim$thetas[,7])

i1 <- ggplot(itthetas, aes(x = beta)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", beta))) +
  xlab(expression(paste(beta))) + ylab("Density") + theme(legend.position = "bottom")

i2 <- ggplot(itthetas, aes(x = gamma)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", gamma))) +
  xlab(expression(paste(gamma))) + ylab("Density") + theme(legend.position = "bottom")

i3 <- ggplot(itthetas, aes(x = delta)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", delta))) +
  xlab(expression(paste(delta))) + ylab("Density") + theme(legend.position = "bottom")

i4 <- ggplot(itthetas, aes(x = lamda0)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\lambda_0$")) +
  xlab(TeX("$\\lambda_0$")) + ylab("Density") + theme(legend.position = "bottom")

i5 <- ggplot(itthetas, aes(x = lamda1)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\lambda_1$")) +
  xlab(TeX("$\\lambda_1$")) + ylab("Density") + theme(legend.position = "bottom")

i6 <- ggplot(itthetas, aes(x = kappa0)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\kappa_0$")) +
  xlab(TeX("$\\kappa_0$")) + ylab("Density") + theme(legend.position = "bottom")

i7 <- ggplot(itthetas, aes(x = kappa1)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\kappa_1$")) +
  xlab(TeX("$\\kappa_1$")) + ylab("Density") + theme(legend.position = "bottom")

itfig <- ggarrange(i1, i2, i3, i4, i5, i6, i7,
                   ncol = 2, nrow = 4)


########### UK PARAMETERS 
ukthetas <- data.frame(ind = 1:1000, beta = uksim$thetas[,1], gamma = uksim$thetas[,2],
                       delta = uksim$thetas[,3], lamda0 = uksim$thetas[,4], lamda1 = uksim$thetas[,5],
                       kappa0 = uksim$thetas[,6], kappa1 = uksim$thetas[,7])

u1 <- ggplot(ukthetas, aes(x = beta)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", beta))) +
  xlab(expression(paste(beta))) + ylab("Density") + theme(legend.position = "bottom")

u2 <- ggplot(ukthetas, aes(x = gamma)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", gamma))) +
  xlab(expression(paste(gamma))) + ylab("Density") + theme(legend.position = "bottom")

u3 <- ggplot(ukthetas, aes(x = delta)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", delta))) +
  xlab(expression(paste(delta))) + ylab("Density") + theme(legend.position = "bottom")

u4 <- ggplot(ukthetas, aes(x = lamda0)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\lambda_0$")) +
  xlab(TeX("$\\lambda_0$")) + ylab("Density") + theme(legend.position = "bottom")

u5 <- ggplot(ukthetas, aes(x = lamda1)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\lambda_1$")) +
  xlab(TeX("$\\lambda_1$")) + ylab("Density") + theme(legend.position = "bottom")

u6 <- ggplot(ukthetas, aes(x = kappa0)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\kappa_0$")) +
  xlab(TeX("$\\kappa_0$")) + ylab("Density") + theme(legend.position = "bottom")

u7 <- ggplot(ukthetas, aes(x = kappa1)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\kappa_1$")) +
  xlab(TeX("$\\kappa_1$")) + ylab("Density") + theme(legend.position = "bottom")

ukfig <- ggarrange(u1, u2, u3, u4, u5, u6, u7,
                   ncol = 2, nrow = 4)








########### GERMANY PARAMETERS 
gemthetas <- data.frame(ind = 1:1000, beta = germanysim$thetas[,1], gamma = germanysim$thetas[,2],
                       delta = germanysim$thetas[,3], lamda0 = germanysim$thetas[,4], lamda1 = germanysim$thetas[,5],
                       kappa0 = germanysim$thetas[,6], kappa1 = germanysim$thetas[,7])
g1 <- ggplot(gemthetas, aes(x = beta)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", beta))) +
  xlab(expression(paste(beta))) + ylab("Density") + theme(legend.position = "bottom")

g2 <- ggplot(gemthetas, aes(x = gamma)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", gamma))) +
  xlab(expression(paste(gamma))) + ylab("Density") + theme(legend.position = "bottom")

g3 <- ggplot(gemthetas, aes(x = delta)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", delta))) +
  xlab(expression(paste(delta))) + ylab("Density") + theme(legend.position = "bottom")

g4 <- ggplot(gemthetas, aes(x = lamda0)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\lambda_0$")) +
  xlab(TeX("$\\lambda_0$")) + ylab("Density") + theme(legend.position = "bottom")

g5 <- ggplot(gemthetas, aes(x = lamda1)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\lambda_1$")) +
  xlab(TeX("$\\lambda_1$")) + ylab("Density") + theme(legend.position = "bottom")

g6 <- ggplot(gemthetas, aes(x = kappa0)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\kappa_0$")) +
  xlab(TeX("$\\kappa_0$")) + ylab("Density") + theme(legend.position = "bottom")

g7 <- ggplot(gemthetas, aes(x = kappa1)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\kappa_1$")) +
  xlab(TeX("$\\kappa_1$")) + ylab("Density") + theme(legend.position = "bottom")

grfig <- ggarrange(g1, g2, g3, g4, g5, g6, g7,
                   ncol = 2, nrow = 4)






########### FRANCE PARAMETERS 
frthetas <- data.frame(ind = 1:1000, beta = francesim$thetas[,1], gamma = francesim$thetas[,2],
                        delta = francesim$thetas[,3], lamda0 = francesim$thetas[,4], lamda1 = francesim$thetas[,5],
                        kappa0 = francesim$thetas[,6], kappa1 = francesim$thetas[,7])
f1 <- ggplot(frthetas, aes(x = beta)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", beta))) +
  xlab(expression(paste(beta))) + ylab("Density") + theme(legend.position = "bottom")

f2 <- ggplot(frthetas, aes(x = gamma)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", gamma))) +
  xlab(expression(paste(gamma))) + ylab("Density") + theme(legend.position = "bottom")

f3 <- ggplot(frthetas, aes(x = delta)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", delta))) +
  xlab(expression(paste(delta))) + ylab("Density") + theme(legend.position = "bottom")

f4 <- ggplot(frthetas, aes(x = lamda0)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\lambda_0$")) +
  xlab(TeX("$\\lambda_0$")) + ylab("Density") + theme(legend.position = "bottom")

f5 <- ggplot(frthetas, aes(x = lamda1)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\lambda_1$")) +
  xlab(TeX("$\\lambda_1$")) + ylab("Density") + theme(legend.position = "bottom")

f6 <- ggplot(frthetas, aes(x = kappa0)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\kappa_0$")) +
  xlab(TeX("$\\kappa_0$")) + ylab("Density") + theme(legend.position = "bottom")

f7 <- ggplot(frthetas, aes(x = kappa1)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\kappa_1$")) +
  xlab(TeX("$\\kappa_1$")) + ylab("Density") + theme(legend.position = "bottom")

frfig <- ggarrange(f1, f2, f3, f4, f5, f6, f7,
                   ncol = 2, nrow = 4)





########### SPAIN PARAMETERS 
spthetas <- data.frame(ind = 1:1000, beta = spainsim$thetas[,1], gamma = spainsim$thetas[,2],
                       delta = spainsim$thetas[,3], lamda0 = spainsim$thetas[,4], lamda1 = spainsim$thetas[,5],
                       kappa0 = spainsim$thetas[,6], kappa1 = spainsim$thetas[,7])
s1 <- ggplot(spthetas, aes(x = beta)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", beta))) +
  xlab(expression(paste(beta))) + ylab("Density") + theme(legend.position = "bottom")

s2 <- ggplot(spthetas, aes(x = gamma)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", gamma))) +
  xlab(expression(paste(gamma))) + ylab("Density") + theme(legend.position = "bottom")

s3 <- ggplot(spthetas, aes(x = delta)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste("Histogram of ", delta))) +
  xlab(expression(paste(delta))) + ylab("Density") + theme(legend.position = "bottom")

s4 <- ggplot(spthetas, aes(x = lamda0)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\lambda_0$")) +
  xlab(TeX("$\\lambda_0$")) + ylab("Density") + theme(legend.position = "bottom")

s5 <- ggplot(spthetas, aes(x = lamda1)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\lambda_1$")) +
  xlab(TeX("$\\lambda_1$")) + ylab("Density") + theme(legend.position = "bottom")

s6 <- ggplot(spthetas, aes(x = kappa0)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\kappa_0$")) +
  xlab(TeX("$\\kappa_0$")) + ylab("Density") + theme(legend.position = "bottom")

s7 <- ggplot(spthetas, aes(x = kappa1)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="lightblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(TeX(" Histogram of $\\kappa_1$")) +
  xlab(TeX("$\\kappa_1$")) + ylab("Density") + theme(legend.position = "bottom")

spfig <- ggarrange(s1, s2, s3, s4, s5, s6, s7,
                   ncol = 2, nrow = 4)


#### Reproduction Numbers 

# Italy
itrep2 <- data.frame(time = 1:T, rp = rowMeans(italysim$Rep.no))
itrepfig2 <- ggplot() +
  geom_line(data = itrep2,  aes(x = time, y=rp), col = 1)  + 
  ggtitle("Italy Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 

# UK
ukrep2 <- data.frame(time = 1:T, rp = rowMeans(uksim$Rep.no))
ukrepfig2 <- ggplot() +
  geom_line(data = ukrep2,  aes(x = time, y=rp), col = 1)  + 
  ggtitle("UK Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 


# Germany
germanyrep2 <- data.frame(time = 1:T, rp = rowMeans(germanysim$Rep.no))
germanyrepfig2 <- ggplot() +
  geom_line(data = germanyrep2,  aes(x = time, y=rp), col = 1)  + 
  ggtitle("Germany Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 


# France
francerep2 <- data.frame(time = 1:T, rp = rowMeans(francesim$Rep.no))
francerepfig2 <- ggplot() +
  geom_line(data = francerep2,  aes(x = time, y=rp), col = 1)  + 
  ggtitle("France Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN") 



# Spain
spainrep2 <- data.frame(time = 1:T, rp = rowMeans(spainsim$Rep.no))
spainrepfig2 <- ggplot() +
  geom_line(data = spainrep2,  aes(x = time, y=rp), col = 1)  + 
  ggtitle("Spain Basic Reproduction Number") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("BRN")


reps2 <- ggarrange(itrepfig2, ukrepfig2, germanyrepfig2, 
                  francerepfig2, spainrepfig2, 
                   ncol = 2, nrow = 4)







###### PLOTS OF PMCMC PARTICLES
### ITALY
italypmmc <- data.frame(time = 1:(T+1), mean1 = rowMeans(italysim$X[,,1]), mean2 = rowMeans(italysim$X[,,2]),
                         mean3 = rowMeans(italysim$X[,,3]), mean4 = rowMeans(italysim$X[,,4]),
                         mean5 = rowMeans(italysim$X[,,5]), mean6 = rowMeans(italysim$X[,,6]),
                         mean7 = rowMeans(italysim$X[,,7]))

itpmcmc1 <- ggplot() +
  geom_line(data =italypmmc,  aes(x = time, y=mean1), col = 1)  + 
  ggtitle("Italy Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

itpmcmc2 <- ggplot() +
  geom_line(data = italypmmc, aes(x = time, y=mean2), col = 1)  +
  ggtitle("Italy Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

itpmcmc3 <- ggplot() +
  geom_line(data = italypmmc, aes(x = time, y=mean3), col = 1)  +
  ggtitle("Italy Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

itpmcmc4 <- ggplot() +
  geom_line(data = italypmmc, aes(x = time, y=mean4), col = 1)  +
  geom_point(data =italymeans2,  aes(x = time, y= truecases), col = 2)  +
  ggtitle("Italy Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

itpmcmc5 <- ggplot() +
  geom_line(data = italypmmc, aes(x = time, y=mean5), col = 1)  +
  geom_point(data =italymeans2,  aes(x = time, y= truerec), col = 2)  +
  ggtitle("Italy Recovered Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

itpmcmc6 <- ggplot() +
  geom_line(data = italypmmc, aes(x = time, y=mean6), col = 1)  +
  geom_point(data =italymeans2,  aes(x = time, y= trueds), col = 2)  +
  ggtitle("Italy Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


itpmcmc7 <- ggplot() +
  geom_line(data = italypmmc, aes(x = time, y=mean7), col = 1)  +
  ggtitle("Italy Insusceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

pmcmcit <- ggarrange(itpmcmc1, itpmcmc2, itpmcmc3, itpmcmc4, itpmcmc5, itpmcmc6, itpmcmc7,
                    ncol = 2, nrow = 4)





##################### UK
ukpmcmc <- data.frame(time = 1:(T+1), mean1 = rowMeans(uksim$X[,,1]), mean2 = rowMeans(uksim$X[,,2]),
                      mean3 = rowMeans(uksim$X[,,3]), mean4 = rowMeans(uksim$X[,,4]),
                      mean5 = rowMeans(uksim$X[,,5]), mean6 = rowMeans(uksim$X[,,6]),
                      mean7 = rowMeans(uksim$X[,,7]))

ukpmcmc1 <- ggplot() +
  geom_line(data =ukpmcmc,  aes(x = time, y=mean1), col = 1)  + 
  ggtitle("UK Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

ukpmcmc2 <- ggplot() +
  geom_line(data =ukpmcmc,  aes(x = time, y=mean2), col = 1)  + 
  ggtitle("UK Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

ukpmcmc3 <- ggplot() +
  geom_line(data =ukpmcmc,  aes(x = time, y=mean3), col = 1)  + 
  ggtitle("UK Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

ukpmcmc4 <- ggplot() +
  geom_line(data =ukpmcmc,  aes(x = time, y=mean4), col = 1)  + 
  ggtitle("UK Quarantined Population") +
  geom_point(data =ukmeans2,  aes(x = time, y= truecases), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

ukpmcmc5 <- ggplot() +
  geom_line(data =ukpmcmc,  aes(x = time, y=mean5), col = 1)  + 
  geom_point(data =ukmeans2,  aes(x = time, y= truerec), col = 2)  + 
  ggtitle("UK Recovered Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

ukpmcmc6 <- ggplot() +
  geom_line(data =ukpmcmc,  aes(x = time, y=mean6), col = 1)  + 
  geom_point(data =ukmeans2,  aes(x = time, y= trueds), col = 2)  + 
  ggtitle("UK Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

ukpmcmc7 <- ggplot() +
  geom_line(data =ukpmcmc,  aes(x = time, y=mean7), col = 1)  + 
  ggtitle("UK Insusceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


figukpmcmc <- ggarrange(ukpmcmc1, ukpmcmc2, ukpmcmc3, ukpmcmc4, ukpmcmc5, ukpmcmc6, ukpmcmc7,
                     ncol = 2, nrow = 4)




############# GERMANY
grpmcmc <- data.frame(time = 1:(T+1), mean1 = rowMeans(germanysim$X[,,1]), mean2 = rowMeans(germanysim$X[,,2]),
                      mean3 = rowMeans(germanysim$X[,,3]), mean4 = rowMeans(germanysim$X[,,4]),
                      mean5 = rowMeans(germanysim$X[,,5]), mean6 = rowMeans(germanysim$X[,,6]),
                      mean7 = rowMeans(germanysim$X[,,7]))


grpmcmc1 <- ggplot() +
  geom_line(data =grpmcmc,  aes(x = time, y=mean1, colour="M1"), col = 1)  + 
  ggtitle("Germany Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") +  theme(legend.position = "bottom")

grpmcmc2 <-ggplot() +
  geom_line(data =grpmcmc,  aes(x = time, y=mean2), col = 1)  + 
  ggtitle("Germany Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

grpmcmc3 <-ggplot() +
  geom_line(data =grpmcmc,  aes(x = time, y=mean3), col = 1)  + 
  ggtitle("Germany Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

grpmcmc4 <-ggplot() +
  geom_line(data =grpmcmc,  aes(x = time, y=mean4), col = 1)  + 
  geom_point(data =grmeans2,  aes(x = time, y= truecases), col = 2)  + 
  ggtitle("Germany Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

grpmcmc5 <-ggplot() +
  geom_line(data =grpmcmc,  aes(x = time, y=mean5), col = 1)  + 
  ggtitle("Germany Recovered Population") +
  geom_point(data =grmeans2,  aes(x = time, y= truerec), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

grpmcmc6 <-ggplot() +
  geom_line(data =grpmcmc,  aes(x = time, y=mean6), col = 1)  + 
  geom_point(data =grmeans2,  aes(x = time, y= trueds), col = 2)  + 
  ggtitle("Germany Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


grpmcmc7 <-ggplot() +
  geom_line(data =grpmcmc,  aes(x = time, y=mean7), col = 1)  + 
  ggtitle("Germany Insusceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

figgrpmcmc <- ggarrange(grpmcmc1, grpmcmc2, grpmcmc3, grpmcmc4, grpmcmc5, grpmcmc6, grpmcmc7,
                     ncol = 2, nrow = 4)





############# FRANCE
frpmcmc <- data.frame(time = 1:(T+1), mean1 = rowMeans(francesim$X[,,1]), mean2 = rowMeans(francesim$X[,,2]),
                      mean3 = rowMeans(francesim$X[,,3]), mean4 = rowMeans(francesim$X[,,4]),
                      mean5 = rowMeans(francesim$X[,,5]), mean6 = rowMeans(francesim$X[,,6]),
                      mean7 = rowMeans(francesim$X[,,7]))
frpmcmc1  <- ggplot() +
  geom_line(data =frpmcmc,  aes(x = time, y=mean1), col = 1)  + 
  ggtitle("France Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

frpmcmc2 <- ggplot() +
  geom_line(data =frpmcmc,  aes(x = time, y=mean2), col = 1)  + 
  ggtitle("France Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

frpmcmc3 <- ggplot() +
  geom_line(data =frpmcmc,  aes(x = time, y=mean3), col = 1)  + 
  ggtitle("France Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

frpmcmc4 <- ggplot() +
  geom_line(data =frpmcmc,  aes(x = time, y=mean4), col = 1)  + 
  geom_point(data =frmeans2,  aes(x = time, y=truecases), col = 2)  + 
  ggtitle("France Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

frpmcmc5 <- ggplot() +
  geom_line(data =frpmcmc,  aes(x = time, y=mean5), col = 1)  + 
  ggtitle("France Recovered Population") +
  geom_point(data =frmeans2,  aes(x = time, y=truerec), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

frpmcmc6 <- ggplot() +
  geom_line(data =frpmcmc,  aes(x = time, y=mean6), col = 1)  + 
  ggtitle("France Deaths") +
  geom_point(data =frmeans2,  aes(x = time, y=trueds), col = 2)  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

frpmcmc7 <- ggplot() +
  geom_line(data =frpmcmc,  aes(x = time, y=mean7), col = 1)  + 
  ggtitle("France Insusceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)")

figfrpmcmc <- ggarrange(frpmcmc1, frpmcmc2, frpmcmc3, frpmcmc4, frpmcmc5, frpmcmc6, frpmcmc7,
                     ncol = 2, nrow = 4)





############# FRANCE
sppmcmc <- data.frame(time = 1:(T+1), mean1 = rowMeans(spainsim$X[,,1]), mean2 = rowMeans(spainsim$X[,,2]),
                      mean3 = rowMeans(spainsim$X[,,3]), mean4 = rowMeans(spainsim$X[,,4]),
                      mean5 = rowMeans(spainsim$X[,,5]), mean6 = rowMeans(spainsim$X[,,6]),
                      mean7 = rowMeans(spainsim$X[,,7]))

sppmcmc1 <- ggplot() +
  geom_line(data =sppmcmc,  aes(x = time, y=mean1), col = 1)  + 
  ggtitle("Spain Susceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

sppmcmc2 <- ggplot() +
  geom_line(data =sppmcmc,  aes(x = time, y=mean2), col = 1)  + 
  ggtitle("Spain Exposed Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

sppmcmc3 <- ggplot() +
  geom_line(data =sppmcmc,  aes(x = time, y=mean3), col = 1)  + 
  ggtitle("Spain Infected Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

sppmcmc4 <- ggplot() +
  geom_line(data =sppmcmc,  aes(x = time, y=mean4), col = 1)  + 
  geom_point(data =spmeans2,  aes(x = time, y=truecases), col = 2)  + 
  ggtitle("Spain Quarantined Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

sppmcmc5 <- ggplot() +
  geom_line(data =sppmcmc,  aes(x = time, y=mean5), col = 1)  + 
  geom_point(data =spmeans2,  aes(x = time, y=truerec), col = 2)  + 
  ggtitle("Spain Recovered Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


sppmcmc6 <- ggplot() +
  geom_line(data =sppmcmc,  aes(x = time, y=mean6), col = 1)  + 
  geom_point(data =spmeans2,  aes(x = time, y=trueds), col = 2)  + 
  ggtitle("Spain Deaths") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 

sppmcmc7 <- ggplot() +
  geom_line(data =sppmcmc,  aes(x = time, y=mean7), col = 1)  + 
  ggtitle("Spain Insusceptible Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("E(X)") 


figsppmcmc <- ggarrange(sppmcmc1, sppmcmc2, sppmcmc3, sppmcmc4, sppmcmc5, sppmcmc6, sppmcmc7,
                     ncol = 2, nrow = 4)

##### PLOTS OF TRUE TRAJECTORIES 
truetr_it <- ggplot(data = Italy,aes(x = Day)) +
  geom_line(aes(y=Cumulative.Confirmed, colour = "Confirmed Cases"), size = 1)  + 
  geom_line(aes(y=Cumulative.Recovered, colour = "Recovered Cases"), size = 1)  + 
  geom_line(aes(y=Cumulative.Deaths, colour = "Deaths"), size = 1)  + 
  ggtitle("Italy True Observations") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("X") 



truetr_uk <- ggplot(data = UK,aes(x = Day)) +
  geom_line(aes(y=Cumulative.Confirmed, colour = "Confirmed Cases"), size = 1)  + 
  geom_line(aes(y=Recovered, colour = "Recovered Cases"), size = 1)  + 
  geom_line(aes(y=Cumulative.Deaths, colour = "Deaths"), size = 1)  + 
  ggtitle("UK True Observations") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("X") 


truetr_gm <- ggplot(data = Germany, aes(x = Day)) +
  geom_line(aes(y=Cumulative.Confirmed, colour = "Confirmed Cases"), size = 1)  + 
  geom_line(aes(y=Cumulative.Recovered, colour = "Recovered Cases"), size = 1)  + 
  geom_line(aes(y=Cumulative.Deaths, colour = "Deaths"), size = 1)  + 
  ggtitle("Germany True Observations") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("X") 

truetr_fr <- ggplot(data = France, aes(x = Day)) +
  geom_line(aes(y=Cumulative.Confirmed, colour = "Confirmed Cases"), size = 1)  + 
  geom_line(aes(y=Cumulative.Recovered, colour = "Recovered Cases"), size = 1)  + 
  geom_line(aes(y=Cumulative.Deaths, colour = "Deaths"), size = 1)  + 
  ggtitle("France True Observations") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("X") 


truetr_sp <- ggplot(data = Spain, aes(x = Day)) +
  geom_line(aes(y=Cumulative.Confirmed, colour = "Confirmed Cases"), size = 1)  + 
  geom_line(aes(y=Cumulative.Recovered, colour = "Recovered Cases"), size = 1)  + 
  geom_line(aes(y=Cumulative.Deaths, colour = "Deaths"), size = 1)  + 
  ggtitle("Spain True Observations") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) +
  xlab("time") + ylab("X") 

trajectories <- ggarrange(truetr_it, truetr_uk, truetr_gm, truetr_fr, truetr_sp, 
                        ncol = 2, nrow = 3)

