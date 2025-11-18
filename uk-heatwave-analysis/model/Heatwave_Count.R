##################### Data #####################################
library(mvtnorm)
all_simulators_data <- read.csv("all_simulators_data.csv")
observations <- read.csv("observation.csv")

simulators <- subset(all_simulators_data, Year >= 1981 & Year <= 2000 & Region_ID == 1)
obs <- subset(observations, Year >= 1981 & Year <= 2000 & Region_ID == 1)
######## Def ###########
OneWay_Heatwave<- function(){
  for (i in 1:num_GCMs) {   
    for (t in 1:20) {      
      y[i,t] ~ dnorm(theta[i], tau[i])  
    }
    
    # theta_i
    theta[i] <- theta_0 + delta[i]
    tau[i] ~ dgamma(0.001, 0.001) 
    sigma2[i] <- 1/tau[i]
    sigma[i] <- 1/sqrt(tau[i]) 
    # delta 
    delta[i] ~ dnorm(omega, tau_C[i])      
    tau_C[i] ~ dgamma(0.001, 0.001)     
    # C_i
    C[i] <- 1/tau_C[i]
  }

  
  # theta_0 
  theta_0 ~ dnorm(mu_0, tau_0)
  theta_0_hat ~ dnorm(theta_0, tau_0_hat)
  tau_0_hat ~ dgamma(0.001,0.001)
  sigma2_0 <- 1/tau_0_hat
  sigma_0 <- 1/sqrt(tau_0_hat)
  
  # omega
  omega ~ dnorm(0, tau_lambda)
  tau_lambda ~ dgamma(0.001, 0.001)
  lambda <- 1/tau_lambda
}

write.model(OneWay_Heatwave,"OneWay_HeatwaveCount.txt")


############# Bayes data #################
GCM.names <- unique(simulators$GCM)
num_GCMs <- length(GCM.names)

y <- matrix(NA, nrow = num_GCMs, ncol = 20)  

for (i in 1:num_GCMs) {
  gcm_name <- GCM.names[i]
  gcm_data <- subset(simulators, GCM == gcm_name)
  gcm_data <- gcm_data[order(gcm_data$Year), ]
  y[i, ] <- gcm_data$Heatwave_Count[1:20]
}
dim(y)  # 6 20 

BayesModelData <- list(
  num_GCMs = num_GCMs,
  y = y,
  mu_0 = mean(obs$Heatwave_Count),
  tau_0 = 1e-3 ,
  theta_0_hat = mean(obs$Heatwave_Count)
)

bugs.data(BayesModelData, data.file = "Data.txt") 


####### Initial values ########
BayesInits <- list(
  list(
    tau = rep(0.001, num_GCMs),
    theta_0 = 0.001,
    tau_0_hat = 0.001,
    delta = rep(0, num_GCMs),
    tau_C = rep(0.01, num_GCMs),
    omega = 0,  
    tau_lambda = 0.01
  )
)


######## Bugs #########
OneWay_GCM_Model <- bugs(
  data = BayesModelData, 
  n.chains=1,
  inits = BayesInits,   
  n.burnin=0,
  n.iter=1000,
  n.thin = 1,
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_HeatwaveCount.txt", 
  debug = FALSE
)

OneWay_GCM_Model 

if (dev.cur()==1) x11(width=50, height=30) 
plot(as.mcmc.list(OneWay_GCM_Model),smooth=FALSE)


######### More chains ##################
BayesInits <- list(
  list(
    tau = rep(1, num_GCMs),
    theta_0 = 0.0248, 
    tau_0_hat = 1,
    delta = rep(0.01, num_GCMs),
    tau_C = rep(0.1, num_GCMs),
    omega = 0.01,
    tau_lambda = 0.1
  ),
  list(
    tau = rep(1, num_GCMs),
    theta_0 = 0.025,
    tau_0_hat = 1,
    delta = rep(0.015, num_GCMs),
    tau_C = rep(0.1, num_GCMs),
    omega = 0.015,
    tau_lambda = 0.1
  ),
  list(
    tau = rep(1, num_GCMs),
    theta_0 = 0.0245,
    tau_0_hat = 1,
    delta = rep(0.005, num_GCMs),
    tau_C = rep(0.1, num_GCMs),
    omega = 0.005,
    tau_lambda = 0.1
  ),
  list(
    tau = rep(1, num_GCMs),
    theta_0 = 0.024,
    tau_0_hat = 1,
    delta = rep(0.02, num_GCMs),
    tau_C = rep(0.1, num_GCMs),
    omega = 0.02,
    tau_lambda = 0.1
  )
)



OneWay_Model1a <- bugs(    ### 4 chains, 100 burn-in, 10000 iterations 
  data = BayesModelData, 
  inits = BayesInits,   
  parameters.to.save = c("theta_0","sigma","sigma_0"),   
  model.file = "OneWay_HeatwaveCount.txt", 
  n.chains = 4,  
  n.burnin = 100,       
  n.iter = 10000, 
  n.thin = 1, 
  debug = FALSE
)

OneWay_Model1a
plot(as.mcmc.list(OneWay_Model1a),smooth=FALSE)
plot(acfplot(as.mcmc.list(OneWay_Model1a)))




OneWay_Heatwave1b <- bugs(    ### 4 chains, 100 burn-in, 10000 iterations, 10 thins
  data = BayesModelData, 
  inits = BayesInits,   
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_HeatwaveCount.txt", 
  n.chains = 4,  
  n.burnin = 1000,       
  n.iter = 10000, 
  n.thin = 24, 
  debug = FALSE
)

OneWay_Heatwave1b
plot(as.mcmc.list(OneWay_Heatwave1b),smooth=FALSE)
plot(acfplot(as.mcmc.list(OneWay_Heatwave1b)))


OneWay_Heatwave1c <- bugs(    ### 4 chains, 100 burn-in, 10000 iterations, 10 thins
  data = BayesModelData, 
  inits = BayesInits,   
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_HeatwaveCount.txt", 
  n.chains = 4,  
  n.burnin = 100,       
  n.iter = 10000, 
  n.thin = 1, 
  debug = FALSE
)

OneWay_Heatwave1c

########################## Half Cauchy #################################
OneWay_Heatwave_half <- function() {
  for (i in 1:num_GCMs) {   
    for (t in 1:20) {      
      y[i,t] ~ dnorm(theta[i], tau[i])  
    }
    
    # theta_i
    theta[i] <- theta_0 + delta[i]
    tau[i] <- pow(sigma[i], -2)
    d[i] ~ dt(0,10,1)		# Scaled Cauchy
    sigma[i] <- abs(d[i])
    
    # delta_i
    delta[i] ~ dnorm(omega, tau_C[i])      
    tau_C[i] <- pow(C[i], -1)  
    
    # C_i
    t_C[i] ~ dt(0,10,1)
    C[i] <- abs(t_C[i])
  }
  
  # theta_0 
  theta_0 ~ dnorm(mu_0, tau_0)
  theta_0_hat ~ dnorm(theta_0, tau_0_hat)
  tau_0_hat <- pow(sigma_0, -2)
  t_0 ~ dt(0,10,1)
  sigma_0 <- abs(t_0)
  
  # omega
  omega ~ dnorm(0, tau_lambda)
  t_lambda ~ dt(0,10,1)
  lambda <- abs(t_lambda)
  tau_lambda <- pow(lambda, -1)  
}


write.model(OneWay_Heatwave_half,"OneWay_Heatwave_half.txt")

BayesModelData <- list(
  num_GCMs = num_GCMs,
  y = y,
  mu_0 = mean(obs$Heatwave_Count),
  tau_0 = 1e-3 ,
  theta_0_hat = mean(obs$Heatwave_Count)
)

bugs.data(BayesModelData, data.file = "Data.txt") 


####### Initial values 
BayesInits <- list(
  list(
    d = rep(0.001, num_GCMs),
    theta_0 = 0.001,
    t_0 = 0.001,
    delta = rep(0, num_GCMs),
    t_C = rep(0.01, num_GCMs),
    omega = 0,  
    t_lambda = 0.01
  )
)



OneWay_GCM_Model <- bugs(
  data = BayesModelData, 
  n.chains=1,
  inits = BayesInits,   
  n.burnin=0,
  n.iter=1000,
  n.thin = 1,
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_Heatwave_half.txt", 
  debug = TRUE
)

OneWay_GCM_Model 



############## More chains 
BayesInits <- list(
  list(
    d = rep(0.001, num_GCMs),
    theta_0 = 0.001,
    t_0 = 0.001,
    delta = rep(0, num_GCMs),
    t_C = rep(0.01, num_GCMs),
    omega = 0,  
    t_lambda = 0.01
  ),
  list(
    d = rep(0.002, num_GCMs),
    theta_0 = 0.002,
    t_0 = 0.002,
    delta = rep(0.001, num_GCMs),
    t_C = rep(0.02, num_GCMs),
    omega = 0.001,  
    t_lambda = 0.02
  ),
  list(
    d = rep(0.003, num_GCMs),
    theta_0 = 0.003,
    t_0 = 0.003,
    delta = rep(-0.001, num_GCMs),
    t_C = rep(0.015, num_GCMs),
    omega = -0.001,  
    t_lambda = 0.015
  ),
  list(
    d = rep(0.004, num_GCMs),
    theta_0 = 0.004,
    t_0 = 0.004,
    delta = rep(0.002, num_GCMs),
    t_C = rep(0.025, num_GCMs),
    omega = 0.002,  
    t_lambda = 0.025
  )
)



OneWay_half1a <- bugs(
  data = BayesModelData, 
  n.chains=4,
  inits = BayesInits,   
  n.burnin=0,
  n.iter=1000,
  n.thin = 1,
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_Heatwave_half.txt", 
  debug = FALSE
)

OneWay_half1a

OneWay_half1b <- bugs(
  data = BayesModelData, 
  n.chains = 4,
  inits = BayesInits,   
  n.burnin = 10,
  n.iter=10000,
  n.thin = 1,
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_Heatwave_half.txt", 
  debug = FALSE
)

OneWay_half1b






