##################### Data #####################################
library(mvtnorm)
region_ID <- 2
all_simulators_data <- process_simulator_data(region_ID)
all_simulators_data_nozeros <- filter_nonzero(all_simulators_data)
observations <- read.csv("observation.csv")
observations_nozeros <- filter_nonzero(observations)

simulators_historic_DI <- subset(all_simulators_data_nozeros, Year >= 1981 & Year <= 2000 & Region_ID == region_ID)
simulators_future_DI <- subset(all_simulators_data_nozeros, Year >= 2026 & Year <= 2045 & Region_ID == region_ID)
obs <- subset(observations_nozeros, Year >= 1981 & Year <= 2000 & Region_ID == region_ID)
######## Def ##########
TwoWay_DI <- function(){
  for (i in 1:num_GCMs) {   
    for (j in 1:num_RCMs) {  
      for (r in 1:nruns[i, j, 1]) {  
        y[i, j, r, 1:2] ~ dmnorm(theta[i, j, 1:2], Tau_h[i, j, 1:2, 1:2]) 
      }
      for (r in 1:nruns[i, j, 2]) { 
        y[i, j, r, 3:4] ~ dmnorm(theta[i, j, 3:4], Tau_f[i, j, 1:2, 1:2]) 
      }
      
      # theta_ij
      for (k in 1:4) {
        theta[i, j, k] <- theta_0[k] + delta[i, k] + zeta[j, k]
      }
      
      # Tau_ij
      Tau_h[i, j, 1:2, 1:2] ~ dwish(V_h[i, j, 1:2, 1:2], 3)
      Tau_f[i, j, 1:2, 1:2] ~ dwish(V_f[i, j, 1:2, 1:2], 3)
      sigma2_h[i, j, 1:2, 1:2] <- inverse(Tau_h[i, j, 1:2, 1:2])
      for (k in 1:2) {
        sigma_h[i, j, k] <- sqrt(sigma2_h[i, j, k, k])
      }
      sigma2_f[i, j, 1:2, 1:2] <- inverse(Tau_f[i, j, 1:2, 1:2])
      for (k in 1:2) {
        sigma_f[i, j, k] <- sqrt(sigma2_f[i, j, k, k])
      }
    }
    
    # delta
    delta[i, 1:4] ~ dmnorm(omega[1:4], Tau_C[i, 1:4, 1:4]) 
    # C
    Tau_C[i, 1:4, 1:4] ~ dwish(V_C[i,1:4, 1:4], 5)
    C[i, 1:4, 1:4] <- inverse(Tau_C[i, ,]) 
  }
  
  # zeta
  for (j in 1:num_RCMs) {
    zeta[j, 1:4] ~ dmnorm(eta[1:4], Tau_D[j, 1:4, 1:4])
    Tau_D[j, 1:4, 1:4] ~ dwish(V_D[j, 1:4, 1:4], 5)
    D[j, 1:4, 1:4] <- inverse(Tau_D[j, ,])
  }
  
  # theta_0 
  theta_0[1:4] ~ dmnorm(mu_0[1:4], Tau_0[1:4, 1:4])
  # theta_0_hat
  theta_0_hat[1:4] ~ dmnorm(theta_0[1:4], Tau_0_hat[1:4, 1:4])
  Tau_0_hat[1:4,1:4] ~ dwish(V_0[1:4, 1:4], 5)
  sigma2_0[1:4,1:4] <- inverse(Tau_0_hat[1:4,1:4])
  for (j in 1:4) {
    sigma_0[j] <- sqrt(sigma2_0[j, j])
  }
  
  # omega
  for (i in 1:4) {
    zero[i] <- 0
  }
  omega[1:4] ~ dmnorm(zero[], Tau_lambda[1:4, 1:4])
  # lambda
  Tau_lambda[1:4, 1:4] ~ dwish(V_lambda[1:4, 1:4], 5)
  lambda[1:4,1:4] <- inverse(Tau_lambda[1:4,1:4])
  
  # eta 
  eta[1:4] ~ dmnorm(zero[], Tau_chi[1:4, 1:4])
  # chi
  Tau_chi[1:4, 1:4] ~ dwish(V_chi[1:4, 1:4], 5)
  chi[1:4,1:4] <- inverse(Tau_chi[1:4,1:4])
}

library(R2OpenBUGS); library(coda); library(nlme); library(lme4)
write.model(TwoWay_DI,"TwoWay_DI.txt")

########## Data ###########
GCM.names <- unique(all_simulators_data$GCM)
num_GCMs <- length(GCM.names)

RCM.names <- unique(all_simulators_data$RCM)
num_RCMs <- length(RCM.names)

runs_table_historic <- table(simulators_historic_DI$GCM, simulators_historic_DI$RCM)
runs_table_future <- table(simulators_future_DI$GCM,simulators_future_DI$RCM)


y <- array(NA, dim = c(num_GCMs, num_RCMs, 20, 4))
for (gcm in 1:num_GCMs){
  for (rcm in 1:num_RCMs){
    wanted.rows_historic <- subset(simulators_historic_DI, GCM == GCM.names[gcm] & RCM == RCM.names[rcm])
    wanted.rows_future <- subset(simulators_future_DI, GCM == GCM.names[gcm] & RCM == RCM.names[rcm])
    num_runs_historic <- nrow(wanted.rows_historic)
    num_runs_future <- nrow(wanted.rows_future)
    if (num_runs_historic > 0) {
      y[gcm, rcm, 1:num_runs_historic, 1] <- wanted.rows_historic$Mean_Duration
      y[gcm, rcm, 1:num_runs_historic, 2] <- wanted.rows_historic$Mean_Intensity
    }
    if (num_runs_future > 0) {
      y[gcm, rcm, 1:num_runs_future, 3] <- wanted.rows_future$Mean_Duration
      y[gcm, rcm, 1:num_runs_future, 4] <- wanted.rows_future$Mean_Intensity
    }
  }
}

dim(y) #6 10 20  4

nruns <- array(NA, dim = c(num_GCMs, num_RCMs,2))
nruns[, , 1] <- runs_table_historic
nruns[, , 2] <- runs_table_future

dim(nruns)  # 6 10  2


BayesModelData <- list(
  num_GCMs = num_GCMs,
  num_RCMs = num_RCMs,  
  nruns = nruns,  # pigs
  y = y, 
  V_h = array(rep(diag(c(3,30)), each = num_GCMs * num_RCMs), dim = c(num_GCMs, num_RCMs, 2, 2)), 
  V_f = array(rep(diag(c(3,30)), each = num_GCMs * num_RCMs), dim = c(num_GCMs, num_RCMs, 2, 2)), 
  V_C = array(rep(diag(c(3,30,3,30)), each = num_GCMs), dim = c(num_GCMs, 4, 4)),  
  V_D = array(rep(diag(c(3,30,3,30)), each = num_RCMs), dim = c(num_RCMs, 4, 4)),  
  V_lambda = diag(c(3,30,3,30)), 
  V_chi = diag(c(3,30,3,30)),
  mu_0 = c(A$theta_0_hat[1:2],A$E_theta_hat_future[1:2]),  
  Tau_0 = solve(diag(c(3,30,3,30))),  
  theta_0_hat = c(mean(obs$Mean_Duration), mean(obs$Mean_Intensity),A$E_theta_hat_future[1:2]), 
  V_0 = diag(c(3,30,3,30))               
)


####### Initial values ########
BayesInits <- list(
  list(
    Tau_f = array(rep(solve(diag(c(3,30))), each = num_GCMs * num_RCMs), 
                dim = c(num_GCMs, num_RCMs, 2, 2)),  
    Tau_h = array(rep(solve(diag(c(3,30))), each = num_GCMs * num_RCMs), 
                  dim = c(num_GCMs, num_RCMs, 2, 2)),
    theta_0 = c(3, 30,3,30),  
    Tau_0_hat = solve(diag(c(3,30,3,30))),  
    
    # Delta
    delta = matrix(0, num_GCMs, 4),  
    Tau_C = array(rep(solve(diag(c(3,30,3,30))), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
    
    # Zeta
    zeta = matrix(0, num_RCMs, 4),  
    Tau_D = array(rep(solve(diag(c(3,30,3,30))), each = num_RCMs), dim = c(num_RCMs, 4, 4)), 
    
    # Omega & Eta
    omega = c(0, 0,0,0),  
    Tau_lambda = solve(diag(c(3,30,3,30))),
    eta = c(0, 0,0,0),  
    Tau_chi = solve(diag(c(3,30,3,30)))
  )
)


TwoWay_DI <- bugs(
  data = BayesModelData, 
  inits = BayesInits,   
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "TwoWay_DI.txt", 
  n.chains = 1,  
  n.burnin = 0,       
  n.iter = 1000, 
  n.thin = 1,
  debug = TRUE
)

TwoWay_DI

############### More chains ########################
BayesInits <- list(
  list(
    Tau_f = array(rep(solve(diag(c(3,30))), each = num_GCMs * num_RCMs), 
                   dim = c(num_GCMs, num_RCMs, 2, 2)),  
    Tau_h = array(rep(solve(diag(c(3,30))), each = num_GCMs * num_RCMs), 
                  dim = c(num_GCMs, num_RCMs, 2, 2)),  
    theta_0 = c(3, 30),  
    Tau_0_hat = solve(diag(c(3,30))),  
    
    # Delta
    delta = matrix(rnorm(num_GCMs * 2, mean = 0, sd = 0.1), num_GCMs, 2),  
    Tau_C = array(rep(solve(diag(c(3,30))), each = num_GCMs), dim = c(num_GCMs, 2, 2)), 
    
    # Zeta
    zeta = matrix(rnorm(num_RCMs * 2, mean = 0, sd = 0.1), num_RCMs, 2),  
    Tau_D = array(rep(solve(diag(c(3,30))), each = num_RCMs), dim = c(num_RCMs, 2, 2)), 
    
    # Omega & Eta
    omega = c(0, 0 ,0,0),  
    Tau_lambda = solve(diag(c(3,30))),
    eta = c(0, 0,0,0),  
    Tau_chi = solve(diag(c(3,30)))
  ),
  
  list(
    Tau_f = array(rep(solve(diag(c(3,30) * 1.1)), each = num_GCMs * num_RCMs), 
                  dim = c(num_GCMs, num_RCMs, 2, 2)),  
    Tau_h = array(rep(solve(diag(c(3,30) * 1.1)), each = num_GCMs * num_RCMs), 
                  dim = c(num_GCMs, num_RCMs, 2, 2)),
    theta_0 = c(3.2, 29.8),  
    Tau_0_hat = solve(diag(c(3,30) * 1.1)),  
    
    delta = matrix(rnorm(num_GCMs * 2, mean = 0.05, sd = 0.1), num_GCMs, 2),  
    Tau_C = array(rep(solve(diag(c(3,30) * 1.1)), each = num_GCMs), dim = c(num_GCMs, 2, 2)), 
    
    zeta = matrix(rnorm(num_RCMs * 2, mean = 0.05, sd = 0.1), num_RCMs, 2),  
    Tau_D = array(rep(solve(diag(c(3,30) * 1.1)), each = num_RCMs), dim = c(num_RCMs, 2, 2)), 
    
    omega = c(0.05, -0.05,0.05,-0.05),  
    Tau_lambda = solve(diag(c(3,30) * 1.1)),
    eta = c(0.05, -0.05,0.05,-0.05),  
    Tau_chi = solve(diag(c(3,30) * 1.1))
  ),
  
  list(
    Tau_f = array(rep(solve(diag(c(3,30) * 0.9)), each = num_GCMs * num_RCMs), 
                  dim = c(num_GCMs, num_RCMs, 2, 2)),  
    Tau_h = array(rep(solve(diag(c(3,30) * 0.9)), each = num_GCMs * num_RCMs), 
                  dim = c(num_GCMs, num_RCMs, 2, 2)),
    theta_0 = c(2.8, 30.2),  
    Tau_0_hat = solve(diag(c(3,30) * 0.9)),  
    
    delta = matrix(rnorm(num_GCMs * 2, mean = -0.05, sd = 0.1), num_GCMs, 2),  
    Tau_C = array(rep(solve(diag(c(3,30) * 0.9)), each = num_GCMs), dim = c(num_GCMs, 2, 2)), 
    
    zeta = matrix(rnorm(num_RCMs * 2, mean = -0.05, sd = 0.1), num_RCMs, 2),  
    Tau_D = array(rep(solve(diag(c(3,30) * 0.9)), each = num_RCMs), dim = c(num_RCMs, 2, 2)), 
    
    omega = c(0.05, -0.05,0.05,-0.05),  
    Tau_lambda = solve(diag(c(3,30) * 0.9)),
    eta = c(0.05, -0.05,0.05,-0.05),  
    Tau_chi = solve(diag(c(3,30) * 0.9))
  ),
  
  list(
    Tau_f = array(rep(solve(diag(c(3,30) * 1.2)), each = num_GCMs * num_RCMs), 
                  dim = c(num_GCMs, num_RCMs, 2, 2)),  
    Tau_h = array(rep(solve(diag(c(3,30) * 1.2)), each = num_GCMs * num_RCMs), 
                  dim = c(num_GCMs, num_RCMs, 2, 2)),
    theta_0 = c(3.5, 29.5),  
    Tau_0_hat = solve(diag(c(3,30) * 1.2)),  
    
    delta = matrix(rnorm(num_GCMs * 2, mean = 0.1, sd = 0.1), num_GCMs, 2),  
    Tau_C = array(rep(solve(diag(c(3,30) * 1.2)), each = num_GCMs), dim = c(num_GCMs, 2, 2)), 
    
    zeta = matrix(rnorm(num_RCMs * 2, mean = 0.1, sd = 0.1), num_RCMs, 2),  
    Tau_D = array(rep(solve(diag(c(3,30) * 1.2)), each = num_RCMs), dim = c(num_RCMs, 2, 2)), 
    
    omega = c(0.1, -0.1,0.1,-0.1),  
    Tau_lambda = solve(diag(c(3,30) * 1.2)),
    eta = c(0.1, -0.1,0.1,-0.1),  
    Tau_chi = solve(diag(c(3,30) * 1.2))
  )
)




TwoWay_Model1a <- bugs(
  data = BayesModelData, 
  inits = BayesInits,   
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "TwoWay_DI.txt", 
  n.chains = 4,  
  n.burnin = 100,       
  n.iter = 10000, 
  debug = FALSE
)
TwoWay_Model1a   # 全是1 DIC = 1771

TwoWay_Model1b <- bugs(
  data = BayesModelData, 
  inits = BayesInits,   
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "TwoWay_DI.txt", 
  n.chains = 4,  
  n.burnin = 500,       
  n.iter = 10000, 
  n.thin = 4,
  debug = FALSE
)
TwoWay_Model1b





plot(as.mcmc.list(TwoWay_Model1b),smooth=FALSE)
plot(acfplot(as.mcmc.list(TwoWay_Model1b)))
plot(TwoWay_Model1b)



