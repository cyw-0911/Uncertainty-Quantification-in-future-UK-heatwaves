##################### Data #####################################
library(mvtnorm)
all_simulators_data <- read.csv("all_simulators_data.csv")
all_simulators_data_nozeros <- filter_nonzero(all_simulators_data)
observations <- read.csv("observation.csv")
observations_nozeros <- filter_nonzero(observations)

simulators_historic_DI <- subset(all_simulators_data_nozeros, Year >= 1981 & Year <= 2000 & Region_ID == 1)
simulators_future_DI <- subset(all_simulators_data_nozeros, Year >= 2026 & Year <= 2045 & Region_ID == 1)
obs <- subset(observations_nozeros, Year >= 1981 & Year <= 2000 & Region_ID == 1)
######## Def ##########
OneWay_DI <- function(){
  for (i in 1:num_GCMs) {   
    for (j in 1:nruns[i]) {      
      y[i, j, 1:4] ~ dmnorm(theta[i, 1:4], Tau[i, 1:4, 1:4])  
    }
    
    # theta_i
    for (j in 1:4) {
      theta[i, j] <- theta_0[j] + delta[i, j]
    }
    
    # Tau
    Tau[i, 1:4, 1:4] ~ dwish(V[i, 1:4, 1:4], 5)
    sigma2[i, 1:4, 1:4] <- inverse(Tau[i, 1:4, 1:4])
    for (j in 1:4){
      sigma[i, j] <- sqrt(sigma2[i, j, j])
    }
    
    # delta
    delta[i, 1:4] ~ dmnorm(omega[1:4], Tau_C[i, 1:4, 1:4]) 
    
    # C
    Tau_C[i, 1:4, 1:4] ~ dwish(V_C[i, 1:4, 1:4], 5)
    C[i, 1:4, 1:4] <- inverse(Tau_C[i, ,]) 
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
  lambda[1:4, 1:4] <- inverse(Tau_lambda[1:4, 1:4])
}


write.model(OneWay_DI,"OneWay_DI.txt")


########## Data ###########
GCM.names <- unique(simulators_future_DI$GCM)
num_GCMs <- length(GCM.names)

runs_table_historic <- table(simulators_historic_DI$GCM)
runs_table_future <- table(simulators_future_DI$GCM)

y <- array(NA, dim = c(num_GCMs,max(runs_table_historic,runs_table_future),4))
row.names(y) <- GCM.names

for (model in 1:num_GCMs){
  historic_data <- simulators_historic_DI$GCM == GCM.names[model]
  future_data <- simulators_future_DI$GCM == GCM.names[model]
  if (runs_table_historic[model] > 0 & runs_table_future[model] > 0){
    y[model, 1:runs_table_historic[model],1] <- simulators_historic_DI$Mean_Duration[historic_data]
    y[model, 1:runs_table_historic[model],2] <- simulators_historic_DI$Mean_Intensity[historic_data]
    y[model, 1:runs_table_future[model],3] <- simulators_future_DI$Mean_Duration[future_data]
    y[model, 1:runs_table_future[model],4] <- simulators_future_DI$Mean_Intensity[future_data]
  }
}

dim(y)  # 6 81  4


BayesModelData <- list(
  num_GCMs = num_GCMs,
  nruns = c(as.integer(runs_table_historic),as.integer(runs_table_future)),  
  # nruns_future = as.integer(runs_table_future),
  y = y,                           
  V = array(rep(diag(c(3,30,3,30)), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
  mu_0 = c(A$theta_0_hat[1:2],A$E_theta_hat_future[1:2]), 
  Tau_0 = solve(diag(c(3,30,3,30))), 
  theta_0_hat = c(mean(obs$Mean_Duration), mean(obs$Mean_Intensity),A$E_theta_hat_future[1:2]), 
  V_0 = diag(c(3,30,3,30)),               
  V_C = array(rep(diag(c(3,30,3,30)), each=num_GCMs), dim = c(num_GCMs, 4, 4)),  
  V_lambda = diag(c(3,30,3,30)) 
)

####### Initial values ########
BayesInits <- list(
  list(
    Tau = array(rep(solve(diag(c(3,30,3,30))), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
    theta_0 = c(3, 30,3,30), 
    Tau_0_hat = solve(diag(c(3,30,3,30))),
    delta = matrix(0, num_GCMs, 4), 
    Tau_C = array(rep(solve(diag(c(3,30,3,30))), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
    omega = c(0, 0,0,0), 
    Tau_lambda = solve(diag(c(3,30,3,30)))
  )
)


OneWay_DI <- bugs(
  data = BayesModelData, 
  inits = BayesInits,   
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_DI.txt", 
  n.chains = 1,  
  n.burnin = 1000,       
  n.iter = 10000, 
  debug = TRUE
)

OneWay_HC


######### More chains ############
BayesInits <- list(
  list(
    Tau = array(rep(solve(diag(c(3,30,3,30))), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
    theta_0 = c(3, 30, 3, 30), 
    Tau_0_hat = solve(diag(c(3,30,3,30))),
    delta = matrix(0, num_GCMs, 4), 
    Tau_C = array(rep(solve(diag(c(3,30,3,30))), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
    omega = c(0, 0, 0, 0), 
    Tau_lambda = solve(diag(c(3,30,3,30)))
  ),
  list(
    Tau = array(rep(solve(diag(c(3.5, 29, 3.5, 29))), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
    theta_0 = c(3.5, 29, 3.5, 29), 
    Tau_0_hat = solve(diag(c(3.5, 29, 3.5, 29))),
    delta = matrix(0.1, num_GCMs, 4), 
    Tau_C = array(rep(solve(diag(c(3.5, 29, 3.5, 29))), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
    omega = c(0.1, 0.1, 0.1, 0.1), 
    Tau_lambda = solve(diag(c(3.5, 29, 3.5, 29)))
  ),
  list(
    Tau = array(rep(solve(diag(c(2.8, 31, 2.8, 31))), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
    theta_0 = c(2.8, 31, 2.8, 31), 
    Tau_0_hat = solve(diag(c(2.8, 31, 2.8, 31))),
    delta = matrix(-0.1, num_GCMs, 4), 
    Tau_C = array(rep(solve(diag(c(2.8, 31, 2.8, 31))), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
    omega = c(-0.1, -0.1, -0.1, -0.1), 
    Tau_lambda = solve(diag(c(2.8, 31, 2.8, 31)))
  ),
  list(
    Tau = array(rep(solve(diag(c(3.2, 28.5, 3.2, 28.5))), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
    theta_0 = c(3.2, 28.5, 3.2, 28.5), 
    Tau_0_hat = solve(diag(c(3.2, 28.5, 3.2, 28.5))),
    delta = matrix(0.05, num_GCMs, 4), 
    Tau_C = array(rep(solve(diag(c(3.2, 28.5, 3.2, 28.5))), each = num_GCMs), dim = c(num_GCMs, 4, 4)), 
    omega = c(0.05, 0.05, 0.05, 0.05), 
    Tau_lambda = solve(diag(c(3.2, 28.5, 3.2, 28.5)))
  )
)



OneWay_DI1a <- bugs(
  data = BayesModelData, 
  inits = BayesInits,   
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_DI.txt", 
  n.chains = 4,  
  n.burnin = 0,       
  n.iter = 10000, 
  debug = FALSE
)

OneWay_DI1a


OneWay_DI1b <- bugs(
  data = BayesModelData, 
  inits = BayesInits,   
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_DI.txt", 
  n.chains = 4,  
  n.burnin = 1000,       
  n.iter = 10000, 
  n.thin = 1,
  debug = FALSE
)

OneWay_DI1b




OneWay_DI1c <- bugs(
  data = BayesModelData, 
  inits = BayesInits,   
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_DI.txt", 
  n.chains = 4,  
  n.burnin = 1000,       
  n.iter = 10000, 
  n.thin = 3,
  debug = FALSE
) # 3

OneWay_DI1c





OneWay_DI1d <- bugs(
  data = BayesModelData, 
  inits = BayesInits,   
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_DI.txt", 
  n.chains = 4,  
  n.burnin = 1000,       
  n.iter = 10000, 
  n.thin = 5,
  debug = FALSE
) # 3

OneWay_DI1d





