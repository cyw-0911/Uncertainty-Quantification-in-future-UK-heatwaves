##################### Data #####################################
library(mvtnorm)
all_simulators_data <- read.csv("all_simulators_data.csv")
region_ID = 2
all_simulators_data <- process_simulator_data(region_ID)
observations <- read.csv("observation.csv")

simulators_historic <- subset(all_simulators_data, Year >= 1981 & Year <= 2000 & Region_ID == region_ID)
simulators_future <- subset(all_simulators_data, Year >= 2026 & Year <= 2045 & Region_ID == region_ID) 
obs <- subset(observations, Year >= 1981 & Year <= 2000 & Region_ID == region_ID)
######## Def ###########
OneWay_Heatwave <- function(){
  for (i in 1:num_GCMs) {   
    for (j in 1:20) {      
      y[i, j, 1] ~ dnorm(theta[i, 1], tau[i, 1])  # Historical
      y[i, j, 2] ~ dnorm(theta[i, 2], tau[i, 2])  # Future
    }
    
    # theta_i 
    for (j in 1:2) {
      theta[i, j] <- theta_0[j] + delta[i, j]
    }
    
    for (j in 1:2) {
      tau[i, j] ~ dgamma(0.01, 0.01)  
      sigma2[i, j] <- 1/tau[i, j]
      sigma[i, j] <- sqrt(sigma2[i, j])
    }
    
    # delta
    for (j in 1:2) {
      delta[i, j] ~ dnorm(omega[j], tau_C[i, j])      
      tau_C[i, j] ~ dgamma(0.01, 0.01) 
      C[i, j] <- 1/tau_C[i, j]
    }
  }
  
  # theta_0 
  for (j in 1:2) {
    theta_0[j] ~ dnorm(mu_0[j], tau_0[j]) 
    
    theta_0_hat[j] ~ dnorm(theta_0[j], tau_0_hat[j])  
    tau_0_hat[j] ~ dgamma(0.01, 0.01)  
    sigma2_0[j] <- 1/tau_0_hat[j]
    sigma_0[j] <- sqrt(sigma2_0[j])
    
    # omega
    omega[j] ~ dnorm(0, tau_lambda[j])  
    tau_lambda[j] ~ dgamma(0.01, 0.01)  
    lambda[j] <- 1/tau_lambda[j]
  }
}
write.model(OneWay_Heatwave,"OneWay_HeatwaveCount.txt")



############# Bayes data #################
GCM.names <- unique(simulators_historic$GCM)
num_GCMs <- length(GCM.names)

y <- array(NA, dim = c(num_GCMs,20,2))  
for (i in 1:num_GCMs) {
  gcm_name <- GCM.names[i]
  historic_data <- subset(simulators_historic, GCM == gcm_name)
  historic_mean <- tapply(historic_data$Heatwave_Count, historic_data$Year, mean)  
  y[i, ,1] <- historic_mean
  
  
  future_data <- subset(simulators_future, GCM == gcm_name)
  future_mean <- tapply(future_data$Heatwave_Count, future_data$Year, mean)  
  y[i, ,2] <- future_mean
}
dim(y)  # 6 20 2

BayesModelData <- list(
  num_GCMs = num_GCMs,
  y = y,
  mu_0 = c(mean(obs$Heatwave_Count),mean(simulators_future$Heatwave_Count)),
  tau_0 = c(0.01,0.01),
  theta_0_hat = c(mean(obs$Heatwave_Count),mean(simulators_future$Heatwave_Count))
)

BayesInits <- list(
  list(
    tau = matrix(rep(0.01 , num_GCMs * 2), nrow = num_GCMs),
    theta_0 = c(mean(obs$Heatwave_Count),mean(simulators_future$Heatwave_Count)),
    tau_0_hat = c(0.01,0.01),
    delta = matrix(rep(0 , num_GCMs * 2), nrow = num_GCMs),
    tau_C = matrix(rep(0.01 , num_GCMs * 2), nrow = num_GCMs),
    omega = c(0,0),  
    tau_lambda = c(0.01,0.01)
  )
)

######## Bugs #########
OneWay_HC <- bugs(
  data = BayesModelData, 
  n.chains= 1,
  inits = BayesInits,   
  n.burnin= 0,
  n.iter= 1000,
  n.thin = 1,
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_HeatwaveCount.txt", 
  debug = FALSE
)

OneWay_HC


######### More chains ############
BayesInits <- list(
  list(
    tau = matrix(runif(num_GCMs * 2, 0.005, 0.02), nrow = num_GCMs),
    theta_0 = c(mean(obs$Heatwave_Count) + runif(1, -0.05, 0.05),
                mean(simulators_future$Heatwave_Count) + runif(1, -0.05, 0.05)),
    tau_0_hat = runif(2, 0.005, 0.02),
    delta = matrix(runif(num_GCMs * 2, -0.05, 0.05), nrow = num_GCMs),
    tau_C = matrix(runif(num_GCMs * 2, 0.005, 0.02), nrow = num_GCMs),
    omega = runif(2, -0.05, 0.05),  
    tau_lambda = runif(2, 0.005, 0.02)
  ),
  list(
    tau = matrix(runif(num_GCMs * 2, 0.01, 0.03), nrow = num_GCMs),
    theta_0 = c(mean(obs$Heatwave_Count) + runif(1, -0.1, 0.1),
                mean(simulators_future$Heatwave_Count) + runif(1, -0.1, 0.1)),
    tau_0_hat = runif(2, 0.01, 0.03),
    delta = matrix(runif(num_GCMs * 2, -0.1, 0.1), nrow = num_GCMs),
    tau_C = matrix(runif(num_GCMs * 2, 0.01, 0.03), nrow = num_GCMs),
    omega = runif(2, -0.1, 0.1),  
    tau_lambda = runif(2, 0.01, 0.03)
  ),
  list(
    tau = matrix(runif(num_GCMs * 2, 0.002, 0.015), nrow = num_GCMs),
    theta_0 = c(mean(obs$Heatwave_Count) + runif(1, -0.02, 0.02),
                mean(simulators_future$Heatwave_Count) + runif(1, -0.02, 0.02)),
    tau_0_hat = runif(2, 0.002, 0.015),
    delta = matrix(runif(num_GCMs * 2, -0.02, 0.02), nrow = num_GCMs),
    tau_C = matrix(runif(num_GCMs * 2, 0.002, 0.015), nrow = num_GCMs),
    omega = runif(2, -0.02, 0.02),  
    tau_lambda = runif(2, 0.002, 0.015)
  ),
  list(
    tau = matrix(runif(num_GCMs * 2, 0.007, 0.025), nrow = num_GCMs),
    theta_0 = c(mean(obs$Heatwave_Count) + runif(1, -0.08, 0.08),
                mean(simulators_future$Heatwave_Count) + runif(1, -0.08, 0.08)),
    tau_0_hat = runif(2, 0.007, 0.025),
    delta = matrix(runif(num_GCMs * 2, -0.08, 0.08), nrow = num_GCMs),
    tau_C = matrix(runif(num_GCMs * 2, 0.007, 0.025), nrow = num_GCMs),
    omega = runif(2, -0.08, 0.08),  
    tau_lambda = runif(2, 0.007, 0.025)
  )
)


OneWay_HC1a <- bugs(
  data = BayesModelData, 
  n.chains= 4,
  inits = BayesInits,   
  n.burnin= 100,
  n.iter= 10000,
  n.thin = 1,
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_HeatwaveCount.txt", 
  debug = FALSE
)

OneWay_HC1a

OneWay_HC1b <- bugs(
  data = BayesModelData, 
  n.chains= 4,
  inits = BayesInits,   
  n.burnin= 500,
  n.iter= 10000,
  n.thin = 19,
  parameters.to.save = c("theta_0","sigma_0"),   
  model.file = "OneWay_HeatwaveCount.txt", 
  debug = TRUE
) # 19 

if (dev.cur()==1) x11(width=9, height=7) # Open large-ish graphics window 

OneWay_HC1b

plot(as.mcmc.list(OneWay_HC1b),smooth=FALSE)
plot(acfplot(as.mcmc.list(OneWay_HC1b)))
png("acfplot.png", width = 900, height = 700)
plot(OneWay_HC1b)

library(coda)

# 把你的MCMC结果拿出来
mcmc_list <- as.mcmc.list(OneWay_HC1b)

# （可选）重命名参数
varnames(mcmc_list) <- gsub("theta_0\\[1\\]", "hc[1]", varnames(mcmc_list))
varnames(mcmc_list) <- gsub("theta_0\\[2\\]", "hc[2]", varnames(mcmc_list))
varnames(mcmc_list) <- gsub("sigma_0\\[1\\]", "sigma0[1]", varnames(mcmc_list))
varnames(mcmc_list) <- gsub("sigma_0\\[2\\]", "sigma0[2]", varnames(mcmc_list))
# deviance 通常名字不会变

# 只提取需要的参数
mcmc_sel <- mcmc_list[, c("deviance", "sigma0[1]", "sigma0[2]", "hc[1]", "hc[2]")]

# 画 traceplot 和 density plot
par(mfrow=c(5,2))  # 每个参数2个图（trace + density），一列10格
for (param in varnames(mcmc_sel)) {
  traceplot(mcmc_sel[,param], main=paste("Traceplot:", param), smooth=FALSE)
  densplot(mcmc_sel[,param], main=paste("Density:", param))
}
par(mfrow=c(1,1))  # 画完恢复默认













































































































