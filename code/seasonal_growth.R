library(tidyverse)
library(nimble)
library(MCMCvis)
library(parallel)

# read in data
data <- readRDS("data/model_data/growth_data.rds")

# nimble code
seasonal_growth <- nimbleCode({
  # prior distributions
  k ~ dunif(0, 2) # growth rate
  A ~ dunif(0, 4) # amplitude of seasonal variation
  ts ~ dunif(-1, 0) # time between t=0 and start of growth oscillation
  t0 ~ dunif(-1, 1) # age organism has 0 size
  tau_w ~ dunif(0, 100) # process error sd
  tau_y ~ dunif(0, 100) # year random effect sd
  xinf ~ dunif(70, 140) # asymptotic size
  
  for(i in 1:nsizes){
    W[i] ~ dnorm(W_hat[i], tau_w)
    W_hat[i] <- xinf * (1 - exp(-k * (age[i] - t0) - S_t[i] + S_t0)) + 
      ranef[year[i]]
    S_t[i] <- (A * k / (2 * pi)) * sin(2 * pi * (age[i] - ts))
  }
  S_t0 <- (A * k / (2 * pi)) * sin(2 * pi * (t0 - ts))
  
  for(y in 1:nyears){
    ranef[y] ~ dnorm(0, tau_y)
  }
  
})


# bundle up data and constants
constants <- list(nsizes = length(data$CW), # number of crabs
                  year = data$year_index, # year index
                  nyears = length(unique(data$year_index)), # number of years
                  pi = pi # pi
)

model_data <- list(W = data$CW, # size captured
                   age = data$age # ages
)

# set initial values
inits_season <- list (xinf = 100,
                      k = 0.6,
                      A = 0.7,
                      ts = -0.6,
                      t0 = 0,
                      tau_w = 0.05,
                      tau_y = 0.05,
                      ranef = rep(0, length(unique(data$year_index)))
                      )

# run MCMC chains in parallel
cl <- makeCluster(4)

set.seed(10120)

clusterExport(cl, c("seasonal_growth", "inits_season", 
                    "model_data", "constants"))

# Create a function with all the needed code
out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  # build the model
  model_season <- nimbleModel(seasonal_growth, constants, model_data,
                              inits_season)
  
  # build the MCMC
  mcmcConf_season <- configureMCMC(model_season,
                                   monitors = c("xinf", "k", "A", "ts",
                                                "t0", "tau_w", "tau_y"
                                   ))
  model_mcmc_season <- buildMCMC(mcmcConf_season)
  
  # compile the model and MCMC
  # model
  cmodel_season <- compileNimble(model_season)
  # MCMC
  cmodel_mcmc_season <- compileNimble(model_mcmc_season, project = model_season)
  
  cmodel_mcmc_season$run(100000, thin = 10,
                         reset = FALSE)
  
  return(as.mcmc(as.matrix(cmodel_mcmc_season$mvSamples)))
})

stopCluster(cl)

# save samples
saveRDS(out, "data/posterior_samples/savedsamples_growth.rds")
