library(nimble)
library(MCMCvis)
library(parallel)
library(expm)

# import data
counts <- readRDS('data/model_data/counts.rds')
ncap <- readRDS('data/model_data/n_cap.rds')
index_frac <- readRDS('data/model_data/index_frac.rds')
index <- readRDS('data/model_data/index.rds')
totalt <- readRDS('data/model_data/totalt.rds')
totalo <- readRDS('data/model_data/totalo.rds')
soak_days <- readRDS('data/model_data/soak_days.rds')
m_index <- readRDS('data/model_data/m_index.rds')
f_index <- readRDS('data/model_data/f_index.rds')
s_index <- readRDS('data/model_data/s_index.rds')
recruit_intro <- readRDS('data/model_data/recruit_intro.rds')
y <- readRDS('data/model_data/y.rds')
b <- readRDS('data/model_data/b.rds')

# get mark-recapture data
n1_mc <- readRDS('data/model_data/n1_mc.rds')
m2_mc <- readRDS('data/model_data/m2_mc.rds')
totalo_mc <- readRDS('data/model_data/totalo_mc.rds')
soak_days_mc <- readRDS('data/model_data/soak_days_mc.rds')
f_index_mc <- readRDS('data/model_data/f_index_mc.rds')
s_index_mc <- readRDS('data/model_data/s_index_mc.rds')
m_index_mc <- readRDS('data/model_data/m_index_mc.rds')
mc_index <- readRDS('data/model_data/mc_index.rds')

# Write model code 
code_inter_annual <- nimbleCode({
  
  ## Prior distributions
  
  # size selectivity params.
  trapm_pmax ~ dunif(0,0.1) # minnow max. probability of capture
  trapm_xmax ~ dunif(35,60) # minnow max. size of capture
  trapm_sigma ~ dunif(3,10) # minnow sigma of gaussian size selectivity curve
  trapf_pmax ~ dunif(0,0.1) # fukui max. probability of capture
  trapf_k ~ dunif(0.1,1.5) # fukui k of logistic size selectivity curve
  trapf_midpoint ~ dunif(30,100) # fukui midpoint of logistic size selectivity curve
  traps_pmax ~ dunif(0,0.1) # shrimp max. probability of capture
  traps_k ~ dunif(0.1,1.5) # shrimp k of logistic size selectivity curve
  traps_midpoint ~ dunif(30,100) # shrimp midpoint of logistic size selectivity curve
  
  # IPM - natural mortality
  nmort_shape_s ~ dunif(0,1000) # gamma distributions shape for instantaneous probability of mortality (natural)
  nmort_rate_s ~ dunif(0,1000) # gamma distributions rate for instantaneous probability of mortality (natural)
  nmort_size ~ dunif(0,10000) # size-dependent overwinter natural mortality, shared across all sites
  wnmort_shape_s ~ dunif(0,1000) # gamma distributions shape for instantaneous probability of overwinter mortality (natural)
  wnmort_rate_s ~ dunif(0,1000) # gamma distributions rate for instantaneous probability of overwinter mortality (natural)
  
  # IPM - growth
  growth_xinf ~ dnorm(95.4,sd=5.52) # asymptotic size -- informative
  growth_k ~ dnorm(0.67,sd=0.1)  # growth rate -- informative
  growth_C ~ dnorm(0.79,sd=0.11)  # amplitude of growth oscillations -- informative
  growth_ts ~ dnorm(-0.642,sd=0.027)  # inflection point of growth oscillations -- informative
  growth_sd ~ dunif(0.01,4) # growth error
  
  # expected density at initial population size
  init_lmean_adult ~ dunif(3.25,4.5)
  init_mean_recruit ~ dunif(1, 25)
  
  init_sd_r ~ dunif(0.01,20)
  init_lsd_adult ~ dunif(0.1,1)
  
  
  N_mu_recruit ~ dunif(1,1000000)
  N_sd_recruit ~ dunif(0,10000)
  
  N_adult ~ dunif(1,1000000)
  for(m in 1:nsite){
    N_recruit[m] ~ T(dnorm(N_mu_recruit, sd = N_sd_recruit),0,Inf)
  }
  
  # mark recapture data
  hazard_mc[1:totalo_mc,1:nsize] <- calc_hazard(totalo_mc,nsize,
                                                trapf_pmax,trapf_k,trapf_midpoint,
                                                traps_pmax,traps_k,traps_midpoint,
                                                trapm_pmax,trapm_xmax,trapm_sigma,
                                                f_index_mc[1:totalo_mc],
                                                s_index_mc[1:totalo_mc],
                                                m_index_mc[1:totalo_mc],
                                                soak_days_mc[1:totalo_mc],
                                                y[1:nsize])
  
  # p(capture at all)
  p_mc[1:nsize] <- calc_prob(totalo_mc,nsize,
                             hazard_mc[1:totalo_mc,1:nsize])
  
  # project n1
  n1_project[1:nsize] <- gkernel(growth_xinf,growth_k,
                                 growth_sd,growth_C,growth_ts,
                                 mc_index[1],mc_index[2],
                                 nsize,pi,y[1:nsize],
                                 lower[1:nsize],upper[1:nsize]) %*% n1[1:nsize]
  
  for(k in 1:nsize){
    # n2[k] ~ dbinom(size = round(N_mc[k]), prob = p_mc[k])
    m2[k] ~ dbinom(size = round(n1_project[k]*exp(-(nmortality_s[5]+nmort_size/y[k]^2))), 
                   prob = p_mc[k])
  }
  
  # calculate proportions in each size class -- recruits
  prop_recruit[1:nsize] <- init_sizes_gamma(init_mean_recruit,
                                            init_sd_r,lower[1:nsize],upper[1:nsize])
  
  # get initial sizes of recruits
  for(i in 1:nsite){
    # initial pop. size in each size class, k (just recruits)
    N_init_recruit[i,1:nsize] <- prop_recruit[1:nsize]*N_recruit[i]
    
  }
  
  # get initial size of first year of adults
  # calculate proportions in each size class of adults - site 1
  prop_adult[1:nsize] <- init_sizes_lnorm(init_lmean_adult,init_lsd_adult,
                                          lower[1:nsize],upper[1:nsize])
  
  # initial pop. size in each size class, k - site 1
  N_init[1:nsize] <- prop_adult[1:nsize]*N_adult
  
  # project to first observed time period - site 1
  N[1,1,1:nsize] <- gkernel(growth_xinf,growth_k,
                            growth_sd,growth_C,growth_ts,
                            0,index_frac[1,1],
                            nsize,pi,y[1:nsize],
                            lower[1:nsize],upper[1:nsize]) %*% N_init[1:nsize]
  
  # project to first observed time period - site 2:4
  for(i in 2:nsite){
    # project to first observed time period
    N[1,i,1:nsize] <- gkernel(growth_xinf,growth_k,
                              growth_sd,growth_C,growth_ts,
                              0,index_frac[1,i],
                              nsize,pi,y[1:nsize],
                              lower[1:nsize],upper[1:nsize]) %*% N_overwinter[i-1,1:nsize]
  }
  
  # get observations for t = 1
  for(i in 1:nsite){
    
    # calculate hazard rate
    hazard[1,1:totalo[1,i],i,1:nsize] <- calc_hazard(totalo[1,i],nsize,
                                                     trapf_pmax,trapf_k,trapf_midpoint,
                                                     traps_pmax,traps_k,traps_midpoint,
                                                     trapm_pmax,trapm_xmax,trapm_sigma,
                                                     f_index[1,1:totalo[1,i],i],
                                                     s_index[1,1:totalo[1,i],i],
                                                     m_index[1,1:totalo[1,i],i],
                                                     soak_days[1,1:totalo[1,i],i],
                                                     y[1:nsize])
    
    # Conditional multinomial cell probabilities
    pic_trap[1,1:totalo[1,i],i,1:nsize] <- calc_cond_prob(totalo[1,i],nsize,
                                                          hazard[1,1:totalo[1,i],i,1:nsize])
    # p(capture at all)
    p[1,i,1:nsize] <- calc_prob(totalo[1,i],nsize,
                                hazard[1,1:totalo[1,i],i,1:nsize])
    
    for(k in 1:nsize){
      
      # multinomial distribution -- trap level, observed sample size
      counts[1,1:totalo[1,i],i,k] ~ dmulti(pic_trap[1,1:totalo[1,i],i,k], n_cap[1,i,k])
      
      # binomial distribution -- observed sample size
      n_cap[1,i,k] ~ dbinom(size = round(N[1,i,k]), prob = p[1,i,k])
      
    }
    
  }
  
  # times t = 2+
  for(i in 1:nsite){
    for(t in 2:totalt[i]){
      
      # calculate hazard rate
      hazard[t,1:totalo[t,i],i,1:nsize] <- calc_hazard(totalo[t,i],nsize,
                                                       trapf_pmax,trapf_k,trapf_midpoint,
                                                       traps_pmax,traps_k,traps_midpoint,
                                                       trapm_pmax,trapm_xmax,trapm_sigma,
                                                       f_index[t,1:totalo[t,i],i],
                                                       s_index[t,1:totalo[t,i],i],
                                                       m_index[t,1:totalo[t,i],i],
                                                       soak_days[t,1:totalo[t,i],i],
                                                       y[1:nsize])
      
      
      # Conditional multinomial cell probabilities
      pic_trap[t,1:totalo[t,i],i,1:nsize] <- calc_cond_prob(totalo[t,i],nsize,
                                                            hazard[t,1:totalo[t,i],i,1:nsize])
      # p(capture at all)
      p[t,i,1:nsize] <- calc_prob(totalo[t,i],nsize,
                                  hazard[t,1:totalo[t,i],i,1:nsize])
      
      for(k in 1:nsize){
        
        # multinomial distribution -- trap level, observed sample size
        counts[t,1:totalo[t,i],i,k] ~ dmulti(pic_trap[t,1:totalo[t,i],i,k], n_cap[t,i,k])
        
        # binomial distribution -- observed sample size
        n_cap[t,i,k] ~ dbinom(size = round(N[t,i,k]), prob = p[t,i,k])
        
        # multiply next_N by nmortality (with process error)
        N[t,i,k] <- next_N[t-1,i,k]*exp(-(nmortality_s[i]+nmort_size/y[k]^2)) + recruit_intro[t,i]*N_init_recruit[i,k]
        
      }
      
      
      # project to next time period
      next_N[t-1,i,1:nsize] <- gkernel(growth_xinf,growth_k,
                                       growth_sd,growth_C,growth_ts,
                                       index_frac[t-1,i],index_frac[t,i],
                                       nsize,pi,y[1:nsize],
                                       lower[1:nsize],upper[1:nsize]) %*% (N[t-1,i,1:nsize] - n_cap[t-1,i,1:nsize])
    }
    
  }
  
  # project all sites to the end of size-independent natural mortality
  for(i in 1:nsite_short){
    for(t in additionalt[site_short[i],1]:additionalt[site_short[i],2]){
      
      # project to next time period
      next_N[t,site_short[i],1:nsize] <- gkernel(growth_xinf,growth_k,
                                                 growth_sd,growth_C,growth_ts,
                                                 index_frac[t,site_short[i]],
                                                 index_frac[t+1,site_short[i]],
                                                 nsize,pi,y[1:nsize],
                                                 lower[1:nsize],
                                                 upper[1:nsize]) %*% N[t,site_short[i],1:nsize]

      for(k in 1:nsize){
        # multiply next_N by nmortality (with process error)
        N[t+1,site_short[i],k] <- next_N[t,site_short[i],k]*exp(-(nmortality_s[site_short[i]]+nmort_size/y[k]^2))
      }
    }
  }
  
  # gamma distributed logmortality process error
  for(m in 1:(nsite+1)){
    nmortality_s[m] ~ dgamma(nmort_shape_s, nmort_rate_s)  
  }
  for(m in 1:(nsite-1)){
    wnmortality_s[m] ~ dgamma(wnmort_shape_s, wnmort_rate_s)
  }
  
  # project to new year with size-dependent natural mortality
  for(i in 1:(nsite-1)){
    wgrowth_N[i,1:nsize] <- gkernel(growth_xinf,growth_k,
                                    growth_sd,growth_C,growth_ts,
                                    max_index_frac,1,
                                    nsize,pi,y[1:nsize],
                                    lower[1:nsize],
                                    upper[1:nsize]) %*% N[additionalt[i,2]+1,i,1:nsize]
    
    for(k in 1:nsize){
      N_overwinter[i,k] ~ dbinom(size = round(wgrowth_N[i,k]), 
                                 prob = exp(-(wnmortality_s[i]*sum(wgrowth_N[i,1:nsize])/y[k]^2)))
    }
  }
  
  
})


# data structure for additional time for size-independent mortality
additionalt <- cbind(as.vector(totalt),
                     totalt+max(index,
                                na.rm=TRUE) - apply(index, 
                                                    2, max, na.rm = TRUE)-1)

# bundle up data and constants
constants <- list(nsite = dim(counts)[3], # number of sites
                  totalt = totalt, # number of time periods for each site
                  totalo = totalo, # trap obs, site i, time t
                  additionalt = additionalt,
                  site_short = which(apply(index, 2, 
                                           max, na.rm = TRUE) < max(index,
                                                                    na.rm=TRUE)),
                  nsite_short = length(which(apply(index, 
                                                   2, max, 
                                                   na.rm = TRUE) < max(index,
                                                                       na.rm=TRUE))),
                  nsize = length(y), # number of sizes in IPM mesh
                  upper = b[2:length(b)], # upper bounds in IPM mesh
                  lower = b[1:length(y)], # lower bounds in IPM mesh
                  f_index = f_index, # binary indicator of fukui traps at site, i, time, t, obs, j
                  m_index = m_index, # binary indicator of minnow traps at site, i, time, t, obs, j
                  s_index = s_index, # binary indicator of shrimp traps at site, i, time, t, obs, j
                  y = y, # midpoints in IPM mesh
                  index_frac = index_frac,
                  max_index_frac = max(index_frac,na.rm=TRUE),
                  soak_days = soak_days,
                  recruit_intro = recruit_intro,
                  totalo_mc = totalo_mc,
                  f_index_mc = f_index_mc,
                  s_index_mc = s_index_mc,
                  m_index_mc = m_index_mc,
                  soak_days_mc = soak_days_mc,
                  mc_index = mc_index,
                  pi = pi
)

data <- list(n_cap = ncap, # total captured within each time point
             counts = counts, # crab counts
             n1 = n1_mc,
             m2 = m2_mc
) 


# create some initial values
# log nmortality
nmortality_init <- rbeta(dim(counts)[3]+1,0.1,1)
wnmortality_init <- rbeta(dim(counts)[3]-1,2,75)
# N
N_overwinter <- matrix(NA,nrow=dim(counts)[3]-1,ncol=length(y))
mean_adult <- c(log(75),log(78.5),log(82))
mean_recruit <- c(log(55),log(55),log(52))
N_adult <- c(500,400,100)
N_recruit <- c(300,75,400)
adult_sd <- 0.08
recruit_sd <- 0.1
for(i in 1:(dim(counts)[3]-1)){
  prob_r <- dlnorm(y,mean_recruit[i],recruit_sd)/sum(dlnorm(y,mean_recruit[i],recruit_sd))
  prob_a <- dlnorm(y,mean_adult[i],adult_sd)/sum(dlnorm(y,mean_adult[i],adult_sd))
  N_overwinter[i,] <- round(prob_a*N_adult[i]+prob_r*N_recruit[i])
}

# initial values
inits <- function(){
  list(trapm_pmax = runif(1,0.0001,0.0008),
       trapm_xmax = 44,
       trapm_sigma = 6.67,
       trapf_pmax = runif(1,0.0001,0.0008),
       trapf_k = 0.2,
       trapf_midpoint = 45,
       traps_pmax = runif(1,0.001,0.005),
       traps_k = 0.2,
       traps_midpoint = 45,
       nmort_shape_s = 2,
       nmort_rate_s = 1,
       nmort_size = 0.1,
       wnmort_shape_s = 2,
       wnmort_rate_s = 1,
       growth_k = 1,
       growth_xinf = 85,
       growth_C = 0.79,
       growth_ts = -0.64,
       growth_sd = 2.5,
       init_sd_r = 1,
       init_mean_recruit = 18,
       init_lmean_adult = 4,
       init_lsd_adult = 0.2,
       N_adult = 1800,
       N_recruit = c(1000,100,1000,100),
       N_mu_recruit = 500,
       N_sd_recruit = 100,
       nmortality_s = rep(0.001,5),
       wnmortality_s = c(0.0308, 0.00669, 0.0346),
      # wnmortality_s = wnmortality_init,
       N_overwinter = N_overwinter
  )
}



#################
#parallelization#
#################

cl <- makeCluster(4)

set.seed(10120)

clusterExport(cl, c("code_inter_annual", "inits", "data", "constants"))
for (j in seq_along(cl)) {
  set.seed(j)
  init <- replicate(15, inits(), simplify = FALSE)
  clusterExport(cl[j], "init")
}

# Create a function with all the needed code
out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  library(expm)
  
  #write function for size selectivity of trap type m
  # gaussian function (drawn from Jorgensen et al. 2009)
  size_sel_norm <- nimbleFunction(
    #input and output types
    run = function(pmax = double(0), #scalar
                   xmax = double(0), # scalar
                   sigma = double(0), # scalar
                   y = double(1)) #vector 
    {
      returnType(double(1))  # vector return
      
      # Conditional multinomial cell probabilities
      vector <- pmax * exp(-(y-xmax)^2/(2*sigma^2))
      
      return(vector) 
    }#,
    #buildDerivs = TRUE
  )
  assign('size_sel_norm', size_sel_norm, envir = .GlobalEnv)
  
  #write function for size selectivity of trap type f and s
  # logistic function
  size_sel_log <- nimbleFunction(
    #input and output types
    run = function(pmax = double(0), #scalar
                   k = double(0), # scalar
                   midpoint = double(0), # scalar
                   y = double(1)) #vector 
    {
      returnType(double(1))  # vector return
      
      # Conditional multinomial cell probabilities
      
      vector <- pmax/(1+exp(-k*(y-midpoint)))
      
      return(vector)
    }#,
    #buildDerivs = TRUE
  )
  assign('size_sel_log', size_sel_log, envir = .GlobalEnv)
  
  
  #write function for calculating trap hazard rates
  calc_hazard <- nimbleFunction(
    #input and output types
    run = function(nobs = double(0), # scalar
                   nsize = double(0), # scalar
                   #maxo = double(0), # scalar
                   trapf_pmax = double(0), # scalar
                   trapf_k = double(0), # scalar
                   trapf_midpoint = double(0), # scalar 
                   traps_pmax = double(0), # scalar
                   traps_k = double(0), # scalar
                   traps_midpoint = double(0), # scalar 
                   trapm_pmax = double(0), # scalar
                   trapm_xmax = double(0), # scalar
                   trapm_sigma = double(0), # scalar 
                   obs_ref_f = double(1),  # 1-d array
                   obs_ref_s = double(1), # 1-d array
                   obs_ref_m = double(1), # 1-d array
                   soak_days = double(1), # 1-d array
                   y = double(1)) # vector 
    {
      returnType(double(2))  # 2-d array return
      
      #create empty array
      array <- array(init=FALSE,dim=c(nobs,nsize))
      #loop through sites, visits, observations
      for(j in 1:nobs){
        #P(capture)
        array[j,1:nsize] <- size_sel_log(pmax=trapf_pmax,k=trapf_k,midpoint=trapf_midpoint,y=y)*obs_ref_f[j]*soak_days[j] +
          size_sel_log(pmax=traps_pmax,k=traps_k,midpoint=traps_midpoint,y=y)*obs_ref_s[j]*soak_days[j] +
          size_sel_norm(pmax=trapm_pmax,xmax=trapm_xmax,sigma=trapm_sigma,y=y)*obs_ref_m[j]*soak_days[j] 
      }
      
      return(array)
    }#,
    #buildDerivs = TRUE
  )
  assign('calc_hazard', calc_hazard, envir = .GlobalEnv)
  
  #write function for calculating conditional probabilities
  calc_cond_prob <- nimbleFunction(
    #input and output types
    run = function(nobs = double(0), # scalar
                   nsize = double(0), #scalar
                   #maxo = double(0), # scalar
                   hazard = double(2)) # 2-d array 
    {
      returnType(double(2))  # 2-d array return
      
      # Conditional multinomial cell probabilities
      #create empty array
      array <- array(init=FALSE,dim=c(nobs,nsize))
      #loop through sizes
      for(k in 1:nsize){
        #P(capture in trap j | captured at all)
        array[1:nobs,k] <- hazard[1:nobs,k] / sum(hazard[1:nobs,k])
      }
      # fill empty space with zeros
      # if(maxo > nobs){
      #   array[(nobs+1):maxo,1:nsize] <- 0
      # }
      return(array)
    }#,
    #buildDerivs = TRUE
  )
  assign('calc_cond_prob', calc_cond_prob, envir = .GlobalEnv)
  
  
  #write function for calculating capture probabilities
  calc_prob <- nimbleFunction(
    #input and output types
    run = function(nobs = double(0), # scalar
                   nsize = double(0), #scalar
                   hazard = double(2)) # 2-d array 
    {
      returnType(double(1))  # 1-d array return
      
      #create empty array
      p <- rep(NA,nsize)
      #loop through sizes
      for(k in 1:nsize){
        #1-exp(-P(captured at all))
        p[k] <- 1-exp(-sum(hazard[1:nobs,k]))
      }
      return(p)
    }#,
    #buildDerivs = TRUE
  )
  assign('calc_prob', calc_prob, envir = .GlobalEnv)
  
  #write function for growth kernel -- biweekly
  gkernel <- nimbleFunction(
    #input and output types
    run = function(growth_xinf = double(0), #scalar
                   growth_k = double(0), #scalar
                   growth_sd = double(0), #scalar
                   growth_C = double(0), # scalar
                   growth_ts = double(0), # scalar
                   t1 = double(0), # scalar
                   t2 = double(0), # scalar
                   nsize = double(0), #scalar
                   pi = double(0), #scalar
                   y = double(1), #vector
                   lower = double(1), #vector
                   upper = double(1)) #vector 
    {
      returnType(double(2))  # array return
      
      #create empty array
      array <- matrix(NA,ncol=nsize,nrow=nsize)
      
      if(!is.na(t2) & t2 == 0){
        array <- diag(nsize)
      } else {
        # season adjusted params
        S_t <- (growth_C*growth_k/(2*pi))*sin(2*pi*(t2-(1+growth_ts)))
        S_t0 <- (growth_C*growth_k/(2*pi))*sin(2*pi*(t1-(1+growth_ts)))
        
        #p(y'|y)
        for(i in 1:nsize){
          increment <- (growth_xinf-y[i])*(1-exp(-growth_k*(t2-t1)-S_t+S_t0))
          mean <- y[i] + increment
          array[1:nsize,i] <- (pnorm(upper,mean,sd=growth_sd)-pnorm(lower,mean,sd=growth_sd)) 
        }
        
        # normalize
        for(i in 1:nsize){
          array[,i] <- array[,i]/sum(array[,i])
        }
      }
      
      return(array)
    }
  )
  assign('gkernel', gkernel, envir = .GlobalEnv)
  
  #write function for initial sizes -- gamma
  init_sizes_gamma <- nimbleFunction(
    #input and output types
    run = function(init_mean = double(0),
                   init_sd = double(0),
                   lower = double(1), #vector
                   upper = double(1)) #vector 
    {
      returnType(double(1))  # vector return
      
      var <- init_sd^2
      shape <- init_mean^2/var
      rate <- init_mean/var
      out <- pgamma(q=upper,shape=shape,rate=rate)-pgamma(q=lower,shape=shape,rate=rate)
      
      return(out)
    }
  )
  assign('init_sizes_gamma', init_sizes_gamma, envir = .GlobalEnv)
  
  #write function for initial sizes -- lognormal
  init_sizes_lnorm <- nimbleFunction(
    #input and output types
    run = function(init_lmean = double(0),
                   init_lsd = double(0),
                   lower = double(1), #vector
                   upper = double(1)) #vector 
    {
      returnType(double(1))  # vector return
      
      out <- plnorm(q=upper,meanlog=init_lmean,sdlog=init_lsd)-plnorm(q=lower,meanlog=init_lmean,sdlog=init_lsd)
      
      return(out)
    }
  )
  assign('init_sizes_lnorm', init_sizes_lnorm, envir = .GlobalEnv)
  
  loglik <- Inf
  index <- 1
  while(is.infinite(loglik)){
    myModel <- nimbleModel(code = code_inter_annual,
                           data = data,
                           constants = constants,
                           inits = init[[index]])
      loglik <- myModel$calculate()
    index <- index + 1
  }  
  
  
  #build the MCMC
  mcmcConf_myModel <- configureMCMC(myModel,
                                    monitors=c("trapm_pmax","trapm_xmax","trapm_sigma",
                                               "trapf_pmax","trapf_k","trapf_midpoint",
                                               "traps_pmax","traps_k","traps_midpoint",
                                               "nmort_shape_s","nmort_rate_s",
                                               "growth_k", "growth_xinf","growth_C",
                                               "growth_ts","growth_sd","init_sd_r",
                                               "init_mean_recruit","init_lmean_adult",
                                               "init_lsd_adult","N_mu_recruit","N_sd_recruit",
                                               "N_recruit","N_adult","nmortality_s","wnmortality_s",
                                               "nmort_size","nmort_shape_s","nmort_rate_s",
                                               "wnmort_shape_s","wnmort_rate_s",
                                               "N_overwinter"
                                    ),
                                    useConjugacy=FALSE,
                                    enableWAIC=TRUE)
  
  # add block sampler for nmort params
  mcmcConf_myModel$removeSamplers(c('nmort_shape_s','nmort_rate_s'))	
  mcmcConf_myModel$addSampler(c('nmort_shape_s','nmort_rate_s'),
                                type='RW_block')
  mcmcConf_myModel$removeSamplers(c('wnmort_shape_s','wnmort_rate_s'))	
  mcmcConf_myModel$addSampler(c('wnmort_shape_s','wnmort_rate_s'),
                              type='RW_block')
  
  #build MCMC
  myMCMC <- buildMCMC(mcmcConf_myModel)
  
  #compile the model and MCMC
  CmyModel <- compileNimble(myModel)
  #compile the MCMC
  cmodel_mcmc <- compileNimble(myMCMC,project=myModel)
  
  cmodel_mcmc$run(100000,thin=10,
                  reset=FALSE)
  
  return(as.mcmc(as.matrix(cmodel_mcmc$mvSamples)))
})
# save samples
saveRDS(out,'20241216_savedsamples.rds')

stopCluster(cl)



