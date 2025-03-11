library(tidyverse)
library(bayestestR)
library(patchwork)
library(cowplot)


##################################
# simulate final population size #
##################################

# read in samples
out <- readRDS("data/posterior_samples/savedsamples_IPM.rds")

# subset samples
lower <- 2000
upper <- 10001
out_sub <- list(
  out[[1]][lower:upper, ], out[[2]][lower:upper, ],
  out[[3]][lower:upper, ], out[[4]][lower:upper, ]
)

# function for getting sample
get_params <- function(samples, index) {
  samples <- rbind(out_sub[[1]], out_sub[[2]],
                   out_sub[[3]], out_sub[[4]])

  # set parameters
  params <- list(
    h_M_max = as.numeric(samples[index, "h_M_max"]),
    h_M_A = as.numeric(samples[index, "h_M_A"]),
    h_M_sigma = as.numeric(samples[index, "h_M_sigma"]),
    h_F_max = as.numeric(samples[index, "h_F_max"]),
    h_F_k = as.numeric(samples[index, "h_F_k"]),
    h_F_0 = as.numeric(samples[index, "h_F_0"]),
    h_S_max = as.numeric(samples[index, "h_S_max"]),
    h_S_k = as.numeric(samples[index, "h_S_k"]),
    h_S_0 = as.numeric(samples[index, "h_S_0"]),
    beta_alpha = as.numeric(samples[index, "beta_alpha"]),
    beta_theta = as.numeric(samples[index, "beta_theta"]),
    alpha = as.numeric(samples[index, "alpha"]),
    alpha_o_alpha = as.numeric(samples[index, "alpha_o_alpha"]),
    alpha_o_theta = as.numeric(samples[index, "alpha_o_theta"]),
    gk = as.numeric(samples[index, "gk"]),
    yinf = as.numeric(samples[index, "yinf"]),
    C = as.numeric(samples[index, "C"]),
    ts = as.numeric(samples[index, "ts"]),
    sigma_G = as.numeric(samples[index, "sigma_G"]),
    sigma_R_y = as.numeric(samples[index, "sigma_R_y"]),
    mu_R_y = as.numeric(samples[index, "mu_R_y"]),
    mu_A_y = as.numeric(samples[index, "mu_A_y"]),
    sigma_A_y = as.numeric(samples[index, "sigma_A_y"]),
    lambda_A = as.numeric(samples[index, "lambda_A"]),
    mu_R_lambda = as.numeric(samples[index, "mu_R_lambda"]),
    sigma_R_lambda = as.numeric(samples[index, "sigma_R_lambda"])
  )

  return(params)
}

# create parameters for IPM mesh
min_size <- 0
max_size <- 110
n <- 22 # number of cells in the discretized kernel
b <- min_size + c(0:n) * (max_size - min_size) / n # boundary points
y <- 0.5 * (b[1:n] + b[2:(n + 1)]) # mesh points (midpoints of the cells)
h <- y[2] - y[1] # width of the cells
recruit_index <- 6

# set constant values
index_frac <- readRDS("data/model_data/index_frac.rds")
t <- index_frac[, 2]

constant <- list(
  n_size = n,
  upper = b[2:length(b)],
  lower = b[1:length(y)],
  t_index = t,
  y = y,
  ntime = length(t),
  recruit_intro = recruit_index,
  n_years = 50,
  n_years_burn = 5
)

####################
# helper functions #
####################

# size selective hazard rate of trap type minnow - bell-shaped curve
size_sel_norm <- function(pmax, xmax, sigma, y) {

  vector <- pmax * exp(-(y - xmax) ^ 2 / (2 * sigma ^ 2))

  return(vector)
}

# size selective hazard rate of trap type fukui and shrimp - logistic
size_sel_log <- function(pmax, k, midpoint, y) {
  
  vector <- pmax / (1 + exp(-k * (y - midpoint)))

  return(vector)
}

# calculate trap hazard rate of obs j, based on trap type and soak days
calc_hazard <- function(params, constant, nobs,
                        obs_ref_f, obs_ref_s, obs_ref_m) {
  
  # create empty array
  array <- array(NA, dim = c(nobs, constant$n_size))
  # loop through observations
  for (j in 1:nobs) {
    # P(capture)
    array[j, ] <- size_sel_log(pmax = params$h_F_max,
                               k = params$h_F_k,
                               midpoint = params$h_F_0,
                               y = constant$y) * obs_ref_f[j] +
      size_sel_log(pmax = params$h_S_max,
                   k = params$h_S_k,
                   midpoint = params$h_S_0,
                   y = constant$y) * obs_ref_s[j] +
      size_sel_norm(pmax = params$h_M_max,
                    xmax = params$h_M_A,
                    sigma = params$h_M_sigma,
                    y = constant$y) * obs_ref_m[j]
  }
  return(array)

}

# calculate conditional probability of capture
calc_cond_prob <- function(nobs, n_size, hazard) {

  # create empty array
  array <- array(NA, dim = c(nobs, n_size))

  # loop through sizes
  for (k in 1:n_size) {
    # P(capture in trap j | captured at all)
    array[, k] <- hazard[, k] / sum(hazard[, k])
  }
  return(array)
}

# calculate total capture probability
calc_prob <- function(nobs, n_size, hazard) {

  # create empty array
  p <- rep(NA, n_size)

  # loop through sizes
  for (k in 1:n_size) {
    # 1 - exp(-P(captured at all))
    p[k] <- 1 - exp(-sum(hazard[1:nobs, k]))
  }
  return(p)
}

# function for seasonal growth and mortality kernel -- biweekly time step
get_kernel <- function(params, t1, t2,
                       n_size, y, lower, upper, S) {
  #create empty array
  array <- matrix(NA, ncol = n_size, nrow = n_size)

  # season adjusted params
  S_t <- (params$C * params$gk / (2 * pi)) *
    sin(2 * pi * (t2 - (1 + params$ts)))
  S_t0 <- (params$C * params$gk / (2 * pi)) *
    sin(2 * pi * (t1 - (1 + params$ts)))

  #p(y"|y)
  for (i in 1:n_size) {
    increment <- (params$yinf - y[i]) *
      (1 - exp(-params$gk * (t2 - t1) - S_t + S_t0))
    mean <- y[i] + increment
    array[1:n_size, i] <- (pnorm(upper, mean,
                                 sd = params$sigma_G) -
                             pnorm(lower, mean, sd = params$sigma_G))
  }

  # normalize and apply natural survival
  for (i in 1:n_size) {
    array[, i] <- array[, i] / sum(array[, i]) * S
  }

  return(array)

}

# natural (non-winter) survival
survival <- function(params, beta, y, t1, t2) {
  
  # number of biweeks
  deltat <- round((t2 - t1) * (52.1429 / 2))
    
  # get survival rate
  out <- exp(-deltat * (beta + params$alpha / y ^ 2))
    
  return(out)
}

# density- and size-dependent overwinter survival
overwinter_survival <- function(alpha_o, N_sum, y) {
  
  # get probability of survival
  out <- exp(-(alpha_o * N_sum / y ^ 2))
    
  return(out)
}


# gamma distribution for initial size distribution of recruits
init_sizes_gamma <- function(params, lower, upper) {
  
  var <- params$sigma_R_y ^ 2
  shape <- params$mu_R_y ^ 2 / var
  rate <- params$mu_R_y / var
  
  out <- pgamma(q = upper, shape = shape, rate = rate) -
    pgamma(q = lower, shape = shape, rate = rate)

  return(out)
}

# lognormal distribution for initial size distribution of recruits
init_sizes_lnorm <- function(params, lower, upper) {
  
  out <- plnorm(q = upper, meanlog = params$mu_A_y,
                sdlog = params$sigma_A_y) -
    plnorm(q = lower, meanlog = params$mu_A_y,
           sdlog = params$sigma_A_y)

  return(out)
}

# function to simulate one year
simulate_year <- function(params, constant, nobs, N, N_recruit, beta,
                          alpha_o, year, n_years,
                          obs_ref_f, obs_ref_s, obs_ref_m) {

  # check N_recruit
  N_recruit <- ifelse(as.numeric(N_recruit) < 0, 0, as.numeric(N_recruit))

  # initialize recruits
  init_recruit <- init_sizes_gamma(params, constant$lower,
                                   constant$upper) * N_recruit

  # calculate hazard rate
  hazard <- calc_hazard(params, constant, nobs,
                        obs_ref_f, obs_ref_s, obs_ref_m)
  # p(capture at all)
  p <- calc_prob(nobs, constant$n_size, hazard)

  # create vector for captured
  ncap <- rep(NA, constant$n_size)
  for (t in 1:(length(constant$t_index) - 1)) {


    # remove crabs
    for (k in 1:constant$n_size) {

      # number captured
      ncap[k] <- rbinom(1, size = ifelse(round(as.numeric(N[k])) > 0,
                                         round(as.numeric(N[k])), 0),
                        prob = p[k])
    }

    # get kernel
    S <- survival(params, beta, constant$y, 0, 14 / 365)
    kernel <- get_kernel(params, constant$t_index[t], constant$t_index[t + 1],
                         constant$n_size, constant$y,
                         constant$lower, constant$upper, S)

    # project to next time period
    if (t == constant$recruit_intro) {
      N <- kernel %*% (
        as.numeric(N) - ncap
      ) + init_recruit
    } else {
      N <- kernel %*% (
        as.numeric(N)  - ncap
      )
    }

  }

  # remove crabs in last time point
  for (k in 1:constant$n_size) {

    # number captured
    ncap[k] <- rbinom(1, size = ifelse(round(N[k]) > 0, round(N[k]), 0),
                      prob = p[k])
  }
  N <- N - ncap

  # overwinter growth
  wgrowth_N <- get_kernel(params, max(constant$t_index), 1,
                          constant$n_size, constant$y,
                          constant$lower, constant$upper, 
                          rep(1, length(constant$y))) %*% N

  N_overwinter <- rep(NA, constant$n_size)

  wgrowth_N_sum <- ifelse(round(sum(wgrowth_N)) > 0,
                          round(sum(wgrowth_N)), 0)

  # get overwinter survival
  S_o <- overwinter_survival(alpha_o, wgrowth_N_sum, constant$y)
  for (k in 1:constant$n_size) {
    new_N <- ifelse(round(wgrowth_N[k]) > 0,
                    round(wgrowth_N[k]), 0)
    N_overwinter[k] <- rbinom(n = 1, size = new_N, prob = S_o[k])
  }
  
  return(N_overwinter)

}

# function to simulate one year - no effort
simulate_year_noeffort <- function(params, constant, N, N_recruit, beta,
                                   alpha_o, year, n_years) {

  # check N_recruit
  N_recruit <- ifelse(as.numeric(N_recruit) < 0, 0, as.numeric(N_recruit))

  # initialize recruits
  init_recruit <- init_sizes_gamma(params, constant$lower,
                                   constant$upper) * N_recruit

  for (t in 1:(length(constant$t_index) - 1)) {
    
    # get kernel
    S <- survival(params, beta, constant$y, 0, 14 / 365)
    kernel <- get_kernel(params, constant$t_index[t], constant$t_index[t + 1],
                         constant$n_size, constant$y,
                         constant$lower, constant$upper, S)
    
    # project to next time period
    if (t == constant$recruit_intro) {
      N <- kernel %*% as.numeric(N) + init_recruit
    } else {
      N <- kernel %*% as.numeric(N)
    }
  }

  # overwinter growth
  wgrowth_N <- get_kernel(params, max(constant$t_index), 1,
                          constant$n_size, constant$y,
                          constant$lower, constant$upper, 
                          rep(1, length(constant$y))) %*% N

  N_overwinter <- rep(NA, constant$n_size)

  wgrowth_N_sum <- ifelse(round(sum(wgrowth_N)) > 0,
                          round(sum(wgrowth_N)), 0)

  # get overwinter survival
  S_o <- overwinter_survival(alpha_o, wgrowth_N_sum, constant$y)
  for (k in 1:constant$n_size) {
    new_N <- ifelse(round(wgrowth_N[k]) > 0,
                    round(wgrowth_N[k]), 0)
    N_overwinter[k] <- rbinom(n = 1, size = new_N, prob = S_o[k])
  }
  
  return(N_overwinter)
}

# function to simulate multiple years
simulate_interannual <- function(params, constant, nobs,
                                 n_years, n_years_burn, N, N_recruit,
                                 beta, alpha_o,
                                 obs_ref_f, obs_ref_s, obs_ref_m) {

  out <- as.data.frame(matrix(NA, nrow = n_years, ncol = length(N)))

  for (i in 1:n_years) {
    out[i, ] <- simulate_year(params, constant, nobs, N, N_recruit[i],
                              beta[i], alpha_o[i], i, n_years,
                              obs_ref_f, obs_ref_s, obs_ref_m)
  }

  # subset to averaged years
  sub <- out[(n_years_burn + 1):n_years, ]
  return(colMeans(sub))
}
simulate_interannual_noeffort <- function(params, constant, n_years,
                                          n_years_burn, N, N_recruit,
                                          beta, alpha_o) {

  out <- as.data.frame(matrix(NA, nrow = n_years, ncol = length(N)))

  for (i in 1:n_years) {
    out[i, ] <- simulate_year_noeffort(params, constant, N, N_recruit[i],
                                       beta[i], alpha_o[i], i, n_years)
  }
  # subset to averaged years
  sub <- out[(n_years_burn + 1):n_years, ]
  return(colMeans(sub))
}


############
# simulate #
############

effort <- c(2, 8, 40, 60, 100, 200)
prop <- list(
  list(prop_s = 1,
       prop_m = 0,
       prop_f = 0),
  list(prop_s = 0,
       prop_m = 1,
       prop_f = 0),
  list(prop_s = 0,
       prop_m = 0,
       prop_f = 1)
)

out <- matrix(NA, nrow = 0, ncol = constant$n_size + 5)
colnames(out) <- c(1:constant$n_size,
                   "nobs", "prop_s", "prop_m", "prop_f", "iter")

# get indices of samples to use
set.seed(123)
index <- sample(dim(out_sub[[1]])[1] * length(out_sub), 1000)

# get starting N and N_recruit
N <- as.data.frame(matrix(nrow = length(index), ncol = length(constant$y)))
N_recruit <- as.data.frame(matrix(nrow = length(index),
                                  ncol = constant$n_years))
for (i in 1:length(index)) {
  params <- get_params(out_sub, i)
  N[i, ] <- init_sizes_lnorm(params, constant$lower,
                             constant$upper) * params$lambda_A
  N_recruit[i, ] <- rnorm(constant$n_years, params$mu_R_lambda,
                          params$sigma_R_lambda)
}

# get natural mortality params
beta <- as.data.frame(matrix(nrow = length(index), ncol = constant$n_years))
alpha_o <- as.data.frame(matrix(nrow = length(index), ncol = constant$n_years))

for (i in 1:length(index)) {
  params <- get_params(out_sub, i)
  beta[i, ] <- rgamma(constant$n_years, shape = params$beta_alpha,
                      rate = params$beta_theta)
  alpha_o[i, ] <- rgamma(constant$n_years, shape = params$alpha_o_alpha,
                         rate = params$alpha_o_theta)
}

for (e in effort) {
  for (i in 1:length(prop)) {

    # varying effort
    nobs <- e
    prop_s <- prop[[i]]$prop_s
    prop_m <- prop[[i]]$prop_m
    prop_f <- prop[[i]]$prop_f

    # effort
    if (prop_s == 1) {
      obs_ref_s <- rep(1, nobs)
      obs_ref_m <- rep(0, nobs)
      obs_ref_f <- rep(0, nobs)
    } else if (prop_m == 1) {
      obs_ref_s <- rep(0, nobs)
      obs_ref_m <- rep(1, nobs)
      obs_ref_f <- rep(0, nobs)
    } else if (prop_f == 1) {
      obs_ref_s <- rep(0, nobs)
      obs_ref_m <- rep(0, nobs)
      obs_ref_f <- rep(1, nobs)
    }

    for (k in 1:length(index)) {
      out <- rbind(out,
                   c(simulate_interannual(get_params(out_sub, index[k]),
                                          constant, nobs, constant$n_years,
                                          constant$n_years_burn,
                                          N[k, ], N_recruit[k, ],
                                          as.numeric(beta[k, ]),
                                          as.numeric(alpha_o[k, ]),
                                          obs_ref_f, obs_ref_s, obs_ref_m),
                     nobs * constant$ntime, prop_s, prop_f, prop_m, k))

    }
  }
  print(paste0("finished effort: ", e))
}
saveRDS(out, "data/simulations/simulations_effort_20250311.rds")

# simulate no effort
out_noeffort <- matrix(NA, nrow = length(index), ncol = constant$n_size)
colnames(out_noeffort) <- 1:constant$n_size
for (i in 1:length(index)) {
  out_noeffort[i, ] <- simulate_interannual_noeffort(
    get_params(out_sub, index[i]), constant, constant$n_years,
    constant$n_years_burn, N[i, ], N_recruit[i, ], as.numeric(beta[i, ]),
    as.numeric(alpha_o[i, ])
  )

}
saveRDS(out_noeffort, "data/simulations/simulations_noeffort_20250311.rds")
