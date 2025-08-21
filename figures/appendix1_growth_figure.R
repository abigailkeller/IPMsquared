library(tidyverse)
library(bayestestR)

# read in data
x <- readRDS("data/model_data/x.rds")
n_size <- length(x)
b <- readRDS("data/model_data/b.rds")
upper <- b[2:length(b)]
lower <- b[1:length(x)]
data <- readRDS("data/model_data/growth_data.rds")

# read in samples
out <- readRDS("data/posterior_samples/savedsamples_IPM.rds")
# subset samples
low <- 2000
up <- 10001
out_sub <- list(
  out[[1]][low:up, ], out[[2]][low:up, ],
  out[[3]][low:up, ], out[[4]][low:up, ]
)

summary <- MCMCsummary(out_sub)

# biweek <- c(59, 76, 91, 106, 121, 137, 152, 167, 182, 198, 213, 
#             229, 244, 259, 274, 290, 305, 320, 335)
# biweek <- 59 + 0:18 * 14
# biweek_full <- c(biweek, max(biweek) + 0:6 * 14)
# index_full <- index_frac / 14 * 365 + 3
# 
# julian_days <- c(91, 182, 274)
# index <- c(3, 9, 16)
D1 <- c(0.08767123, 0.08767123 + 0.25, 0.08767123 + 0.5, 0.08767123 + 0.75)
D2 <- D1 + 14 / 365

# function for seasonal growth kernel -- biweekly time step
get_kernel <- nimbleFunction (
  # input and output types
  run = function(xinf = double(0), k = double(0),
                 sigma_G = double(0), A = double(0),
                 ds = double(0), t1 = double(0), t2 = double(0),
                 n_size = double(0), pi = double(0), x = double(1),
                 lower = double(1), upper = double(1), S = double(1))
  {
    returnType(double(2))
    
    # create empty array
    array <- matrix(NA, ncol = n_size, nrow = n_size)
    
    if(!is.na(t2) & t2 == 0) {
      array <- diag(n_size)
    } else {
      # season adjusted params
      S_t <- (A * k / (2 * pi)) * sin(2 * pi * (t2 - (1 + ds)))
      S_t0 <- (A * k / (2 * pi)) * sin(2 * pi * (t1 - (1 + ds)))
      
      # p(y"|y)
      for (i in 1:n_size) {
        increment <- (xinf - x[i]) *
          (1 - exp(-k * (t2 - t1) - S_t + S_t0))
        mean <- x[i] + increment
        array[1:n_size, i] <- (pnorm(upper, mean, sd = sigma_G) -
                                 pnorm(lower, mean, sd = sigma_G))
      }
      
      # normalize and apply natural mortality
      for (i in 1:n_size) {
        array[, i] <- array[, i] / sum(array[, i]) * S
      }
    }
    return(array)
  }
)

# function for seasonal growth kernel -- biweekly time step
get_kernel_nomort <- nimbleFunction (
  # input and output types
  run = function(xinf = double(0), k = double(0),
                 sigma_G = double(0), A = double(0),
                 ds = double(0), t1 = double(0), t2 = double(0),
                 n_size = double(0), pi = double(0), x = double(1),
                 lower = double(1), upper = double(1))
  {
    returnType(double(2))
    
    # create empty array
    array <- matrix(NA, ncol = n_size, nrow = n_size)
    
    if(!is.na(t2) & t2 == 0) {
      array <- diag(n_size)
    } else {
      # season adjusted params
      S_t <- (A * k / (2 * pi)) * sin(2 * pi * (t2 - (1 + ds)))
      S_t0 <- (A * k / (2 * pi)) * sin(2 * pi * (t1 - (1 + ds)))
      
      # p(y"|y)
      for (i in 1:n_size) {
        increment <- (xinf - x[i]) *
          (1 - exp(-k * (t2 - t1) - S_t + S_t0))
        mean <- x[i] + increment
        array[1:n_size, i] <- (pnorm(upper, mean, sd = sigma_G) -
                                 pnorm(lower, mean, sd = sigma_G))
      }
      
      # normalize and apply natural mortality
      for (i in 1:n_size) {
        array[, i] <- array[, i] / sum(array[, i])
      }
    }
    return(array)
  }
)

# D4
kernel_D4_no_mort <- as.data.frame(get_kernel_nomort(xinf = summary["xinf", "50%"], 
                                                     k = summary["gk", "50%"], 
                                                     sigma_G = summary["sigma_G", "50%"], 
                                                     A = summary["A", "50%"], 
                                                     ds = summary["ds", "50%"], 
                                                     t1 = D1[4], t2 = D2[4],
                                                     n_size = n_size, pi = pi, 
                                                     x = x, lower = lower, upper = upper)) %>% 
  mutate(xprime = x)
colnames(kernel_D4_no_mort) <- c(x, "xprime")

kernel_D4_no_mort <- kernel_D4_no_mort %>% 
  pivot_longer(cols = -xprime,
               values_to = "prob",
               names_to = "x")

plot_D4 <- ggplot() +
  geom_tile(data = kernel_D4_no_mort,
            aes(x = as.numeric(x), y = xprime, fill = prob)) +
  scale_fill_viridis(limits = c(0, max(kernel_D3_no_mort$prob))) +
  geom_line(aes(x = seq(min(lower), max(upper), 1), 
                y = seq(min(lower), max(upper), 1)),
            linetype = "dashed", color = "firebrick") +
  labs(x = "x", y = "x'", fill = "probability") +
  ggtitle("January 1st") +
  theme_minimal()

# D1
kernel_D1_no_mort <- as.data.frame(get_kernel_nomort(xinf = summary["xinf", "50%"], 
                                       k = summary["gk", "50%"], 
                                       sigma_G = summary["sigma_G", "50%"], 
                                       A = summary["A", "50%"], 
                                       ds = summary["ds", "50%"], 
                                       t1 = D1[1], t2 = D2[1],
                                       n_size = n_size, pi = pi, 
                                       x = x, lower = lower, upper = upper)) %>% 
  mutate(xprime = x)
colnames(kernel_D1_no_mort) <- c(x, "xprime")

kernel_D1_no_mort <- kernel_D1_no_mort %>% 
  pivot_longer(cols = -xprime,
               values_to = "prob",
               names_to = "x")

plot_D1 <- ggplot() +
  geom_tile(data = kernel_D1_no_mort,
            aes(x = as.numeric(x), y = xprime, fill = prob)) +
  scale_fill_viridis(limits = c(0, max(kernel_D3_no_mort$prob))) +
  geom_line(aes(x = seq(min(lower), max(upper), 1), 
                y = seq(min(lower), max(upper), 1)),
            linetype = "dashed", color = "firebrick") +
  labs(x = "x", y = "x'", fill = "probability") +
  ggtitle("April 1st") +
  theme_minimal()

# D2
kernel_D2_no_mort <- as.data.frame(get_kernel_nomort(xinf = summary["xinf", "50%"], 
                                                     k = summary["gk", "50%"], 
                                                     sigma_G = summary["sigma_G", "50%"], 
                                                     A = summary["A", "50%"], 
                                                     ds = summary["ds", "50%"], 
                                                     t1 = D1[2], t2 = D2[2],
                                                     n_size = n_size, pi = pi, 
                                                     x = x, lower = lower, upper = upper)) %>% 
  mutate(xprime = x)
colnames(kernel_D2_no_mort) <- c(x, "xprime")

kernel_D2_no_mort <- kernel_D2_no_mort %>% 
  pivot_longer(cols = -xprime,
               values_to = "prob",
               names_to = "x")

plot_D2 <- ggplot() +
  geom_tile(data = kernel_D2_no_mort,
            aes(x = as.numeric(x), y = xprime, fill = prob)) +
  scale_fill_viridis(limits = c(0, max(kernel_D3_no_mort$prob))) +
  geom_line(aes(x = seq(min(lower), max(upper), 1), 
                y = seq(min(lower), max(upper), 1)),
            linetype = "dashed", color = "firebrick") +
  labs(x = "x", y = "x'", fill = "probability") +
  ggtitle("July 1st") +
  theme_minimal()

# D3
kernel_D3_no_mort <- as.data.frame(get_kernel_nomort(xinf = summary["xinf", "50%"], 
                                                     k = summary["gk", "50%"], 
                                                     sigma_G = summary["sigma_G", "50%"], 
                                                     A = summary["A", "50%"], 
                                                     ds = summary["ds", "50%"], 
                                                     t1 = D1[3], t2 = D2[3],
                                                     n_size = n_size, pi = pi, 
                                                     x = x, lower = lower, upper = upper)) %>% 
  mutate(xprime = x)
colnames(kernel_D3_no_mort) <- c(x, "xprime")

kernel_D3_no_mort <- kernel_D3_no_mort %>% 
  pivot_longer(cols = -xprime,
               values_to = "prob",
               names_to = "x")

plot_D3 <- ggplot() +
  geom_tile(data = kernel_D3_no_mort,
            aes(x = as.numeric(x), y = xprime, fill = prob)) +
  scale_fill_viridis(limits = c(0, max(kernel_D3_no_mort$prob))) +
  geom_line(aes(x = seq(min(lower), max(upper), 1), 
                y = seq(min(lower), max(upper), 1)),
            linetype = "dashed", color = "firebrick") +
  labs(x = "x", y = "x'", fill = "probability") +
  ggtitle("October 1st") +
  theme_minimal()

plot_D2 + plot_D4 + plot_layout(nrow = 1)

survival <- nimbleFunction (
  
  run = function(alpha = double(0), beta = double(0), x = double(1),
                 t1 = double(0), t2 = double(0))
  {
    returnType(double(1))
    
    # number of biweeks
    deltat <- round((t2 - t1) * (52.1429 / 2))
    
    # get survival rate
    out <- exp(-deltat * (beta + alpha / x ^ 2))
    
    return(out)
  }
)

# D2 - mort
kernel_D2_mort <- as.data.frame(get_kernel(xinf = summary["xinf", "50%"], 
                                                     k = summary["gk", "50%"], 
                                                     sigma_G = summary["sigma_G", "50%"], 
                                                     A = summary["A", "50%"], 
                                                     ds = summary["ds", "50%"], 
                                                     t1 = D1[2], t2 = D2[2],
                                                     n_size = n_size, pi = pi, 
                                                     x = x, lower = lower, upper = upper,,
                                survival(summary["alpha", "50%"], 
                                         summary["beta", "50%"], 
                                         x = x, t1 = D1[2], t2 = D2[2]))) %>% 
  mutate(xprime = x)
colnames(kernel_D2_no_mort) <- c(x, "xprime")

kernel_D2_no_mort <- kernel_D2_no_mort %>% 
  pivot_longer(cols = -xprime,
               values_to = "prob",
               names_to = "x")

plot_D2 <- ggplot() +
  geom_tile(data = kernel_D2_no_mort,
            aes(x = as.numeric(x), y = xprime, fill = prob)) +
  scale_fill_viridis(limits = c(0, max(kernel_D3_no_mort$prob))) +
  geom_line(aes(x = seq(min(lower), max(upper), 1), 
                y = seq(min(lower), max(upper), 1)),
            linetype = "dashed", color = "firebrick") +
  labs(x = "x", y = "x'", fill = "probability") +
  ggtitle("July 1st") +
  theme_minimal()



t <- seq(from = 0.3, to = 3.59, by = 0.01)

# function to get length
get_length <- function(c, k, ts, t0, xinf) {
  s_t <- (c * k / (2 * pi)) * sin(2 * pi * (t - ts))
  s_t0 <- (c * k / (2 * pi)) * sin(2 * pi * (t0 - ts))
  length <- xinf * (1 - exp(-k * (t - t0) - s_t + s_t0))
  return(length)
}

# function for getting samples
get_samples <- function(samples, param) {
  samples <- c(out_sub[[1]][, param], out_sub[[2]][, param],
               out_sub[[3]][, param], out_sub[[4]][, param])
  return(samples)
}

get_hdi_low <- function(input, ci) {
  out <- as.numeric(hdi(input, ci)[2])
  return(out)
}
get_hdi_high <- function(input, ci) {
  out <- as.numeric(hdi(input, ci)[3])
  return(out)
}

# get params
c <- get_samples(out_sub, "A")
k <- get_samples(out_sub, "gk")
ts <- get_samples(out_sub, "ds")
t0 <- get_samples(out_sub, "t0")
xinf <- get_samples(out_sub, "xinf")

# get size for each posterior sample
growth_len <- mapply(get_length, c, k, ts, t0, xinf)

# get summaries
growth_summary <- as.data.frame(matrix(NA, nrow = 3, ncol = length(t)))
growth_summary[1, ] <- apply(growth_len, 1, median)
growth_summary[2, ] <- apply(growth_len, 1,
                             function(row) get_hdi_low(row, 0.95))
growth_summary[3, ] <- apply(growth_len, 1,
                             function(row) get_hdi_high(row, 0.95))
colnames(growth_summary) <- t
growth_summary$stat <- c("median", "low", "high")
growth_summary_long <- growth_summary %>%
  pivot_longer(cols = -stat,
               names_to = "t",
               values_to = "value") %>%
  pivot_wider(values_from = "value",
              names_from = "stat")

# plot sizes
growth_plot <- ggplot() +
  geom_point(data = data,
             aes(x = age,
                 y = CW),
             size = 2, pch = 21, alpha = 0.8) +
  geom_ribbon(data = growth_summary_long,
              aes(x = as.numeric(t),
                  ymin = low, ymax = high),
              alpha = 0.4, fill = "indianred2") +
  geom_line(data = growth_summary_long,
            aes(x = as.numeric(t), y = median),
            linewidth = 1, color = "indianred2") +
  scale_x_continuous(limits = c(0.2, 3.6)) +
  labs(x = "age (years)",
       y = "size (mm)") +
  theme_minimal() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
ggsave("figures/appendix_growth_plot.png", growth_plot,
       dpi = 400, width = 5, height = 4)
