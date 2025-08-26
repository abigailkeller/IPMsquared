library(tidyverse)
library(bayestestR)

# read in data
x <- readRDS("data/model_data/x.rds")[1:20]
n_size <- length(x)
b <- readRDS("data/model_data/b.rds")[1:21]
upper <- b[2:length(b)]
lower <- b[1:length(x)]

# read in samples
out <- readRDS("data/posterior_samples/savedsamples_IPM.rds")
summary <- MCMCsummary(out)

# time points
D1 <- c(0, 0.25, 0.5, 0.75)
D2 <- D1 + 14 / 365

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

# get max fill
fill_max <- max(c(kernel_D1_no_mort$prob, kernel_D2_no_mort$prob, 
                  kernel_D3_no_mort$prob, kernel_D4_no_mort$prob))

plot_D1 <- ggplot() +
  geom_tile(data = kernel_D1_no_mort,
            aes(x = as.numeric(x), y = xprime, fill = prob)) +
  scale_fill_viridis(limits = c(0, fill_max)) +
  geom_line(aes(x = seq(min(lower), max(upper), 1), 
                y = seq(min(lower), max(upper), 1)),
            linetype = "dashed", color = "firebrick") +
  labs(x = "x", y = "x'", fill = "probability") +
  ggtitle("April 1st") +
  theme_minimal() +
  theme(legend.position = "None",
        text = element_text(family = "Arial"))

plot_D2 <- ggplot() +
  geom_tile(data = kernel_D2_no_mort,
            aes(x = as.numeric(x), y = xprime, fill = prob)) +
  scale_fill_viridis(limits = c(0, fill_max)) +
  geom_line(aes(x = seq(min(lower), max(upper), 1), 
                y = seq(min(lower), max(upper), 1)),
            linetype = "dashed", color = "firebrick") +
  labs(x = "x", y = "x'", fill = "probability") +
  ggtitle("July 1st") +
  theme_minimal() +
  theme(legend.position = "None",
        text = element_text(family = "Arial"))

plot_D3 <- ggplot() +
  geom_tile(data = kernel_D3_no_mort,
            aes(x = as.numeric(x), y = xprime, fill = prob)) +
  scale_fill_viridis(limits = c(0, fill_max)) +
  geom_line(aes(x = seq(min(lower), max(upper), 1), 
                y = seq(min(lower), max(upper), 1)),
            linetype = "dashed", color = "firebrick") +
  labs(x = "x", y = "x'", fill = "probability") +
  ggtitle("October 1st") +
  theme_minimal() +
  theme(legend.position = "None",
        text = element_text(family = "Arial"))

plot_D4 <- ggplot() +
  geom_tile(data = kernel_D4_no_mort,
            aes(x = as.numeric(x), y = xprime, fill = prob)) +
  scale_fill_viridis(limits = c(0, fill_max)) +
  geom_line(aes(x = seq(min(lower), max(upper), 1), 
                y = seq(min(lower), max(upper), 1)),
            linetype = "dashed", color = "firebrick") +
  labs(x = "x", y = "x'", fill = "probability") +
  ggtitle("January 1st") +
  theme_minimal() +
  theme(legend.position = "None",
        text = element_text(family = "Arial"))

plot_w_legend <- ggplot() +
  geom_tile(data = kernel_D4_no_mort,
            aes(x = as.numeric(x), y = xprime, fill = prob)) +
  scale_fill_viridis(limits = c(0, fill_max)) +
  geom_line(aes(x = seq(min(lower), max(upper), 1), 
                y = seq(min(lower), max(upper), 1)),
            linetype = "dashed", color = "firebrick") +
  labs(x = "x", y = "x'", fill = "probability") +
  ggtitle("January 1st") +
  theme_minimal() +
  theme(text = element_text(family = "Arial"))

legend <- get_legend(plot_w_legend)

layout <- c(
  "AABB#
  AABBE
  CCDDE
  CCDD#"
)

final_plot <- plot_D1 + plot_D2 + plot_D3 + plot_D4 + legend + 
  plot_layout(design = layout)


ggsave("figures/appendix_growth_plot.png", final_plot,
       dpi = 400, width = 6, height = 5)
