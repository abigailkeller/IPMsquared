library(tidyverse)


############
# get data #
############

# read in time series data
counts <- readRDS("data/model_data/counts_w_NA.rds")
ncap <- readRDS("data/model_data/n_cap_w_NA.rds")

# read in mark-recapture data
m2_mc <- readRDS("data/model_data/m2_mc.rds")

#################
# deviance plot #
#################

# read in data
deviance_yrep <- readRDS("data/posterior_predictive_check/deviance_yrep.rds")
deviance_y <- readRDS("data/posterior_predictive_check/deviance_y.rds")

# join data
joined_data <- as.data.frame(cbind(deviance_yrep, deviance_y)) %>% 
  pivot_longer(cols = c("deviance_yrep", "deviance_y"),
               names_to = "type", values_to = "deviance")

ppp_plot <- ggplot(data = joined_data) +
  geom_histogram(aes(x = deviance, fill = type), alpha = 0.3) +
  scale_fill_manual(labels = c("y", bquote(y^rep)), 
                    values = c("red", "blue")) +
  labs(x = "deviance", y = "count", fill = "") +
  theme_minimal()
ggsave("figures/ppp_deviance_plot.png", dpi = 440,
       width = 5, height = 4)

pvalue_deviance <- mean(deviance_y > deviance_yrep)

#######################
# zero inflation plot #
#######################

# get total count
zero_count_total <- as.data.frame(matrix(0, ncol = 3, nrow = 6404))
colnames(zero_count_total) <- c("iter", "total", "zero")
zero_count_total$iter <- 1:6404
for(i in 1:14) {
  file <- paste0("posterior_predictive_check/ppSamples_count",i,".rds")
  input <- readRDS(file)
  total <- ncol(input)
  zero <- apply(input, 1, function(x) sum(x == 0))
  zero_count_total$total <- zero_count_total$total + total
  zero_count_total$zero <- zero_count_total$zero + zero
}

zero_count_total <- zero_count_total %>% 
  mutate(proportion = zero/total)

data_prop_total <- (
  sum(counts[,,,] == 0, na.rm = TRUE) /
    sum(!is.na(counts[,,,] == 0))
)

pvalue_zero <- sum(data_prop_total > zero_count_total$prop) / 
  length(zero_count_total$prop)

total_zeros_plot <- ggplot() +
  geom_histogram(aes(x = zero_count_total$prop)) +
  geom_vline(aes(xintercept = data_prop_total), color = "red") +
  labs(x = "proportion of zeros", y = "count") +
  theme_minimal()

ggsave("figures/ppp_zeros.png", total_zeros_plot, 
       dpi = 600, height = 4, width = 4)
