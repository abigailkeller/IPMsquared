library(tidyverse)


############
# get data #
############

# read in time series data
counts <- readRDS("data/model_data/counts_w_NA.rds")

#################
#################
# deviance plot #
#################
#################

##############
# Roche Cove #
##############

# read in data
deviance_Drep <- readRDS("data/posterior_predictive_check/deviance_Drep.rds")
deviance_D <- readRDS("data/posterior_predictive_check/deviance_D.rds")

# join data
joined_data <- as.data.frame(cbind(deviance_Drep, deviance_D)) %>% 
  pivot_longer(cols = c("deviance_Drep", "deviance_D"),
               names_to = "type", values_to = "deviance")

ppp_plot <- ggplot(data = joined_data) +
  geom_histogram(aes(x = deviance, fill = type), alpha = 0.3) +
  scale_fill_manual(labels = c("D", bquote(D^rep)), 
                    values = c("red", "blue")) +
  labs(x = "deviance", y = "count", fill = "") +
  theme_minimal()
ggsave("figures/ppp_deviance_plot.png", ppp_plot, dpi = 440,
       width = 5, height = 4)

pvalue_deviance <- mean(deviance_D > deviance_Drep)

###################
# Seadrift Lagoon #
###################

# read in data
deviance_Drep_SL <- readRDS("data/posterior_predictive_check/seadrift/deviance_Drep.rds")
deviance_D_SL <- readRDS("data/posterior_predictive_check/seadrift/deviance_D.rds")

# join data
joined_data_SL <- as.data.frame(cbind(deviance_Drep_SL, deviance_D_SL)) %>% 
  pivot_longer(cols = c("deviance_Drep_SL", "deviance_D_SL"),
               names_to = "type", values_to = "deviance")

ppp_plot_SL <- ggplot(data = joined_data_SL) +
  geom_histogram(aes(x = deviance, fill = type), alpha = 0.3) +
  scale_fill_manual(labels = c("D", bquote(D^rep)), 
                    values = c("red", "blue")) +
  labs(x = "deviance", y = "count", fill = "") +
  theme_minimal()
ggsave("figures/ppp_deviance_plot_SL.png", ppp_plot_SL, dpi = 440,
       width = 5, height = 4)

pvalue_deviance_SL <- mean(deviance_D_SL > deviance_Drep_SL)


#######################
#######################
# zero inflation plot #
#######################
#######################

##############
# Roche Cove #
##############

# get total count
zero_count_total <- as.data.frame(matrix(0, ncol = 3, nrow = 6404))
colnames(zero_count_total) <- c("iter", "total", "zero")
zero_count_total$iter <- 1:6404
for(i in 1:14) {
  file <- paste0("data/posterior_predictive_check/PPS_",i,".rds")
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

###################
# Seadrift Lagoon #
###################

# get total count
zero_count_total_SL <- as.data.frame(matrix(0, ncol = 3, nrow = 6404))
colnames(zero_count_total_SL) <- c("iter", "total", "zero")
zero_count_total_SL$iter <- 1:6404
for(i in 1:14) {
  file <- paste0("data/posterior_predictive_check/seadrift/PPS_",i,".rds")
  input <- readRDS(file)
  total <- ncol(input)
  zero <- apply(input, 1, function(x) sum(x == 0))
  zero_count_total_SL$total <- zero_count_total$total + total
  zero_count_total_SL$zero <- zero_count_total$zero + zero
}

zero_count_total_SL <- zero_count_total_SL %>% 
  mutate(proportion = zero/total)

data_prop_total_SL <- (
  sum(counts[,,,] == 0, na.rm = TRUE) /
    sum(!is.na(counts[,,,] == 0))
)

pvalue_zero_SL <- sum(data_prop_total_SL > zero_count_total_SL$prop) / 
  length(zero_count_total_SL$prop)

total_zeros_plot_SL <- ggplot() +
  geom_histogram(aes(x = zero_count_total_SL$prop)) +
  geom_vline(aes(xintercept = data_prop_total), color = "red") +
  labs(x = "proportion of zeros", y = "count") +
  theme_minimal()

ggsave("figures/ppp_zeros_SL.png", total_zeros_plot_SL, 
       dpi = 600, height = 4, width = 4)
