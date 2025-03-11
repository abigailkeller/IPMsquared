library(tidyverse)
library(bayestestR)
library(patchwork)
library(cowplot)

# read in samples
out <- readRDS("data/posterior_samples/savedsamples_IPM.rds")

# subset samples
lower <- 2000
upper <- 10001
out_sub <- list(
  out[[1]][lower:upper, ], out[[2]][lower:upper, ],
  out[[3]][lower:upper, ], out[[4]][lower:upper, ]
)

# function for getting samples
get_samples <- function(samples, param) {
  samples <- c(samples[[1]][, param], samples[[2]][, param],
               samples[[3]][, param], samples[[4]][, param])
  return(samples)
}

# get params
trapf_k <- get_samples(out_sub, "h_F_k")
trapf_midpoint <- get_samples(out_sub, "h_F_0")
trapf_pmax <- get_samples(out_sub, "h_F_max")
trapm_pmax <- get_samples(out_sub, "h_M_max")
trapm_sigma <- get_samples(out_sub, "h_M_sigma")
trapm_xmax <- get_samples(out_sub, "h_M_A")
traps_k <- get_samples(out_sub, "h_S_k")
traps_midpoint <- get_samples(out_sub, "h_S_0")
traps_pmax <- get_samples(out_sub, "h_S_max")

# make plot of size selectivity
min_size <- 0
max_size <- 110
y <- seq(min_size, max_size, 1)
selective_summaries <- data.frame(
  type = c("minnow", "fukui", "shrimp",
           "minnow", "fukui", "shrimp",
           "minnow", "fukui", "shrimp"),
  stat = c(rep("median", 3), rep("lower_ci", 3), rep("upper_ci", 3))
)
selective_summaries <- cbind(
  selective_summaries,
  as.data.frame(matrix(NA, nrow = dim(selective_summaries), ncol = length(y)))
)
colnames(selective_summaries)[3:113] <- y
size_sel_norm <- function(pmax, xmax, sigma) {
  vector <- pmax * exp(-(y - xmax) ^ 2 / (2 * sigma ^ 2))
  return(1 - exp(-vector))
}
size_sel_log <- function(pmax, k, midpoint) {
  vector <- pmax / (1 + exp(-k * (y - midpoint)))
  return(1 - exp(-vector))
}
get_hdi_low <- function(input, ci) {
  out <- as.numeric(hdi(input, ci)[2])
  return(out)
}
get_hdi_high <- function(input, ci) {
  out <- as.numeric(hdi(input, ci)[3])
  return(out)
}
minnow_select <- mapply(size_sel_norm, trapm_pmax,
                        trapm_xmax, trapm_sigma)
fukui_select <- mapply(size_sel_log, trapf_pmax,
                       trapf_k, trapf_midpoint)
shrimp_select <- mapply(size_sel_log, traps_pmax,
                        traps_k, traps_midpoint)
selective_summaries[1, 3:113] <- apply(minnow_select, 1, median)
selective_summaries[2, 3:113] <- apply(fukui_select, 1, median)
selective_summaries[3, 3:113] <- apply(shrimp_select, 1, median)
selective_summaries[4, 3:113] <- apply(minnow_select, 1,
                                       function(row) get_hdi_low(row, 0.95))
selective_summaries[5, 3:113] <- apply(fukui_select, 1,
                                       function(row) get_hdi_low(row, 0.95))
selective_summaries[6, 3:113] <- apply(shrimp_select, 1,
                                       function(row) get_hdi_low(row, 0.95))
selective_summaries[7, 3:113] <- apply(minnow_select, 1,
                                       function(row) get_hdi_high(row, 0.95))
selective_summaries[8, 3:113] <- apply(fukui_select, 1,
                                       function(row) get_hdi_high(row, 0.95))
selective_summaries[9, 3:113] <- apply(shrimp_select, 1,
                                       function(row) get_hdi_high(row, 0.95))

# convert from wide to long
selective_summaries_long <- selective_summaries %>%
  pivot_longer(cols = ! c("type", "stat"),
               values_to = "p",
               names_to = "size") %>%
  pivot_wider(values_from = "p",
              names_from = "stat")
selective_summaries_long$type <- factor(selective_summaries_long$type,
                                        levels = c("fukui", "shrimp", "minnow"))

##### size selectivity plot
sizesel_plot <- ggplot(data = selective_summaries_long) +
  geom_ribbon(aes(x = as.numeric(size),
                  ymin = lower_ci, ymax = upper_ci,
                  fill = as.factor(type)), alpha = 0.4) +
  geom_line(aes(x = as.numeric(size), y = median,
                color = as.factor(type)),
            linewidth = 1) +
  coord_trans(y = "sqrt") +
  scale_fill_manual(values = c("violet", "goldenrod", "darkviolet"),
                    guide = "none") +
  scale_color_manual(values = c("violet", "goldenrod", "darkviolet"),
                     labels = c("Fukui", "Shrimp", "Minnow")) +
  scale_y_continuous(breaks = c(0.0001, 0.001, 0.002, 0.003, 0.004)) +
  labs(x = "size (mm)", y = "probability of capture", color = "trap type") +
  theme_minimal() +
  theme(text = element_text(size = 14))

ggsave("figures/Figure4_sizesel.png", sizesel_plot,
       dpi = 400, width = 5, height = 4,
       bg = "white")
