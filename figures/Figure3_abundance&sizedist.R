library(tidyverse)
library(bayestestR)
library(patchwork)
library(cowplot)
library(viridis)

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
n_adult_1 <- get_samples(out_sub, "N_adult")
n_recruit_1 <- get_samples(out_sub, "N_recruit[1]")
n_recruit_2 <- get_samples(out_sub, "N_recruit[2]")
n_recruit_3 <- get_samples(out_sub, "N_recruit[3]")
n_recruit_4 <- get_samples(out_sub, "N_recruit[4]")
init_mean_recruit <- get_samples(out_sub, "init_mean_recruit")
init_sd_r <- get_samples(out_sub, "init_sd_r")
init_lsd_adult <- get_samples(out_sub, "init_lsd_adult")
init_lmean_adult <- get_samples(out_sub, "init_lmean_adult")
adult_dist_2 <- as.data.frame(matrix(NA, ncol = 22,
                                     nrow = dim(out_sub[[1]])[1] * 4))
adult_dist_3 <- as.data.frame(matrix(NA, ncol = 22,
                                     nrow = dim(out_sub[[1]])[1] * 4))
adult_dist_4 <- as.data.frame(matrix(NA, ncol = 22,
                                     nrow = dim(out_sub[[1]])[1] * 4))
for (i in 1:22) {
  adult_dist_2[, i] <- get_samples(out_sub, paste0("N_overwinter[1, ", i, "]"))
  adult_dist_3[, i] <- get_samples(out_sub, paste0("N_overwinter[2, ", i, "]"))
  adult_dist_4[, i] <- get_samples(out_sub, paste0("N_overwinter[3, ", i, "]"))
}
n_adult2 <- rowSums(adult_dist_2)
n_adult3 <- rowSums(adult_dist_3)
n_adult4 <- rowSums(adult_dist_4)
# make plots of total N
n_summaries <- data.frame(
  param = c("n_adult_1", "n_adult2", "n_adult3", "n_adult4",
            "n_recruit_1", "n_recruit_2", "n_recruit_3", "n_recruit_4"),
  type = c("adult", "adult", "adult", "adult", "recruit", "recruit",
           "recruit", "recruit"),
  year = c("2020", "2021", "2022", "2023",
           "2020", "2021", "2022", "2023"),
  median = c(median(n_adult_1), median(n_adult2),
             median(n_adult3), median(n_adult4),
             median(n_recruit_1), median(n_recruit_2),
             median(n_recruit_3), median(n_recruit_4)),
  lower_CI = c(as.numeric(hdi(n_adult_1)[2]), as.numeric(hdi(n_adult2)[2]),
               as.numeric(hdi(n_adult3)[2]), as.numeric(hdi(n_adult4)[2]),
               as.numeric(hdi(n_recruit_1)[2]), as.numeric(hdi(n_recruit_2)[2]),
               as.numeric(hdi(n_recruit_3)[2]),
               as.numeric(hdi(n_recruit_4)[2])),
  upper_CI = c(as.numeric(hdi(n_adult_1)[3]), as.numeric(hdi(n_adult2)[3]),
               as.numeric(hdi(n_adult3)[3]), as.numeric(hdi(n_adult4)[3]),
               as.numeric(hdi(n_recruit_1)[3]), as.numeric(hdi(n_recruit_2)[3]),
               as.numeric(hdi(n_recruit_3)[3]), as.numeric(hdi(n_recruit_4)[3]))
)
# make plots of N distribution - recruit
min_size <- 0
max_size <- 110
y_recruit <- seq(min_size, max_size, 1)
n <- 22
b <- min_size + c(0:n) * (max_size - min_size) / n
y <- 0.5 * (b[1:n] + b[2:(n + 1)])
init_mean_recruit <- get_samples(out_sub, "init_mean_recruit")
init_sd_r <- get_samples(out_sub, "init_sd_r")
init_lsd_adult <- get_samples(out_sub, "init_lsd_adult")
init_lmean_adult_1 <- get_samples(out_sub, "init_lmean_adult")
dist_summaries_recruit <- data.frame(
  param = c("n_recruit_1_dist", "n_recruit_2_dist", "n_recruit_3_dist",
            "n_recruit_4_dist", "n_recruit_1_dist", "n_recruit_2_dist",
            "n_recruit_3_dist", "n_recruit_4_dist", "n_recruit_1_dist",
            "n_recruit_2_dist", "n_recruit_3_dist", "n_recruit_4_dist"),
  year = c("2020", "2021", "2022", "2023",
           "2020", "2021", "2022", "2023",
           "2020", "2021", "2022", "2023"),
  stat = c(rep("median", 4), rep("lower_ci", 4), rep("upper_ci", 4))
)
dist_summaries_recruit <- cbind(
  dist_summaries_recruit,
  as.data.frame(matrix(NA, nrow = dim(dist_summaries_recruit)[1],
                       ncol = length(y)))
)
colnames(dist_summaries_recruit)[4:25] <- y
dist_summaries_adult <- data.frame(
  param = c("n_adult_1_dist", "n_adult2_dist", "n_adult3_dist", "n_adult4_dist",
            "n_adult_1_dist", "n_adult2_dist", "n_adult3_dist", "n_adult4_dist",
            "n_adult_1_dist", "n_adult2_dist", "n_adult3_dist",
            "n_adult4_dist"),
  year = c("2020", "2021", "2022", "2023",
           "2020", "2021", "2022", "2023",
           "2020", "2021", "2022", "2023"),
  stat = c(rep("median", 4), rep("lower_ci", 4), rep("upper_ci", 4))
)
dist_summaries_adult <- cbind(
  dist_summaries_adult, as.data.frame(matrix(NA,
                                             nrow = dim(dist_summaries_adult),
                                             ncol = length(y)))
)
colnames(dist_summaries_adult)[4:25] <- y
get_dist_log <- function(mean, sd, n) {
  out <- dlnorm(y, mean, sd) / sum(dlnorm(y, mean, sd)) * n
  return(out)
}
get_dist <- function(mean, sd, n) {
  out <- dnorm(y, mean, sd) / sum(dnorm(y, mean, sd)) * n
  return(out)
}
get_hdi_low <- function(input, ci) {
  out <- as.numeric(hdi(input, ci)[2])
  return(out)
}
get_hdi_high <- function(input, ci) {
  out <- as.numeric(hdi(input, ci)[3])
  return(out)
}
n_adult_1_dist <- mapply(get_dist_log, init_lmean_adult_1,
                         init_lsd_adult, n_adult_1)
n_recruit_1_dist <- mapply(get_dist, init_mean_recruit,
                           init_sd_r, n_recruit_1)
n_recruit_2_dist <- mapply(get_dist, init_mean_recruit,
                           init_sd_r, n_recruit_2)
n_recruit_3_dist <- mapply(get_dist, init_mean_recruit,
                           init_sd_r, n_recruit_3)
n_recruit_4_dist <- mapply(get_dist, init_mean_recruit,
                           init_sd_r, n_recruit_4)
dist_summaries_adult[1, 4:25] <- apply(n_adult_1_dist, 1, median)
dist_summaries_adult[2, 4:25] <- apply(adult_dist_2, 2, median)
dist_summaries_adult[3, 4:25] <- apply(adult_dist_3, 2, median)
dist_summaries_adult[4, 4:25] <- apply(adult_dist_4, 2, median)
dist_summaries_adult[5, 4:25] <- apply(n_adult_1_dist, 1,
                                       function(row) get_hdi_low(row, 0.95))
dist_summaries_adult[6, 4:25] <- apply(adult_dist_2, 2,
                                       function(row) get_hdi_low(row, 0.95))
dist_summaries_adult[7, 4:25] <- apply(adult_dist_3, 2,
                                       function(row) get_hdi_low(row, 0.95))
dist_summaries_adult[8, 4:25] <- apply(adult_dist_4, 2,
                                       function(row) get_hdi_low(row, 0.95))
dist_summaries_adult[9, 4:25] <- apply(n_adult_1_dist, 1,
                                       function(row) get_hdi_high(row, 0.95))
dist_summaries_adult[10, 4:25] <- apply(adult_dist_2, 2,
                                        function(row) get_hdi_high(row, 0.95))
dist_summaries_adult[11, 4:25] <- apply(adult_dist_3, 2,
                                        function(row) get_hdi_high(row, 0.95))
dist_summaries_adult[12, 4:25] <- apply(adult_dist_4, 2,
                                        function(row) get_hdi_high(row, 0.95))
dist_summaries_recruit[1, 4:25] <- apply(n_recruit_1_dist, 1,
                                         function(x) median(x, na.rm = TRUE))
dist_summaries_recruit[2, 4:25] <- apply(n_recruit_2_dist, 1,
                                         function(x) median(x, na.rm = TRUE))
dist_summaries_recruit[3, 4:25] <- apply(n_recruit_3_dist, 1,
                                         function(x) median(x, na.rm = TRUE))
dist_summaries_recruit[4, 4:25] <- apply(n_recruit_4_dist, 1,
                                         function(x) median(x, na.rm = TRUE))
dist_summaries_recruit[5, 4:25] <- apply(n_recruit_1_dist, 1,
                                         function(row) get_hdi_low(row, 0.5))
dist_summaries_recruit[6, 4:25] <- apply(n_recruit_2_dist, 1,
                                         function(row) get_hdi_low(row, 0.5))
dist_summaries_recruit[7, 4:25] <- apply(n_recruit_3_dist, 1,
                                         function(row) get_hdi_low(row, 0.5))
dist_summaries_recruit[8, 4:25] <- apply(n_recruit_4_dist, 1,
                                         function(row) get_hdi_low(row, 0.5))
dist_summaries_recruit[9, 4:25] <- apply(n_recruit_1_dist, 1,
                                         function(row) get_hdi_high(row, 0.5))
dist_summaries_recruit[10, 4:25] <- apply(n_recruit_2_dist, 1,
                                          function(row) get_hdi_high(row, 0.5))
dist_summaries_recruit[11, 4:25] <- apply(n_recruit_3_dist, 1,
                                          function(row) get_hdi_high(row, 0.5))
dist_summaries_recruit[12, 4:25] <- apply(n_recruit_4_dist, 1,
                                          function(row) get_hdi_high(row, 0.5))
# convert from wide to long
dist_summaries_adult_long <- dist_summaries_adult %>%
  pivot_longer(cols = ! c("param", "year", "stat"),
               values_to = "N",
               names_to = "size") %>%
  pivot_wider(values_from = "N",
              names_from = "stat")
dist_summaries_recruit_long <- dist_summaries_recruit %>%
  pivot_longer(cols = ! c("param", "year", "stat"),
               values_to = "N",
               names_to = "size") %>%
  pivot_wider(values_from = "N",
              names_from = "stat")
##### total abundance plots
adult_nplot <- ggplot(data = n_summaries[n_summaries$type == "adult", ]) +
  geom_point(aes(x = param, y = median,
                 color = as.factor(year))) +
  geom_errorbar(aes(x = param,
                    ymin = lower_CI,
                    ymax = upper_CI,
                    color = as.factor(year)),
                width = 0.2,
                linewidth = 1) +
  scale_color_manual(values = viridis_pal()(4)) +
  scale_y_continuous(breaks = c(200, 300, 400)) +
  labs(x = "year", y = expression(N[total]), color = "year") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "None")
recruit_nplot <- ggplot(data = n_summaries[n_summaries$type == "recruit", ]) +
  geom_point(aes(x = param, y = log(median),
                 color = as.factor(year))) +
  geom_errorbar(aes(x = param,
                    ymin = log(lower_CI),
                    ymax = log(upper_CI),
                    color = as.factor(year)),
                width = 0.2,
                linewidth = 1) +
  scale_y_continuous(breaks = c(log(10), log(100), log(1000)),
                     labels = c("10", "100", "1000")) +
  scale_color_manual(values = viridis_pal()(4)) +
  labs(x = "year", y = expression(N[total]), color = "year") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "None")
plot_w_legend <- ggplot(data = n_summaries[n_summaries$type == "recruit", ]) +
  geom_point(aes(x = param, y = log(median),
                 color = as.factor(year))) +
  geom_errorbar(aes(x = param,
                    ymin = log(lower_CI),
                    ymax = log(upper_CI),
                    color = as.factor(year)),
                width = 0.2,
                linewidth = 1) +
  scale_color_manual(values = viridis_pal()(4)) +
  scale_y_continuous(breaks = c(log(20), log(200)),
                     labels = c("20", "200")) +
  labs(x = "year", y = expression(N[total]), color = "year") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
legend <- get_legend(plot_w_legend)
######## size abundance plots
adult_sizeplot <- ggplot(data = dist_summaries_adult_long) +
  geom_ribbon(aes(x = as.numeric(size),
                  ymin = lower_ci, ymax = upper_ci,
                  fill = as.factor(year)), alpha = 0.4) +
  geom_line(aes(x = as.numeric(size), y = median,
                color = as.factor(year)),
            linewidth = 1) +
  scale_x_continuous(limits = c(2.5, 107.5)) +
  scale_color_manual(values = viridis_pal()(4)) +
  scale_fill_manual(values = viridis_pal()(4)) +
  labs(x = "size (mm)", y = expression(N[size]), color = "year") +
  ggtitle("A. adult abundance and size distribution") +
  theme_minimal() +
  theme(legend.position = "None")
recruit_sizeplot <- ggplot(data = dist_summaries_recruit_long) +
  geom_ribbon(aes(x = as.numeric(size),
                  ymin = lower_ci, ymax = upper_ci,
                  fill = as.factor(year)), alpha = 0.4) +
  geom_line(aes(x = as.numeric(size), y = median,
                color = as.factor(year)),
            linewidth = 1) +
  scale_x_continuous(limits = c(2.5, 30)) +
  scale_color_manual(values = viridis_pal()(4)) +
  scale_fill_manual(values = viridis_pal()(4)) +
  ggtitle("B. recruit abundance and size distribution") +
  labs(x = "size (mm)", y = expression(N[size]), color = "year") +
  theme_minimal() +
  theme(legend.position = "None")
layout <- "
AACC#
AACCE
BBDDE
BBDD#
"
final_plot <- adult_sizeplot + recruit_sizeplot + adult_nplot +
  recruit_nplot + legend + plot_layout(design = layout)
ggsave("figures/Figure3_abundance_sizedist.png",
       dpi = 400, width = 6, height = 5)
