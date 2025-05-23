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
min_size <- 0
max_size <- 110
n <- 22
b <- min_size + c(0:n) * (max_size - min_size) / n
x <- 0.5 * (b[1:n] + b[2:(n + 1)])

# function for getting samples
get_samples <- function(samples, param) {
  samples <- c(samples[[1]][, param], samples[[2]][, param],
               samples[[3]][, param], samples[[4]][, param])
  return(samples)
}
get_survival <- function(size_ind, size_dep, deltat) {
  surv <- exp(-deltat * (size_ind + size_dep / x ^ 2))
  return(surv)
}
get_o_survival <- function(size_dep, n) {
  surv <- exp(-(n * size_dep / x ^ 2))
  return(surv)
}
get_hdi_low <- function(input, ci) {
  out <- as.numeric(hdi(input, ci)[2])
  return(out)
}
get_hdi_high <- function(input, ci) {
  out <- as.numeric(hdi(input, ci)[3])
  return(out)
}

# get natural and overwinter survival
nmort_size <- get_samples(out_sub, "alpha")
survival <- list()
o_survival <- list()
for (i in seq_len(3)) {
  survival[[i]] <- mapply(get_survival,
                          get_samples(out_sub, paste0("beta[", i, "]")),
                          nmort_size, deltat = 14)
  o_survival[[i]] <- mapply(get_o_survival,
                            get_samples(out_sub,
                                        paste0("alpha_o[", i, "]")),
                            get_samples(out_sub,
                                        paste0("wgrowth_N_sum[", i, "]")))
}

# get natural survival summaries
nsurv_summaries <- data.frame(
  year = c("2020", "2021", "2022",
           "2020", "2021", "2022",
           "2020", "2021", "2022"),
  stat = c(rep("median", 3), rep("lower_ci", 3), rep("upper_ci", 3))
)
nsurv_summaries <- cbind(nsurv_summaries,
                         as.data.frame(matrix(NA,
                                              nrow = dim(nsurv_summaries)[1],
                                              ncol = length(x))))
colnames(nsurv_summaries)[3:24] <- x

for (i in seq_len(3)) {
  nsurv_summaries[i, 3:24] <- apply(survival[[i]], 1, median)
  nsurv_summaries[i + 3, 3:24] <- apply(survival[[i]], 1,
                                        function(row) get_hdi_low(row, 0.95))
  nsurv_summaries[i + 6, 3:24] <- apply(survival[[i]], 1,
                                        function(row) get_hdi_high(row, 0.95))
}

# get overwinter survival summaries
wsurv_summaries <- data.frame(
  year = c("2020", "2021", "2022",
           "2020", "2021", "2022",
           "2020", "2021", "2022"),
  stat = c(rep("median", 3), rep("lower_ci", 3), rep("upper_ci", 3))
)
wsurv_summaries <- cbind(wsurv_summaries,
                         as.data.frame(matrix(NA,
                                              nrow = dim(wsurv_summaries)[1],
                                              ncol = length(x))))
colnames(wsurv_summaries)[3:24] <- x

for (i in seq_len(3)) {
  wsurv_summaries[i, 3:24] <- apply(o_survival[[i]], 1, median)
  wsurv_summaries[i + 3, 3:24] <- apply(o_survival[[i]], 1,
                                        function(row) get_hdi_low(row, 0.95))
  wsurv_summaries[i + 6, 3:24] <- apply(o_survival[[i]], 1,
                                        function(row) get_hdi_high(row, 0.95))
}

# convert from wide to long
nsurv_summaries_long <- nsurv_summaries %>%
  pivot_longer(cols = ! c("year", "stat"),
               values_to = "p",
               names_to = "size") %>%
  pivot_wider(values_from = "p",
              names_from = "stat")
wsurv_summaries_long <- wsurv_summaries %>%
  pivot_longer(cols = ! c("year", "stat"),
               values_to = "p",
               names_to = "size") %>%
  pivot_wider(values_from = "p",
              names_from = "stat")

# plot
viridis_colors <- viridis_pal()(4)
nsurv_plot <- ggplot(data = nsurv_summaries_long) +
  geom_ribbon(aes(x = as.numeric(size),
                  ymin = lower_ci, ymax = upper_ci,
                  fill = as.factor(year)), alpha = 0.4,
              show.legend = c(fill = FALSE)) +
  geom_line(aes(x = as.numeric(size), y = median,
                color = as.factor(year)),
            linewidth = 1, show.legend = c(fill = FALSE)) +
  labs(x = "", y = "natural survival", color = "year") +
  ggtitle("A. natural survival") +
  scale_color_manual(values = c(viridis_colors[1],
                                viridis_colors[2],
                                viridis_colors[3])) +
  scale_fill_manual(values = c(viridis_colors[1],
                               viridis_colors[2],
                               viridis_colors[3])) +
  theme_minimal() +
  theme(legend.title = element_text(0.5))
wsurv_plot <- ggplot(data = wsurv_summaries_long) +
  geom_ribbon(aes(x = as.numeric(size),
                  ymin = lower_ci, ymax = upper_ci,
                  fill = as.factor(year)), alpha = 0.4,
              show.legend = c(fill = FALSE)) +
  geom_line(aes(x = as.numeric(size), y = median,
                color = as.factor(year)),
            linewidth = 1, show.legend = c(fill = FALSE)) +
  labs(x = "size (mm)", y = "overwinter survival", color = "year") +
  ggtitle("B. overwinter natural survival") +
  scale_color_manual(labels = c("2020 to 2021",
                                "2021 to 2022",
                                "2022 to 2023"),
                     values = c(viridis_colors[1],
                                viridis_colors[2],
                                viridis_colors[3])) +
  scale_fill_manual(labels = c("2020 to 2021",
                               "2021 to 2022",
                               "2022 to 2023"),
                    values = c(viridis_colors[1],
                               viridis_colors[2],
                               viridis_colors[3])) +
  theme_minimal() +
  theme(legend.title = element_text(0.5))

final_plot <- nsurv_plot + wsurv_plot + plot_layout(ncol = 1)

ggsave("figures/Figure5_survival.png", final_plot,
       dpi = 400, width = 4, height = 6)
