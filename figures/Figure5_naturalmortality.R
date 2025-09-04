library(tidyverse)
library(bayestestR)
library(patchwork)
library(cowplot)
library(viridis)

# read in samples
out <- readRDS("data/posterior_samples/savedsamples_IPM.rds")

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
get_o_survival <- function(alpha_o, n, eps) {
  surv <- exp(-(n * alpha_o / x ^ 2 + eps))
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
survival <- mapply(get_survival,
                   get_samples(out, "beta"),
                   get_samples(out, "alpha"), deltat = 14)
o_survival <- list()
for (i in seq_len(3)) {
  o_survival[[i]] <- mapply(get_o_survival,
                            get_samples(out, "alpha_o"),
                            get_samples(out,
                                        paste0("wgrowth_N_sum[", i, "]")),
                            get_samples(out,
                                        paste0("eps_y[", i, "]")))
}

# get natural survival summaries
nsurv_summaries <- data.frame(
  stat = c("median", "lower_ci", "upper_ci")
)
nsurv_summaries <- cbind(nsurv_summaries,
                         as.data.frame(matrix(NA,
                                              nrow = dim(nsurv_summaries)[1],
                                              ncol = length(x))))
colnames(nsurv_summaries)[2:23] <- x

nsurv_summaries[1, 2:23] <- apply(survival, 1, median)
nsurv_summaries[2, 2:23] <- apply(survival, 1,
                                  function(row) get_hdi_low(row, 0.95))
nsurv_summaries[3, 2:23] <- apply(survival, 1,
                                  function(row) get_hdi_high(row, 0.95))


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
  pivot_longer(cols = ! "stat",
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
                  ymin = lower_ci, ymax = upper_ci),
                  alpha = 0.4) +
  geom_line(aes(x = as.numeric(size), y = median),
            linewidth = 1) +
  labs(x = "size (mm)", y = "natural survival", color = "year") +
  ggtitle("A. April - October") +
  theme_minimal() +
  theme(legend.title = element_text(0.5),
        text = element_text(family = "Arial"))
wsurv_plot <- ggplot(data = wsurv_summaries_long) +
  geom_ribbon(aes(x = as.numeric(size),
                  ymin = lower_ci, ymax = upper_ci,
                  fill = as.factor(year)), alpha = 0.4,
              show.legend = c(fill = FALSE)) +
  geom_line(aes(x = as.numeric(size), y = median,
                color = as.factor(year)),
            linewidth = 1, show.legend = c(fill = FALSE)) +
  labs(x = "size (mm)", y = "overwinter survival", color = "year") +
  ggtitle("B. November - March") +
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
  theme(legend.title = element_text(0.5),
        text = element_text(family = "Arial"))

final_plot <- nsurv_plot + wsurv_plot + plot_layout(nrow = 1)

ggsave("figures/Figure5_survival.png", final_plot,
       dpi = 400, width = 8, height = 3)
ggsave("figures/Figure5_survival.pdf", final_plot,
       device = cairo_pdf,
       width = 8, height = 3)
