library(tidyverse)
library(bayestestR)
library(patchwork)
library(cowplot)


# read in simulation results
out <- readRDS("data/simulations/simulations_effort.rds")
out_noeffort <- readRDS("data/simulations/simulations_noeffort.rds")

# calculate difference
difference <- as.data.frame(matrix(NA, nrow = dim(out)[1],
                                   ncol = length(colnames(out)) - 1))
colnames(difference) <- colnames(out)[1:(length(colnames(out)) - 1)]
for (i in seq_len(nrow(out))) {
  difference[i, 1:22] <- (
    out_noeffort[as.numeric(out[i, "iter"]), ] -
      out[i, 1:22]
  )
  difference[i,
             23:(length(colnames(out)) - 1)] <- out[i,
                                                    23:(length(colnames(out)) -
                                                          1)]
}

# calculate proportion removed
proportion <- as.data.frame(matrix(NA, nrow = dim(out)[1],
                                   ncol = length(colnames(out)) - 1))
colnames(proportion) <- colnames(out)[1:(length(colnames(out)) - 1)]
for (i in seq_len(nrow((out)))) {
  proportion[i, 1:22] <- (
    (out[i, 1:22] + 0.1) /
      (out_noeffort[as.numeric(out[i, "iter"]), ] + 0.1)
  )
  proportion[i,
             23:(length(colnames(out)) - 1)] <- out[i,
                                                    23:(length(colnames(out)) -
                                                          1)]
}

effort_summary <- as.data.frame(proportion) %>%
  pivot_longer(cols = -c("nobs", "prop_s", "prop_m", "prop_f"),
               names_to = "size",
               values_to = "p") %>%
  group_by(size, nobs, prop_s, prop_m, prop_f) %>%
  summarise(median = median(p),
            upper_sd = median + sd(p),
            lower_sd = ifelse(median - sd(p) > 0,
                              median - sd(p), 0)) %>%
  mutate(type = case_when(prop_s == 1 ~ "Shrimp",
                          prop_f == 1 ~ "Fukui",
                          prop_m == 1 ~ "Minnow"))

# create parameters for IPM mesh
min_size <- 0
max_size <- 110
n <- 22 # number of cells in the discretized kernel
b <- min_size + c(0:n) * (max_size - min_size) / n # boundary points
y <- 0.5 * (b[1:n] + b[2:(n + 1)])

plot_28 <- ggplot() +
  geom_ribbon(data = effort_summary[effort_summary$nobs == 28, ],
              aes(x = as.numeric(size),
                  ymin = lower_sd, ymax = upper_sd,
                  fill = type),
              alpha = 0.2) +
  geom_line(data = effort_summary[effort_summary$nobs == 28, ],
            aes(x = as.numeric(size), y = median,
                color = type),
            linewidth = 1) +
  scale_fill_manual(values = c("violet", "darkviolet", "goldenrod"),
                    guide = "none") +
  scale_color_manual(values = c("violet", "darkviolet", "goldenrod")) +
  scale_x_continuous(breaks = c(0, 5.5, 10.5, 15.5, 20.5),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, 5.2)) +
  labs(x = "size (mm)", y = "N(effort = 0):N(effort = 28)",
       color = "trap type") +
  theme_minimal() +
  theme(legend.position = "None")

plot_112 <- ggplot() +
  geom_ribbon(data = effort_summary[effort_summary$nobs == 112, ],
              aes(x = as.numeric(size),
                  ymin = lower_sd, ymax = upper_sd,
                  fill = type),
              alpha = 0.2) +
  geom_line(data = effort_summary[effort_summary$nobs == 112, ],
            aes(x = as.numeric(size), y = median,
                color = type),
            linewidth = 1) +
  scale_fill_manual(values = c("violet", "darkviolet", "goldenrod"),
                    guide = "none") +
  scale_color_manual(values = c("violet", "darkviolet", "goldenrod")) +
  scale_x_continuous(breaks = c(0, 5.5, 10.5, 15.5, 20.5),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, 5.2)) +
  labs(x = "size (mm)", y = "N(effort = 0):N(effort = 112)",
       color = "trap type") +
  theme_minimal() +
  theme(legend.position = "None",
        plot.title = element_text(size = 12)) +
  ggtitle("A. effort = 112")

plot_560 <- ggplot() +
  geom_ribbon(data = effort_summary[effort_summary$nobs == 560, ],
              aes(x = as.numeric(size),
                  ymin = lower_sd, ymax = upper_sd,
                  fill = type),
              alpha = 0.2) +
  geom_line(data = effort_summary[effort_summary$nobs == 560, ],
            aes(x = as.numeric(size), y = median,
                color = type),
            linewidth = 1) +
  scale_fill_manual(values = c("violet", "darkviolet", "goldenrod"),
                    guide = "none") +
  scale_color_manual(values = c("violet", "darkviolet", "goldenrod")) +
  scale_x_continuous(breaks = c(0, 5.5, 10.5, 15.5, 20.5),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, 5.2)) +
  labs(x = "size (mm)", y = "N(effort = 0):N(effort = 560)",
       color = "trap type") +
  theme_minimal() +
  theme(legend.position = "None",
        plot.title = element_text(size = 12)) +
  ggtitle("C. effort = 560")

plot_840 <- ggplot() +
  geom_ribbon(data = effort_summary[effort_summary$nobs == 840, ],
              aes(x = as.numeric(size),
                  ymin = lower_sd, ymax = upper_sd,
                  fill = type),
              alpha = 0.2) +
  geom_line(data = effort_summary[effort_summary$nobs == 840, ],
            aes(x = as.numeric(size), y = median,
                color = type),
            linewidth = 1) +
  scale_fill_manual(values = c("violet", "darkviolet", "goldenrod"),
                    guide = "none") +
  scale_color_manual(values = c("violet", "darkviolet", "goldenrod")) +
  scale_x_continuous(breaks = c(0, 5.5, 10.5, 15.5, 20.5),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, 5.2)) +
  labs(x = "size (mm)", y = "N(effort = 0):N(effort = 840)",
       color = "trap type") +
  theme_minimal() +
  theme(legend.position = "None")

plot_1400 <- ggplot() +
  geom_ribbon(data = effort_summary[effort_summary$nobs == 1400, ],
              aes(x = as.numeric(size),
                  ymin = lower_sd, ymax = upper_sd,
                  fill = type),
              alpha = 0.2) +
  geom_line(data = effort_summary[effort_summary$nobs == 1400, ],
            aes(x = as.numeric(size), y = median,
                color = type),
            linewidth = 1) +
  scale_fill_manual(values = c("violet", "darkviolet", "goldenrod"),
                    guide = "none") +
  scale_color_manual(values = c("violet", "darkviolet", "goldenrod")) +
  scale_x_continuous(breaks = c(0, 5.5, 10.5, 15.5, 20.5),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, 5.2)) +
  labs(x = "size (mm)", y = "N(effort = 0):N(effort = 1400)",
       color = "trap type") +
  theme_minimal() +
  theme(legend.position = "None",
        plot.title = element_text(size = 12)) +
  ggtitle("D. effort = 1400")

plot_2800 <- ggplot() +
  geom_ribbon(data = effort_summary[effort_summary$nobs == 2800, ],
              aes(x = as.numeric(size),
                  ymin = lower_sd, ymax = upper_sd,
                  fill = type),
              alpha = 0.2) +
  geom_line(data = effort_summary[effort_summary$nobs == 2800, ],
            aes(x = as.numeric(size), y = median,
                color = type),
            linewidth = 1) +
  scale_fill_manual(values = c("violet", "darkviolet", "goldenrod"),
                    guide = "none") +
  scale_color_manual(values = c("violet", "darkviolet", "goldenrod")) +
  scale_x_continuous(breaks = c(0, 5.5, 10.5, 15.5, 20.5),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, 5.2)) +
  labs(x = "size (mm)", y = "N(effort = 0):N(effort = 2800)",
       color = "trap type") +
  theme_minimal() +
  theme(legend.position = "None",
        plot.title = element_text(size = 12)) +
  ggtitle("D. effort = 2800")

plot_w_legend <- ggplot() +
  geom_ribbon(data = effort_summary[effort_summary$nobs == 1400, ],
              aes(x = as.numeric(size),
                  ymin = lower_sd2, ymax = upper_sd2,
                  fill = type),
              alpha = 0.2) +
  geom_line(data = effort_summary[effort_summary$nobs == 1400, ],
            aes(x = as.numeric(size), y = median,
                color = type),
            linewidth = 1) +
  scale_fill_manual(values = c("violet", "darkviolet", "goldenrod"),
                    guide = "none") +
  scale_color_manual(values = c("violet", "darkviolet", "goldenrod")) +
  scale_x_continuous(breaks = c(0, 5.5, 10.5, 15.5, 20.5),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(x = "size (mm)", y = expression(N[size]),
       color = "trap type") +
  theme_minimal()
legend <- get_legend(plot_w_legend)

final_plot <- plot_112 + plot_560 + plot_2800 + legend +
  plot_layout(nrow = 1)

ggsave("figures/supplemental_IPM_simulations_diff.png", dpi = 400,
       width = 7, height = 3)
