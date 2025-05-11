library(tidyverse)
library(bayestestR)
library(patchwork)
library(cowplot)


# read in simulation results
out <- readRDS("data/simulations/simulations_effort.rds")
out_noeffort <- readRDS("data/simulations/simulations_noeffort.rds")

out_no_effort_summary <- as.data.frame(out_noeffort) %>%
  mutate(iter = seq_len(nrow(out_noeffort))) %>%
  pivot_longer(cols = ! iter,
               names_to = "size",
               values_to = "N") %>%
  group_by(size) %>%
  summarise(median = median(N),
            upper_sd = median + sd(N),
            lower_sd = ifelse(median - sd(N) > 0,
                              median - sd(N), 0))

effort_summary <- as.data.frame(out) %>%
  select(-iter) %>%
  pivot_longer(cols = -c("nobs", "prop_s", "prop_m", "prop_f"),
               names_to = "size",
               values_to = "N") %>%
  group_by(size, nobs, prop_s, prop_m, prop_f) %>%
  summarise(median = median(N),
            upper_sd = median + sd(N),
            lower_sd = ifelse(median - sd(N) > 0,
                              median - sd(N), 0)) %>%
  mutate(type = case_when(prop_s == 1 ~ "Shrimp",
                          prop_f == 1 ~ "Fukui",
                          prop_m == 1 ~ "Minnow"))

# create parameters for IPM mesh
min_size <- 0
max_size <- 110
n <- 22 # number of cells in the discretized kernel
b <- min_size + c(0:n) * (max_size - min_size) / n # boundary points
x <- 0.5 * (b[1:n] + b[2:(n + 1)])

# no effort plot
plot_0 <- ggplot(data = out_no_effort_summary) +
  geom_ribbon(aes(x = as.numeric(size),
                  ymin = lower_sd, ymax = upper_sd),
              alpha = 0.2) +
  geom_line(aes(x = as.numeric(size), y = median),
            linewidth = 1) +
  scale_x_continuous(breaks = c(0, 5.5, 10.5, 15.5, 20.5),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(x = "size (mm)", y = expression(N(x)^E)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12)) +
  ggtitle("A. effort = 0")

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
  scale_y_continuous(limits = c(0, 50)) +
  labs(x = "size (mm)", y = NULL,
       color = "trap type") +
  theme_minimal() +
  theme(legend.position = "None",
        plot.title = element_text(size = 12)) +
  ggtitle("B. effort = 112")

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
  scale_y_continuous(limits = c(0, 50)) +
  labs(x = "size (mm)", y = NULL,
       color = "trap type") +
  theme_minimal() +
  theme(legend.position = "None",
        plot.title = element_text(size = 12)) +
  ggtitle("C. effort = 560")

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
  scale_y_continuous(limits = c(0, 50)) +
  labs(x = "size (mm)", y = NULL,
       color = "trap type") +
  theme_minimal() +
  theme(legend.position = "None",
        plot.title = element_text(size = 12)) +
  ggtitle("D. effort = 2800")

plot_w_legend <- ggplot() +
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
  scale_y_continuous(limits = c(0, 70)) +
  labs(x = "size (mm)", y = expression(N[size]),
       color = "trap type") +
  theme_minimal()
legend <- get_legend(plot_w_legend)

final_plot <- plot_0 + plot_112 + plot_560 + plot_2800 + legend +
  plot_layout(nrow = 1)
ggsave("figures/Figure6_IPM_simulations.png", dpi = 400,
       width = 9, height = 3)
