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

# subset to 112, 560, 2800
effort_summary_sub <- effort_summary[effort_summary$nobs %in% 
                                       c(112, 560, 2800), ]

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
  scale_y_continuous(limits = c(0, max(out_no_effort_summary$upper_sd))) +
  labs(x = "size (mm)", y = expression(N^E * (x))) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12),
        text = element_text(family = "Arial")) +
  ggtitle("A. None")

plot_fukui <- ggplot() +
  geom_ribbon(data = effort_summary_sub[effort_summary_sub$type == "Fukui", ],
              aes(x = as.numeric(size),
                  ymin = lower_sd, ymax = upper_sd,
                  fill = as.factor(nobs)),
              alpha = 0.2) +
  geom_line(data = effort_summary_sub[effort_summary_sub$type == "Fukui", ],
            aes(x = as.numeric(size), y = median,
                color = as.factor(nobs)),
            linewidth = 1) +
  scale_fill_manual(values = c("paleturquoise4", "dodgerblue", "blue4"),
                    guide = "none") +
  scale_color_manual(values = c("paleturquoise4", "dodgerblue", "blue4")) +
  scale_x_continuous(breaks = c(0, 5.5, 10.5, 15.5, 20.5),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, max(out_no_effort_summary$upper_sd))) +
  labs(x = "size (mm)", y = NULL,
       color = "removal effort\n(traps)") +
  theme_minimal() +
  theme(legend.position = "None",
        plot.title = element_text(size = 12),
        text = element_text(family = "Arial")) +
  ggtitle("B. Fukui traps")

plot_minnow <- ggplot() +
  geom_ribbon(data = effort_summary_sub[effort_summary_sub$type == "Minnow", ],
              aes(x = as.numeric(size),
                  ymin = lower_sd, ymax = upper_sd,
                  fill = as.factor(nobs)),
              alpha = 0.2) +
  geom_line(data = effort_summary_sub[effort_summary_sub$type == "Minnow", ],
            aes(x = as.numeric(size), y = median,
                color = as.factor(nobs)),
            linewidth = 1) +
  scale_fill_manual(values = c("paleturquoise4", "dodgerblue", "blue4"),
                    guide = "none") +
  scale_color_manual(values = c("paleturquoise4", "dodgerblue", "blue4")) +
  scale_x_continuous(breaks = c(0, 5.5, 10.5, 15.5, 20.5),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, max(out_no_effort_summary$upper_sd))) +
  labs(x = "size (mm)", y = NULL,
       color = "removal effort\n(traps)") +
  theme_minimal() +
  theme(legend.position = "None",
        plot.title = element_text(size = 12),
        text = element_text(family = "Arial")) +
  ggtitle("C. Minnow traps")

plot_shrimp <- ggplot() +
  geom_ribbon(data = effort_summary_sub[effort_summary_sub$type == "Shrimp", ],
              aes(x = as.numeric(size),
                  ymin = lower_sd, ymax = upper_sd,
                  fill = as.factor(nobs)),
              alpha = 0.2) +
  geom_line(data = effort_summary_sub[effort_summary_sub$type == "Shrimp", ],
            aes(x = as.numeric(size), y = median,
                color = as.factor(nobs)),
            linewidth = 1) +
  scale_fill_manual(values = c("paleturquoise4", "dodgerblue", "blue4"),
                    guide = "none") +
  scale_color_manual(values = c("paleturquoise4", "dodgerblue", "blue4")) +
  scale_x_continuous(breaks = c(0, 5.5, 10.5, 15.5, 20.5),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, max(out_no_effort_summary$upper_sd))) +
  labs(x = "size (mm)", y = NULL,
       color = "removal effort\n(traps)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12),
        legend.title = element_text(hjust = 0.5),
        text = element_text(family = "Arial")) +
  ggtitle("D. Shrimp traps")

final_plot <- plot_0 + plot_fukui + plot_minnow + plot_shrimp +
  plot_layout(nrow = 1)
ggsave("figures/Figure6_IPM_simulations.png", dpi = 400,
       width = 9, height = 3)
