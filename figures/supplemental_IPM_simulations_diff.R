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
x <- 0.5 * (b[1:n] + b[2:(n + 1)])

# subset to 112, 560, 2800
effort_summary_sub <- effort_summary[effort_summary$nobs %in% 
                                       c(112, 560, 2800), ]

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
  scale_y_continuous(limits = c(0, max(effort_summary_sub$upper_sd))) +
  labs(x = "size (mm)", 
       y = expression(N^E *(x)[paste("effort=0")] ~ ":" 
                      ~ N^E * (x)[paste("effort=Z")]),
       color = "Z (number\nof traps)") +
  theme_minimal() +
  theme(legend.position = "None",
        plot.title = element_text(size = 12),
        text = element_text(family = "Arial")) +
  ggtitle("A. Fukui traps")

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
  scale_y_continuous(limits = c(0, max(effort_summary_sub$upper_sd))) +
  labs(x = "size (mm)", 
       y = expression(N^E *(x)[paste("effort=0")] ~ ":" 
                      ~ N^E * (x)[paste("effort=Z")]),
       color = "Z (number\nof traps)") +
  theme_minimal() +
  theme(legend.position = "None",
        plot.title = element_text(size = 12),
        text = element_text(family = "Arial")) +
  ggtitle("B. Minnow traps")

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
  scale_y_continuous(limits = c(0, max(effort_summary_sub$upper_sd))) +
  labs(x = "size (mm)", 
       y = expression(N^E *(x)[paste("effort=0")] ~ ":" 
                      ~ N^E * (x)[paste("effort=Z")]),
       color = "Z (number\nof traps)") +
  theme_minimal() +
  theme(legend.title = element_text(hjust = 0.5),
        plot.title = element_text(size = 12),
        text = element_text(family = "Arial")) +
  ggtitle("C. Shrimp traps")

final_plot <- plot_fukui + plot_minnow + plot_shrimp +
  plot_layout(nrow = 1)

ggsave("figures/supplemental_IPM_simulations_diff.png", dpi = 400,
       width = 8, height = 3.5)
