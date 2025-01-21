library(tidyverse)
library(bayestestR)
library(patchwork)
library(cowplot)
library(viridis)

# read in samples
out <- readRDS('data/posterior_samples/savedsamples_IPM.rds')

# subset samples
lower <- 2000
upper <- 10001
out_sub <- list(
  out[[1]][lower:upper,],out[[2]][lower:upper,],
  out[[3]][lower:upper,],out[[4]][lower:upper,]
)

min.size <- 0
max.size <- 110
n <- 22 
b <- min.size+c(0:n)*(max.size-min.size)/n
y <- 0.5*(b[1:n]+b[2:(n+1)])

# function for getting samples
get_samples <- function(samples, param){
  samples <- c(samples[[1]][,param],samples[[2]][,param],
               samples[[3]][,param],samples[[4]][,param])
  return(samples)
}

get_survival <- function(size_ind,size_dep,deltat){
  surv <- exp(-(size_ind+size_dep/y^2))
  return(surv)
}
get_o_survival <- function(size_dep,N){
  surv <- exp(-(N*size_dep/y^2))
  return(surv)
}

get_hdi_low <- function(input, ci){
  out <- as.numeric(hdi(input,ci)[2])
  return(out)
}
get_hdi_high <- function(input, ci){
  out <- as.numeric(hdi(input,ci)[3])
  return(out)
}

wnmortality_s1 <- get_samples(out_sub,'wnmortality_s[1]')
wnmortality_s2 <- get_samples(out_sub,'wnmortality_s[2]')
wnmortality_s3 <- get_samples(out_sub,'wnmortality_s[3]')
wnmort_size <- get_samples(out_sub,'wnmort_size')

wsurvival1 <- mapply(get_survival, wnmortality_s1,
                     wnmort_size)
wsurvival2 <- mapply(get_survival, wnmortality_s2,
                     wnmort_size)
wsurvival3 <- mapply(get_survival, wnmortality_s3,
                     wnmort_size)

# get summaries
wmort_summaries <- data.frame(
  year = c('2020','2021','2022',
           '2020','2021','2022',
           '2020','2021','2022'),
  stat = c(rep('median',3),rep('lower_ci',3),rep('upper_ci',3))
)
wmort_summaries <- cbind(wmort_summaries,
                         as.data.frame(matrix(NA,
                                              nrow=dim(wmort_summaries)[1],
                                              ncol=length(y))))
colnames(wmort_summaries)[3:24] <- y

wmort_summaries[1,3:24] <- apply(wsurvival1, 1, median)
wmort_summaries[2,3:24] <- apply(wsurvival2, 1, median)
wmort_summaries[3,3:24] <- apply(wsurvival3, 1, median)
wmort_summaries[4,3:24] <- apply(wsurvival1, 1, function(row) get_hdi_low(row, 0.95))
wmort_summaries[5,3:24] <- apply(wsurvival2, 1, function(row) get_hdi_low(row, 0.95))
wmort_summaries[6,3:24] <- apply(wsurvival3, 1, function(row) get_hdi_low(row, 0.95))
wmort_summaries[7,3:24] <- apply(wsurvival1, 1, function(row) get_hdi_high(row, 0.95))
wmort_summaries[8,3:24] <- apply(wsurvival2, 1, function(row) get_hdi_high(row, 0.95))
wmort_summaries[9,3:24] <- apply(wsurvival3, 1, function(row) get_hdi_high(row, 0.95))

# convert from wide to long 
wmort_summaries_long <- wmort_summaries %>% 
  pivot_longer(cols=!c('year','stat'),
               values_to='p',
               names_to='size') %>% 
  pivot_wider(values_from='p',
              names_from='stat')

viridis_colors <- viridis_pal()(4)

final_plot <- ggplot(data = wmort_summaries_long) +
  geom_ribbon(aes(x=as.numeric(size),
                  ymin=lower_ci, ymax=upper_ci,
                  fill=as.factor(year)),alpha=0.4,
              show.legend = c(fill = FALSE))+
  geom_line(aes(x=as.numeric(size), y=median,
                color=as.factor(year)),
            linewidth=1,show.legend = c(fill = FALSE))+
  labs(x = "size (mm)", y = "overwinter survival", color='year') +
  scale_color_manual(labels=c('2020 → 2021',
                              '2021 → 2022',
                              '2022 → 2023'),
                     values=c(viridis_colors[1],
                              viridis_colors[2],
                              viridis_colors[3]))+
  scale_fill_manual(labels=c('2020 → 2021',
                              '2021 → 2022',
                              '2022 → 2023'),
                     values=c(viridis_colors[1],
                              viridis_colors[2],
                              viridis_colors[3]))+
  theme_minimal()
ggsave('figures/fig_wmort.png',final_plot,dpi=400,width=4,height=3)
