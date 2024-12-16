library(tidyverse)
library(bayestestR)

# read in data
data <- readRDS('data/model_data/growth_data.rds')

# read in samples
out <- readRDS('data/posterior_samples/savedsamples_growth.rds')
out_sub <- list(out[[1]][100:10001,],out[[2]][100:10001,],
                out[[3]][100:10001,],out[[4]][100:10001,])

t <- seq(from=0.3,to=3.59,by=0.01)
# function to get length
get_L <- function(C,k,ts,t0,size_inf){
  S_t <- (C*k/(2*pi))*sin(2*pi*(t-ts))
  S_t0 <- (C*k/(2*pi))*sin(2*pi*(t0-ts))
  L <- size_inf*(1-exp(-k*(t-t0)-S_t+S_t0))
  return(L)
}

# function for getting samples
get_samples <- function(samples, param){
  samples <- c(out_sub[[1]][,param],out_sub[[2]][,param],
               out_sub[[3]][,param],out_sub[[4]][,param])
  return(samples)
}

get_hdi_low <- function(input, ci){
  out <- as.numeric(hdi(input,ci)[2])
  return(out)
}
get_hdi_high <- function(input, ci){
  out <- as.numeric(hdi(input,ci)[3])
  return(out)
}

# get params
C <- get_samples(out_sub,'C')
k <- get_samples(out_sub,'k')
ts <- get_samples(out_sub,'ts')
t0 <- get_samples(out_sub,'t0')
size_inf <- get_samples(out_sub,'size_inf')

# get size for each posterior sample
growth_L <- mapply(get_L, C, k, ts, t0, size_inf)

# get summaries
growth_summary <- as.data.frame(matrix(NA,nrow=3,ncol=length(t)))
growth_summary[1,] <- apply(growth_L, 1, median)
growth_summary[2,] <- apply(growth_L, 1, function(row) get_hdi_low(row, 0.95))
growth_summary[3,] <- apply(growth_L, 1, function(row) get_hdi_high(row, 0.95))
colnames(growth_summary) <- t
growth_summary$stat <- c('median','low','high')
growth_summary_long <- growth_summary %>% 
  pivot_longer(cols = -stat,
               names_to = 't',
               values_to = 'value') %>% 
  pivot_wider(values_from='value',
              names_from='stat')

# plot sizes
growth_plot <- ggplot()+
  geom_point(data=data,
             aes(x=age,
                 y=CW),
             size=2,pch=21,alpha=0.8)+
  geom_ribbon(data=growth_summary_long,
              aes(x=as.numeric(t),
                  ymin=low, ymax=high),
              alpha=0.4,fill='indianred2')+
  geom_line(data=growth_summary_long,
            aes(x=as.numeric(t), y=median),
            linewidth=1,color='indianred2')+
  scale_x_continuous(limits=c(0.2,3.6))+
  labs(x='age (years)',
       y='size (mm)')+
  theme_minimal()+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))
ggsave('figures/appendix_growth_plot.png',growth_plot,dpi=400,width=5,height=4)
