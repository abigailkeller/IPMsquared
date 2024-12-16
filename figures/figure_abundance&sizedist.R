library(tidyverse)
library(bayestestR)
library(patchwork)
library(cowplot)
library(viridis)

# read in samples
out <- readRDS('IPM/20241030b_savedsamples.rds')

# subset samples
lower <- 2000
upper <- 10001
out_sub <- list(
  out[[1]][lower:upper,],out[[2]][lower:upper,],
  out[[3]][lower:upper,],out[[4]][lower:upper,]
)

# function for getting samples
get_samples <- function(samples, param){
  samples <- c(samples[[1]][,param],samples[[2]][,param],
               samples[[3]][,param],samples[[4]][,param])
  return(samples)
}

# get params
N_adult1 <- get_samples(out_sub,'N_adult')
# N_adult2 <- get_samples(out_sub,'N_adult[2]')
# N_adult3 <- get_samples(out_sub,'N_adult[3]')
# N_adult4 <- get_samples(out_sub,'N_adult[4]')
N_recruit1 <- get_samples(out_sub,'N_recruit[1]')
N_recruit2 <- get_samples(out_sub,'N_recruit[2]')
N_recruit3 <- get_samples(out_sub,'N_recruit[3]')
N_recruit4 <- get_samples(out_sub,'N_recruit[4]')
init_mean_recruit <- get_samples(out_sub,'init_mean_recruit')
init_sd_r <- get_samples(out_sub,'init_sd_r')
init_lsd_adult <- get_samples(out_sub,'init_lsd_adult')
init_lmean_adult <- get_samples(out_sub,'init_lmean_adult')
# init_lmean_adult2 <- get_samples(out_sub,'init_lmean_adult[2]')
# init_lmean_adult3 <- get_samples(out_sub,'init_lmean_adult[3]')
# init_lmean_adult4 <- get_samples(out_sub,'init_lmean_adult[4]')
adult_dist2 <- as.data.frame(matrix(NA,ncol=22,nrow=dim(out_sub[[1]])[1]*4))
adult_dist3 <- as.data.frame(matrix(NA,ncol=22,nrow=dim(out_sub[[1]])[1]*4))
adult_dist4 <- as.data.frame(matrix(NA,ncol=22,nrow=dim(out_sub[[1]])[1]*4))
for(i in 1:22){
  adult_dist2[,i] <- get_samples(out_sub, paste0('N_overwinter[1, ',i,']'))
  adult_dist3[,i] <- get_samples(out_sub, paste0('N_overwinter[2, ',i,']'))
  adult_dist4[,i] <- get_samples(out_sub, paste0('N_overwinter[3, ',i,']'))
}
N_adult2 <- rowSums(adult_dist2)
N_adult3 <- rowSums(adult_dist3)
N_adult4 <- rowSums(adult_dist4)


# make plots of total N
N_summaries <- data.frame(
  param = c('N_adult1','N_adult2','N_adult3','N_adult4',
            'N_recruit1','N_recruit2','N_recruit3','N_recruit4'),
  type = c('adult','adult','adult','adult','recruit','recruit','recruit','recruit'),
  year = c('2020','2021','2022','2023',
           '2020','2021','2022','2023'),
  median = c(median(N_adult1),median(N_adult2),median(N_adult3),median(N_adult4),
             median(N_recruit1),median(N_recruit2),
             median(N_recruit3),median(N_recruit4)),
  lower_CI = c(as.numeric(eti(N_adult1)[2]),as.numeric(eti(N_adult2)[2]),
               as.numeric(eti(N_adult3)[2]),as.numeric(eti(N_adult4)[2]),
               as.numeric(eti(N_recruit1)[2]),as.numeric(eti(N_recruit2)[2]),
               as.numeric(eti(N_recruit3)[2]),as.numeric(eti(N_recruit4)[2])),
  upper_CI = c(as.numeric(eti(N_adult1)[3]),as.numeric(eti(N_adult2)[3]),
               as.numeric(eti(N_adult3)[3]),as.numeric(eti(N_adult4)[3]),
               as.numeric(eti(N_recruit1)[3]),as.numeric(eti(N_recruit2)[3]),
               as.numeric(eti(N_recruit3)[3]),as.numeric(eti(N_recruit4)[3]))
)

# make plots of N distribution - recruit
min.size <- 0
max.size <- 110
y_recruit <- seq(min.size,max.size,1)
n <- 22 
b <- min.size+c(0:n)*(max.size-min.size)/n
y <- 0.5*(b[1:n]+b[2:(n+1)])

init_mean_recruit <- get_samples(out_sub,'init_mean_recruit')
init_sd_r <- get_samples(out_sub,'init_sd_r')
init_lsd_adult <- get_samples(out_sub,'init_lsd_adult')
init_lmean_adult1 <- get_samples(out_sub,'init_lmean_adult')
# init_lmean_adult2 <- get_samples(out_sub,'init_lmean_adult[2]')
# init_lmean_adult3 <- get_samples(out_sub,'init_lmean_adult[3]')
# init_lmean_adult4 <- get_samples(out_sub,'init_lmean_adult[4]')


dist_summaries_recruit <- data.frame(
  param = c('N_recruit1_dist','N_recruit2_dist','N_recruit3_dist','N_recruit4_dist',
            'N_recruit1_dist','N_recruit2_dist','N_recruit3_dist','N_recruit4_dist',
            'N_recruit1_dist','N_recruit2_dist','N_recruit3_dist','N_recruit4_dist'),
  year = c('2020','2021','2022','2023',
           '2020','2021','2022','2023',
           '2020','2021','2022','2023'),
  stat = c(rep('median',4),rep('lower_ci',4),rep('upper_ci',4))
) 

dist_summaries_recruit <- cbind(dist_summaries_recruit, 
                                as.data.frame(matrix(NA,
                                                     nrow=dim(dist_summaries_recruit)[1],
                                                     ncol=length(y))))
colnames(dist_summaries_recruit)[4:25] <- y

dist_summaries_adult <- data.frame(
  param = c('N_adult1_dist','N_adult2_dist','N_adult3_dist','N_adult4_dist',
            'N_adult1_dist','N_adult2_dist','N_adult3_dist','N_adult4_dist',
            'N_adult1_dist','N_adult2_dist','N_adult3_dist','N_adult4_dist'),
  year = c('2020','2021','2022','2023',
           '2020','2021','2022','2023',
           '2020','2021','2022','2023'),
  stat = c(rep('median',4),rep('lower_ci',4),rep('upper_ci',4))
) 

dist_summaries_adult <- cbind(dist_summaries_adult,
                              as.data.frame(matrix(NA,
                                                   nrow=dim(dist_summaries_adult),
                                                   ncol=length(y))))
colnames(dist_summaries_adult)[4:25] <- y

get_dist_log <- function(mean,sd,N){
  out <- dlnorm(y,mean,sd)/sum(dlnorm(y,mean,sd))*N
  return(out)
}
get_dist <- function(mean,sd,N){
  out <- dnorm(y,mean,sd)/sum(dnorm(y,mean,sd))*N
  return(out)
}
get_eti_low <- function(input, ci){
  out <- as.numeric(eti(input,ci)[2])
  return(out)
}
get_eti_high <- function(input, ci){
  out <- as.numeric(eti(input,ci)[3])
  return(out)
}

N_adult1_dist <- mapply(get_dist_log, init_lmean_adult1,
                        init_lsd_adult, N_adult1)
# N_adult2_dist <- mapply(get_dist_log, init_lmean_adult2,
#                         init_lsd_adult, N_adult2)
# N_adult3_dist <- mapply(get_dist_log, init_lmean_adult3,
#                         init_lsd_adult, N_adult3)
# N_adult4_dist <- mapply(get_dist_log, init_lmean_adult4,
#                         init_lsd_adult, N_adult4)
N_recruit1_dist <- mapply(get_dist, init_mean_recruit,
                          init_sd_r, N_recruit1)
N_recruit2_dist <- mapply(get_dist, init_mean_recruit,
                          init_sd_r, N_recruit2)
N_recruit3_dist <- mapply(get_dist, init_mean_recruit,
                          init_sd_r, N_recruit3)
N_recruit4_dist <- mapply(get_dist, init_mean_recruit,
                          init_sd_r, N_recruit4)

dist_summaries_adult[1,4:25] <- apply(N_adult1_dist, 1, median)
dist_summaries_adult[2,4:25] <- apply(adult_dist2, 2, median)
dist_summaries_adult[3,4:25] <- apply(adult_dist3, 2, median)
dist_summaries_adult[4,4:25] <- apply(adult_dist4, 2, median)
dist_summaries_adult[5,4:25] <- apply(N_adult1_dist, 1, function(row) get_eti_low(row, 0.95))
dist_summaries_adult[6,4:25] <- apply(adult_dist2, 2, function(row) get_eti_low(row, 0.95))
dist_summaries_adult[7,4:25] <- apply(adult_dist3, 2, function(row) get_eti_low(row, 0.95))
dist_summaries_adult[8,4:25] <- apply(adult_dist4, 2, function(row) get_eti_low(row, 0.95))
dist_summaries_adult[9,4:25] <- apply(N_adult1_dist, 1, function(row) get_eti_high(row, 0.95))
dist_summaries_adult[10,4:25] <- apply(adult_dist2, 2, function(row) get_eti_high(row, 0.95))
dist_summaries_adult[11,4:25] <- apply(adult_dist3, 2, function(row) get_eti_high(row, 0.95))
dist_summaries_adult[12,4:25] <- apply(adult_dist4, 2, function(row) get_eti_high(row, 0.95))

dist_summaries_recruit[1,4:25] <- apply(N_recruit1_dist, 1, function(x) median(x, na.rm = TRUE))
dist_summaries_recruit[2,4:25] <- apply(N_recruit2_dist, 1, function(x) median(x, na.rm = TRUE))
dist_summaries_recruit[3,4:25] <- apply(N_recruit3_dist, 1, function(x) median(x, na.rm = TRUE))
dist_summaries_recruit[4,4:25] <- apply(N_recruit4_dist, 1, function(x) median(x, na.rm = TRUE))
dist_summaries_recruit[5,4:25] <- apply(N_recruit1_dist, 1, function(row) get_eti_low(row, 0.5))
dist_summaries_recruit[6,4:25] <- apply(N_recruit2_dist, 1, function(row) get_eti_low(row, 0.5))
dist_summaries_recruit[7,4:25] <- apply(N_recruit3_dist, 1, function(row) get_eti_low(row, 0.5))
dist_summaries_recruit[8,4:25] <- apply(N_recruit4_dist, 1, function(row) get_eti_low(row, 0.5))
dist_summaries_recruit[9,4:25] <- apply(N_recruit1_dist, 1, function(row) get_eti_high(row, 0.5))
dist_summaries_recruit[10,4:25] <- apply(N_recruit2_dist, 1, function(row) get_eti_high(row, 0.5))
dist_summaries_recruit[11,4:25] <- apply(N_recruit3_dist, 1, function(row) get_eti_high(row, 0.5))
dist_summaries_recruit[12,4:25] <- apply(N_recruit4_dist, 1, function(row) get_eti_high(row, 0.5))


# convert from wide to long 
dist_summaries_adult_long <- dist_summaries_adult %>% 
  pivot_longer(cols=!c('param','year','stat'),
               values_to='N',
               names_to='size') %>% 
  pivot_wider(values_from='N',
              names_from='stat')
dist_summaries_recruit_long <- dist_summaries_recruit %>% 
  pivot_longer(cols=!c('param','year','stat'),
               values_to='N',
               names_to='size') %>% 
  pivot_wider(values_from='N',
              names_from='stat')

##### total abundance plots
adult_Nplot <- ggplot(data=N_summaries[N_summaries$type=='adult',])+
  geom_point(aes(x=param, y=median,
                 color=as.factor(year)))+
  geom_errorbar(aes(x=param,
                    ymin = lower_CI, 
                    ymax = upper_CI,
                    color=as.factor(year)),
                width = 0.2,
                linewidth=1)+
  scale_color_manual(values=viridis_pal()(4))+
  labs(x='year',y=expression(N[total]),color='year')+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        legend.position='None')
recruit_Nplot <- ggplot(data=N_summaries[N_summaries$type=='recruit',])+
  geom_point(aes(x=param, y=log(median),
                 color=as.factor(year)))+
  geom_errorbar(aes(x=param,
                    ymin = log(lower_CI), 
                    ymax = log(upper_CI),
                    color=as.factor(year)),
                width = 0.2,
                linewidth=1)+
  scale_y_continuous(breaks=c(log(20),log(200)),
                     labels=c('20','200'))+
  scale_color_manual(values=viridis_pal()(4))+
  labs(x='year',y=expression(N[total]),color='year')+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        legend.position='None')

plot_w_legend <- ggplot(data=N_summaries[N_summaries$type=='recruit',])+
  geom_point(aes(x=param, y=log(median),
                 color=as.factor(year)))+
  geom_errorbar(aes(x=param,
                    ymin = log(lower_CI), 
                    ymax = log(upper_CI),
                    color=as.factor(year)),
                width = 0.2,
                linewidth=1)+
  scale_color_manual(values=viridis_pal()(4))+
  scale_y_continuous(breaks=c(log(20),log(200)),
                     labels=c('20','200'))+
  labs(x='year',y=expression(N[total]),color='year')+
  theme_minimal()+
  theme(axis.text.x=element_blank())
legend <- get_legend(plot_w_legend)


######## size abundance plots
adult_sizeplot <- ggplot(data=dist_summaries_adult_long)+
  geom_ribbon(aes(x=as.numeric(size),
                  ymin=lower_ci, ymax=upper_ci,
                  fill=as.factor(year)),alpha=0.4)+
  geom_line(aes(x=as.numeric(size), y=median,
                color=as.factor(year)),
            linewidth=1)+
  scale_x_continuous(limits=c(2.5,107.5))+
  scale_color_manual(values=viridis_pal()(4))+
  scale_fill_manual(values=viridis_pal()(4))+
  labs(x='size (mm)',y=expression(N[size]),color='year')+
  ggtitle('A. adult abundance and size distribution')+
  theme_minimal()+
  theme(legend.position = 'None')
recruit_sizeplot <- ggplot(data=dist_summaries_recruit_long)+
  geom_ribbon(aes(x=as.numeric(size),
                  ymin=lower_ci, ymax=upper_ci,
                  fill=as.factor(year)),alpha=0.4)+
  geom_line(aes(x=as.numeric(size), y=median,
                color=as.factor(year)),
            linewidth=1)+
  scale_x_continuous(limits=c(2.5,30))+
  scale_color_manual(values=viridis_pal()(4))+
  scale_fill_manual(values=viridis_pal()(4))+
  ggtitle('B. recruit abundance and size distribution')+
  labs(x='size (mm)',y=expression(N[size]),color='year')+
  theme_minimal()+
  theme(legend.position = 'None')

layout <- '
AACC#
AACCE
BBDDE
BBDD#
'

final_plot <- adult_sizeplot+recruit_sizeplot+adult_Nplot+recruit_Nplot+legend+plot_layout(design=layout)
ggsave('figures/fig_abundance_sizedist.png',dpi=400,width=6,height=5)
