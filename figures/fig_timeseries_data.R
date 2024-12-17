library(tidyverse)
library(cowplot)
library(patchwork)

# read in data
DH2020_catch <- readRDS('data/catch_data/DH2020_catch.rds')
DH2021_catch <- readRDS('data/catch_data/DH2021_catch.rds')
DH2022_catch <- readRDS('data/catch_data/DH2022_catch.rds')
DH2023_catch <- readRDS('data/catch_data/DH2023_catch.rds')

# plots
plot_2020 <- ggplot()+
  geom_point(data=DH2020_catch,
             aes(x=Date, y=as.numeric(Size_mm),
                 color=Trap_Type))+
  scale_color_manual(values=c('violet','darkviolet'))+
  labs(x='date',y='size (mm)',color='trap type')+
  scale_x_date(limits=c(as.Date("2020-04-14"),as.Date("2020-10-26")),
               breaks=c(as.Date('2020-05-01'),as.Date('2020-07-01'),
                        as.Date('2020-09-01'),as.Date('2020-11-01')),
               labels=c('May','July','Sept.','Nov.'))+
  scale_y_continuous(limits=c(30,95))+
  ggtitle('i. 2020')+
  theme_minimal()+
  theme(legend.position='None',
        plot.title=element_text(size=12))

plot_2021 <- ggplot()+
  geom_point(data=DH2021_catch,
             aes(x=Date_Observed, y=as.numeric(CW_mm),
                 color=Trap_Type))+
  scale_color_manual(values=c('violet','darkviolet','goldenrod'))+
  labs(x='date',y='',color='trap type')+
  scale_x_date(limits=c(as.Date("2021-04-14"),as.Date("2021-10-26")),
               breaks=c(as.Date('2021-05-01'),as.Date('2021-07-01'),
                        as.Date('2021-09-01'),as.Date('2021-11-01')),
               labels=c('May','July','Sept.','Nov.'))+
  scale_y_continuous(limits=c(30,95))+
  ggtitle('ii. 2021')+
  theme_minimal()+
  theme(legend.position='None',
        plot.title=element_text(size=12)) 

plot_2022 <- ggplot()+
  geom_point(data=DH2022_catch,
             aes(x=Date_Observed, y=as.numeric(CW_mm),
                 color=Trap_Type))+
  scale_color_manual(values=c('violet','darkviolet','goldenrod'))+
  labs(x='date',y='',color='trap type')+
  scale_x_date(limits=c(as.Date("2022-04-14"),as.Date("2022-10-26")),
               breaks=c(as.Date('2022-05-01'),as.Date('2022-07-01'),
                        as.Date('2022-09-01'),as.Date('2022-11-01')),
               labels=c('May','July','Sept.','Nov.'))+
  scale_y_continuous(limits=c(30,95))+
  ggtitle('iii. 2022')+
  theme_minimal()+
  theme(legend.position='None',
        plot.title=element_text(size=12))

plot_2023 <- ggplot()+
  geom_point(data=DH2023_catch,
             aes(x=Date_Observed, y=as.numeric(CW_mm),
                 color=Trap_Type))+
  scale_color_manual(values=c('violet','darkviolet','goldenrod'))+
  scale_x_date(limits=c(as.Date("2023-04-14"),as.Date("2023-10-26")),
               breaks=c(as.Date('2023-05-01'),as.Date('2023-07-01'),
                        as.Date('2023-09-01'),as.Date('2023-11-01')),
               labels=c('May','July','Sept.','Nov.'))+
  labs(x='date',y='',color='trap type')+
  ggtitle('iv. 2023')+
  theme_minimal()+
  theme(legend.position='None',
        plot.title=element_text(size=12))
plot_2022_w_legend <- ggplot()+
  geom_point(data=DH2022_catch,
             aes(x=Date_Observed, y=as.numeric(CW_mm),
                 color=Trap_Type))+
  scale_color_manual(values=c('violet','darkviolet','goldenrod'))+
  scale_x_date(limits=c(as.Date("2022-04-14"),as.Date("2022-10-26")),
               breaks=c(as.Date('2022-05-01'),as.Date('2022-07-01'),
                        as.Date('2022-09-01'),as.Date('2022-11-01')),
               labels=c('May','July','Sept.','Nov.'))+
  scale_y_continuous(limits=c(30,95))+
  labs(x='date',y='size (mm)',color='trap type')+
  ggtitle('Drayton Harbor 2022')+
  theme_minimal()+
  theme(plot.title=element_text(size=12))
legend <- get_legend(plot_2022_w_legend)

layout <- '
AAABBBCCCDDD##
AAABBBCCCDDD#E
AAABBBCCCDDD##
'

final_plot <- plot_2020+plot_2021+plot_2022+plot_2023+legend+plot_layout(design=layout)
ggsave('figures/fig_timeseries_data.svg',final_plot,
       dpi=400,width=10,height=3)
