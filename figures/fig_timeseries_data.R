library(tidyverse)
library(cowplot)
library(patchwork)

biweek <- c(59,76,91,106,121,137,152,167,182,198,213,229,244,259,274,290,305,320,335)

#######################
# Drayton Harbor 2020 #
#######################

# import catch data
DH2020_catch <- read.csv('IPM/data/csv/DraytonHarbor/EGC.DH.2020_catch.csv')

# convert date
DH2020_catch$Date <- as.Date(DH2020_catch$Date,'%m/%d/%Y')

# rename site names in catch data not present in effort data
DH2020_catch$Site_Name[DH2020_catch$Site_Name=='Prospecting - California Creek - Bridge'] <- 'Prospecting - California Creek'

# add Julian day
DH2020_catch$julian_day <- as.POSIXlt(DH2020_catch$Date,format='%m/%d/%Y')$yday

# keep only shrimp, fukui and minnow traps
DH2020_catch <- DH2020_catch[c(DH2020_catch$Trap_Type=='Fukui'|
                                 DH2020_catch$Trap_Type=='Shrimp'|
                                 DH2020_catch$Trap_Type=='Minnow'),]

# consolidate subsites
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Prospecting - Willse's"] <- "Willsei's"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Prospecting - Willse's/DNR East"] <- "Willsei's"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Prospecting - Willse's/DNR West"] <- "Willsei's"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Prospecting - Willse's/DNR West"] <- "Willsei's"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Prospecting - California Creek"] <- "Core - California Creek"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Prospecting - California Creek "] <- "Core - California Creek"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Index - CC - Mouth"] <- "Core - California Creek"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "California Creek - Mouth Index"] <- "Core - California Creek"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Index - CC - Graveyard "] <- "CC - Graveyard"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Index - CC - Graveyard"] <- "CC - Graveyard"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Prospecting - CC Graveyard"] <- "CC - Graveyard"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Prospecting - Pillars"] <- "Core - Pillars"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Prospecting - Jetty"] <- "Core - Jetty"
DH2020_catch$Site_Name[DH2020_catch$Site_Name == "Prospecting - No Name Creek"] <- "Core - No Name Creek"

# keep only subsites that are present in 90% of biweeks
DH2020_catch <- DH2020_catch[DH2020_catch$Site_Name%in%c('Core - California Creek',
                                                         'Core - Jetty','Core - No Name Creek',
                                                         'Core - Pillars'),]

#######################
# Drayton Harbor 2021 #
#######################

# import catch data
DH2021_catch <- read.csv('IPM/data/csv/DraytonHarbor/DraytonHarbor_2021_catch.csv')

# convert date
DH2021_catch$Date_Observed <- as.Date(DH2021_catch$Date_Observed,'%m/%d/%Y')

# remove traps without lat/lon
DH2021_catch <- DH2021_catch[!is.na(DH2021_catch$Latitude),]

# standard shrimp trap type
DH2021_catch$Trap_Type[DH2021_catch$Trap_Type == 'Shrimp - 1/2"'] <- 'Shrimp'

# add Julian day
DH2021_catch$julian_day <- as.POSIXlt(DH2021_catch$Date_Observed,format='%m/%d/%Y')$yday

# remove data from 'Central'
DH2021_catch <- DH2021_catch[!DH2021_catch$Site_Name=='Central',]

# remove grid survey traps in effort
DH2021_gridsurvey <- DH2021_catch[c(DH2021_catch$julian_day==235|
                                      DH2021_catch$julian_day==236),]
DH2021_catch <- DH2021_catch[!c(DH2021_catch$julian_day==235|
                                  DH2021_catch$julian_day==236),]
# add back Dakota Creek grid survey
DH2021_catch <- rbind(DH2021_catch,
                       DH2021_gridsurvey[DH2021_gridsurvey$Site_Name=='Dakota Creek',])

# keep only shrimp, fukui and minnow traps
DH2021_catch <- DH2021_catch[c(DH2021_catch$Trap_Type=='Fukui'|
                                 DH2021_catch$Trap_Type=='Shrimp'|
                                 DH2021_catch$Trap_Type=='Minnow'),]


#######################
# Drayton Harbor 2022 #
#######################

# import catch data
DH2022_catch <- read.csv('IPM/data/csv/DraytonHarbor/DH_catch_20221216.csv')

# convert date
DH2022_catch$Date_Observed <- as.Date(DH2022_catch$Date_Observed,'%m/%d/%Y')

# keep only shrimp, fukui and minnow traps
DH2022_catch <- DH2022_catch[c(DH2022_catch$Trap_Type=='Fukui'|
                                 DH2022_catch$Trap_Type=='Shrimp'|
                                 DH2022_catch$Trap_Type=='Minnow'),]

# add julian day to catch data
DH2022_catch$julian_day <- as.POSIXlt(DH2022_catch$Date_Observed,format='%m/%d/%Y')$yday

# remove early recruits
DH2022_catch <- DH2022_catch[!c(DH2022_catch$CW_mm < 40 & DH2022_catch$julian_day < 196),]
DH2022_catch <- DH2022_catch[!c(DH2022_catch$CW_mm < 45 & DH2022_catch$Month==7),]

# add biweekly index
for(i in 1:length(DH2022_catch$julian_day)){
  index <- tail(biweek[as.numeric(DH2022_catch[i,'julian_day']) > biweek],1)
  DH2022_catch[i,'biweek'] <- which(index == biweek)
}

# remove data from last biweek
max_biweek <- max(DH2022_catch$biweek)
DH2022_catch <- DH2022_catch[DH2022_catch$biweek < max_biweek,]


#######################
# Drayton Harbor 2023 #
#######################

# import catch data
DH2023_catch <- read.csv('IPM/data/csv/DraytonHarbor/DraytonHarbor_2023_catch_edited.csv')

# convert date
DH2023_catch$Date_Observed <- as.Date(DH2023_catch$Date_Observed,'%m/%d/%Y')

# standard shrimp trap type
DH2023_catch$Trap_Type[DH2023_catch$Trap_Type == 'Shrimp \x96 other'] <- 'Shrimp'

# consolidate subsites
DH2023_catch$Site_Name[DH2023_catch$Site_Name == 'East '] <- 'East'
DH2023_catch$Site_Name[DH2023_catch$Site_Name == 'California Creek - Mouth Index'] <- 'California Creek'
DH2023_catch$Site_Name[DH2023_catch$Site_Name == 'California Creek - Graveyard Index'] <- 'California Creek'
DH2023_catch$Site_Name[DH2023_catch$Site_Name == 'California Creek - Bridgeway Index'] <- 'California Creek'
DH2023_catch$Site_Name[DH2023_catch$Site_Name == 'Dakota Creek - I5 Index'] <- 'Dakota Creek'
DH2023_catch$Site_Name[DH2023_catch$Site_Name == 'Dakota Creek - Mouth Index'] <- 'Dakota Creek'

# keep only consistent sites
DH2023_catch <- DH2023_catch[DH2023_catch$Site_Name %in% c('California Creek',
                                                              'Dakota Creek',
                                                              'North','Northeast'),]

# keep only shrimp, fukui and minnow traps
DH2023_catch <- DH2023_catch[c(DH2023_catch$Trap_Type=='Fukui'|
                                 DH2023_catch$Trap_Type=='Shrimp'|
                                 DH2023_catch$Trap_Type=='Minnow'),]



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
  ggtitle('A. 2020')+
  theme_minimal()+
  theme(legend.position='None',
        plot.title=element_text(size=10))

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
  ggtitle('B. 2021')+
  theme_minimal()+
  theme(legend.position='None',
        plot.title=element_text(size=10)) 

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
  ggtitle('C. 2022')+
  theme_minimal()+
  theme(legend.position='None',
        plot.title=element_text(size=10))

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
  ggtitle('D. 2023')+
  theme_minimal()+
  theme(legend.position='None',
        plot.title=element_text(size=10))
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
  theme(plot.title=element_text(size=10))
plot_2021_w_legend <- ggplot()+
  geom_point(data=DH2021_catch,
             aes(x=Date_Observed, y=as.numeric(CW_mm),
                 color=Trap_Type))+
  scale_color_manual(values=c('violet','darkviolet','goldenrod'))+
  scale_x_date(limits=c(as.Date("2021-04-14"),as.Date("2021-10-26")),
               breaks=c(as.Date('2021-05-01'),as.Date('2021-07-01'),
                        as.Date('2021-09-01'),as.Date('2021-11-01')),
               labels=c('May','July','Sept.','Nov.'))+
  scale_y_continuous(limits=c(30,95))+
  labs(x='date',y='size (mm)',color='trap type')+
  ggtitle('Drayton Harbor 2021')+
  theme_minimal()+
  theme(plot.title=element_text(size=10))
legend <- get_legend(plot_2022_w_legend)

layout <- '
AAABBBCCCDDD##
AAABBBCCCDDD#E
AAABBBCCCDDD##
'

final_plot <- plot_2020+plot_2021+plot_2022+plot_2023+legend+plot_layout(design=layout)
ggsave('figures/fig_timeseries_data.png',final_plot,dpi=400,width=10,height=3)
