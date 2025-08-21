library(tidyverse)

# create parameters for IPM mesh
min.size <- 0
max.size <- 110
nsize <- 22 #number of cells in the discretized kernel
b=min.size+c(0:nsize)*(max.size-min.size)/nsize #boundary points (the edges of the cells defining the kernel)
y=0.5*(b[1:nsize]+b[2:(nsize+1)]) # mesh points (midpoints of the cells)

# get size column names
size_colnames <- rep(NA,length(b)-1)
for(i in 1:(length(b)-1)){
  size_colnames[i] <- paste0('s',b[i],'_',b[i+1])
}

biweek <- 59 + 0:18 * 14





#######################
# Drayton Harbor 2020 #
#######################

# import effort data
DH2020_effort <- read.csv('IPM/data/csv/DraytonHarbor/EGC.DH.2020_effort.csv')

# import catch data
DH2020_catch <- read.csv('IPM/data/csv/DraytonHarbor/EGC.DH.2020_catch.csv')

# convert date
DH2020_effort$Date_Retrieved <- as.Date(DH2020_effort$Date_Retrieved,'%m/%d/%Y')
DH2020_effort$Date_Deployed <- as.Date(DH2020_effort$Date_Deployed,'%m/%d/%Y')
DH2020_catch$Date <- as.Date(DH2020_catch$Date,'%m/%d/%Y')

# remove traps that were not retrieved
DH2020_effort <- DH2020_effort[!is.na(DH2020_effort$Date_Retrieved),]

# rename site names in catch data not present in effort data
DH2020_catch$Site_Name[DH2020_catch$Site_Name=='Prospecting - California Creek - Bridge'] <- 'Prospecting - California Creek'

# calculate soak days
DH2020_effort$Soak_days <- as.numeric(DH2020_effort[,'Date_Retrieved']-DH2020_effort[,'Date_Deployed'])

# add Julian day
DH2020_effort$julian_day <- as.POSIXlt(DH2020_effort$Date_Retrieved,format='%m/%d/%Y')$yday

# add biweekly index
for(i in 1:length(DH2020_effort$julian_day)){
  index <- tail(biweek[as.numeric(DH2020_effort[i,'julian_day']) > biweek],1)
  DH2020_effort[i,'biweek'] <- which(index == biweek)
}

# keep only shrimp, fukui and minnow traps
DH2020_effort <- DH2020_effort[c(DH2020_effort$Trap_Type=='Fukui'|
                                   DH2020_effort$Trap_Type=='Shrimp'|
                                   DH2020_effort$Trap_Type=='Minnow'),]

# add dummy variables for trap type and create trap_ID column
DH2020_effort <- DH2020_effort %>%
  mutate(Shrimp = ifelse(Trap_Type=='Shrimp',1,0),
         Fukui = ifelse(Trap_Type=='Fukui',1,0),
         Minnow = ifelse(Trap_Type=='Minnow',1,0)) %>%
  unite('ID',c('Trap_Number','julian_day'),
        sep='_',remove=FALSE)

# consolidate subsites
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Prospecting - Willse's"] <- "Willsei's"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Prospecting - Willse's/DNR East"] <- "Willsei's"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Prospecting - Willse's/DNR West"] <- "Willsei's"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Prospecting - Willse's/DNR West"] <- "Willsei's"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Prospecting - California Creek"] <- "Core - California Creek"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Prospecting - California Creek "] <- "Core - California Creek"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Index - CC - Mouth"] <- "Core - California Creek"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "California Creek - Mouth Index"] <- "Core - California Creek"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Index - CC - Graveyard "] <- "CC - Graveyard"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Index - CC - Graveyard"] <- "CC - Graveyard"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Prospecting - CC Graveyard"] <- "CC - Graveyard"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Prospecting - Pillars"] <- "Core - Pillars"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Prospecting - Jetty"] <- "Core - Jetty"
DH2020_effort$Site_Name[DH2020_effort$Site_Name == "Prospecting - No Name Creek"] <- "Core - No Name Creek"


# keep only subsites that are present in 90% of biweeks
DH2020_effort <- DH2020_effort[DH2020_effort$Site_Name%in%c('Core - California Creek',
                                                            'Core - Jetty','Core - No Name Creek',
                                                            'Core - Pillars'),]

# add unique trap id for biweekly index 
DH2020_effort_unique <- DH2020_effort %>%
  group_by(biweek) %>%
  tally()

# add unique site id for subsites
DH2020_subsite_unique <- DH2020_effort %>%
  group_by(Site_Name) %>%
  summarize(lat=median(Latitude),
            lon=median(Longitude),
            minnow_total=sum(Minnow),
            fukui_total=sum(Fukui),
            shrimp_total=sum(Shrimp)) %>% 
  mutate(subsite_index=1:n())
DH2020_shrimp_subsites <- DH2020_subsite_unique[DH2020_subsite_unique$shrimp_total>0,]$subsite_index
DH2020_minnow_subsites <- DH2020_subsite_unique[DH2020_subsite_unique$minnow_total>0,]$subsite_index
DH2020_fukui_subsites <- DH2020_subsite_unique[DH2020_subsite_unique$fukui_total>0,]$subsite_index

# add unique subsite id to effort
DH2020_effort <- left_join(DH2020_effort,DH2020_subsite_unique,by='Site_Name')

# create new empty dataframe
DH2020_effort2 <- as.data.frame(matrix(NA,nrow=0,ncol=ncol(DH2020_effort)+1))
colnames(DH2020_effort2) <- c(colnames(DH2020_effort),'trap_int')

# create new effort df with trap integer for each trap in each biweekly index of trapping
for(i in 1:length(DH2020_effort_unique$biweek)){
  #make subset of original effort data
  subset <- DH2020_effort[c(DH2020_effort$biweek==as.numeric(DH2020_effort_unique[i,'biweek'])),]
  #find count for unique combination
  count <- DH2020_effort_unique[c(DH2020_effort_unique$biweek==as.numeric(DH2020_effort_unique[i,'biweek'])),'n']
  #add count to df
  subset <- subset %>%
    mutate(trap_int=1:as.integer(count))
  
  #rbind with empty df
  DH2020_effort2 <- rbind(DH2020_effort2,subset)
}

# keep only crabs caught in fukui, minnow, and shrimp traps
DH2020_catch <- DH2020_catch[c(DH2020_catch$Trap_Type=='Fukui'|
                                 DH2020_catch$Trap_Type=='Shrimp'|
                                 DH2020_catch$Trap_Type=='Minnow'),]

# convert catch size to numeric
DH2020_catch$Size_mm <- as.numeric(DH2020_catch$Size_mm)

# add size ID to catch data
for(i in 1:(length(b)-1)){
  DH2020_catch[[size_colnames[i]]] <- ifelse(DH2020_catch$Size_mm>=b[i] & DH2020_catch$Size_mm < b[i+1],1,0)
}

# find the trap window associated with each trap and add the julian day
for(i in 1:length(DH2020_catch$Date)){
  greater <- DH2020_effort2[DH2020_effort2$Date_Deployed<=DH2020_catch[i,'Date'],] 
  all <- greater[greater$Date_Retrieved>=DH2020_catch[i,'Date'],]
  trapnum <- all[all$Trap_Number==DH2020_catch[i,'Trap_Number'],]
  DH2020_catch[i,'julian_day'] <- as.numeric(trapnum[trapnum$Site_Name==DH2020_catch[i,'Site_Name'],]['julian_day'])
}

#add trap ID to catch data
DH2020_catch_total <- DH2020_catch %>%
  unite('ID',c('Trap_Number','julian_day'),
        sep='_',remove=FALSE) %>%
  group_by(ID,julian_day) %>%
  summarise_at(size_colnames,sum)

#add catch to effort data
DH2020_effort2 <- left_join(DH2020_effort2,DH2020_catch_total,by=c('ID')) %>%
  replace(is.na(.),0)

# create counts
counts_DH2020 <- array(data=NA,
                       dim=c(length(unique(DH2020_effort2$biweek)),
                             max(DH2020_effort2$trap_int),
                             length(y)))
for(i in 1:length(unique(DH2020_effort2$biweek))){
  time <- unique(DH2020_effort2$biweek)[i]
  subset <- DH2020_effort2[DH2020_effort2$biweek==time,size_colnames]
  counts_DH2020[i,1:dim(subset)[1],] <- as.matrix(subset)
}

# create ncap
ncap_DH2020 <- DH2020_effort2 %>% 
  group_by(biweek) %>% 
  summarize(across(all_of(size_colnames), ~sum(.), .names = "{.col}"))
ncap_DH2020 <- as.data.frame(ncap_DH2020)
rownames(ncap_DH2020) <- ncap_DH2020[,1]
ncap_DH2020 <- ncap_DH2020[,-1]

# create julian day index
index_DH2020 <- unique(DH2020_effort2$biweek)

# total time
totalt_DH2020 <- length(unique(DH2020_effort2$biweek))

# total obs
totalo_DH2020 <- as.vector(table(DH2020_effort2$biweek))

# create soak days
soak_days_DH2020 <- DH2020_effort2 %>% 
  dplyr::select(biweek,trap_int,Soak_days) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Soak_days)
soak_days_DH2020 <- as.data.frame(soak_days_DH2020)
rownames(soak_days_DH2020) <- soak_days_DH2020[,1]
soak_days_DH2020 <- soak_days_DH2020[,-1]

# trap type binary indicators
# m
m_index_DH2020 <- DH2020_effort2 %>% 
  dplyr::select(biweek,trap_int,Minnow) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Minnow)
m_index_DH2020 <- as.data.frame(m_index_DH2020)
rownames(m_index_DH2020) <- m_index_DH2020[,1]
m_index_DH2020 <- m_index_DH2020[,-1]
# f
f_index_DH2020 <- DH2020_effort2 %>% 
  dplyr::select(biweek,trap_int,Fukui) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Fukui)
f_index_DH2020 <- as.data.frame(f_index_DH2020)
rownames(f_index_DH2020) <- f_index_DH2020[,1]
f_index_DH2020 <- f_index_DH2020[,-1]
# s
s_index_DH2020 <- DH2020_effort2 %>% 
  dplyr::select(biweek,trap_int,Shrimp) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Shrimp)
s_index_DH2020 <- as.data.frame(s_index_DH2020)
rownames(s_index_DH2020) <- s_index_DH2020[,1]
s_index_DH2020 <- s_index_DH2020[,-1]

#######################
# Drayton Harbor 2021 #
#######################

# import effort data
DH2021_effort <- read.csv('IPM/data/csv/DraytonHarbor/DraytonHarbor_2021_effort.csv')

# import catch data
DH2021_catch <- read.csv('IPM/data/csv/DraytonHarbor/DraytonHarbor_2021_catch.csv')

# convert date
DH2021_effort$Date_Retrieved <- as.Date(DH2021_effort$Date_Retrieved,'%m/%d/%Y')
DH2021_effort$Date_Deployed <- as.Date(DH2021_effort$Date_Deployed,'%m/%d/%Y')
DH2021_effort$Date_Checked <- as.Date(DH2021_effort$Date_Checked,'%m/%d/%Y')
DH2021_catch$Date_Observed <- as.Date(DH2021_catch$Date_Observed,'%m/%d/%Y')

# remove traps that were not retrieved
DH2021_effort <- DH2021_effort[!is.na(DH2021_effort$Date_Retrieved),]

# remove spaces in data
DH2021_catch$Trap_Number <- gsub(" ","",DH2021_catch$Trap_Number)
DH2021_effort$Trap_Number <- gsub(" ","",DH2021_effort$Trap_Number)

# remove traps without lat/lon
DH2021_catch <- DH2021_catch[!is.na(DH2021_catch$Latitude),]

# standard shrimp trap type
DH2021_effort$Trap_Type[DH2021_effort$Trap_Type == 'Shrimp - 1/2"'] <- 'Shrimp'
DH2021_catch$Trap_Type[DH2021_catch$Trap_Type == 'Shrimp - 1/2"'] <- 'Shrimp'

# add Julian day
DH2021_effort$julian_day <- as.POSIXlt(DH2021_effort$Date_Checked,format='%m/%d/%Y')$yday

# remove traps that were not checked
DH2021_effort <- DH2021_effort[!is.na(DH2021_effort$julian_day),]

# remove data from 'Central'
DH2021_effort <- DH2021_effort[!DH2021_effort$Site_Name=='Central',]
DH2021_catch <- DH2021_catch[!DH2021_catch$Site_Name=='Central',]

# remove grid survey traps in effort
DH2021_gridsurvey <- DH2021_effort[c(DH2021_effort$julian_day==235|
                                       DH2021_effort$julian_day==236),]
DH2021_effort <- DH2021_effort[!c(DH2021_effort$julian_day==235|
                                    DH2021_effort$julian_day==236),]
# add back Dakota Creek grid survey
DH2021_effort <- rbind(DH2021_effort,
                       DH2021_gridsurvey[DH2021_gridsurvey$Site_Name=='Dakota Creek',])

# add biweekly index
for(i in 1:length(DH2021_effort$julian_day)){
  index <- tail(biweek[as.numeric(DH2021_effort[i,'julian_day']) > biweek],1)
  DH2021_effort[i,'biweek'] <- which(index == biweek)
}

# add dummy variables for trap type and create trap_ID column
DH2021_effort <- DH2021_effort %>%
  mutate(Shrimp = ifelse(Trap_Type=='Shrimp',1,0),
         Fukui = ifelse(Trap_Type=='Fukui',1,0),
         Minnow = ifelse(Trap_Type=='Minnow',1,0)) %>%
  unite('ID',c('Trap_Number','julian_day'),
        sep='_',remove=FALSE)

# keep only core sites
DH2021_effort <- DH2021_effort[DH2021_effort$Site_Name%in%c('California Creek',
                                                            'Dakota Creek',
                                                            'North'),]

# ggplot()+
#   geom_point(data=DH2021_effort,
#              aes(x=as.numeric(Longitude),
#                  y=as.numeric(Latitude),color=Site_Name))

# add unique site id for subsites
DH2021_subsite_unique <- DH2021_effort %>%
  group_by(Site_Name) %>%
  summarize(lat=median(as.numeric(Latitude)),
            lon=median(as.numeric(Longitude)),
            minnow_total=sum(Minnow),
            fukui_total=sum(Fukui),
            shrimp_total=sum(Shrimp)) %>% 
  mutate(subsite_index=1:n())
DH2021_shrimp_subsites <- DH2021_subsite_unique[DH2021_subsite_unique$shrimp_total>0,]$subsite_index
DH2021_minnow_subsites <- DH2021_subsite_unique[DH2021_subsite_unique$minnow_total>0,]$subsite_index
DH2021_fukui_subsites <- DH2021_subsite_unique[DH2021_subsite_unique$fukui_total>0,]$subsite_index

# add unique subsite id to effort
DH2021_effort <- left_join(DH2021_effort,DH2021_subsite_unique,by='Site_Name')

# add unique trap id for biweekly index
DH2021_effort_unique <- DH2021_effort %>%
  group_by(biweek) %>%
  tally()

# create new empty dataframe
DH2021_effort2 <- as.data.frame(matrix(NA,nrow=0,ncol=ncol(DH2021_effort)+1))
colnames(DH2021_effort2) <- c(colnames(DH2021_effort),'trap_int')

# create new effort df with trap integer for each trap in each biweekly index of trapping
for(i in 1:length(DH2021_effort_unique$biweek)){
  #make subset of original effort data
  subset <- DH2021_effort[c(DH2021_effort$biweek==as.numeric(DH2021_effort_unique[i,'biweek'])),]
  #find count for unique combination
  count <- DH2021_effort_unique[c(DH2021_effort_unique$biweek==as.numeric(DH2021_effort_unique[i,'biweek'])),'n']
  #add count to df
  subset <- subset %>%
    mutate(trap_int=1:as.integer(count))
  
  #rbind with empty df
  DH2021_effort2 <- rbind(DH2021_effort2,subset)
}


DH2021_effort2$Soak_days <- 1

# keep only crabs caught in fukui, minnow, and shrimp traps
DH2021_catch <- DH2021_catch[c(DH2021_catch$Trap_Type=='Fukui'|
                                 DH2021_catch$Trap_Type=='Shrimp'|
                                 DH2021_catch$Trap_Type=='Minnow'),]

# convert catch size to numeric
DH2021_catch$CW_mm <- as.numeric(DH2021_catch$CW_mm)

# add size ID to catch data
for(i in 1:(length(b)-1)){
  DH2021_catch[[size_colnames[i]]] <- ifelse(DH2021_catch$CW_mm>=b[i] & DH2021_catch$CW_mm < b[i+1],1,0)
}

# add julian day to catch data
DH2021_catch$julian_day <- as.POSIXlt(DH2021_catch$Date_Observed,format='%m/%d/%Y')$yday

# remove grid survey traps
DH2021_catch_gridsurvey <- DH2021_catch[c(DH2021_catch$julian_day==235|
                                            DH2021_catch$julian_day==236),]
DH2021_catch <- DH2021_catch[!c(DH2021_catch$julian_day==235|
                                  DH2021_catch$julian_day==236),]
# add back Dakota Creek grid DH2021_catch
DH2021_catch <- rbind(DH2021_catch,
                      DH2021_catch_gridsurvey[DH2021_catch_gridsurvey$Site_Name=='Dakota Creek',])


#add trap ID to catch data
DH2021_catch_total <- DH2021_catch %>%
  unite('ID',c('Trap_Number','julian_day'),
        sep='_',remove=FALSE) %>%
  group_by(ID,julian_day,Trap_Type) %>%
  summarise_at(size_colnames,sum)

#add catch to effort data
DH2021_effort2 <- left_join(DH2021_effort2,DH2021_catch_total,by=c('ID','Trap_Type')) %>%
  replace(is.na(.),0)

# create counts
counts_DH2021 <- array(data=NA,
                       dim=c(length(unique(DH2021_effort2$biweek)),
                             max(DH2021_effort2$trap_int),
                             length(y)))
for(i in 1:length(unique(DH2021_effort2$biweek))){
  time <- unique(DH2021_effort2$biweek)[i]
  subset <- DH2021_effort2[DH2021_effort2$biweek==time,size_colnames]
  counts_DH2021[i,1:dim(subset)[1],] <- as.matrix(subset)
}

# create ncap
ncap_DH2021 <- DH2021_effort2 %>% 
  group_by(biweek) %>% 
  summarize(across(all_of(size_colnames), ~sum(.), .names = "{.col}"))
ncap_DH2021 <- as.data.frame(ncap_DH2021)
rownames(ncap_DH2021) <- ncap_DH2021[,1]
ncap_DH2021 <- ncap_DH2021[,-1]

# create julian day index
index_DH2021 <- unique(DH2021_effort2$biweek)

# total time
totalt_DH2021 <- length(unique(DH2021_effort2$biweek))

# total obs
totalo_DH2021 <- as.vector(table(DH2021_effort2$biweek))

# create soak days
soak_days_DH2021 <- DH2021_effort2 %>% 
  dplyr::select(biweek,trap_int,Soak_days) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Soak_days)
soak_days_DH2021 <- as.data.frame(soak_days_DH2021)
rownames(soak_days_DH2021) <- soak_days_DH2021[,1]
soak_days_DH2021 <- soak_days_DH2021[,-1]

# trap type binary indicators
# m
m_index_DH2021 <- DH2021_effort2 %>% 
  dplyr::select(biweek,trap_int,Minnow) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Minnow)
m_index_DH2021 <- as.data.frame(m_index_DH2021)
rownames(m_index_DH2021) <- m_index_DH2021[,1]
m_index_DH2021 <- m_index_DH2021[,-1]
# f
f_index_DH2021 <- DH2021_effort2 %>% 
  dplyr::select(biweek,trap_int,Fukui) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Fukui)
f_index_DH2021 <- as.data.frame(f_index_DH2021)
rownames(f_index_DH2021) <- f_index_DH2021[,1]
f_index_DH2021 <- f_index_DH2021[,-1]
# s
s_index_DH2021 <- DH2021_effort2 %>% 
  dplyr::select(biweek,trap_int,Shrimp) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Shrimp)
s_index_DH2021 <- as.data.frame(s_index_DH2021)
rownames(s_index_DH2021) <- s_index_DH2021[,1]
s_index_DH2021 <- s_index_DH2021[,-1]




#######################
# Drayton Harbor 2022 #
#######################

# import effort data
DH2022_effort <- read.csv('IPM/data/csv/DraytonHarbor/DH_effort_20221216.csv')

# import catch data
DH2022_catch <- read.csv('IPM/data/csv/DraytonHarbor/DH_catch_20221216.csv')

# convert date
DH2022_effort$Date_Retrieved <- as.Date(DH2022_effort$Date_Retrieved,'%m/%d/%Y')
DH2022_effort$Date_Deployed <- as.Date(DH2022_effort$Date_Deployed,'%m/%d/%Y')
DH2022_effort$Date_Checked <- as.Date(DH2022_effort$Date_Checked,'%m/%d/%Y')
DH2022_catch$Date_Observed <- as.Date(DH2022_catch$Date_Observed,'%m/%d/%Y')

# remove traps that were not retrieved
DH2022_effort <- DH2022_effort[!is.na(DH2022_effort$Date_Retrieved),]

# remove spaces in data
DH2022_catch$Trap_Number <- gsub(" ","",DH2022_catch$Trap_Number)
DH2022_effort$Trap_Number <- gsub(" ","",DH2022_effort$Trap_Number)
# clean catch data
DH2022_catch$Trap_Number[DH2022_catch$Trap_Number == '1191-Pink'] <- '1191-R'
DH2022_catch$Trap_Number[DH2022_catch$Trap_Number == '1191-Pi'] <- '1191-R'
DH2022_catch$Trap_Number[DH2022_catch$Trap_Number == '1398-G'] <- '1398-P'
DH2022_catch$Trap_Number[DH2022_catch$Trap_Number == '73-p'] <- '73-P'
for(i in 1:length(DH2022_effort$Date_Checked)){
  if(DH2022_effort[i,'Date_Checked']>DH2022_effort[i,'Date_Retrieved']){
    DH2022_effort[i,'Date_Checked'] <- DH2022_effort[i,'Date_Retrieved']
  }
}

# add Julian day
DH2022_effort$julian_day <- as.POSIXlt(DH2022_effort$Date_Checked,format='%m/%d/%Y')$yday

# add biweekly index
for(i in 1:length(DH2022_effort$julian_day)){
  index <- tail(biweek[as.numeric(DH2022_effort[i,'julian_day']) > biweek],1)
  DH2022_effort[i,'biweek'] <- which(index == biweek)
}

# keep only shrimp, fukui and minnow traps
DH2022_effort <- DH2022_effort[c(DH2022_effort$Trap_Type=='Fukui'|
                                   DH2022_effort$Trap_Type=='Shrimp'|
                                   DH2022_effort$Trap_Type=='Minnow'),]

# remove traps without lat/lon
DH2022_effort <- DH2022_effort[!is.na(DH2022_effort$Latitude),]

# add dummy variables for trap type and create trap_ID column
DH2022_effort <- DH2022_effort %>%
  mutate(Shrimp = ifelse(Trap_Type=='Shrimp',1,0),
         Fukui = ifelse(Trap_Type=='Fukui',1,0),
         Minnow = ifelse(Trap_Type=='Minnow',1,0)) %>%
  unite('ID',c('Trap_Number','julian_day'),
        sep='_',remove=FALSE)

# keep only consistent sites
DH2022_effort <- DH2022_effort[DH2022_effort$Site_Name %in% c('California Creek',
                                                              'Dakota Creek','North'),]

# ggplot()+
#   geom_point(data=DH2022_effort,
#              aes(x=Longitude,y=Latitude,color=Site_Name))

# add unique trap id for biweekly index
DH2022_effort_unique <- DH2022_effort %>%
  group_by(biweek) %>%
  tally()

# add unique site id for subsites
DH2022_subsite_unique <- DH2022_effort %>%
  group_by(Site_Name) %>%
  summarize(lat=median(as.numeric(Latitude)),
            lon=median(as.numeric(Longitude)),
            minnow_total=sum(Minnow),
            fukui_total=sum(Fukui),
            shrimp_total=sum(Shrimp)) %>% 
  mutate(subsite_index=1:n())

DH2022_shrimp_subsites <- DH2022_subsite_unique[DH2022_subsite_unique$shrimp_total>0,]$subsite_index
DH2022_minnow_subsites <- DH2022_subsite_unique[DH2022_subsite_unique$minnow_total>0,]$subsite_index
DH2022_fukui_subsites <- DH2022_subsite_unique[DH2022_subsite_unique$fukui_total>0,]$subsite_index


# add unique subsite id to effort
DH2022_effort <- left_join(DH2022_effort,DH2022_subsite_unique,by='Site_Name')

# create new empty dataframe
DH2022_effort2 <- as.data.frame(matrix(NA,nrow=0,ncol=ncol(DH2022_effort)+1))
colnames(DH2022_effort2) <- c(colnames(DH2022_effort),'trap_int')

# create new effort df with trap integer for each trap in each biweekly index of trapping
for(i in 1:length(DH2022_effort_unique$biweek)){
  #make subset of original effort data
  subset <- DH2022_effort[c(DH2022_effort$biweek==as.numeric(DH2022_effort_unique[i,'biweek'])),]
  #find count for unique combination
  count <- DH2022_effort_unique[c(DH2022_effort_unique$biweek==as.numeric(DH2022_effort_unique[i,'biweek'])),'n']
  #add count to df
  subset <- subset %>%
    mutate(trap_int=1:as.integer(count))
  
  #rbind with empty df
  DH2022_effort2 <- rbind(DH2022_effort2,subset)
}


DH2022_effort2$Soak_days <- 1

# keep only crabs caught in fukui, minnow, and shrimp traps
DH2022_catch <- DH2022_catch[c(DH2022_catch$Trap_Type=='Fukui'|
                                 DH2022_catch$Trap_Type=='Shrimp'|
                                 DH2022_catch$Trap_Type=='Minnow'),]

# convert catch size to numeric
DH2022_catch$CW_mm <- as.numeric(DH2022_catch$CW_mm)

# add size ID to catch data
for(i in 1:(length(b)-1)){
  DH2022_catch[[size_colnames[i]]] <- ifelse(DH2022_catch$CW_mm>=b[i] & DH2022_catch$CW_mm < b[i+1],1,0)
}

# add julian day to catch data
DH2022_catch$julian_day <- as.POSIXlt(DH2022_catch$Date_Observed,format='%m/%d/%Y')$yday

# remove early recruits
DH2022_catch <- DH2022_catch[!c(DH2022_catch$CW_mm < 40 & DH2022_catch$julian_day < 196),]
DH2022_catch <- DH2022_catch[!c(DH2022_catch$CW_mm < 45 & DH2022_catch$Month==7),]

#add trap ID to catch data
DH2022_catch_total <- DH2022_catch %>%
  unite('ID',c('Trap_Number','julian_day'),
        sep='_',remove=FALSE) %>%
  group_by(ID,julian_day) %>%
  summarise_at(size_colnames,sum)

#add catch to effort data
DH2022_effort2 <- left_join(DH2022_effort2,DH2022_catch_total,by=c('ID')) %>%
  replace(is.na(.),0)

# create counts
counts_DH2022 <- array(data=NA,
                       dim=c(length(unique(DH2022_effort2$biweek)),
                             max(DH2022_effort2$trap_int),
                             length(y)))
for(i in 1:length(unique(DH2022_effort2$biweek))){
  time <- unique(DH2022_effort2$biweek)[i]
  subset <- DH2022_effort2[DH2022_effort2$biweek==time,size_colnames]
  counts_DH2022[i,1:dim(subset)[1],] <- as.matrix(subset)
}

# create ncap
ncap_DH2022 <- DH2022_effort2 %>% 
  group_by(biweek) %>% 
  summarize(across(all_of(size_colnames), ~sum(.), .names = "{.col}"))
ncap_DH2022 <- as.data.frame(ncap_DH2022)
rownames(ncap_DH2022) <- ncap_DH2022[,1]
ncap_DH2022 <- ncap_DH2022[,-1]

# create julian day index
index_DH2022 <- unique(DH2022_effort2$biweek)

# total time
totalt_DH2022 <- length(unique(DH2022_effort2$biweek))

# total obs
totalo_DH2022 <- as.vector(table(DH2022_effort2$biweek))

# create soak days
soak_days_DH2022 <- DH2022_effort2 %>% 
  dplyr::select(biweek,trap_int,Soak_days) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Soak_days)
soak_days_DH2022 <- as.data.frame(soak_days_DH2022)
rownames(soak_days_DH2022) <- soak_days_DH2022[,1]
soak_days_DH2022 <- soak_days_DH2022[,-1]

# trap type binary indicators
# m
m_index_DH2022 <- DH2022_effort2 %>% 
  dplyr::select(biweek,trap_int,Minnow) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Minnow)
m_index_DH2022 <- as.data.frame(m_index_DH2022)
rownames(m_index_DH2022) <- m_index_DH2022[,1]
m_index_DH2022 <- m_index_DH2022[,-1]
# f
f_index_DH2022 <- DH2022_effort2 %>% 
  dplyr::select(biweek,trap_int,Fukui) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Fukui)
f_index_DH2022 <- as.data.frame(f_index_DH2022)
rownames(f_index_DH2022) <- f_index_DH2022[,1]
f_index_DH2022 <- f_index_DH2022[,-1]
# s
s_index_DH2022 <- DH2022_effort2 %>% 
  dplyr::select(biweek,trap_int,Shrimp) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Shrimp)
s_index_DH2022 <- as.data.frame(s_index_DH2022)
rownames(s_index_DH2022) <- s_index_DH2022[,1]
s_index_DH2022 <- s_index_DH2022[,-1]


#######################
# Drayton Harbor 2023 #
#######################

# import effort data
DH2023_effort <- read.csv('IPM/data/csv/DraytonHarbor/DraytonHarbor_2023_effort.csv')

# import catch data
DH2023_catch <- read.csv('IPM/data/csv/DraytonHarbor/DraytonHarbor_2023_catch_edited.csv')

# convert date
DH2023_effort$Date_Retrieved <- as.Date(DH2023_effort$Date_Retrieved,'%m/%d/%Y')
DH2023_effort$Date_Deployed <- as.Date(DH2023_effort$Date_Deployed,'%m/%d/%Y')
DH2023_effort$Date_Checked <- as.Date(DH2023_effort$Date_Checked,'%m/%d/%Y')
DH2023_catch$Date_Observed <- as.Date(DH2023_catch$Date_Observed,'%m/%d/%Y')

# remove spaces in data
DH2023_catch$Trap_Number <- gsub(" ","",DH2023_catch$Trap_Number)
DH2023_effort$Trap_Number <- gsub(" ","",DH2023_effort$Trap_Number)

# standard shrimp trap type
DH2023_effort$Trap_Type[DH2023_effort$Trap_Type == 'Shrimp \x96 other'] <- 'Shrimp'

# add Julian day
DH2023_effort$julian_day <- as.POSIXlt(DH2023_effort$Date_Checked,format='%m/%d/%Y')$yday

# add biweekly index
for(i in 1:length(DH2023_effort$julian_day)){
  index <- tail(biweek[as.numeric(DH2023_effort[i,'julian_day']) > biweek],1)
  DH2023_effort[i,'biweek'] <- which(index == biweek)
}


# add dummy variables for trap type and create trap_ID column
DH2023_effort <- DH2023_effort %>%
  mutate(Shrimp = ifelse(Trap_Type=='Shrimp',1,0),
         Fukui = ifelse(Trap_Type=='Fukui',1,0),
         Minnow = ifelse(Trap_Type=='Minnow',1,0)) %>%
  unite('ID',c('Trap_Number','julian_day'),
        sep='_',remove=FALSE)

# consolidate subsites
DH2023_effort$Site_Name[DH2023_effort$Site_Name == 'East '] <- 'East'
DH2023_effort$Site_Name[DH2023_effort$Site_Name == 'California Creek - Mouth Index'] <- 'California Creek'
DH2023_effort$Site_Name[DH2023_effort$Site_Name == 'California Creek - Graveyard Index'] <- 'California Creek'
DH2023_effort$Site_Name[DH2023_effort$Site_Name == 'California Creek - Bridgeway Index'] <- 'California Creek'
DH2023_effort$Site_Name[DH2023_effort$Site_Name == 'Dakota Creek - I5 Index'] <- 'Dakota Creek'
DH2023_effort$Site_Name[DH2023_effort$Site_Name == 'Dakota Creek - Mouth Index'] <- 'Dakota Creek'

# keep only consistent sites
DH2023_effort <- DH2023_effort[DH2023_effort$Site_Name %in% c('California Creek',
                                                              'Dakota Creek',
                                                              'North','Northeast'),]

# ggplot()+
#   geom_point(data=DH2023_effort,
#              aes(x=Longitude,y=Latitude,color=Site_Name))


# add unique trap id for biweekly index
DH2023_effort_unique <- DH2023_effort %>%
  group_by(biweek) %>%
  tally()

# add unique site id for subsites
DH2023_subsite_unique <- DH2023_effort %>%
  group_by(Site_Name) %>%
  summarize(lat=median(as.numeric(Latitude)),
            lon=median(as.numeric(Longitude)),
            minnow_total=sum(Minnow),
            fukui_total=sum(Fukui),
            shrimp_total=sum(Shrimp)) %>% 
  mutate(subsite_index=1:n())

DH2023_shrimp_subsites <- DH2023_subsite_unique[DH2023_subsite_unique$shrimp_total>0,]$subsite_index
DH2023_minnow_subsites <- DH2023_subsite_unique[DH2023_subsite_unique$minnow_total>0,]$subsite_index
DH2023_fukui_subsites <- DH2023_subsite_unique[DH2023_subsite_unique$fukui_total>0,]$subsite_index

# add unique subsite id to effort
DH2023_effort <- left_join(DH2023_effort,DH2023_subsite_unique,by='Site_Name')

# create new empty dataframe
DH2023_effort2 <- as.data.frame(matrix(NA,nrow=0,ncol=ncol(DH2023_effort)+1))
colnames(DH2023_effort2) <- c(colnames(DH2023_effort),'trap_int')

# create new effort df with trap integer for each trap in each biweekly index of trapping
for(i in 1:length(DH2023_effort_unique$biweek)){
  #make subset of original effort data
  subset <- DH2023_effort[c(DH2023_effort$biweek==as.numeric(DH2023_effort_unique[i,'biweek'])),]
  #find count for unique combination
  count <- DH2023_effort_unique[c(DH2023_effort_unique$biweek==as.numeric(DH2023_effort_unique[i,'biweek'])),'n']
  #add count to df
  subset <- subset %>%
    mutate(trap_int=1:as.integer(count))
  
  #rbind with empty df
  DH2023_effort2 <- rbind(DH2023_effort2,subset)
}

DH2023_effort2$Soak_days <- 1

# keep only crabs caught in fukui, minnow, and shrimp traps
DH2023_catch <- DH2023_catch[c(DH2023_catch$Trap_Type=='Fukui'|
                                 DH2023_catch$Trap_Type=='Shrimp'|
                                 DH2023_catch$Trap_Type=='Minnow'),]

# convert catch size to numeric
DH2023_catch$CW_mm <- as.numeric(DH2023_catch$CW_mm)

# add size ID to catch data
for(i in 1:(length(b)-1)){
  DH2023_catch[[size_colnames[i]]] <- ifelse(DH2023_catch$CW_mm>=b[i] & DH2023_catch$CW_mm < b[i+1],1,0)
}

# add julian day to catch data
DH2023_catch$julian_day <- as.POSIXlt(DH2023_catch$Date_Observed,format='%m/%d/%Y')$yday


#add trap ID to catch data
DH2023_catch_total <- DH2023_catch %>%
  unite('ID',c('Trap_Number','julian_day'),
        sep='_',remove=FALSE) %>%
  group_by(ID,julian_day,Trap_Type) %>%
  summarise_at(size_colnames,sum)

#add catch to effort data
DH2023_effort2 <- left_join(DH2023_effort2,DH2023_catch_total,by=c('ID','Trap_Type')) %>%
  replace(is.na(.),0)

# create counts
counts_DH2023 <- array(data=NA,
                       dim=c(length(unique(DH2023_effort2$biweek)),
                             max(DH2023_effort2$trap_int),
                             length(y)))
for(i in 1:length(unique(DH2023_effort2$biweek))){
  time <- unique(DH2023_effort2$biweek)[i]
  subset <- DH2023_effort2[DH2023_effort2$biweek==time,size_colnames]
  counts_DH2023[i,1:dim(subset)[1],] <- as.matrix(subset)
}

# create ncap
ncap_DH2023 <- DH2023_effort2 %>% 
  group_by(biweek) %>% 
  summarize(across(all_of(size_colnames), ~sum(.), .names = "{.col}"))
ncap_DH2023 <- as.data.frame(ncap_DH2023)
rownames(ncap_DH2023) <- ncap_DH2023[,1]
ncap_DH2023 <- ncap_DH2023[,-1]

# create julian day index
index_DH2023 <- unique(DH2023_effort2$biweek)

# total time
totalt_DH2023 <- length(unique(DH2023_effort2$biweek))

# total obs
totalo_DH2023 <- as.vector(table(DH2023_effort2$biweek))

# create soak days
soak_days_DH2023 <- DH2023_effort2 %>% 
  dplyr::select(biweek,trap_int,Soak_days) %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Soak_days)
soak_days_DH2023 <- as.data.frame(soak_days_DH2023)
rownames(soak_days_DH2023) <- soak_days_DH2023[,1]
soak_days_DH2023 <- soak_days_DH2023[,-1]

# trap type binary indicators
# m
m_index_DH2023 <- DH2023_effort2 %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Minnow)
m_index_DH2023 <- as.data.frame(m_index_DH2023)
rownames(m_index_DH2023) <- m_index_DH2023[,1]
m_index_DH2023 <- m_index_DH2023[,-1]
# f
f_index_DH2023 <- DH2023_effort2 %>%
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Fukui)
f_index_DH2023 <- as.data.frame(f_index_DH2023)
rownames(f_index_DH2023) <- f_index_DH2023[,1]
f_index_DH2023 <- f_index_DH2023[,-1]
# s
s_index_DH2023 <- DH2023_effort2 %>% 
  pivot_wider(id_cols=biweek,
              names_from=trap_int,
              values_from=Shrimp)
s_index_DH2023 <- as.data.frame(s_index_DH2023)
rownames(s_index_DH2023) <- s_index_DH2023[,1]
s_index_DH2023 <- s_index_DH2023[,-1]




####################
# final dataframes #
####################

maxt <- max(c(dim(counts_DH2020)[1],dim(counts_DH2021)[1],
              dim(counts_DH2022)[1],dim(counts_DH2023)[1]))
maxo <- max(c(dim(counts_DH2020)[2],dim(counts_DH2021)[2],
              dim(counts_DH2022)[2],dim(counts_DH2023)[2]))

# create counts
count_list <- list(counts_DH2020,counts_DH2021,
                   counts_DH2022,counts_DH2023
                   )

counts <- array(data=NA,dim=c(maxt,maxo,length(count_list),nsize))
for(i in 1:length(count_list)){
  counts[1:dim(count_list[[i]])[1], 1:dim(count_list[[i]])[2], i,] <- count_list[[i]]
}


# create ncap
ncap_list <- list(ncap_DH2020,ncap_DH2021,
                  ncap_DH2022,ncap_DH2023
                  )

ncap <- array(data=NA,dim=c(maxt,length(ncap_list),length(y)))
for(i in 1:length(ncap_list)){
  ncap[1:dim(ncap_list[[i]])[1], i,] <- as.matrix(ncap_list[[i]])
}

# create biweek index
index_list <- list(index_DH2020,index_DH2021,
                   index_DH2022,index_DH2023
                   )

index <- matrix(NA,nrow=maxt,ncol=length(index_list))
for(i in 1:length(index_list)){
  index[1:length(index_list[[i]]),i] <- index_list[[i]]
}

# fill in index for missing weeks
index_full <- index
for(i in 1:dim(index)[2]){
  if(max(index[,i],na.rm=TRUE) < max(index,na.rm=TRUE)){
    add <- (max(index[,i],na.rm = TRUE)+1):max(index,na.rm = TRUE)
    index_full[which(is.na(index[,i]))[1:length(add)],i] <- add
  }
}

# create delta t
delta_t <- matrix(NA,nrow=maxt,ncol=length(index_list))
# first delta t
delta_t[1,] <- index[1,] - 1
for(i in 1:length(index_list)){
  for(j in 1:(length(index[,i][!is.na(index[,i])])-1)){
    delta_t[j+1,i] <- (index[j+1,i]-index[j,i])
  }
}

# total time
totalt <- c(totalt_DH2020,totalt_DH2021,
            totalt_DH2022,totalt_DH2023
            )

# total obs
totalo_list <- list(totalo_DH2020,totalo_DH2021,
                    totalo_DH2022,totalo_DH2023
                    )

totalo <- matrix(NA,nrow=maxt,ncol=length(totalo_list))
for(i in 1:length(totalo_list)){
  totalo[1:length(totalo_list[[i]]),i] <- totalo_list[[i]]
}

# create soak days
soak_days_list <- list(soak_days_DH2020,soak_days_DH2021,
                       soak_days_DH2022,soak_days_DH2023
                       )

soak_days <- array(data=NA,dim=c(maxt,maxo,length(soak_days_list)))
for(i in 1:length(soak_days_list)){
  soak_days[1:dim(soak_days_list[[i]])[1], 1:dim(soak_days_list[[i]])[2], i] <- as.matrix(soak_days_list[[i]])
}

# trap type binary indicators
# m
m_index_list <- list(m_index_DH2020,m_index_DH2021,
                     m_index_DH2022,m_index_DH2023
                     )

m_index <- array(data=NA,dim=c(maxt,maxo,length(m_index_list)))
for(i in 1:length(m_index_list)){
  m_index[1:dim(m_index_list[[i]])[1], 1:dim(m_index_list[[i]])[2], i] <- as.matrix(m_index_list[[i]])
}

# f
f_index_list <- list(f_index_DH2020,f_index_DH2021,
                     f_index_DH2022,f_index_DH2023
                     )

f_index <- array(data=NA,dim=c(maxt,maxo,length(f_index_list)))
for(i in 1:length(f_index_list)){
  f_index[1:dim(f_index_list[[i]])[1], 1:dim(f_index_list[[i]])[2], i] <- as.matrix(f_index_list[[i]])
}

# s
s_index_list <- list(s_index_DH2020,s_index_DH2021,
                     s_index_DH2022,s_index_DH2023
                     )

s_index <- array(data=NA,dim=c(maxt,maxo,length(s_index_list)))
for(i in 1:length(s_index_list)){
  s_index[1:dim(s_index_list[[i]])[1], 1:dim(s_index_list[[i]])[2], i] <- as.matrix(s_index_list[[i]])
}

# subsite index
subsite_index_list <- list(DH2020_effort2$subsite_index,DH2021_effort2$subsite_index,
                           DH2022_effort2$subsite_index,DH2023_effort2$subsite_index)
max_effort <- max(length(DH2020_effort2$subsite_index),length(DH2021_effort2$subsite_index),
                  length(DH2022_effort2$subsite_index),length(DH2023_effort2$subsite_index))
subsite_index <- matrix(NA,nrow=max_effort,ncol=length(subsite_index_list))
for(i in 1:length(subsite_index_list)){
  subsite_index[1:length(subsite_index_list[[i]]),i] <- subsite_index_list[[i]]
}

# nsubsites
shrimp_subsites <- list(DH2020_shrimp_subsites,DH2021_shrimp_subsites,DH2022_shrimp_subsites,
                        DH2023_shrimp_subsites)
minnow_subsites <- list(DH2020_minnow_subsites,DH2021_minnow_subsites,DH2022_minnow_subsites,
                        DH2023_minnow_subsites)
fukui_subsites <- list(DH2020_fukui_subsites,DH2021_fukui_subsites,DH2022_fukui_subsites,
                        DH2023_fukui_subsites)

additionalt <- cbind(as.vector(totalt),
                     totalt + max(index, na.rm = TRUE) - 
                       apply(index, 2, max, na.rm = TRUE) - 1)

# save data
saveRDS(counts, "data/model_data/counts_new.rds")
saveRDS(ncap,"data/model_data/ncap_new.rds")
saveRDS(index,"data/model_data/index_new.rds")
saveRDS(index_full,"data/model_data/index_full_new.rds")
saveRDS((index_full - 3) / 365 * 14, "data/model_data/index_frac_new.rds")
saveRDS(totalt,"data/model_data/totalt_new.rds")
saveRDS(totalo,"data/model_data/totalo_new.rds")
saveRDS(soak_days,"data/model_data/soak_days_new.rds")
saveRDS(m_index,"data/model_data/m_index_new.rds")
saveRDS(f_index,"data/model_data/f_index_new.rds")
saveRDS(s_index,"data/model_data/s_index_new.rds")
saveRDS(additionalt, "data/model_data/additionalt_new.rds")
