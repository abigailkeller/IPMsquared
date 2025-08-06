library(tidyverse)
library(lubridate)

# create parameters for IPM mesh
min.size <- 0
max.size <- 110
nsize <- 22 #number of cells in the discretized kernel
b <- min.size + c(0:nsize) * (max.size - min.size) / nsize #boundary points (the edges of the cells defining the kernel)
y <- 0.5 * (b[1:nsize] + b[2:(nsize + 1)]) # mesh points (midpoints of the cells)

# get size column names
size_colnames <- rep(NA, length(b) - 1)
for(i in 1:(length(b) - 1)){
  size_colnames[i] <- paste0("s", b[i] , "_", b[i + 1])
}

biweek <- c(59, 76, 91, 106, 121, 137, 152, 167, 182, 198, 213,
            229, 244, 259, 274, 290, 305, 320, 335)

# read in mark-recapture data
data_mc <- read.csv("data/DFO_markrecapture.csv") %>%
  mutate(date = make_date(year = Year, month = Month, day = Day),
         julian_day = yday(date))

# add biweekly index
for(i in 1:length(data_mc$julian_day)){
  index <- tail(biweek[as.numeric(data_mc[i, "julian_day"]) > biweek], 1)
  data_mc[i, "biweek"] <- which(index == biweek)
}

# add size ID to catch data
for(i in 1:(length(b)-1)){
  data_mc[[size_colnames[i]]] <- ifelse(data_mc$WidthNotch >= b[i] & 
                                          data_mc$WidthNotch < b[i + 1], 1 , 0)
}

data_mc_recap <- data_mc[data_mc$Recap == 1, ]
data_mc_mark <- data_mc[data_mc$Recap == 0, ]

marked <- matrix(NA, nrow = length(unique(data_mc$biweek)) - 1, ncol = nsize)
caught <- matrix(NA, nrow = length(unique(data_mc$biweek)) - 1, ncol = nsize)

# get number marked and caught
for (t in 1:(length(unique(data_mc$biweek)) - 1)) {
  
  marked[t, ] <- colSums(
    data_mc_mark[data_mc_mark$biweek == unique(data_mc$biweek)[t], 
                 size_colnames]
    )
  
  caught[t, ] <- colSums(
    data_mc_recap[data_mc_recap$biweek == unique(data_mc$biweek)[t + 1],
                  size_colnames])
    
}

# get total traps
total_traps <- data_mc %>% 
  group_by(date, biweek) %>% 
  summarise(count = length(unique(SetNum)) * 6) %>% 
  group_by(biweek) %>% 
  summarise(n = sum(count))
  
saveRDS(total_traps, "data/model_data/roche_mc_totalo.rds")
saveRDS(caught, "data/model_data/roche_mc_catch.rds")
saveRDS(marked, "data/model_data/roche_mc_mark.rds")
saveRDS((total_traps$biweek - 3) / 365 * 14, "data/model_data/mc_index_rc.rds")
saveRDS(matrix(1, nrow = nrow(total_traps),
               ncol = max(total_traps$n)), 
        "data/model_data/soak_days_mc_rc.rds")
saveRDS(matrix(1, nrow = nrow(total_traps),
               ncol = max(total_traps$n)), 
        "data/model_data/f_index_mc_rc.rds")
saveRDS(matrix(0, nrow = nrow(total_traps),
               ncol = max(total_traps$n)), 
        "data/model_data/m_index_mc_rc.rds")
saveRDS(matrix(0, nrow = nrow(total_traps),
               ncol = max(total_traps$n)), 
        "data/model_data/s_index_mc_rc.rds")
