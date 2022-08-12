
# Load prewhitening code and required libraries
source('~/Thesis/9 Climate Run Analysis/general_test_mk_prewhiten.R')
library(readr)
library(trend)
library(Kendall)
library(lubridate)
library(ggplot2)
library(tidyverse)
library(hydroGOF)


#***********************************************************************************************
# Set up basic information
# Doing time series analysis in the Nelson Churchill basin, for wtmp variable (stream temperature)
#***********************************************************************************************

basin <- "NC"
varOfInt <- "wtmp"


#***********************************************************************************************
# Open Observed/Measured Data (from ECCC) and file that connects observed data station to 
# corresponding HYPE subbasin
#***********************************************************************************************
Obs_data <- read.csv("G:/Michele/6 Primary Comparison Statistics/ForMatlab_OBS_wtmp_NC.csv", header = TRUE)
Obs_to_HYPE <- read.delim("G:/Michele/5 Arctic HYPE/OBS_to_AHYPE_station_overlap_NC.txt", header = TRUE)


#***********************************************************************************************
#Set up tables to output data
#***********************************************************************************************
names_cols <- c("OBS", "SUBID")
T1 <- data.frame(matrix(ncol = 2, nrow = nrow(Obs_to_HYPE)))
colnames(T1) <- names_cols
T1$OBS <- Obs_to_HYPE$pointID
T1$SUBID <- Obs_to_HYPE$polygonID

names_cols_HYPE <- c("jday", "wtmp_HYPE")
T_HYPE <- data.frame(matrix(ncol = 2, nrow = 0))

names_cols_obs <- c("jday", "wtmp_obs", "mn", "szn")
T_obs <- data.frame(matrix(ncol = 4, nrow = 0))

Stns <- read.delim("G:/Michele/11 Writing/timeseries/NC_stations.txt", header = FALSE)


#***********************************************************************************************
#run the following for every station in the Nelson-Churchill (or whichever stations identified in Stns)
#***********************************************************************************************

for(i in seq(1,nrow(Stns))){ 
  
  
  #***************************************************************************************
  # Get observed/measured data
  #***************************************************************************************
  
  # Pull out observed data only for correct station
  idx_obs <- which(Obs_data$StationID == Obs_stn)
  wtmp_obs <- data.frame(date = Obs_data$Date[idx_obs], wtmp_obs = Obs_data$Tw[idx_obs])
  dates <- as.Date(wtmp_obs$date, "%m/%d/%Y")
  obs_data_i <- data.frame(date = dates, wtmp_obs = wtmp_obs$wtmp_obs)
  
  #remove data from feb 29
  obs_data_i$ymd <- make_date(year=2001L, month = month(obs_data_i$date), day = day(obs_data_i$date))
  idx_feb29 <- which(is.na(obs_data_i$ymd) == TRUE)
  if(length(idx_feb29) > 0){ 
    obs_data_i <- obs_data_i[(-1*idx_feb29), ]
  }
  
  #add data to output table
  T_obs <- rbind(T_obs, wtmp_obs_jday)

  #***************************************************************************************
  # Get corresponding model data
  #***************************************************************************************
  
  # Get HYPE ID NAME corresponding to the station
  Obs_stn <- Stns$V1[i]
  idx <- match(Obs_stn, as.character(T1$OBS))
  if(is.na(idx) == TRUE){next}
  SUBID <- T1$SUBID[idx]
  
  # Get data
  inData <- read.delim(paste0("G:/Michele/5 Arctic HYPE/Results_NC/", SUBID, ".txt"))
  inData <- inData[-1, ]
  inData[, match(varOfInt, colnames(inData))] <- as.numeric(inData$wtmp)
  inData$wtmp <- ifelse(inData$wtmp < 0, 0, inData$wtmp)
  
  # remove feb 29
  inData$ymd <- make_date(year = 2001L, month = month(inData$DATE), day = day(inData$DATE))
  idx_feb29 <- which(is.na(inData$ymd) == TRUE)
  inData <- inData[(-1*idx_feb29), ]
  inData$jday <- yday(inData$ymd)
  
  # Add data to output table
  wtmp_HYPE_jday <- data.frame(jday = inData$jday, wtmp_HYPE = inData$wtmp)
  wtmp_HYPE <- aggregate(data=wtmp_HYPE_jday, wtmp_HYPE ~ jday, mean)
  
  T_HYPE <- rbind(T_HYPE, wtmp_HYPE)
  
  
}

#***************************************************************************************
# Plot
#***************************************************************************************

colnames(T_HYPE) <- c("jday", "wtmp")
T_HYPE_plot <- aggregate(data=T_HYPE, wtmp ~ jday, mean)

# Mean obs line
T_obs_mean <- aggregate(data=T_obs, wtmp_obs ~ jday, mean)

# Plot
ggplot(T_obs,aes(jday,wtmp_obs)) +geom_point(size = 1) + 
  geom_line(T_HYPE_plot, mapping=aes(jday, wtmp), colour = 'red', size = 2) + 
  geom_line(T_obs_mean, mapping=aes(jday, wtmp_obs), colour = 'blue', size = 1) +
  #theme(panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(color = 'gray87'),
  #panel.grid.minor = element_line(color = 'gray87')) + 
  ylim(0, 30) + xlab("Julian Date") + ylab("Stream Temperature (°C)")


#***************************************************************************************
# Calculate KGE, NSE
#***************************************************************************************

wtmp_match <- T_HYPE_plot$wtmp[T_obs_mean$jday]
KGE <- KGE(wtmp_match, T_obs_mean$wtmp_obs)
NSE <- NSE(wtmp_match, T_obs_mean$wtmp_obs)

#compare
idx_jdays <- match(T_obs_mean$jday, T_HYPE_plot$jday)
wtmp_compare <- data.frame(obs = T_obs_mean$wtmp_obs, hype = T_HYPE_plot$wtmp[idx_jdays])
wtmp_dif <- data.frame(date = T_obs_mean$jday, dif = (wtmp_compare$obs - wtmp_compare$hype))

