

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
# Doing time series analysis in the Nelson Churchill basin, for cmri variable (river ice thickness)
#***********************************************************************************************

basin <- "NC"
varOfInt = "cmri"


#***********************************************************************************************
# Open Observed/Measured Data (from ECCC) and file that connects observed data station to 
# corresponding HYPE subbasin
#***********************************************************************************************
Obs_data <- read.csv("G:/Michele/6 Primary Comparison Statistics/ForMatlab_ICE_NC.csv", header = TRUE)
Obs_to_HYPE <- read.delim("G:/Michele/5 Arctic HYPE/WSC_to_AHYPE_station_overlap_NC.txt", header = TRUE)

Stns <- read.delim("G:/Michele/11 Writing/timeseries/NC_stations.txt", header = FALSE)


#***********************************************************************************************
#Set up tables to output data
#***********************************************************************************************
names_cols <- c("WSC", "SUBID")
T1 <- data.frame(matrix(ncol = 2, nrow = nrow(Obs_to_HYPE)))
colnames(T1) <- names_cols
T1$WSC <- Obs_to_HYPE$pointID
T1$SUBID <- Obs_to_HYPE$polygonID

names_cols_HYPE <- c("jday", "cmri_HYPE")
T_HYPE <- data.frame(matrix(ncol = 2, nrow = 0))

names_cols_obs <- c("jday", "cmri_obs", "mn", "szn")
T_obs <- data.frame(matrix(ncol = 4, nrow = 0))


#***********************************************************************************************
#run the following for every station in the Nelson-Churchill (or whichever stations identified in Stns)
#***********************************************************************************************

for(i in seq(1,nrow(Stns))){ 
  
  
  #***************************************************************************************
  # Get observed/measured data
  #***************************************************************************************
  idx_obs <- which(Obs_data$ID == Obs_stn)
  yrs <- seq(1971, 2014)
  cmris <- t(Obs_data[idx_obs, ])
  cmris <- as.numeric(cmris[(5:48), ])
  cmri_obs <- data.frame(date = yrs, cmri_obs = cmris)
  
  T_obs <- rbind(T_obs, cmri_obs)
  
  #***************************************************************************************
  # Get corresponding model data
  #***************************************************************************************

  # Get HYPE ID NAME corresponding to the station
  Obs_stn <- Stns$V1[i]
  idx <- match(Obs_stn, as.character(T1$WSC))
  if(is.na(idx) == TRUE){next}
  SUBID <- T1$SUBID[idx]
  
  # Get Data
  inData <- read.delim(paste0("G:/Michele/5 Arctic HYPE/River Ice/Results_MK/", SUBID, ".txt"))
  inData <- inData[-1, ]
  cmriData <- data.frame(date = inData$DATE, cmri = inData$cmri)
  
  # Calculate mean ice thickness each year
  T_ice_mean <- data.frame(matrix(ncol = 2, nrow = 0))
  
  for(yr in seq(1971, 2014)){
    print(yr)
    d1 <- make_date(year = yr, month = 07, day = 01)
    d2 <- make_date(year = yr + 1, month = 06, day = 30)
    idx_d1 <- which(cmriData$date == d1)
    idx_d2 <- which(cmriData$date == d2)
    mean_cmri_data <- as.numeric(cmriData$cmri[(idx_d1:idx_d2)])
    idx_0 <- which(mean_cmri_data != 0)
    mean_cmri <- mean(mean_cmri_data[idx_0])
    outData <- data.frame(year = yr, mean_cmri = mean_cmri)
    T_ice_mean <- rbind(T_ice_mean, outData)
  }
  
  # Add data to output table
  T_HYPE <- rbind(T_HYPE, T_ice_mean)
  
}

#***************************************************************************************
# Plot
#***************************************************************************************
colnames(T_HYPE) <- c("year", "cmri")
T_test <- data.frame(yr = as.numeric(T_HYPE$year), ice = as.numeric(T_HYPE$cmri))
T_HYPE_plot <- aggregate(data=T_test, ice ~ yr, mean)

T_obs_mean <- aggregate(data=T_obs, cmri_obs ~ date, mean)

ggplot(T_obs, aes(x = date, y = cmri_obs*100)) + 
  geom_point(data = T_obs, mapping = aes(x = date, y = cmri_obs*100), size = 1, color = "black") + 
  geom_line(data = T_HYPE_plot, mapping=aes(x = yr, y = ice), size = 2, color = 'red') +
  geom_line(T_obs_mean, mapping=aes(date, cmri_obs*100), colour = 'blue', pch = 24, size = 1) + 
  labs(x = "Julian Date", y = "Mean River Ice Thickness (cm)")


#***************************************************************************************
# Calculate KGE, NSE
#***************************************************************************************

KGE <- KGE(T_HYPE_plot$ice, T_obs_mean$cmri_obs*100)
NSE <- NSE(T_HYPE_plot$ice, T_obs_mean$cmri_obs*100)

idx_years <- match(T_obs_mean$date, T_HYPE_plot$yr)
cmri_compare <- data.frame(obs = T_obs_mean$cmri_obs*100, hype = T_HYPE_plot$ice[idx_years])
cmri_dif <- data.frame(date = T_obs_mean$date, dif = (cmri_compare$obs - cmri_compare$hype))
