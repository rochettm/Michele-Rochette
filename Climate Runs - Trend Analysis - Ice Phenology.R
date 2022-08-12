
# Load prewhitening code and required libraries
source('~/Thesis/9 Climate Run Analysis/general_test_mk_prewhiten.R')
library(readr)
library(trend)
library(Kendall)
library(lubridate)


#***********************************************************************************************
# Set up basic information
# Doing time series analysis for Climate run 1 in the Nelson Churchill basin, for ice phenology
#***********************************************************************************************
runOfInt <- "1"
varOfInt <- "cmri"


#***********************************************************************************************
# List of stations in the Nelson Churchill Basin and create path to find files later
#***********************************************************************************************
Stns <- read.delim("G:/Michele/9 Climate Run Analysis/3 R Codes/pmsf_NC.txt", header = FALSE)
runList <- list.files(path = "G:/Michele/9 Climate Run Analysis/1 Climate Run Results",pattern = "NC - 1 - ")


#***********************************************************************************************
#Set up tables to output data
#***********************************************************************************************
names_cols <- c("SUBID", "mean", "d1", "dL", "daysOn", "d1_2020", "dL_2020")
T1 <- data.frame(matrix(ncol = 7, nrow = 0))
names_cols2 <- c("SUBID", "MK_daysOn_pw", "pVal_daysOn_pw", "MK_daysOn", "pVal_daysOn",
                  "MK_d1_pw", "pVal_d1_pw", "MK_d1", "pVal_d1", "MK_dL_pw", "pVal_dL_pw", "MK_dL", "pVal_dL",
                 "MK_avg_pw", "pVal_avg_pw", "MK_avg", "pVal_avg")
T2 <- data.frame(matrix(ncol = 17, nrow = 0))


#***********************************************************************************************
#run the following for every station in the Nelson-Churchill (or whichever stations identified in Stns)
#***********************************************************************************************

for(i in seq(1,nrow(Stns))){   

  SUBID <- Stns[i, ]
  
  # Get Data
  inData <- read_delim(paste0("G:/Michele/9 Climate Run Analysis/1 Climate Run Results/", 
                              runList,"/", SUBID, ".txt"), delim = "\t", col_names = TRUE, col_types = cols())
  inData <- inData[-1, ]
  
  # Prepare Data - sorts it yearly
  inData[, match("DATE", colnames(inData))] <- as.Date(inData$DATE)
  inData[, match(varOfInt, colnames(inData))] <- as.numeric(inData$cmri)
  inData$YEAR <- as.numeric(strftime(inData$DATE, format = "%Y"))
  inData$MONTH <- as.numeric(strftime(inData$DATE, format = "%m"))
  inData$code <- inData$YEAR + inData$MONTH/100
  
  month_data <- aggregate(data=inData, as.numeric(cmri) ~ code, mean)
  names(month_data)[names(month_data) == 'as.numeric(cmri)'] <- 'cmri_avg'
  
  # Get Ice On and Ice Off Dates
  for(j in seq(1,89)){
    # Get Data for Year
    yr1 <- 1980 + j
    yr2 <- yr1 + 1
    
    yr1_rws <- which(inData$YEAR == yr1)
    yr2_rws <- which(inData$YEAR == yr2)
    
    rw1 <- yr1_rws[1] + 181
    if(j%%4 == 0){rw1 <- yr1_rws[1] + 182}
    rw2 <- yr2_rws[1] + 180
    if((j+1)%%4 == 0){rw2 <- yr2_rws[1] + 181}

    cmriData <- inData$cmri[min(yr1_rws):max(yr1_rws)]
    dateData <- inData$DATE[min(yr1_rws):max(yr1_rws)]
    
    #find first and last row with ice
    idx_iceOn <- which(cmriData != 0)
    
    # Calculate Statistics
    avg <- 0
    d1 <- as.Date("01/01/2020", "%m/%d/%y") 
    dL <- as.Date("01/01/2020", "%m/%d/%y")
    daysOn <- 0
    d1_2020 <- 0
    dL_2020 <- 0
    
    # error code
    if(length(idx_iceOn) != 0){
      iceOn <- as.numeric(cmriData[idx_iceOn])
      avg <- mean(iceOn)
      # day1 ice On, last day Ice on, total days with ice
      d1 <- dateData[idx_iceOn[1]]
      dL <- dateData[idx_iceOn[length(iceOn)]]
      daysOn <- as.numeric(dL-d1)
      d1_2020 <- month(d1) + day(d1)/30
      dL_2020 <- month(dL) + day(dL)/30

    }
    
    #add to T1
    outData <- data.frame(col1 = SUBID, col2 = avg, 
                          col3 = d1, col4 = dL, col5 = daysOn,
                          col6 = d1_2020, col7 = dL_2020)
    T1 <- rbind(T1, outData)
  }
  
  #colnames(T1) <- names_cols
  
  #***************************************************************************************
  # Calcs for Mann-Kendall trend test - for days on, first day with ice, last day with ice
  #***************************************************************************************
  
  #Get days on data point for each year
  rw1 <- 89*i - 88
  rw2 <- 89*i
  daysOn_series <- T1$col5[rw1:rw2]
  d1_2020_series <- T1$col6[rw1:rw2]
  dL_2020_series <- T1$col7[rw1:rw2]
  avg_series <- T1$col2[rw1:rw2]
  
  #run MK test
  MK_daysOn_pw <- general_test_mk_prewhiten(daysOn_series, 1, 0.05)
  MK_daysOn <- mk.test(daysOn_series)
  MK_d1_pw <- general_test_mk_prewhiten(d1_2020_series, 1, 0.05)
  MK_d1 <- mk.test(d1_2020_series)
  MK_dL_pw <- general_test_mk_prewhiten(dL_2020_series, 1, 0.05)
  MK_dL <- mk.test(dL_2020_series)
  MK_avg_pw <- general_test_mk_prewhiten(avg_series, 1, 0.05)
  MK_avg <- mk.test(avg_series)
  
  #Add data to output table
  outDataMK <- data.frame(col1 = SUBID, col2 = MK_daysOn_pw$slope, col3 = MK_daysOn_pw$pVal,
                          col4 = MK_daysOn$statistic, col5 = MK_daysOn$p.value,
                          col6 = MK_d1_pw$slope, col7 = MK_d1_pw$pVal,
                          col8 = MK_d1$statistic, col9 = MK_d1$p.value,
                          col10 = MK_dL_pw$slope, col11 = MK_dL_pw$pVal,
                          col12 = MK_dL$statistic, col13 = MK_dL$p.value,
                          col14 = MK_avg_pw$slope, col15 = MK_avg_pw$pVal,
                          col16 = MK_avg$statistic, col17 = MK_avg$p.value)  
                            
  T2 <- rbind(T2, outDataMK)
} 

#***************************************************************************************
# Update column names in tables
#***************************************************************************************
#*
colnames(T1) <- names_cols
colnames(T2) <- names_cols2
  

  
  