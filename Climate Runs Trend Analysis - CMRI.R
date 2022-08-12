
# Load prewhitening code and required libraries
source('G:/Michele/9 Climate Run Analysis/general_test_mk_prewhiten.R')
library(readr)
library(trend)
library(Kendall)
library(lubridate)


#***********************************************************************************************
# Set up basic information
# Doing time series analysis for Climate run 1 in the Nelson Churchill basin, for cmri variable (river ice thickness)
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
names_cols_Td <- c("SUBID", "year1", "ice_on", "ice_off", "ice_days", "avg_wo0")
Td <- data.frame(matrix(ncol = 6, nrow = 0))

names_cols_Tt <- c("SUBID", "MK_8170", "pVal_8170", "MK_8120", "pVal_8120", "MK_2170", 
                  "pVal_2170", "MK_2150", "pVal_2150", "MK_4170", "pVal_4170")   
Tt <- data.frame(matrix(ncol = 11, nrow = 0))



#***********************************************************************************************
#run the following for every station in the Nelson-Churchill (or whichever stations identified in Stns)
#***********************************************************************************************

for(i in seq(1,nrow(Stns))){   

  SUBID <- Stns[i, ]
  
  # Get Data
  inData <- read.delim(paste0("G:/Michele/9 Climate Run Analysis/1 Climate Run Results/", 
                              runList,"/", SUBID, ".txt"))
  inData <- inData[-1, ]
  
  # Prepare Data
  inData[, match("DATE", colnames(inData))] <- as.Date(inData$DATE)
  inData$YEAR <- as.numeric(strftime(inData$DATE, format = "%Y"))
  inData$MONTH <- as.numeric(strftime(inData$DATE, format = "%m"))
  days <- as.numeric(format(inData$DATE, "%j"))
  inData$DAY <- days
  inData$code <- inData$YEAR + as.numeric(inData$DAY)/1000

  #***************************************************************************************
  # Calculate mean ice thickness for every water year (Oct 1 - Sept 30)
  #***************************************************************************************
  for(j in seq(1,89)){
    
    #Get data for Year 1 Oct 1 to Year 2 Sept 30
    leap1 <- 0
    if(j%%4 == 0){leap1 <- 1}
    d1 <- 1980 + j + (274 + leap1)/1000
    leap2 <- 0
    if((j+1)%%4 == 0){leap2 <- 1}
    d2 <- 1981 + j + (273 + leap2)/1000
    oct1 <- which(inData$code == d1)
    sept30 <- which(inData$code == d2)
    
    #Get ice thickness data for water year
    cmriData <- inData$cmri[oct1:sept30]
    dateData <- inData$DATE[oct1:sept30]
    
    #find first and last row with ice
    idx_iceOn <- which(cmriData != 0)
    
    
    # Calculate Statistics
    ice_on_date <- 0
    ice_on <- 0
    ice_on_j <- 0
    ice_off_date <- 0
    ice_off <- 0
    ice_off_j <- 0
    ice_days <- 0
    avg_wo0 <- 0
    avg_w0 <- 0
    
    
    # error code
    if(length(idx_iceOn) != 0){
      ice_on_date <- dateData[idx_iceOn[1]]
	    ice_on_j <- format(ice_on_date, "%j")
      ice_off_date <- dateData[max(idx_iceOn)]
	    ice_off_j <- format(ice_off_date, "%j")
      ice_days <- length(idx_iceOn)
      avg_wo0 <- mean(as.numeric(as.character(cmriData[idx_iceOn])))
      }
    
    
    # Add to output table
    outData1 <- data.frame(col1 = SUBID, col2 = j + 1980, col3 = ice_on_j, col4 = ice_off_j, col5 = ice_days, col6 = avg_wo0)
    Td <- rbind(Td, outData1)
  
    
    }
  
  
  
  #***************************************************************************************
  # Get data for different time periods
  #***************************************************************************************
  
  rw1 <- 1 + 89*(i-1)
  rw2 <- 89 + 89*(i-1)

  year_8170 <- Td$col6[rw1:rw2]
  year_8120 <- Td$col6[rw1:(rw1+39)]
  year_2170 <- Td$col6[(rw1+40):rw2]
  year_2150 <- Td$col6[(rw1+40):(rw1+69)]
  year_4170 <- Td$col6[(rw1+60):rw2]

  #***************************************************************************************
  # Calcs for Mann-Kendall trend test - for different time periods
  #***************************************************************************************
  MK_yr_8170 <- general_test_mk_prewhiten(year_8170, 1, 0.05)
  MK_yr_8120 <- general_test_mk_prewhiten(year_8120, 1, 0.05)
  MK_yr_2170 <- general_test_mk_prewhiten(year_2170, 1, 0.05)
  MK_yr_2150 <- general_test_mk_prewhiten(year_2150, 1, 0.05)
  MK_yr_4170 <- general_test_mk_prewhiten(year_4170, 1, 0.05)
  
  #***************************************************************************************
  # Add to output table
  #***************************************************************************************
  #("SUBID", "MK_8170", "pVal_8170", "MK_8120", "pVal_8120", "MK_2170", "pVal_2170", "MK_2150", "pVal_2150", "MK_4170", "pVal_4170")
  outData2 <- data.frame(col1 = SUBID, col2 = MK_yr_8170$slope, col3 = MK_yr_8170$pVal,
					col4 = MK_yr_8120$slope, col5 = MK_yr_8120$pVal,
					col6 = MK_yr_2170$slope, col7 = MK_yr_2170$pVal,
					col8 = MK_yr_2150$slope, col9 = MK_yr_2150$pVal,
					col10 = MK_yr_4170$slope, col11 = MK_yr_4170$pVal)
  
  Tt <- rbind(Tt, outData2)
  
  
  
}

#***************************************************************************************
# Update column names in tables
#***************************************************************************************

colnames(Td) <- names_cols_Td
colnames(Tt) <- names_cols_Tt
