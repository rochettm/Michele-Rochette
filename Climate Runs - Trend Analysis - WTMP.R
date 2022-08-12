
# Load prewhitening code and required libraries
source('~/Thesis/9 Climate Run Analysis/general_test_mk_prewhiten.R')
library(readr)
library(trend)
library(Kendall)
library(lubridate)


#***********************************************************************************************
# Set up basic information
# Doing trend analysis for climate run 1, in the Nelson Churchill basin, for wtmp variable (stream temperature)
#***********************************************************************************************
runOfInt <- "1"
varOfInt <- "wtmp" 

#***********************************************************************************************
# List of stations in the Nelson Churchill Basin and create path to find files later
#***********************************************************************************************
Stns <- read.delim("G:/Michele/9 Climate Run Analysis/3 R Codes/pmsf_NC.txt", header = FALSE)
runList <- list.files(path = "G:/Michele/9 Climate Run Analysis/1 Climate Run Results",pattern = "NC - 1 - ")

#***********************************************************************************************
#Set up tables to output data
#***********************************************************************************************
names_cols_mn <- c("SUBID", "MK_8119", "pVal_8119", "MK_2070", "pVal_2070", "MK_2050", "pVal_2050", 
                       "MK_4070", "pVal_4070", "MK_8170", "pVal_8170")
T_mn <- data.frame(matrix(ncol = 11, nrow = 0))

names_cols_yr <- c("SUBID", "MK_yr_8119", "pVal_yr_8119", "MK_yr_2070", "pVal_yr_2070", "MK_yr_2050", "pVal_yr_2050", 
                       "MK_yr_4070", "pVal_yr_4070", "MK_yr_8170", "pVal_yr_8170")
T_yr <- data.frame(matrix(ncol = 11, nrow = 0))

names_cols_szns <- c("SUBID", "DJF", "pVal_DJF", "MAM", "pVal+MAM", "JJA", "pVal_JJA", "SON", "pVal_SON")
T_szns <- data.frame(matrix(ncol = 9, nrow = 0))

names_cols_pettitt <- c("SUBID", "DATE", "pVal", "DATE_fut", "pVal_fut", "error")
T_pettitt <- data.frame(matrix(ncol = 6, nrow = 0))


#***********************************************************************************************
#run the following for every station in the Nelson-Churchill (or whichever stations identified in Stns)
#***********************************************************************************************
for(i in seq(1,nrow(Stns))){ 
  

  SUBID <- Stns[i, ]
  
  inData <- read_delim(paste0("G:/Michele/9 Climate Run Analysis/1 Climate Run Results/", 
                              runList,"/", SUBID, ".txt"), delim = "\t", col_names = TRUE, col_types = cols())
  inData <- inData[-1, ]
  
  # Prepare Data - add column for year and month
  inData[, match("DATE", colnames(inData))] <- as.Date(inData$DATE)
  inData[, match(varOfInt, colnames(inData))] <- as.numeric(inData$wtmp)
  inData$YEAR <- as.numeric(strftime(inData$DATE, format = "%Y"))
  inData$MONTH <- as.numeric(strftime(inData$DATE, format = "%m"))
  inData$code <- inData$YEAR + inData$MONTH/100
  
  # Make temps<0 = 0
  inData$wtmp <- ifelse(inData$wtmp < 0, 0, inData$wtmp)
  
  #***************************************************************************************
  # Calcs for MK trend test - for different time periods
  #***************************************************************************************
  #Get monthly data
  month_data <- aggregate(data=inData, as.numeric(wtmp) ~ code, mean)
  names(month_data)[names(month_data) == 'as.numeric(wtmp)'] <- 'wtmp_avg'
  data_8119 <- month_data[(1:468),2]
  data_2070 <- month_data[(469:1080),2]
  data_2050 <- month_data[(469:840),2]
  data_4070 <- month_data[(721:1080),2]
  data_8170 <- month_data[(1:1080), 2]
  
  #Perform trend tests
  MK_8119 <- general_test_mk_prewhiten(data_8119, 1, 0.05)
  MK_2070 <- general_test_mk_prewhiten(data_2070, 1, 0.05)
  MK_2050 <- general_test_mk_prewhiten(data_2050, 1, 0.05)
  MK_4070 <- general_test_mk_prewhiten(data_4070, 1, 0.05)
  MK_8170 <- general_test_mk_prewhiten(data_8170, 1, 0.05)

  
  #Add to Table
  outData1 <- data.frame(col1 = SUBID, col2 = MK_8119$slope, col3 = MK_8119$pVal, 
                        col4 = MK_2070$slope, col5 = MK_2070$pVal,
                        col6 = MK_2050$slope, col7 = MK_2050$pVal,
                        col8 = MK_4070$slope, col9 = MK_4070$pVal,
                        col10 = MK_8170$slope, col11 = MK_8170$pVal)
                         
  T_mn <- rbind(T_mn, outData1)
  
  #Get yearly data
  year_data <- aggregate(data=inData, as.numeric(wtmp) ~ YEAR, mean)
  names(year_data)[names(year_data) == 'as.numeric(wtmp)'] <- 'wtmp_avg'
  year_8119 <- year_data[(1:39),2]
  year_2070 <- year_data[(40:90),2]
  year_2050 <- year_data[(40:70),2]
  year_4070 <- year_data[(60:90),2]
  year_8170 <- year_data[(1:90), 2]
  
  MK_yr_8119 <- general_test_mk_prewhiten(year_8119, 1, 0.05)
  MK_yr_2070 <- general_test_mk_prewhiten(year_2070, 1, 0.05)
  MK_yr_2050 <- general_test_mk_prewhiten(year_2050, 1, 0.05)
  MK_yr_4070 <- general_test_mk_prewhiten(year_4070, 1, 0.05)
  MK_yr_8170 <- general_test_mk_prewhiten(year_8170, 1, 0.05)
  
  outData2 <- data.frame(col1 = SUBID, col2 = MK_yr_8119$slope, col3 = MK_yr_8119$pVal, 
                        col4 = MK_yr_2070$slope, col5 = MK_yr_2070$pVal,
                        col6 = MK_yr_2050$slope, col7 = MK_yr_2050$pVal,
                        col8 = MK_yr_4070$slope, col9 = MK_yr_4070$pVal,
                        col10 = MK_yr_8170$slope, col11 = MK_yr_8170$pVal)
  
  T_yr <- rbind(T_yr, outData2)
  
  
  #***************************************************************************************
  # Calcs for seasonality 
  #***************************************************************************************
  # Calculate monthly averages
  inData$YR_MN <- inData$YEAR + inData$MONTH/100
  month_data <- aggregate(data = inData, as.numeric(wtmp) ~ YR_MN, mean)
  names(month_data)[names(month_data) == 'as.numeric(wtmp)'] <- 'wtmp_month_avg'
  month_data$MONTH <- seq(1,12)
  month_data$YEAR <- round(month_data$YR_MN)
  
  # get one data point for each year for each season per year
  szns <- c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1)
  month_data$szn <- szns
  month_data$szn_mn <- month_data$YEAR + month_data$szn/10
  szn_data <- aggregate(data = month_data, as.numeric(wtmp_month_avg) ~ szn_mn, mean)
  
  # get series
  szn_data$szn <- seq(1,4)
  names(szn_data)[names(szn_data) == 'as.numeric(wtmp_month_avg)'] <- 'szn_avg'
  idx_DJF <- which(szn_data$szn == 1)
  DJF <- szn_data$szn_avg[idx_DJF]
  idx_MAM <- which(szn_data$szn == 2)
  MAM <- szn_data$szn_avg[idx_MAM]
  idx_JJA <- which(szn_data$szn == 3)
  JJA <- szn_data$szn_avg[idx_JJA]
  idx_SON <- which(szn_data$szn == 4)
  SON <- szn_data$szn_avg[idx_SON]
  
  # Put season series through MK tests
  DJF_MK <- general_test_mk_prewhiten(DJF, 1, 0.05)
  MAM_MK <- general_test_mk_prewhiten(MAM, 1, 0.05)
  JJA_MK <- general_test_mk_prewhiten(JJA, 1, 0.05)
  SON_MK <- general_test_mk_prewhiten(SON, 1, 0.05)
  
  # Add to Table
  outData <- data.frame(col1 = SUBID, col2 = DJF_MK$slope, col3 = DJF_MK$pVal, 
                        col4 = MAM_MK$slope, col5 = MAM_MK$pVal, col6 = JJA_MK$slope,
                        col7 = JJA_MK$pVal, col8 = SON_MK$slope, col9 = SON_MK$pVal)
  T_szns <- rbind(T_szns, outData)
  
  
  #***************************************************************************************
  # Calcs for Change Point - use for stations where you want change point analysis
  #***************************************************************************************
  
  # pettitt.test
  #petit <- pettitt.test(inData$wtmp)
  #petit_fut <- pettitt.test(inData$wtmp[14245:32872])
  #K <- petit$estimate
  #K_fut <- petit_fut$estimate
  
  #"SUBID", "DATE", "pVal", "DATE_fut", "pVal_fut"
  
  # Pettitt values into table
  #outData <- data.frame(col1 = SUBID, col2 = inData$DATE[K], col3 = petit$p.value, 
  #                      col4 = inData$DATE[K_fut], col5 = petit_fut$p.value, col6 = " ") 
  #if(nrow(outData)>1){outData <- outData[1, ]
  #outData$col6 = "nrow > 1"}
  
  #final_data <- rbind(final_data, outData)
  
}

print(Sys.time())

#**********************************************************************
#Update column names
colnames(T_mn) <- names_cols_mn
colnames(T_yr) <- names_cols_yr
colnames(T_szns) <- names_cols_szns
colnames(T_pettitt) <- names_cols_pettitt