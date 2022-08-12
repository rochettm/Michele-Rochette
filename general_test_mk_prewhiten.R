# Input description
# timeseries   vector    timeSeriesies of data (can contain NAs, but needs spacing preserved)
# acLag        integer   number of timesteps used for lag
# alpha        double    significance of Mann-Kendall test to be evaluated 

general_test_mk_prewhiten <- function ( timeSeries, acLag, alpha ) {
  # Method for correcting trend tests summarized in pg. 1821 - 1822:
  # Yue et al., 2002, Hydrological Processes, DOI: 10.1002/hyp.1095,
  # "The influence of autocorrelation on the ability to detect trend..." 
   
  # load libraries for Sen's slope and Mann-Kendall testing
  library(trend)
  library(Kendall)
  
  # entries in timesries
  n <- length(timeSeries)
  
  # timeseries mean
  muTS <- mean(timeSeries,
               na.rm = TRUE)
  
  # slope threshold in percent (slope over period / mean over period)
  # slopes shallower than this will be excluded and the Mann-Kendall test will be refused
  slpThr <- 0.01 * mean(timeSeries,
                        na.rm = TRUE) / length(timeSeries)
  
  # critical value of AC coefficient 
  acCrit <- abs(qnorm(alpha / 2)) / sqrt(n - acLag)
  
  # compute AC coefficient for k steps of lag (specified in input)
  ack <- (1 / (n - acLag)) * sum((timeSeries[(1 + acLag):n] - muTS) * (timeSeries[1:(n - acLag)] - muTS)) / (1 / n) / sum((timeSeries - muTS)^2)
  
  # compute Mann-Kendall p-value of non-transformed data
  pValue <- mk.test(timeSeries)$p.value
  # print(paste0("Original p-value: ",pValue))
  
  # compute the raw Sen Thiel slope
  senSlp <- as.numeric(sens.slope(timeSeries)$estimates)
  
  # compute Sen's slope intercept for the same period
  senInt <- median(timeSeries - senSlp * seq(1, n),
                   na.rm = TRUE)
  
  # IF Sen's slope is too small (will corrupt significance test)
  if (abs(senSlp) <= slpThr) {
    # MK p-value is not detectable due to tiny slope, return no-data value
    pValue <- NA
    
  # ELSE sen's slope is not too small to detect significance   
  } else {
    # de-trend data by removing linear equivalent of Thiel-Sen slope
    timeSeriesDT <- timeSeries - senSlp * seq(1, n)
    
    # de-trended timeseries mean
    muDT <- mean(timeSeriesDT,
                 na.rm = TRUE)
    
    # compute autocorreatlion component of timeseries for lag k
    ackDT <- (1 / (n - acLag)) * sum((timeSeriesDT[(1 + acLag):n] - muDT) * (timeSeriesDT[1:(n - acLag)] - muDT)) / (1 / n) / sum((timeSeriesDT - muDT)^2)
    
    # pre-whiten data by removing k-timestep lagged autocorrelation component
    timeSeriesPW <- c(timeSeriesDT[1:acLag],
                      timeSeriesDT[(1 + acLag):n] - ackDT * timeSeriesDT[1:(n - acLag)])
    
    # re-apply the origianl trend to the pre-whitened data (blended timeseries)
    timeSeriesBl <- timeSeriesPW + senSlp * seq(1,
                                                n)
    
    # re-compute Mann-Kendall test on blended timeseries
    pValue <- mk.test(timeSeriesBl)$p.value
  
  # END test if slope is too shallow to trust Mann-Kendall test
  }
  
  # print(paste0("Final    p-value: ",pValue))
  # RETURN p-value, Sen slope, and intercept for Sen slope
  return(list(pVal = pValue, slope = senSlp, interc = senInt))
} 
