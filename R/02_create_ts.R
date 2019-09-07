# Step 2 - turn tidied monthly dataframe of water quality data into time series data (then export .rds file)
# source: provided by the City of Fort Collins
# author: William Raseman

# clear environment
rm(list=ls()) 

create_ts <- function() {
  
  # load packages
  library(tidyverse)  # modern R packages: ggplot2, dplyr, readr, etc.
  library(readxl)     # import data from Excel files
  library(lubridate)  # date-time data manipulation
  library(forecast)   # time series data analysis and visualization
  
  # read water quality dataframe (.rds file)
  clean.path <- "./data/source-water/02_create_ts/"
  file.name <- "sw_cleaned-agg-monthly.rds"
  
  mon_df <- read_rds(path=str_c(clean.path, file.name))
  
  # set period (or frequency) to 12 for monthly data
  m <- 12
  
  # get starting and ending month and year of time series 
  # start.mon <- mon_df$month[1]
  start.mon <- 1  # force starting month to be January
  start.yr  <- mon_df$year[1]
  # end.mon <- mon_df$month[length(mon_df$month)]  
  end.mon <- 12  # force ending month to be December
  end.yr <- mon_df$year[length(mon_df$year)]
  
  # spread dataframe (untidy) to convert to time series
  mon_df <- spread(mon_df, parameter, mean_monthly_value) %>%
    select(-year, -month)
  
  # convert to time series format
  mon_ts <- ts(mon_df, start = c(start.yr, start.mon),
               end = c(end.yr, end.mon),
               frequency = m)
  
  # interpolate between missing values
  int_ts <- mon_ts
  for (i in (1:ncol(mon_df))) {
    int_ts[,i] <- na.interp(int_ts[,i])  # na.interp() ref - https://www.rdocumentation.org/packages/forecast/versions/8.1/topics/na.interp
  }
  
  # visually, the interpolation seems reasonable
  autoplot(mon_ts, facet = TRUE)
  autoplot(int_ts, facet = TRUE)
  
  # finalize time series object
  sw_ts <- int_ts
  
  write_rds(sw_ts, "./data/source-water/02_create_ts/sw_ts.rds")
}

# save function
save("create_ts", file="./lib/create_ts.RData")

# run script
# create_ts() # uncomment to run script
