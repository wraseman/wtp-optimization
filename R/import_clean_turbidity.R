# Import turbidity data

# Step 1 - import, clean, and tidy source water quality data then export .rds file
# author: William Raseman

import_clean_turbidity <- function(show_plots = FALSE) {
  # load packages
  library(tidyverse)  # modern R packages: ggplot2, dplyr, readr, etc.
  library(readxl)     # import data from Excel files
  library(lubridate)  # date-time data manipulation
  library(forecast)   # time series data analysis and visualization
  library(stringr)    # string manipulation
  library(zoo)        # date-time data manipulation
  library(padr)       # fill in records for time points where observations are absent
  
  # read in data
  path = "./data/source-water/01_import_clean/"
  
  ##  There are two sampling locations: one labeled "plant", the other labeled "intake".
  ##  The only difference between these should be any water quality changes that occur due
  ##  to transporting the water from intake to the plant. It should be minimal but these
  ##  datasets are compared to verify.
  
  # read in water quality data #2 (pH and temperature)
  wq3_df <-
    read_excel(
      path = str_c(path, "sw_plant-intake_turbidity.xlsx"),
      range = cell_cols("A:C")
    ) %>%
    mutate(datetime = as.POSIXct(strptime(datetime, "%Y-%m-%d %H:%M"))) %>%
    mutate(date = as.Date(datetime))
  
  # Perform QA/QC on the data
  
  ## figure out what percentage of values are missing for 15-min values
  sapply(wq3_df, function(x)
    sum(is.na(x))) / nrow(wq3_df)
  
  if (show_plots == TRUE) {
    ## plot CDF plots of the data (warning: takes quite a bit of time to plot)
    ggplot(wq3_df, aes(sample = turb_intake)) +
      stat_qq()
    
    ggplot(wq3_df, aes(sample = turb_plant)) +
      stat_qq()
    
    ## looks like there is a strange number of measurements of exactly 100.0 ntu
    ## need to investigate this...
    
    ## look at difference between intake and plant
    ggplot(
      data = mutate(wq3_df, delta_turbidity = turb_intake - turb_plant),
      aes(x = date, y = delta_turbidity)
    ) +
      geom_point() +
      ggtitle("Difference between plant and intake measurements")
  }
  
  
  ## although the difference between the intake and the plant and the repeated 100.0 ntu
  ##  values seem odd, I am not confident enough to say that these data are not representative
  ##  measurements of actual conditions. For that reason, I will not remove any data.
  
  ## max/min values for each column
  summarize_all(wq3_df, funs(max(., na.rm = TRUE)))
  summarize_all(wq3_df, funs(min(., na.rm = TRUE)))
  
  # tidy data
  wq3_tidy <- wq3_df %>% # add missing dates to dataset
    gather(key = param_loc, value = value,-date,-datetime) %>%
    separate(
      col = param_loc,
      into = c("parameter", "samp_loc"),
      sep = "_"
    )
  
  # visualize turbidity at plant and intake
  if (show_plots == TRUE) {
    ggplot(data = wq3_tidy,
           aes(x = date, y = value, color = samp_loc)) +
      geom_point() +
      ggtitle("Turbidity")
  }
  
  # aggregate to daily
  wq3_tidy <- wq3_tidy %>%
    group_by(samp_loc, parameter, date) %>%
    summarize(mean_daily_value = mean(value, na.rm = TRUE))
  
  # visualize daily data
  filt3_df <-
    filter(wq3_tidy, parameter == "turb") # filter the data for each parameter
  
  if (show_plots == TRUE) {
    p <- ggplot(data = filt3_df,
                aes(x = date, y = mean_daily_value, color = samp_loc)) +
      geom_point() +
      ggtitle("Daily values of Turbidity")
  }
  
  ## notice that there are difference in between the "intake" and "plant" locations
  ## for consistency in analysis, just choose to use "plant" data
  
  # aggregate to monthly
  ## aggregate water quality #3 to monthly mean
  wq3_mon <- mutate(
    wq3_tidy,
    year = year(date) %>% as.integer,
    month = month(date) %>% as.integer
  )
  
  ## group, summarize, and create date column
  options(warn = -1) # turn off warnings for zoo::as.yearmon, source: https://github.com/business-science/sweep/issues/5
  wq3_mon <- group_by(wq3_mon, samp_loc, parameter, year, month) %>%
    summarize(mean_monthly_value = mean(mean_daily_value, na.rm = TRUE)) %>%
    mutate(date = zoo::as.yearmon(paste(year, month), "%Y %m"))  # create new date column for visualization
  
  sapply(wq3_mon, function(x)
    sum(is.na(x)))  # count missing values
  
  # plot monthly time series
  filt3_mon <-
    filter(wq3_mon, parameter == "turb") # filter the data for each parameter
  
  if (show_plots == TRUE) {
    p <- ggplot(data = filt3_mon,
                aes(x = date, y = mean_monthly_value, color = samp_loc)) +
      geom_point() +
      geom_line() +
      ggtitle("Monthly values of Turbidity")
    
    print(p)
  }
  
  # in this instance, I think the intake data is of better quality
  wq3_mon <- wq3_mon %>%
    filter(samp_loc == "intake")  %>% # only keep data from "intake" sampling location
    ungroup() %>%
    select(-samp_loc,-date)
  
  options(warn = 0) # warnings back on
  
  return(wq3_mon)
}
