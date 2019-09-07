# Step 1 - import, clean, and tidy source water quality data then export .rds file
# author: William Raseman

# clear environment
rm(list=ls()) 

import_clean <- function(show_plots=FALSE) {
# show_plots <- F
  # load packages
  library(tidyverse)  # modern R packages: ggplot2, dplyr, readr, etc.
  library(readxl)     # import data from Excel files
  library(lubridate)  # date-time data manipulation
  library(forecast)   # time series data analysis and visualization
  library(stringr)    # string manipulation
  library(zoo)        # date-time data manipulation
  library(padr)       # fill in records for time points where observations are absent
  
  # import local functions (addition for wtp-optimize)
  source("import_clean_hardness.R")  
  source("import_clean_turbidity.R")
  
  # read in data
  path = "./data/source-water/01_import_clean/"
  
  # read in water quality data #1 (total organic carbon and alkalinity)
  wq1_df <- read_excel(path = str_c(path, "sw_plant-intake_toc-alk.xlsx"))
  
  ##  There are two sampling locations: one labeled "plant", the other labeled "intake".
  ##  The only difference between these should be any water quality changes that occur due
  ##  to transporting the water from intake to the plant. It should be minimal but these
  ##  datasets are compared to verify.
  
  # visualize the time series for toc for both intake and plant
  param.1 <- c("toc", "alk")
  
  for (i in param.1) {
    filt_df <- filter(wq1_df, parameter==i) # filter the data for each parameter
    
    p <- ggplot(data=filt_df,
                aes(x=date, y=value, color=samp_loc)) +
      geom_point() +
      ggtitle(i)
    if (show_plots==TRUE) print(p)
  }
  
  ##  Hypothesis supported that the "plant" and "intake" samples generally the same (with the
  ##  exception on an outlier in 2013)
  
  ## count how much data there is
  filter(wq1_df, parameter=="toc", samp_loc=="plant") %>% nrow
  filter(wq1_df, parameter=="toc", samp_loc=="intake") %>% nrow
  filter(wq1_df, parameter=="alk", samp_loc=="plant") %>% nrow
  filter(wq1_df, parameter=="alk", samp_loc=="intake") %>% nrow
  
  # tidyr::spread() data and fill in all days for period of record
  ## what is the starting date? what is the ending date?
  start.date1 <- min(wq1_df$date)
  end.date1 <- max(wq1_df$date)
  range.date1 <- seq(from=start.date1, to=end.date1, by="day")
  
  ## spread the data
  ## source: https://groups.google.com/forum/#!topic/manipulatr/oos-1t-e25g
  ## source: https://stackoverflow.com/questions/39053451/using-spread-with-duplicate-identifiers-for-rows
  wq1_pad <- wq1_df %>% 
    select(-units) %>%  # remove units
    unite(param_samp, parameter, samp_loc) %>% 
    distinct(date, param_samp, .keep_all = TRUE) %>%  # first value is kept  %>% # remove duplicates
    spread(param_samp, value) %>% 
    pad(interval="day") 
  
  ## figure out what percentage of values are missing for daily values
  sapply(wq1_pad, function(x) sum(is.na(x)))/nrow(wq1_pad)  
  
  ## troubleshooting: there should be 1302 data points for in-plant TOC. 
  
  wq1_df <- wq1_pad %>% # add missing dates to dataset
    gather(parameter, value, -date) %>% # undo spread
    separate(parameter, into=c("parameter", "samp_loc"), sep = "_")
  
  ## visualize daily values of TOC and alkalinity
  for (i in param.1) { 
    filt_df <- filter(wq1_df, parameter==i) # filter the data for each parameter
    
    p <- ggplot(data=filt_df, 
                aes(x=date, y=value, color=samp_loc)) + 
      geom_point() + 
      geom_line() + 
      ggtitle(paste("Daily values of", i))
    if (show_plots == TRUE) print(p)
  }
  
  # aggregate water quality #1 to monthly mean (keeping intake and plant separate still)
  options(warn=-1) # turn off warnings for zoo::as.yearmon, source: https://github.com/business-science/sweep/issues/5
  wq1_mon <- mutate(wq1_df, year=year(date), month=month(date)) %>% 
    group_by(samp_loc, parameter, year, month) %>%
    summarize(mean_monthly_value = mean(value, na.rm=TRUE)) %>%
    mutate(date = zoo::as.yearmon(paste(year, month), "%Y %m")) %>% # create new date column
    arrange(date)  # order by date
  
  # plot monthly time series
  for (i in param.1) {  
    filt_mon <- filter(wq1_mon, parameter==i) # filter the data for each parameter
    
    p <- ggplot(data=filter(wq1_mon, parameter==i),  
                aes(x=date, y=mean_monthly_value, color=samp_loc)) + 
      geom_point() + 
      geom_line() + 
      ggtitle(paste("Monthly values of", i))
    if (show_plots == TRUE) print(p)
  }
  options(warn=0) # warnings back on
  
  # read in water quality data #2 (pH and temperature)
  
  wq2_df <- read_excel(path = str_c(path, "sw_plant-intake_pH-temp.xlsx"),
                       range = cell_cols("A:E")) %>%
    mutate(datetime=as.POSIXct(strptime(datetime, "%Y-%m-%d %H:%M"))) %>%
    mutate(date=as.Date(datetime))
  
  # Perform QA/QC on the data
  
  ## figure out what percentage of values are missing for 15-min values
  sapply(wq2_df, function(x) sum(is.na(x)))/nrow(wq2_df)  
  
  if (show_plots == TRUE) {
    ## plot Q-Q plots of the data (warning: takes quite a bit of time to plot)
    ggplot(wq2_df, aes(sample=pH_intake)) +
      stat_qq()  ## suspicious values at 14.0 (try keeping only 6 < pH < 10)
    
    ggplot(wq2_df, aes(sample=pH_plant)) +
      stat_qq()  ## suspicious values at 14.0 (try keeping only 6 < pH < 10)
    
    ggplot(wq2_df, aes(sample=temp_intake)) +
      stat_qq()  ## no suspicious values
    
    ggplot(wq2_df, aes(sample=temp_plant)) +
      stat_qq()  # suspicious repeated values at 30 deg C (try eliminating values > 28)
  }
  
  ## based on results, remove suspicious values
  ### pH at intake
  filter(wq2_df, (pH_intake < 6)|(pH_intake > 11))
  wq2_df <- mutate(wq2_df, pH_intake=if_else((pH_intake < 6)|(pH_intake > 11), NA_real_, pH_intake))
  
  
  ### pH at plant
  filter(wq2_df, (pH_plant < 6)|(pH_plant > 11))
  wq2_df <- mutate(wq2_df, pH_plant=if_else((pH_plant < 6)|(pH_plant > 11), NA_real_, pH_plant))
  
  ### temp at intake
  # no suspicious values!
  
  ### temp at plant
  filter(wq2_df, temp_plant > 28)  # all 30.0 deg C. seems like an error
  wq2_df <- mutate(wq2_df, temp_plant=if_else(temp_plant > 28, NA_real_, temp_plant))
  
  ## max/min values for each column
  summarize_all(wq2_df, funs(max(., na.rm=TRUE)))
  summarize_all(wq2_df, funs(min(., na.rm=TRUE)))
  
  # tidy data
  wq2_tidy <- wq2_df %>% # add missing dates to dataset
    gather(key=param_loc, value=value, -date, -datetime) %>%
    separate(col=param_loc, into=c("parameter", "samp_loc"), sep = "_")
  
  # aggregate to daily
  wq2_tidy<- wq2_tidy %>%
    group_by(samp_loc, parameter, date) %>%
    summarize(mean_daily_value = mean(value, na.rm=TRUE))
  
  # visualize daily data
  param.2 <- c("temp", "pH")
  
  for (i in param.2) { 
    filt2_df <- filter(wq2_tidy, parameter==i) # filter the data for each parameter
    
    p <- ggplot(data=filt2_df, 
                aes(x=date, y=mean_daily_value, color=samp_loc)) + 
      geom_point() + 
      ggtitle(paste("Daily values of", i))
    if (show_plots == TRUE) print(p)
  }
  ## notice that there are difference in between the "intake" and "plant" locations 
  ## for consistency in analysis, just choose to use "plant" data
  
  # aggregate to monthly
  ## aggregate water quality #2 to monthly mean
  wq2_mon <- mutate(wq2_tidy, year=year(date) %>% as.integer, 
                    month=month(date) %>% as.integer)
  
  
  ## group, summarize, and create date column
  options(warn=-1) # turn off warnings for zoo::as.yearmon, source: https://github.com/business-science/sweep/issues/5
  wq2_mon <- group_by(wq2_mon, samp_loc, parameter, year, month) %>%
    summarize(mean_monthly_value = mean(mean_daily_value, na.rm = TRUE)) %>% 
    mutate(date=zoo::as.yearmon(paste(year, month), "%Y %m"))  # create new date column for visualization

  sapply(wq2_mon, function(x) sum(is.na(x)))  # count missing values
  
  # plot monthly time series
  for (i in param.2) { 
    filt2_mon <- filter(wq2_mon, parameter==i) # filter the data for each parameter
    
    p <- ggplot(data=filt2_mon, 
                aes(x=date, y=mean_monthly_value, color=samp_loc)) + 
      geom_point() + 
      geom_line() +
      ggtitle(paste("Monthly values of", i))
    if (show_plots == TRUE) print(p)
  }
  options(warn=0) # warnings back on
  
  # save cleaned monthly dateframe for all water quality data (both wq_df #1 and #2)
  wq12_mon <- rbind(wq1_mon, wq2_mon) %>% 
    filter(samp_loc == "plant")  %>% # only keep data from "plant" sampling location
    ungroup() %>%
    select(-samp_loc, -date) %>%
    mutate(year = year %>% as.integer, month = month %>% as.integer)
  
  # addition for wtp-optimize: import hardness and turbidity
  wq3_mon <- import_clean_turbidity()
  wq4_mon <- import_clean_hardness() %>% 
    ungroup() %>%
    mutate(year = year %>% as.integer, month = month %>% as.integer)
    
  wq1234_mon <- rbind(wq12_mon, wq3_mon, wq4_mon)  # append all water quality parameters together
  
  clean.path <- "./data/source-water/02_create_ts/"
  file.name <- "sw_cleaned-agg-monthly.rds"
  
  write_rds(wq1234_mon, str_c(clean.path, file.name))
}

# save function
save("import_clean", file="./lib/import_clean.RData")

# run script
# import_clean()  # uncomment to run script



