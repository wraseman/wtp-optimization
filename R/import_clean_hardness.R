# Import hardness data

# Step 1 - import, clean, and tidy source water quality data then export .rds file
# author: William Raseman

import_clean_hardness <- function(show_plots = FALSE) {
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
  
  # read in water quality data #1 (total organic carbon and alkalinity)
  wq4_df <-
    read_excel(path = str_c(path, "sw_plant-intake_hardness.xlsx")) %>%
    mutate(datetime = as.POSIXct(strptime(datetime, "%Y-%m-%d %H:%M"))) %>%
    mutate(date = as.Date(datetime))
  
  # visualize the time series for hardness
  if (show_plots == TRUE) {
    p <- ggplot(data = wq4_df,
                aes(x = date, y = hardness)) +
      geom_point() +
      ggtitle("Hardness (mg/L CaCO3)")
    
    print(p)
  }
  
  # # tidyr::spread() data and fill in all days for period of record
  # ## what is the starting date? what is the ending date?
  start.date4 <- min(wq4_df$date)
  end.date4 <- max(wq4_df$date)
  range.date4 <- seq(from = start.date4, to = end.date4, by = "month")
  
  # aggregate water quality #1 to monthly mean
  options(warn = -1) # turn off warnings for zoo::as.yearmon, source: https://github.com/business-science/sweep/issues/5
  wq4_mon <- mutate(wq4_df, year = year(date), month = month(date)) %>%
    group_by(year, month) %>%
    summarize(mean_monthly_value = mean(hardness, na.rm = TRUE)) %>%
    mutate(date = zoo::as.yearmon(paste(year, month), "%Y %m")) %>% # create new date column
    arrange(date)  # order by date
  
  # plot monthly time series
  p <- ggplot(data = wq4_mon,
              aes(x = date, y = mean_monthly_value)) +
    geom_point() +
    geom_line() +
    ggtitle("Monthly values of Hardness")
  
  if (show_plots == TRUE)
    print(p)
  
  options(warn = 0) # warnings back on
  
  # missing values: 3
  count.missing <- length(range.date4) - nrow(wq4_mon)
  
  # match format of other monthly water quality dataframes
  wq4_mon <- mutate(wq4_mon, parameter="hard") %>%
    select(parameter, year, month, mean_monthly_value)
  
  return(wq4_mon)
}
