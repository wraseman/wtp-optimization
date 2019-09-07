# Step 3 - plot complete, monthly time series and create Figure 1 for Raseman et al. (2019)
# author: William Raseman

# clear environment
rm(list=ls()) 

visualize_ts <- function(additional_plots=FALSE) {
  
  # load packages
  library(tidyverse)  # modern R packages: ggplot2, dplyr, readr, etc.
  library(GGally)     # scatterplot matrices
  library(forecast)   # time series analysis and visualization  
  library(gridExtra)  # arrange multiple ggplot2 plots on a single plot
  library(seasonal)   # time series decomposition 
  
  # load user-defined libraries
  source("./lib/time-series-sim_lib.R")  # time series simulation library
  
  # read in source water time series data
  path <- "./data/source-water/02_create_ts/sw_ts.rds"
  wq.ts <- wq.fullnames.ts <- readr::read_rds(path)
  
  # Figure 1 - visualize water quality time series for source water
  # colnames(wq.fullnames.ts) <- c("Alkalinity (mg/L)", 
  #                                "pH", "Temperature (°C)", "TOC (mg/L)", 
  #                                "Hardness (mg/L)", "Turbidity (ntu)")  
  colnames(wq.fullnames.ts) <- c("Alk (mg/L)", "Hard (mg/L)",
                                 "pH", "Temp (°C)", "TOC (mg/L)",
                                 "Turb (ntu)")
  ## note: not sure if degree symbol will work correctly with all operating systems. 
  
  # source: https://danieljhocking.wordpress.com/2013/03/12/high-resolution-figures-in-r/
  tiff(filename = "./figures/figure-1.tiff", 
       height = 12, width = 17, units = 'cm', 
       compression = "lzw", res = 600)
  autoplot(wq.fullnames.ts, facet = TRUE, ylab="") %>% print  # need print statement to save within a function
  dev.off()

  ## note: must post-process this figure to add shaded regions associated with interpolated values
  ##       because I could not find a way to do this in a reproducible manner. 
  
  # show additional plots about time series
  if (additional_plots==TRUE) {
    
    ## create user-defined function
    create_seasonalplot_list <- function(data, FUN) {
      
      # inputs:
      #   data - time series data
      #   FUN  - ggplot-style time series plotting function (e.g. ggsubseriesplot(), 
      #          ggseasonplot())
      
      plot_list <- list()
      
      if (identical(FUN,ggseasonplot) | identical(FUN,ggsubseriesplot)) {
        for (i in 1:ncol(data)) {
          local({
            i <- i
            p <- FUN(data[,i], year.labels = TRUE, year.labels.left = TRUE) +
              ggtitle(colnames(data)[i])
            plot_list[[i]] <<- p  # add each plot into plot list
          })
        }
      }
      
      return(plot_list)
    }
    
    ## scatterplot matrix
    ggpairs(wq.ts %>% data.frame())
    
    ## seasonal series plots
    ggseason_plots <- create_seasonalplot_list(data = wq.ts, FUN=ggseasonplot)
    ggseason_plots
    
    ## seasonal subseries plots
    ggsubseries_plots <- create_seasonalplot_list(data = wq.ts, FUN=ggsubseriesplot)
    ggsubseries_plots
    
    ## x11 decomposition
    ggx11_plots <- list()
    
    for (i in 1:ncol(wq.ts)) {
      local({
        i <- i
        p <- wq.ts[,i] %>% seas(x11="") %>%
          autoplot() +
          ggtitle(colnames(wq.ts)[i])
        ggx11_plots[[i]] <<- p  # add each plot into plot list
      })
    }
    
    ggx11_plots
  }
}

# save function
save("visualize_ts", file="./lib/visualize_ts.RData")

# run script
visualize_ts(additional_plots = TRUE) # uncomment to run script
