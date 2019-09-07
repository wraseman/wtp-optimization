# Step 5 - Visualize simulation results and generate figures for manuscript 
# author: William Raseman

# clear environment
rm(list=ls()) 

visualize_statistics <- function(innov, data.type) {
  
  # inputs
  #   innov: random innovations on or off (TRUE/FALSE)
  #   data.type: specify source water ("sw") or mine meteorology ("mine") dataset
  
  # load packages
  library(tidyverse)  # modern R packages: ggplot2, dplyr, readr, etc.
  library(forecast)   # time series analysis and visualization  
  library(gridExtra)  # arrange multiple ggplot2 plots on a single plot
  library(zoo)        # rolling mean calculations
  
  # load user-defined libraries and functions
  source("./lib/time-series-sim_lib.R")  # time series simulation library
  
  fig.resolution = 300  # 300 dpi required for color art
  
  # read in time series data
  if (data.type == "sw") {
    
    ## source water quality
    path <- "./data/source-water/02_create_ts/sw_ts.rds"
    ts.data <- readr::read_rds(path)  
    
  } else if (data.type == "mine"){
    
    ## precipitation and temperature from mine
    ts.data <- read_rds(path="./data/mine/04_simulate_kNN/mine_ts-data_with-noise.RData")
  }
  
  ## simulated data 
  var.names <- colnames(ts.data)
  
  if (data.type == "sw") {
    # var.full.names <- c("Alkalinity (mg/L)", "pH", "Temperature (°C)", "TOC (mg/L)")
    var.full.names <- c("Alk", "Hard", "pH", "Temp", "TOC", "Turb")
  } else if (data.type == "mine") {
    var.full.names <- c("Temperature (°C)", "Precipitation (mm)")
  }
  
  nvars <- length(var.names)
  model.type <- read.path <- c()
  sim.list <- vector(mode="list", length=nvars)
  
  ## read in data for each variable
  for (i in 1:nvars) {
    model.type <- str_c("kNN", str_c("innov", "-", innov), sep="_")
    read.path <- str_c("./data/source-water/04_simulate_kNN/", model.type, 
                       "_", 
                       data.type,
                       "_",
                       var.names[i],
                       ".rds")  # name simulation based on model, data type, and variable 
    sim.list[[i]] <- read_rds(read.path)
  }
  
  n.sims <- dim(sim.list[[1]])[1]/length(ts.data[,1])  # get number of simulations from multivariate data
  
  # Figure 2 and S1 - visualize boxplots of observed and simulated data
  
  p.list <- vector(mode="list", length=nvars)
  
  for (i in 1:nvars) {
    p.list[[i]] <- viz_obs_sim(ts.data[,i], sim.list[[i]]) + 
      ylab(var.full.names[i]) + 
      xlab("Month") + 
      theme(legend.position="none")
  }
  
  ## save figure
  if (data.type == "sw") {
    if (innov == TRUE) {
      tiff(filename = "./figures/figure-2.tiff", 
           height = 12, width = 17, units = 'cm', 
           compression = "lzw", res = fig.resolution)
      grid.arrange(grobs = p.list)
      dev.off()
    } else if (innov == FALSE) {
      tiff(filename = "./figures/figure-S1.tiff", 
           height = 12, width = 17, units = 'cm', 
           compression = "lzw", res = fig.resolution)
      grid.arrange(grobs = p.list)
      dev.off()
    }
  }

  
  # Figure 3, S2, and S3 - visualize mean, standard deviation, minimum, and maximum of Total Organic Carbon
  
  if (data.type == "sw" && innov == TRUE) {
    tiff(filename = "./figures/figure-3.tiff", 
         height = 12, width = 17, units = 'cm', 
         compression = "lzw", res = fig.resolution)
    
    viz_ts_sample_stats(ts.data[,"toc"], sim.list[[5]], title=var.full.names[5])
    dev.off()
  } else if (data.type == "mine") {
    
    if (innov == TRUE) {
      tiff(filename = "./figures/figure-S2.tiff", 
           height = 12, width = 17, units = 'cm', 
           compression = "lzw", res = fig.resolution)
    } else if (innov == FALSE) {
      tiff(filename = "./figures/figure-S3.tiff", 
           height = 12, width = 17, units = 'cm', 
           compression = "lzw", res = fig.resolution)
    }
    
    viz_ts_sample_stats(ts.data[,"temp_C_plus_noise"], sim.list[[1]], title=var.full.names[1])
    dev.off()
  }

  
  # Figure 4 - visualize pairwise correlation between all variables
  
  # TODO: need to relabel these graphs
  
  if (data.type == "sw" && innov == TRUE) {
    tiff(filename = "./figures/figure-4.tiff", 
         height = 12, width = 17, units = 'cm', 
         compression = "lzw", res = fig.resolution)
    viz_pair_corr(y=ts.data, data=sim.list, data.type=data.type)
    dev.off()
  } 

  # Figure 5 - visualize lag-1 correlation for each variable
  
  # TODO: need to redo these graphs (currently, not accurate)
  
  if (data.type == "sw" && innov == TRUE) {
    tiff(filename = "./figures/figure-5.tiff", 
         height = 12, width = 17, units = 'cm', 
         compression = "lzw", res = fig.resolution)
    viz_ts_lag1(y=ts.data, data=sim.list, data.type=data.type, var.names=var.full.names)
    dev.off()
  } else if (data.type == "mine"  && innov == TRUE) {
    tiff(filename = "./figures/figure-S4.tiff", 
         height = 12, width = 17, units = 'cm', 
         compression = "lzw", res = fig.resolution)
    viz_ts_lag1(y=ts.data, data=sim.list, data.type=data.type, var.names=var.full.names)
    dev.off()
  }

  
  # Figure 6 - plot maximum running annual average
  
  ## for source water data, plot TOC
  if (data.type == "sw" && innov == TRUE) {
    ### get just toc data and make sure data is in proper order 
    toc.df <- arrange(sim.list[[5]], sim, year, month) %>%
      transform(value = as.numeric(value), 
                month = as.integer(month), 
                year = as.integer(year), 
                sim = as.integer(sim))
    
    ### determine maximum running annual average for each simulation
    toc.df2 <- toc.df %>%
      group_by(sim) %>%
      mutate(run_avg = rollmean(x = value, k=12, align = "right", fill = NA)) %>%
      summarize(max_run_avg = max(run_avg, na.rm = TRUE))
    
    mean.max_run_avg <- summarize(toc.df2, mean_max_run_avg = mean(max_run_avg, na.rm = TRUE))
    
    ggplot(data = toc.df2) + 
      geom_histogram(aes(max_run_avg))
    
    ### determine maximum running annual average for observed data
    toc.df3 <- ts.data[,"toc"] %>% 
      as.data.frame %>%
      mutate(run_avg = rollmean(x = x, k=12, align = "right", fill = NA)) %>%
      summarize(max_run_avg = max(run_avg, na.rm = TRUE))
    
    # set x axis limits
    if (min(toc.df2$max_run_avg) < toc.df3$max_run_avg) {
      x.min <- min(toc.df2$max_run_avg) %>% floor
    } else {
      x.min <- toc.df3$max_run_avg %>% floor
    }
    
    if (max(toc.df2$max_run_avg) > toc.df3$max_run_avg) {
      x.max <- max(toc.df2$max_run_avg) %>% ceiling
    } else {
      x.max <- toc.df3$max_run_avg %>% ceiling
    }
    
    p.raa <- ggplot() + 
      geom_density(data = toc.df2, aes(max_run_avg)) +
      geom_vline(data=mean.max_run_avg, aes(xintercept=mean_max_run_avg), size=1, color="black") + # mean of simulated data
      geom_vline(data = toc.df3, aes(xintercept=max_run_avg), size=1, color="#FF9999") +  # observed value
      xlim(x.min, x.max) +
      ylab("Density") +
      xlab("Maximum Running Annual Average TOC (mg/L)")
    
    tiff(filename = "./figures/figure-6.tiff", 
         height = 12, width = 17, units = 'cm', 
         compression = "lzw", res = fig.resolution)
    print(p.raa)
    dev.off()
    
    # calculate percentile of observed value compared to simulated
    ecdf(toc.df2$max_run_avg)(toc.df3$max_run_avg)
    
  } else if (data.type == "mine" && innov == TRUE) {
    
    ## plot maximum running annual average for temperature (mine data)
    
    ### get just temp data and make sure data is in proper order 
    temp.df <- arrange(sim.list[[1]], sim, year, month) %>%
      transform(value = as.numeric(value), 
                month = as.integer(month), 
                year = as.integer(year), 
                sim = as.integer(sim))
    
    ### determine maximum running annual average for each simulation
    temp.df2 <- temp.df %>%
      group_by(sim) %>%
      mutate(run_avg = rollmean(x = value, k=12, align = "right", fill = NA)) %>%
      summarize(max_run_avg = max(run_avg, na.rm = TRUE)) 
    
    mean.max_run_avg <- summarize(temp.df2, mean_max_run_avg = mean(max_run_avg, na.rm = TRUE))
    
    ggplot(data = temp.df2) + 
      geom_histogram(aes(max_run_avg))
    
    ### determine maximum running annual average for observed data
    temp.df3 <- ts.data[,"temp_C_plus_noise"] %>% 
      as.data.frame %>%
      mutate(run_avg = rollmean(x = x, k=12, align = "right", fill = NA)) %>%
      summarize(max_run_avg = max(run_avg, na.rm = TRUE))
    
    # set x axis limits
    if (min(temp.df2$max_run_avg) < temp.df3$max_run_avg) {
      x.min <- min(temp.df2$max_run_avg) %>% floor
    } else {
      x.min <- temp.df3$max_run_avg %>% floor
    }
    
    if (max(temp.df2$max_run_avg) > temp.df3$max_run_avg) {
      x.max <- max(temp.df2$max_run_avg) %>% ceiling
    } else {
      x.max <- temp.df3$max_run_avg %>% ceiling
    }
    
    p.raa <- ggplot() + 
      geom_density(data = temp.df2, aes(max_run_avg)) +
      geom_vline(data=mean.max_run_avg, aes(xintercept=mean_max_run_avg), size=1, color="black") + # mean of simulated data
      geom_vline(data = temp.df3, aes(xintercept=max_run_avg), size=1, color="#FF9999") +  # observed value
      xlim(x.min, x.max) +
      ylab("Density") +
      xlab("Maximum Running Annual Average (°C)")
    
    tiff(filename = "./figures/figure-S5.tiff", 
         height = 12, width = 17, units = 'cm', 
         compression = "lzw", res = fig.resolution)
    print(p.raa)
    dev.off()
    
    # calculate percentile of observed value compared to simulated
    # temp.df2$percentile <- ecdf(temp.df3$max_run_avg)
  }
}

# save function
save("visualize_statistics", file="./lib/visualize_statistics.RData")

# run script
visualize_statistics(innov=TRUE, data.type="sw")

