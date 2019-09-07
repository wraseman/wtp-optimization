# Time series simulation library
# purpose: functions for the simulation and visualization of time series data
# author: William Raseman

# load packages
library(tidyverse)  # modern R packages: ggplot2, dplyr, readr, etc.
library(lubridate)  # easy date handling

# make sure dplyr::select is used instead of MASS::select
select <- dplyr::select

read_sw_ts <- function(sw.id, param = 'all') {
  # purpose: read in source water quality time series from .rds file
  # input:
  #   - sw.id: source water id - 1, 2, 3, 4, or 5
  #   - param: water quality parameter of interest: 'alk', 'pH', 'toc', 'temp', 'uv254', 'suva' or 'all'
  # output: time series data file

  path <- "./data/source-water/03_complete-monthly_timeseries/"  # define relative path
  file.ts <- str_c(path, "sw", sw.id, "-ts.rds")  # create relative file path with file included

  if (param=='all') {
  y <- readr::read_rds(file.ts)
  } else {
    y <- readr::read_rds(file.ts)[,param]
  }

  return(y)
}

add_monthly_stats_cols <- function(y, data, stats, distr='norm') {
  # purpose: add columns of monthly statistics to time series (monthly) data frame
  # input: 
  #   - data: dataframe of monthly data [y | year | month]
  #   - stats: datafrme of monthly statistics [month | distr. stats #1 | distr. stats #2 ...]
  #   - distr: type of statistical distribution (e.g. 'norm', 'lnorm', 'beta', 'gamma', 'pois', etc.)
  # output: datafrome of monthly with monthly statistic columns added
  
  freq <- 12  # monthly time series frequency
  
  if (distr=='norm') {
    for (i in 1:freq) {  # statistics
      tmp.df <- filter(data, month == i) %>%
        mutate(Mean = stats[i,"Mean"] %>% unlist) %>% 
        mutate(Sd = stats[i,"Sd"] %>% unlist)
      if (i == 1) {
        new.df <- tmp.df 
      } else {
        new.df <- rbind(new.df, tmp.df)
      }
    }
  } else {
    stop("Statistical distribution entered into add_monthly_stats_cols() currently not supported. Try other distribution like 'norm'.")
  }
  
  data <- group_by(new.df, year, month) %>% 
    as.tibble %>% 
    arrange(year, month)  # arrange data in chronological order
  
}

yearqrt_vect <- function(y) {
  # purpose: create a vector of years and quarters based on time series data
  # input:
  #   y - time series (class:ts, multivariate or univariate)
  # output: 
  # year.qrt - vector of years and quarters
  
  # select only the first column (in case a multivariate time series is input)
  if (is.null(ncol(y)) == FALSE) {
    y <- y[,1]
  }
  
  # initilize necessary dimensions and vectors
  freq = frequency(y)
  sim.length = length(y)
  
  # create vector of years and quarter for period of record
  year.qrt <- matrix(nrow=length(y), ncol=2)
  colnames(year.qrt) <- c("year", "quarter")

  count <- 0
  for (i in start(y)[1]:end(y)[1]) {
    for (j in start(y)[2]:end(y)[2]) {
      count <- count + 1
      year.qrt[count,] <- c(i,j)
    }
  }
  
  return(year.qrt)
}


yearmon_vect <- function(y) {
  # purpose: create a vector of years and months based on time series data
  # input:
  #   y - time series (class:ts, multivariate or univariate)
  # output: 
  # year.mon - vector of years and months
  
  # select only the first column (in case a multivariate time series is input)
  if (is.null(ncol(y)) == FALSE) {
    y <- y[,1]
  }
  
  # initilize necessary dimensions and vectors
  freq = frequency(y)
  sim.length = length(y)
  
  # create vector of years and month for period of record
  year.mon <- matrix(nrow=length(y), ncol=2)
  colnames(year.mon) <- c("year", "month")
  
  # create sequence of from start to end of time series data
  seq.date <- seq(as.Date(str_c(start(y)[1],start(y)[2],"1", sep="/")), 
                  as.Date(str_c(end(y)[1],end(y)[2],"1", sep="/")), 
                  by="month")
  
  for (i in 1:length(seq.date)) {
    year.mon[i,] <- c(year(seq.date[i]), month(seq.date[i]))
  }
  
  return(year.mon)
}

create_seasonalplot_list <- function(data, FUN) {
  
  # create list of seasonal plots
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

ts_sim <- function(y, model, n, future = TRUE, bootstrap = FALSE, xreg = NULL) {
  
  # purpose: visualize simulated and observed statistics for univariate time series
  # inputs:
  #   y - univariate time series (class: ts)
  #   model  - time series model (e.g. "Arima")
  #   n      - number of simulations
  #   future - produce sample paths that are future to and conditional on the data in 'model'
  #   bootstrap - do simulation using resampled errors rather than normally distributed errors
  #   xreg   - externel regressors for model (default: NULL)
  
  # initilize necessary dimensions and vectors
  n = n.sim
  
  # create vector of years and month for period of record
  year.mon <- yearmon_vect(y)
  
  # initilize simulation data frame
  sim.df <- data.frame(value=numeric(0),
                       month=integer(0),
                       year=integer(0),
                       sim=integer(0))
  
  sim.data <- matrix(data=NA, nrow=length(y), ncol=n.sim)
  
  for (j in 1:n.sim) {

    if (any(class(model)[1:length(class(model))]=="tslm")) {
      sim.data[,j] <- simulate(bestfit, future=future, bootstrap=bootstrap) %>% unlist 
    } else {
      if (is.null(model$xreg)) {
        
        sim.data[,j] <- simulate(bestfit, future=future, bootstrap=bootstrap) 
      } else {
        sim.data[,j] <- simulate(bestfit, future=future, bootstrap=bootstrap, xreg=xreg)
      }
    }
    temp.df <- data.frame(value=sim.data[,j],
                          month=year.mon[,"month"],
                          year=year.mon[,"year"]-min(year.mon[,"year"])+1,
                          sim=j)
    sim.df <- rbind(sim.df, temp.df)
  }
  
  return(sim.df)
}

combine_obs_sim <- function(y, data) {
  # purpose: create one tidy dataframe of observed and simulated data
  # input: 
  #   y - univariate time series of observed data (class: ts)
  #   data - simulated data from model of obs. time series with the following column names:
  #       value | month | year (starting at 1) | simulation [if frequency(y) == 12]
  #       value | quarter | year (starting at 1) | simulation [if frequency(y) == 4]
  #       note: this data frame is also the output of ts_sim() also defined in tssim_lib.R
  # output: combined dataframe
  
  # remove any groups
  data <- as.tibble(data) %>% ungroup
  
  if (frequency(y) == 12) {
    # create vector of years and month for period of record
    year.mon <- yearmon_vect(y)
    
    obs.df <- y %>% data.frame %>%
      mutate(month=year.mon[,"month"], 
             obs_sim = rep("obs", length(y)))
    
  } else if (frequency(y) == 4) {
    # create vector of years and quarters for period of record
    year.qrt <- yearqrt_vect(y)
    
    obs.df <- y %>% data.frame %>%
      mutate(quarter=year.qrt[,"quarter"], 
             obs_sim = rep("obs", length(y)))
  } else {
    stop('Only supported timesteps are month and quarter.')
  }

  colnames(obs.df)[1] <- "value"
  
  data <- mutate(data, obs_sim= rep("sim", nrow(data)))
  
  full.df <- select(data, -year, -sim) %>% rbind(obs.df)
  
  # Convert the variables obs_sim and timestep (month or quarter) from a numeric to a factor variable
  full.df$obs_sim <- as.factor(full.df$obs_sim)
  
  if (frequency(y) == 12) full.df$month <- as.factor(full.df$month)
  if (frequency(y) == 4) full.df$quarter <- as.factor(full.df$quarter)
  
  return(full.df)
}

viz_obs_sim <- function(y, data, title="") {
  
  # purpose: create plot of simulated and observed statistics for univariate time series
  # inputs:
  #   y - univariate time series of observed data (class: ts)
  #   data - simulated data from model of obs. time series with the following column names:
  #       value | month | year (starting at 1) | simulation [if frequency(y) == 12]
  #   title - resulting plot title
  # return: plot
  # reference: http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization
  
  # remove any groups
  data <- as.tibble(data) %>% ungroup
  
  if (frequency(y) == 12) {
    
    full.df <- combine_obs_sim(y, data)
    
    p <- ggplot(aes(x=month, y=value, fill=obs_sim), data=full.df) +
      geom_boxplot(outlier.size=0.5) +
      ggtitle(label=title) +
      xlab('month') +
      scale_fill_manual(values=c("#FF9999", "#FFFFFF"))
    
  } else {
    
    stop('Only supported timesteps are month.')
    
  }
  
  return(p)
}

viz_ts_sample_stats <- function(y, data, title="") {
  
  # purpose: create plot of simulated and observed statistics (mean, standard deviation,
  #           minimum, maximum) for univariate time series
  # inputs:
  #   y - univariate time series of observed data (class: ts)
  #   data - simulated data from model of obs. time series with the following column names:
  #       value | month | year (starting at 1) | simulation [if frequency(y) == 12]
  #   title - resulting plot title
  # return: nothing
  
  # remove any groups
  data <- as.tibble(data) %>% ungroup 
  
  freq <- frequency(y)
  
  # arrange in proper order
  if (freq == 12) {
    
    data <- arrange(data, sim, year, month) %>%
      transform(value = as.numeric(value), 
                month = as.integer(month), 
                year = as.integer(year), 
                sim = as.integer(sim))
    
    # create vector of years and month for period of record
    year.mon <- yearmon_vect(y)
    
    obs.df <- y %>% as.numeric %>% data.frame %>%
      mutate(month=year.mon[,"month"],
             year=year.mon[,"year"])
    colnames(obs.df)[1] <- "value" 
    
  } else {
    
    stop('Only supported frequency is monthly.')
    
  }
  
  # plot simulation and observation statistics
  
  stats.sim <- summarize(group_by(data, month, sim),
                         mean=mean(value), median=median(value),
                         sd=sd(value), min=min(value), max=max(value))
  
  stats.sim$month <- as.factor(stats.sim$month)
  
  # calculate statistics for observed data
  stats.obs <- summarize(group_by(obs.df , month),
                         mean=mean(value), median=median(value),
                         sd=sd(value), min=min(value),
                         max=max(value))
  
  stats.obs$month <- as.factor(stats.obs$month)
  x.label <- "Month"
  
  outlier.size <- 0.5
  
  # mean
  p.mean <- ggplot() +
    geom_boxplot(mapping=aes(x=month %>% as.factor, y=mean), data=stats.sim,
                 outlier.size = outlier.size) +
    geom_point(mapping=aes(x=month, y=mean), data=stats.obs, col="#FF9999") +
    xlab(x.label) +
    ylab("Mean")
  
  # sd
  p.sd <- ggplot() +
    geom_boxplot(mapping=aes(x=month, y=sd), data=stats.sim,
                 outlier.size = outlier.size) + 
    geom_point(mapping=aes(x=month, y=sd), data=stats.obs, col="#FF9999") +
    xlab(x.label) +
    ylab("Standard Deviation")
  
  # min
  p.min <- ggplot() +
    geom_boxplot(mapping=aes(x=month, y=min),  data=stats.sim,
                 outlier.size = outlier.size) +
    geom_point(mapping=aes(x=month, y=min), data=stats.obs, col='#FF9999') +
    xlab(x.label) +
    ylab("Minimum")
  
  # max
  p.max <- ggplot() +
    geom_boxplot(mapping=aes(x=month, y=max), data=stats.sim,
                 outlier.size = outlier.size) +
    geom_point(mapping=aes(x=month, y=max), data=stats.obs, col='#FF9999') +
    xlab(x.label) +
    ylab("Maximum")
  
  grid.arrange(p.mean, p.sd, p.min, p.max, top = title)
}

viz_pair_corr <- function(y=y, data, title="", nsims, data.type) {
  # purpose: create plot of simulated and observed pairwise correlation for multivariate time series
  # inputs:
  #   y - multivariate time series of observed data (class: ts)
  #   data - simulated data from model of obs. time series with the following column names:
  #       value | month | year (starting at 1) | simulation [if frequency(y) == 12]
  #   title - resulting plot title
  #   data.type - either "sw" (source water) or "mine" (mine data)
  # return: nothing
  
  sim.lag1 <- list()
  for (i in 1:length(data)) {
    # remove any groups and tidy up formatting
    sim.lag1[[i]] <- as.tibble(data[[i]]) %>% ungroup
  }
  
  nvars <- ncol(y)
  nsims <- max(data[[1]]$sim)
  
  var.names <- colnames(y)
  
  pair.comb <- combn(x=1:nvars, m=2)  # matrix of pairwise combinations
  pair.cor.sim <- matrix(nrow=nsims, ncol=ncol(pair.comb)) %>%
    as.data.frame  # pairwise correlation between two variables
  pair.cor.obs <- matrix(nrow=1, ncol=ncol(pair.comb)) %>% as.data.frame
  pair.names <- vector()
  
  for (i in 1:ncol(pair.comb)) {
    
    var1 <- pair.comb[1,i]
    var2 <- pair.comb[2,i] 
    var1.name <- var.names[var1]
    var2.name <- var.names[var2]
    
    pair.names[i] <- str_c(var1.name, "_" , var2.name, sep="")
    
    # create a dataframe in which observed is appended to the last column and then there is 
    
    for (k in 1:nsims) {
      # get simulated values of first variable
      sim.var1 <- filter(sim.lag1[[var1]], sim==k) %>% 
        select(value) %>% 
        unlist
      # get simulated values of secondvariable
      sim.var2 <- filter(sim.lag1[[var2]], sim==k) %>% 
        select(value) %>% 
        unlist
      # calculate 
      pair.cor.sim[k,i] <- cor(sim.var1, sim.var2) 
      
      # filter(sim.lag1[[var1]], sim==k)
    }
    
    obs.var1 <- y[,var1]
    obs.var2 <- y[,var2]
    pair.cor.obs[i] <- cor(obs.var1, obs.var2)
    
  }  # END pairwise combination loop
  
  colnames(pair.cor.sim) <- colnames(pair.cor.obs) <- pair.names
  
  # "tidy" data 
  pair.cor.sim <- gather(pair.cor.sim)
  pair.cor.obs <- gather(pair.cor.obs)
  
  
  # pairwise correlation
  p.paircor <- ggplot() +
    geom_boxplot(mapping=aes(x=key, y=value), data=pair.cor.sim) +
    geom_point(mapping=aes(x=key, y=value), data=pair.cor.obs, col="#FF9999") +
    ylab("Correlation") + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle=-75, hjust=0, size=12),
          axis.ticks.x=element_blank())
  
  if (data.type == "sw") {
    print(p.paircor + scale_x_discrete(labels=c("Alkalinity, pH", "Alkalinity, Temperature",
                                          "Alkalinity, TOC", "pH, Temperature", 
                                          "pH, TOC", "Temperature, TOC")) 
    )
  } else if (data.type == "mine") {
    print(p.paircor + scale_x_discrete(labels=c("Temperature, Precipitation")))
  }

}

viz_ts_lag1 <- function(y, data, data.type, title="", var.names) {

  # purpose: create plot of simulated and observed lag-1 for multivariate time series
  #   y - multivariate time series of observed data (class: ts)
  #   data - simulated data from model of obs. time series with the following column names:
  #       value | month | year (starting at 1) | simulation [if frequency(y) == 12]
  #       value | quarter | year (starting at 1) | simulation [if frequency(y) == 4]
  #       note: this data frame is also the output of ts_sim() also defined in tssim_lib.R
  #   data.type - either "sw" (source water) or "mine" (mine data)
  #   title - resulting plot title
  #   var.names - names of each variable
  # return: nothing

  p.lag1 <- list()
  sim.lag1 <- list()
  for (i in 1:length(data)) {
    # remove any groups and tidy up formatting
    sim.lag1[[i]] <- as.tibble(data[[i]]) %>% ungroup
  }
  
  for (i in 1:length(sim.lag1)) {
    sim.lag1[[i]] <- arrange(sim.lag1[[i]], sim, year, month) %>%
      transform(value = as.numeric(value),
                month = as.integer(month),
                year = as.integer(year),
                sim = as.integer(sim))

    # create vector of years and month for period of record
    year.mon <- yearmon_vect(y[,i])

    obs.df <- y[,i] %>% as.numeric %>% data.frame %>%
      mutate(month=year.mon[,"month"],
             year=year.mon[,"year"])
    colnames(obs.df)[1] <- "value"

    # add lag variable
    sim.lag1[[i]] <- mutate(sim.lag1[[i]], lag1 = lag(value, n=1))  # add lag variable
    obs.df <- mutate(obs.df, lag1 = lag(value, n=1))  # add lag variable

    stats.sim <- summarize(group_by(sim.lag1[[i]], month, sim),
                           lag1_cor=cor(value,lag1,use="na.or.complete"))

    stats.sim$month <- as.factor(stats.sim$month)

    # calculate statistics for observed data
    stats.obs <- summarize(group_by(obs.df, month), lag1_cor=cor(value,lag1,use="na.or.complete"))

    stats.obs$month <- as.factor(stats.obs$month)

    # lag-1 correlation
    p.lag1[[i]] <- ggplot() +
      geom_boxplot(mapping=aes(x=month %>% as.factor, y=lag1_cor),
                   outlier.size=0.5, data=stats.sim) +
      geom_point(mapping=aes(x=month, y=lag1_cor), data=stats.obs, col="#FF9999") +
      xlab("Month") + ylab("Lag-1 Correlation") +
      ggtitle(var.names[i])
  }

  if (data.type == "sw") {
    grid.arrange(p.lag1[[1]], p.lag1[[2]], p.lag1[[3]], p.lag1[[4]])
  } else if (data.type == "mine") {
    grid.arrange(p.lag1[[1]], p.lag1[[2]])
  }
}


# viz_obs_sim_stats <- function(y, data, title="") {
#   
#   # purpose: create plot of simulated and observed statistics for univariate time series
#   # inputs:
#   #   y - univariate time series of observed data (class: ts)
#   #   data - simulated data from model of obs. time series with the following column names:
#   #       value | month | year (starting at 1) | simulation [if frequency(y) == 12]
#   #       value | quarter | year (starting at 1) | simulation [if frequency(y) == 4]
#   #       note: this data frame is also the output of ts_sim() also defined in tssim_lib.R
#   #   title - resulting plot title
#   # return: nothing
#   
#   # remove any groups
#   data <- as.tibble(data) %>% ungroup 
#   
#   freq <- frequency(y)
#   
#   # arrange in proper order
#   if (freq == 12) {
#     
#     data <- arrange(data, sim, year, month) %>%
#       transform(value = as.numeric(value), 
#                 month = as.integer(month), 
#                 year = as.integer(year), 
#                 sim = as.integer(sim))
#     
#     # create vector of years and month for period of record
#     year.mon <- yearmon_vect(y)
#     
#     obs.df <- y %>% as.numeric %>% data.frame %>%
#       mutate(month=year.mon[,"month"],
#              year=year.mon[,"year"])
#     colnames(obs.df)[1] <- "value" 
#     
#   } else if (freq == 4) {
#     
#     data <- arrange(data, sim, year, quarter) %>%
#       transform(value = as.numeric(value), 
#                 quarter = as.integer(quarter), 
#                 year = as.integer(year), 
#                 sim = as.integer(sim))
#     
#     # create vector of years and quarter for period of record
#     year.qrt <- yearqrt_vect(y)
#     
#     obs.df <- y %>% as.numeric %>% data.frame %>%
#       mutate(quarter=year.qrt[,"quarter"],
#              year=year.qrt[",year"])
#     colnames(obs.df)[1] <- "value" 
#     
#   } else {
#     
#     stop('Only supported timesteps are month and quarter.')
#     
#   }
# 
#   # add columns for important lags 
#   data <- data %>% 
#     mutate(lag1 = lag(value, n=1), 
#            lag2 = lag(value, n=2), 
#            lagfreq1 = lag(value, n=freq*1),
#            lagfreq2 = lag(value, n=freq*2))
# 
#   obs.df <- obs.df %>%  mutate(lag1 = lag(value, n=1), 
#                                lag2 = lag(value, n=2), 
#                                lagfreq1 = lag(value, n=freq*1),
#                                lagfreq2 = lag(value, n=freq*2)) 
#   
#   # plot simulation and observation statistics
# 
#   stats.sim <- summarize(group_by(data, month, sim),
#                          mean=mean(value), median=median(value),
#                          sd=sd(value), min=min(value), max=max(value),
#                          lag1_cor=cor(value,lag1,use="na.or.complete"), 
#                          lag2_cor=cor(value,lag2,use="na.or.complete"), 
#                          lagfreq1_cor=cor(value,lagfreq1,use="na.or.complete"), 
#                          lagfreq2_cor=cor(value,lagfreq2,use="na.or.complete"))
#   
#   stats.sim$month <- as.factor(stats.sim$month)
#   
#   # calculate statistics for observed data
#   stats.obs <- summarize(group_by(obs.df , month),
#                          mean=mean(value), median=median(value),
#                          sd=sd(value), min=min(value),
#                          max=max(value),  
#                          lag1_cor=cor(value,lag1,use="na.or.complete"), 
#                          lag2_cor=cor(value,lag2,use="na.or.complete"), 
#                          lagfreq1_cor=cor(value,lagfreq1,use="na.or.complete"), 
#                          lagfreq2_cor=cor(value,lagfreq2,use="na.or.complete"))
#   
#   stats.obs$month <- as.factor(stats.obs$month)
#   x.label <- "month"
#   
#   # mean
#   p.mean <- ggplot() +
#     geom_boxplot(mapping=aes(x=month %>% as.factor, y=mean), data=stats.sim) +
#     geom_point(mapping=aes(x=month, y=mean), data=stats.obs, col="#FF9999") +
#     xlab(x.label)
#   
#   # sd
#   p.sd <- ggplot() +
#     geom_boxplot(mapping=aes(x=month, y=sd), data=stats.sim) +
#     geom_point(mapping=aes(x=month, y=sd), data=stats.obs, col="#FF9999") +
#     xlab(x.label)
#   
#   # min
#   p.min <- ggplot() +
#     geom_boxplot(mapping=aes(x=month, y=min), data=stats.sim) +
#     geom_point(mapping=aes(x=month, y=min), data=stats.obs, col='#FF9999') +
#     xlab(x.label)
#   
#   # max
#   p.max <- ggplot() +
#     geom_boxplot(mapping=aes(x=month, y=max), data=stats.sim) +
#     geom_point(mapping=aes(x=month, y=max), data=stats.obs, col='#FF9999') +
#     xlab(x.label)
#   
#   grid.arrange(p.mean, p.sd, p.min, p.max, top = title)
#   
#   # lag-1 correlation
#   p.lag1 <- ggplot() +
#     geom_boxplot(mapping=aes(x=month %>% as.factor, y=lag1_cor), data=stats.sim) +
#     geom_point(mapping=aes(x=month, y=lag1_cor), data=stats.obs, col="#FF9999") +
#     xlab(x.label) + ylab("lag-1 correlation")
#   
#   # lag-2 correlation
#   p.lag2 <- ggplot() +
#     geom_boxplot(mapping=aes(x=month, y=lag2_cor), data=stats.sim) +
#     geom_point(mapping=aes(x=month, y=lag2_cor), data=stats.obs, col="#FF9999") +
#     xlab(x.label) + ylab("lag-2 correlation")
#   
#   # lag-(frequency of data) correlation [e.g. for monthly data, lag-12]
#   p.lagfreq1 <- ggplot() +
#     geom_boxplot(mapping=aes(x=month, y=lagfreq1_cor), data=stats.sim) +
#     geom_point(mapping=aes(x=month, y=lagfreq1_cor), data=stats.obs, col='#FF9999') +
#     xlab(x.label) + ylab(str_c("lag-", freq, " correlation"))
#   
#   # lag-(2*frequency of data) correlation [e.g. for monthly data, lag-24]
#   p.lagfreq2 <- ggplot() +
#     geom_boxplot(mapping=aes(x=month, y=lagfreq2_cor), data=stats.sim) +
#     geom_point(mapping=aes(x=month, y=lagfreq2_cor), data=stats.obs, col='#FF9999') +
#     xlab(x.label) + ylab(str_c("lag-", freq*2, " correlation"))
#   
#   grid.arrange(p.lag1, p.lag2, p.lagfreq1, p.lagfreq2, top = title)
# }

ggacf_matrix <- function(data) {
  
  # purpose: plot ggstyle autocorrelation and cross correlation matrix
  # input: multivariate timeseries (class: ts)
  # output: none (but prints plot)
  
  ccf_list <- list()
  count <- 1
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      local({
        if (i == j) {
          p <- ggAcf(x = data[,i]) +
            ggtitle(colnames(data)[i]) +
            theme(text = element_text(size=10))
        } else if (j > i) {
          p <- ggplot() + theme_classic()
        }
        else {
          p <- ggCcf(x = data[,i], y = data[,j]) +
            ggtitle(str_c(colnames(data)[i]," & ",colnames(data)[j])) +
            theme(text = element_text(size=10))
        }
        ccf_list[[count]] <<- p
      })
      count <- count + 1
    }
  }
  
  n.wqvar <- ncol(data)
  grid.arrange(grobs = ccf_list, nrow=n.wqvar, ncol=n.wqvar) %>% print
}