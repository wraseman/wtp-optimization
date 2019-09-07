# figures used to create the Table of Contents image

# clear environment
rm(list=ls()) 

# load packages
library(tidyverse)  # modern R packages: ggplot2, dplyr, readr, etc.
library(GGally)     # scatterplot matrices
library(forecast)   # time series analysis and visualization  
library(gridExtra)  # arrange multiple ggplot2 plots on a single plot
library(seasonal)   # time series decomposition 

# load user-defined libraries
source("./lib/time-series-sim_lib.R")  # time series simulation library

# define path for data
sw.id <- 3
wq.ts <- read_sw_ts(sw.id)

# observed data (suva)
ggseasonplot(wq.ts[,"suva"], col="red", main="") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank()
        )

# "simulated" data #1 (mirror this)
ggseasonplot(wq.ts[,"alk"], col="darkgreen", main="") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank()
        )

# "simulated" data #2 (squeeze this)
ggseasonplot(wq.ts[,"suva"], col="blue", main="") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank()
  )