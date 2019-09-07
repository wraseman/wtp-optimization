# Post-process water quality ensembles for input to WTP Model optimization
# author: William Raseman 

## Tasks: aggregate data to quarterly frequency, add water quality parameters relevant 
##        to water treatment, and create tidy format .csv file with all water quality 
##        parameters included. 

# clear environment
rm(list=ls()) 

# load libraries
library(tidyverse)

# read in .rds files for each simulated water quality parameters:
#  pH, temperature, total organic carbon, turbidity, alkalinity, and hardness

## source water quality
var.names <- c("alk", "hard", "pH", "temp", "toc", "turb")

data.type <- "sw"
innov <- TRUE
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
  
  # create water quality data frame with all the water quality parameters
  if (i == 1) {
    wq_df <- read_rds(read.path) %>%
      mutate(parameter = var.names[i])
  } else {
    temp_df <- read_rds(read.path) %>%
      mutate(parameter = var.names[i])
    
    wq_df <- rbind(wq_df, temp_df)
  }
}

# aggregate to quarterly frequency
wq_df <- as.tibble(wq_df) %>% # make data frame a tibble so it is easier to view
  mutate(month = as.integer(month), year = as.integer(year)) %>% # make month and year integer type
  mutate(quarter = ceiling(month/3) %>% as.integer) %>%  # make quarter
  group_by(year, quarter, sim, parameter) %>% 
  summarize(value = mean(value, na.rm=TRUE))

# "untidy" data with tidyr::spread() into a wide format for easier manipulation
wq_wide_df <- spread(wq_df, key = parameter, value = value)

# water quality parameters to add based on constant concentration assumptions
#  or relationships with simulated parameters

wq_wide_df <- wq_wide_df %>%
  mutate(bromide = 0.015) %>%  # bromide: assume 0.015 mg/L for WTP Model calibration for Fort Collins
  mutate(ammonia = 0.01) %>%  # ammonia (nh3): assume 0.01 mg/L for WTP Model calibration for Fort Collins
  mutate(calcium = hard / 3.81) %>%  # calcium: hardness / 3.81 based on linear regression (R^2 = 0.93 after removal of one outlier)
  mutate(uv254 = toc * 0.03) ## uv254: assume SUVA of 3.0 based on 2001 data reported in Billica et al. (2008)
                             ##        "Design of a Collaborative Water Quality Monitoring Program for the Upper
                             ##         Cache la Poudre River"
## note - we chose not to use uv254 data in "Influent pH temperature turbidity & UVA.xlsx" provided
##        by Fort Collins SCADA because the data appeared to be unreliable (out of the range of 
##        expected values) 

# visualize data as a quality check
wq_wide_df <- ungroup(wq_wide_df) %>%
  arrange(sim, year, quarter) %>%  # arrange based on simulation number, year, and quarter
  select(sim, year, quarter, alk, ammonia, bromide, calcium, hard, pH, temp, toc, turb, uv254)

for (i in 4:ncol(wq_wide_df)) {
  p <- ggplot(wq_wide_df, aes_string(x = colnames(wq_wide_df)[i])) +
    geom_histogram(aes(color = as.factor(quarter), fill = as.factor(quarter)),
                         alpha = 0.4, position = "identity") 
  print(p)
}

# create .csv file with all water quality parameters included
path_csv <- "../in/monte_carlo/influent-wq-data.csv"
# header_comment <- paste0("#", paste(colnames(wq_wide_df), collapse=","), collapse = "")
# write_csv(x = header_comment, path = path_csv, append=FALSE)
write_csv(x = wq_wide_df, path = path_csv, col_names = TRUE)
