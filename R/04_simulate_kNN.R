# Step 4 - Fit a multivariate K-NN lag-1 model, as described in "Nearest Neighbor Resampling Algorithm" 
#   section of "Nearest neighbor bootstrap for generating influent time series for water treatment"
# author: Billy Raseman

# clear environment
rm(list=ls())

simulate_kNN <- function(nsims=100, innov=TRUE, threshold, 
                         thresh.alpha=0.06, data.type,
                         standardize=TRUE) {

  # load packages
  library(abind)
  library(tidyverse)
  
  # load user-defined functions
  source("./lib/time-series-sim_lib.R")  # time series simulation library
  
  # read in time series data
  if (data.type == "sw") {
    
    ## source water quality
    path <- "./data/source-water/02_create_ts/sw_ts.rds"
    ts.data <- readr::read_rds(path)  
  
  } else if (data.type == "mine"){
    
    ## precipitation and temperature from mine
    mine.df <- read_csv(file = "./data/mine/observed/monthly-temp+precip_coppermine.csv")
    
    ## add noise to avoid singular matrix (due to duplicate observations).
    ## alternatively, could just calculate the pseudoinverse of the covariance
    ## matrix for the Mahalanobis calculation.
    noise.temp <- runif(n=nrow(mine.df), min=1e-10, max=1e-8)  # if having issues, can increase minimum and maximum
    noise.precip <- runif(n=nrow(mine.df),  min=1e-10, max=1e-8)
    mine.df <- mutate(mine.df, temp_C_plus_noise=temp_C+noise.temp,
                       precip_mm_plus_noise=precip_mm+noise.precip)
    
    ## turn dataframe into time series dataset
    ts.data <- select(mine.df, -month_year, -temp_C, -precip_mm) %>%
      ts(frequency=12, start=c(1933, 1))
    write_rds(ts.data, path="./data/mine/04_simulate_kNN/mine_ts-data_with-noise.RData")
    
  }
  var.names <- colnames(ts.data)
  
  ## specifications for simulation
  nvars <- ncol(ts.data)  # number of variables to simulate
  nmonths <- 12  # number of months in a year
  nyrs <- dim(ts.data)[1]/nmonths  # number of years on record
  
  ## lambda threshold calculation (for Step 5 - Add random innovations...)
  z.score <- qnorm(thresh.alpha) 
  
  ## create 3-dimensional matrix of data 
  ### where 1st dim: variable
  ###       2nd dim: year
  ###       3rd dim: month
  ### note: this differs than the matrix discussed in the paper. This is because it is easier to 
  ###   implement the algorithm in 3-dimensions and easier to discuss conceptually as a 2-dimensional matrix.
  X <- stand.X <- array(0,dim=c(nvars,nyrs,nmonths))

  monthly.mean <- matrix(nrow=nvars, ncol=nmonths)
  monthly.sd <- matrix(nrow=nvars, ncol=nmonths)
  
  ## convert time series to 3-dimensional matrix format
  for (i in 1:nvars) {
    for (j in 1:nyrs) {
      for (k in 1:nmonths) {
        X[i,j,k] <- ts.data[,i][k+(j-1)*nmonths]
      }
      
      if (standardize == TRUE) {
        # calculate monthly statistics
        monthly.mean[i,] <- apply(X[i,,], 2, mean)
        monthly.sd[i,]   <- apply(X[i,,], 2, sd) 
        
        # standardize monthly
        stand.X[i,,] <- t( (t(X[i,,]) - monthly.mean[i,]) / monthly.sd[i,])  # standardize data based on monthly statistics. can also use built in function scale()
      }
    }
  }
  
  if (standardize == TRUE) {
    X <- stand.X  # set standardized values as X
  }
  
  # choose 'k' nearest-neighbors based on number of years on record
  K = sqrt(nyrs) %>% round  # use heuristic for choosing 'k' discussed in Lall and Sharma (1996)
  ## chose captial 'K' as variable name to not interfere with iteration variable 'k'
  nyrs1=nyrs-1  
  year = 1:nyrs  # vector of years
  
  ## create matrix for current simulation
  x.sim.mat <- array(0,dim=c(nvars,nyrs,nmonths))  
  ### 1st dimension is the variable
  ### 2nd dimension is rows - i.e., the years
  ### 3rd is the columns - i.e, months
  
  ## define kernel: the weighting metric to do the K-NN resampling (step 3 in Lall and Sharma, 1996)
  W=1:K
  W=1/W
  W=W/sum(W)  
  W=cumsum(W)
  
  ## initilize variables for simulation
  sim.df <- as.tibble(matrix(nrow=nsims*nvars*nyrs*nmonths, ncol=5)) 
  colnames(sim.df) <- c("value", "month", "year", "sim", "var")
  
  # initialize variables for random innovations
  if (innov == TRUE) {
    z.next <- vector()  # random variates (length = number of variables) drawn from normal distribution of mean zero and sd of one
      ## estimates of nonparametric distribution fit to k nearest neighbors
    sigma.cond  <- vector()  # conditional standard deviations (length = number of variables) 
    lambda <- vector()  # bandwidth (function of number of samples)
    lambda.prime <- vector()  # acceptable value of lambda (10b in Sharif and Burn, 2007)
    x.tilde <- vector()  # simulated value from basic kNN (i.e., without random innovations)
    x.tilde.prime <- vector()  # simulate value after random innovations added
    rand.innov <- vector()  # vector of random innovations
  }
  
  for(i in 1:nsims){
  
    x.sim = array(0,dim=c((nyrs*nmonths),nvars))  # simulated values from current simulation, i
  
    for(j in 1:(nyrs*nmonths)){
      
      if (j == 1) {
        
        ## for first timestep, use a randomly sampled year on record
        sample.num <- runif(n=1,min=1,max=nyrs) %>% round  # runif() randomly sample from uniform distribution
        x.sim[1,] <- as.matrix(X[,sample.num,1])  # get January data for randomly selected year 
        
      } else {
        
        ## get the month of simulation
        imon = j %% 12  
        if(imon == 0) {
          imon = 12
        }
  
        # Step 1 - Define a feature vector
        ## note: here, we choose a lag-1 dependence structure for the four water quality variables
        D.i <- x.sim[j-1,]  # current feature vector
   
        
        if(imon == 1) {
          imon1 = 12  # if January is the current month (imon) the lag-1 month (imon1) is December
          D.t = rbind(t(X[,1:nyrs1,imon1]))  # neighbors to current feature vector, D_t
          N = nyrs1
        } else {
          imon1 = imon-1  # any month other than January, the lag-1 month (imon1) is just 1 less than the current month (imon)
          D.t = rbind(t(X[,,imon1]))  # neighbors to current feature vector, D_t
          N = nyrs
        }
        
        # Step 2 - Find nearest neighbors
        
        distance <- mahalanobis(x = D.t, center = D.i, cov = cov(rbind(t(D.i), D.t)))  # distance from feature vector and each neighbor
        
        # Step 3 - Rank nearest neighbors and select k neighbors
        
        ## rank nearest neighbors
        ordered.distance = order(distance)  # nearest neighbors from feature vector by year on record
        
        ## define k nearest neighbors
        kNN.index <- ordered.distance[1:K]  # k nearest neighbors (by year on record)
        D.t.kNN <- t(X[,kNN.index,imon])
        
        # Step 4 - Choose successor
        
        ## define discrete kernel for resampling
        rand.samp=runif(1,0,1)  # sample from kernel
        xy=c(rand.samp,W)
        chosen.neighbor=rank(xy)[1]  # neighbor simulation from k nearest neighbors (1 being the nearest neighbor, 2 the second, and so on)
        chosen.year=ordered.distance[chosen.neighbor]  # year on record that chosen neighbor corresponds to
        if(imon == 1) chosen.year=chosen.year+1  # if the simulated month is January, the first year cannot be sampled from
        
        x.tilde = t(X[,chosen.year,imon])  # chosen successor (simulated value without random innovations)
  
        # Step 5 - Add random innovations (i.e., errors) to successor
        
        if (innov == TRUE) {
          
          non.negative.check <- rep(-Inf, times=length(x.tilde))  # check that simulations are non-negative
  
          # 5c. repeat steps until all variables are non-negative
          while (any(non.negative.check < 0)) {
            
            # 5a. generate random variate, calc conditional standard deviation, and bandwidth
            for (k in 1:nvars) {
              current.var <- D.t.kNN[,k]
              z.next[k] <- rnorm(n=1, mean=0, sd=1)  # z_(k,t+1): where k is the variable and t is the current timestep
              sigma.cond[k] <- sd(current.var)  # sigma_k : conditional standard deviation, where k is the variable  
              sigma.kern <- 1  # standard deviation of kernel (in this case, 1)
              
              if (sigma.cond[k]==0) {
                
                stop("Standard deviation of k-Nearest Neighbors is zero. Cannot calculate bandwidth. Likely there are duplicate observations.")
              }
              
              ## calculate bandwidth
              lambda[k] <- bw.nrd0(current.var)  # rule-of-thumb estimation of a Gaussian kernel density estimator from Silverman (1986)
              
              if (threshold[k] == TRUE) {
                
                lambda.prime[k] <- x.tilde[k]/(z.score*sigma.cond[k])
                
                # 5b. account for bounded (non-negative) variables
                if (lambda[k] > lambda.prime[k]) {
                  chosen.lambda <- lambda.prime[k]
                } else {
                  chosen.lambda <- lambda[k]
                }
                
              } else {
                chosen.lambda <- lambda[k]
              }
              
              ## modify successor values with a smoothed boostrap with variance correction (Silverman, 1986)
              rand.innov[k] <- chosen.lambda*sigma.cond[k]*z.next[k]/sqrt(1+(chosen.lambda^2*sigma.kern)/sigma.cond[k]^2)
              
              x.tilde.prime[k] <- x.tilde[k] + rand.innov[k]  # choose successor as x-tilde' for each variable
    
              ## unstandaridize simulation (if applicable)
              if (standardize == TRUE) {
                non.negative.check[k] <- x.tilde.prime[k]*monthly.sd[k,imon] + monthly.mean[k,imon]  # transform back to unstandardized (original space) value to check if non-negative
              } else {
                non.negative.check[k] <- x.tilde.prime[k]
              }
              
              x.sim[j,k] <- x.tilde.prime[k]  # save simulated value
            }
          }
        } else {
          x.sim[j,] <- x.tilde  # save simulated value
        }
      }
  
      ## create dataframe to store simulations
      for (k in 1:nvars) {
        x.sim.mat[k,,] = matrix(x.sim[,k],ncol=12,byrow=T) 

        if (standardize == TRUE) {
          x.sim.mat[k,,] = t(t(x.sim.mat[k,,]) * monthly.sd[k,] + monthly.mean[k,])  # unstandardize data
        }

        temp.mat <- cbind(year, x.sim.mat[k,,]) %>% as.tibble
        colnames(temp.mat) <- c("year", "1", "2", "3", "4", "5", "6",
                                "7", "8", "9", "10", "11", "12")

        temp.df <- temp.mat %>%
          gather(month, value, '1':'12') %>%
          transform(month = as.numeric(month))

        var <- rep(k, nyrs*nmonths)

        if (k==1) {
          x.sim.df <- cbind(temp.df, var)
        } else {
          this.var.df <- cbind(temp.df, var)
          x.sim.df <- rbind(x.sim.df, this.var.df)
        }
      }
    }  # end 'for(j in 1:(nyrs*nmonths)){'
    
    sim <- rep(i, nyrs*nmonths)  # current simulation
    x.sim.df <- cbind(x.sim.df, sim)
    
    if (i == 1) {
      sim.df <- x.sim.df
    } else {
      sim.df <- rbind(sim.df, x.sim.df)
    }
  
  }  # end 'for(i in 1:nsims){'

  sim.df <- select(sim.df, value, month, year, sim, var)
  
  if (data.type == "sw") {
    data.dir <- "source-water"
  } else if (data.type == "mine") {
    data.dir <- "mine"
  }
  
  for (i in 1:nvars) {
    model.type <- str_c("kNN", str_c("innov", "-", innov), sep="_")
    write.path <- str_c("./data/source-water/04_simulate_kNN/", model.type, 
                        "_", 
                        data.type,
                        "_",
                        var.names[i],
                        ".rds")  # name simulation based on model, data type, and variable
    write_rds(filter(sim.df, var==i) %>% 
                select(value, month, year, sim), 
              write.path)
  }
}

# save function
save("simulate_kNN", file="./lib/simulate_kNN.RData")

# run script
simulate_kNN(nsims=5, innov=TRUE, threshold=c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE), data.type="sw", standardize=TRUE)  # uncomment to run script
# simulate_kNN(nsims=5, innov=TRUE, threshold=c(TRUE, TRUE), data.type="mine", standardize=TRUE)  # uncomment to run script

