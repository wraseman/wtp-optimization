README for the observed and simulated data for the "mine" dataset
author: William Raseman

./observed

	Observed data can be found at the following path: ./observed/monthly-temp+precip_coppermine.csv

		Reference for observed data: 
		Hipel, Keith W., and A. Ian McLeod. Time Series Modelling of Water Resources and 
		Environmental Systems. Vol. 45. Elsevier, 1994.
		
		Links for observed data:
			1. Precipitation: https://datamarket.com/data/set/22n8/monthly-rain-coppermine-mm-1933-1976#!ds=22n8&display=line
			2. Temperature: https://datamarket.com/data/set/22t9/monthly-temperature-coppermine-celsius-1933-1976#!ds=22t9&display=line

./04_simulate_kNN
			
	Data imported or created by ../../04_simulate_kNN.R are contained within ./04_simulate_kNN. After running
	all scripts, this will include:

		mine_ts-data_with-noise.RData
			This .RData file contains the observed data converted to a time series class. Moreover, 
			some noise has been added to this dataset to avoid calculating the inverse of a singular
			matrix for the Mahalanobis distance calculation. 
			
		kNN_innov-*_mine_*.rds
			These files are dataframes of simulated data for various kNN runs. The "innov" aspect of the 
			file name refers to whether or not random innovations were included in the algorithm. 
			A distinct dataframe (and .rds file) is created for each simulated variable. 
