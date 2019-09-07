# Stochastic influent water quality generation
Code and data for implementing a *k* nearest neighbor (k-NN) bootstrap resampling approach for generating influent time series for water treatment.

## Dependencies
All dependencies are freely and openly available:

- [R (version 3.5.0)](https://cran.r-project.org/src/base/R-3/)
- [RStudio](https://www.rstudio.com/products/rstudio/download/#download)
- R packages: all R packages contained in the .R files must be installed before running the scripts.

## Running the code
- Download or clone this GitHub repository. If you've downloaded the repo, unzip the directory.
- Navigate to the repository, and open the .Rproj file.
- Open `run_all_scripts.R` in RStudio and click "Source".
- Wait for simulations to run: it may take several hours.

To reduce computation time, you can edit the number of simulations (default is 2500) by altering `nsims` before running `run_all_scripts.R`

## Individual scripts
There are five different scripts that make up the analysis in this repository:

- `01_import_clean.R`: import and clean observed water quality data
- `02_create_ts.R`: interpolate between missing data points and create complete time series dataset
- `03_visualize_ts.R`: plot complete, monthly time series
- `04_simulate_kNN.R`: generate synthetic influent water quality data using *k*-NN resampling algorithm
- `05_visualize_statistics.R`: visualize statistics of both observed and simulated datasets

Each script creates a function that is saved to `./lib` and is loaded be loaded by `run_all_scripts.R`. If any changes are made to the above scripts, they need to be run and reloaded by `run_all_scripts.R` to redo the analysis.

## Data
Two datasets are included in the analysis. The first is a water quality dataset of the Cache la Poudre River from the City of Fort Collins Utility. This dataset has been cleaned (as described in `01_import_clean.R`) and missing values have been interpolated (as described in `02_create_ts.R`). The second dataset is not a water quality dataset, rather it is a record of temperature and precipitation, but is used as a reference because it is a long multivariate dataset. 
