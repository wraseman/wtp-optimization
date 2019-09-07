/* wtp_optimize.h */ 

#ifndef WTP_OPTIMIZE_H
#define WTP_OPTIMIZE_H 1

#include "wtp.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>

// Function declarations
// std::string read_json(const char *filename);
// std::string get_wtp_filename(std::string treatment_type);

void validate_optarg_runmode(std::string runmode);
void validate_optarg_int(char *optarg, char opt);
void validate_optarg_file(char *optarg, char opt);
void display_usage_help();
void save_cli_args(std::string filepath_cli_args, int argc, char** argv);
void copy_file(std::string source_path, std::string dest_path);

// Macros for operational parameters
#define ALKALINITY_SETPT 0
#define PH_SETPT 1
#define DBP_SAFETY_FACTOR 2

// Macros for influent water quality parameters
#define INF_ALKALINITY 0
#define INF_AMMONIA 1
#define INF_BROMIDE 2
#define INF_CALCIUM_HARDNESS 3
#define INF_TOTAL_HARDNESS 4
#define INF_PH 5
#define INF_TEMPERATURE 6
#define INF_TOTAL_ORGANIC_CARBON 7
#define INF_TURBIDITY 8
#define INF_UV254 9

/* Macros for optimization problem formulation */
#define N_YEARS 11                               // number of years on record
#define N_QUARTERS_PER_YEAR 4                    // number of quarters in a year
#define N_TIMESTEPS (N_YEARS*N_QUARTERS_PER_YEAR) // number of timesteps in a single Monte Carlo simulation
#define N_VAR_TYPES 3
#define N_VARS (N_QUARTERS_PER_YEAR*N_VAR_TYPES)
#define N_OBJS 5
#define N_CONSTS 3

#define MIN_ALK 0.0     // minimum value for alkalinity setpoint #1
#define MAX_ALK 150.0   // maximum value for alkalinity setpoint #1
#define MIN_PH 6.0      // minimum value for pH setpoint #1
#define MAX_PH 9.5      // maximum value for pH setpoint #1
#define MIN_DBP_SF 0.02 // minimum value for disinfection byproduct safety factor
#define MAX_DBP_SF 0.75 // maximum value for disinfection byproduct safety factor

// Macros for booleans
#define FALSE 0
#define TRUE 1

// Declare global variables
extern std::string current_datetime; /* String which contains the date and time at the start of running the program */
extern struct SimOptParameters simopt_params;  /* Optimize (i.e., simulation-optimization) paramaters */

/* Purpose: this header contains function declarations and structure definitions to implement the various 
* run modes in the wtp-optimize project. Including:
* - single simulation
* - single simulation (with automatic chemical dosing)
* - multi-simulation 
* - multi-simulation (with automatic dosing)
* - simulation-optimization 
* - validation
*/

/* Structure definitions */ 
// // run configuration data
// struct RunConfig {
//     const char * treatment_type;
//     int run_mode;
//     union {
//         struct SingleSim;
//         struct MultiSim;
//     };
// };

// struct SingleSimParameters
// {
//     bool run_model;
//     bool list_process_train;
//     bool thm_and_disinfection;
// };

struct SingleSimParameters
{
    const char* operations_filename;
    const char* influent_filename;
};

// struct TimeSeriesParameters
// { // water quality time series parameters
//     bool ts_flag;
//     int n_years;
//     int timesteps_per_year;
// };

// struct MonteCarloParameters
// { // water quality monte carlo parameters
//     bool mc_flag;
//     int n_samples;
// };

// struct SimOptParameters
// { // simulation-optimization parameters
//     struct MonteCarloParameters mc;
//     struct TimeSeriesParameters ts;
//     int num_func_evals;
//     int problem;
// };

struct SimOptParameters
{ // simulation-optimization parameters
    int num_func_evals;
    int num_wq_scenarios;
};

struct OperationalParameters
{   // operational parameters
    double alk_setpt;  // alkalinity setpoint at front of plant
    double pH_setpt;   // pH setpoint at front of plant
    double DBPsf;      // disinfection byproduct safety factor at end of treatment system
};

struct InfluentParameters
{   // influent water quality parameters
    double alkalinity;  // alkalinity in mg/L as CaCO3
    double nh3;  // ammonia in mg/L
    double bromide;  // bromide in mg/L 
    double calcium;  // calcium in mg/L as CaCO3
    double hardness;  // total hardness in mg/L as CaCO3
    double pH;  // pH
    double temp;  // temperature in degrees Celcius
    double toc;  // total organic carbon in mg/L as C
    double ntu;  // turbidity in ntu
    double uv254;  // ultraviolet absorbance at 254 nm
};

/* Function declarations */
// read_montecarlo.cpp
// void read_montecarlo(double ***samples, int n_samples, int n_timesteps, int n_params, FILE *fin); /* Read in Monte Carlo samples of influent water quality */
std::vector<double> read_single_sim(std::string filename, bool header);
std::vector< std::vector<double> > read_montecarlo(std::string filename, bool header, int nsims);

// wtp_optimize.cpp
// void single_sim_mode(struct ProcessTrain *train, struct SingleSimParameters params);
// void single_sim(struct ProcessTrain *train,   // Process train control structure    
//                      std::string current_datetime, // Current date and time (string)
//                      FILE *fout,                   // Output file     
//                      struct SingleSimParameters params                  
// );
int single_sim_mode(struct ProcessTrain *train, /* Process train control structure    */
                     struct OperationalParameters *operations,  // operational parameters
                     struct InfluentParameters *influent_wq,
                     FILE *fout /* Output file */
);
// void multi_sim_mode(struct ProcessTrain *train, struct MultiSimParameters params);
// void multi_sim_autodose(struct ProcessTrain *train, struct MultiSimAutoDose params);
// void sim_opt_mode(struct ProcessTrain *train, struct SimOptParameters params);
void sim_opt_mode(struct ProcessTrain *train);
// void validation_mode(struct ProcessTrain *train, struct ValidationParameters params);

// wtp_problem.cpp
void wtp(double* vars, double* objs, double* consts);  // Water Treatment Plant Model problem definition
void validate_wq_nonnegative (double value, const char* name);  // validate that water quality parameter is non-negative
void validate_wq_bounds(double value, const char *name, double minimum, double maximum);  // validate that water quality parameter is between some specified bounds
void validate_vars_bounds(double value, const char *name, double minimum, double maximum);  // validate that the decision variable is between some specified bounds
void validate_sim_year_quarter(int expected, int actual, const char* name, int column, const char* filename);  // validate that the simulation number, year, and quarter are being read in correctly

#endif