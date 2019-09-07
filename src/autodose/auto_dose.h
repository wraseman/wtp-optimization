/* auto_dose.h */

#ifndef AUTO_DOSE_H
#define AUTO_DOSE_H 1

#include "wtp/wtp.h"
#include <string>

/* Automatic chemical dosing functions
 * The following functions and constants are used to find chemical doses which meet water quality setpoints
 */

/* TODO add comments to auto_dose.h */
double mod_dose_check_target(struct ProcessTrain *train, int dose_unit, double dose, int dose_location,
                             char target_param, int target_unit, double target, int target_location, FILE *fout);
void auto_dose(struct ProcessTrain *train,
               double pH_setpt_1, double alk_setpt_1, double cl2_setpt, double pH_setpt_2, double DBP_safety_factor, 
               FILE *fres);
// void meet_TOC_req(struct ProcessTrain *train, double cl2_setpt, double pH_setpt_2);
void write_dose_target(int dose_unit, int dose_location, char target_param, int target_unit, int target_location, double target, FILE *fout);
typedef double (*RF_FUNC_PTR)(struct ProcessTrain *, int, double, int, char, int, double, int, FILE *);
// a function pointer passed into rootfind_and_mod_dose() to evaluate the impact of changing a chemical dose
double rootfind_and_mod_dose(RF_FUNC_PTR func, struct ProcessTrain *train, int dose_unit, int dose_location, char target_param,
                             int target_unit, int target_location, double target, double x_lo, double x_up, FILE *fout);

/* Functions in extrema.cpp used to find max or mins of an array */
double min(double array[], int n); // finds minimum value of an array of doubles of length n
double max(double array[], int n); // finds maximum value of an array of doubles of length n
double find_min_nonneg(double array[], int n);
int find_extrema_index(double array[], int n, int maximize);
double update_xextrema_nonneg(double f_extrema, double f_0, double f_1, double x_extrema, double x_0, double x_1, int maximum);

// Global variable to keep track of date and time of run
// extern std::string current_datetime; /* String which contains the date and time at the start of running the program */

#endif 