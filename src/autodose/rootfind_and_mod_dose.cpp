/* rootfind_and_mod_dose.c */

#include "wtp.h"
#include "auto_dose.h"

/* root-finding for the WTP Model
 * 
 * Purpose: rootfind_and_mod_dose() solves for the root of f(x) (i.e. the x-value which yields f(x) = 0).
 *  Where f(x) is the difference between a given water quality parameter somewhere within the treatment plant 
 *  and the target water quality of that parameter somewhere within the treatment plant. 
 *  Where x is a chemical dose which controls that water quality parameter. Once an x-value (chemical dose) 
 *  is found which yields f(x) = 0, that chemical dosage is altered in the treatment train. 
 * 
 *  To account for multiple roots or negative roots (negative roots would translate into a non-physical
 *  negative chemical dose), rootfind_and_mod_dose() will implement the minimum non-negative root that is found. 
 *  For situations in which there are no roots, the non-negative x-value which yields the minimum |f(x)| found during
 *  the search will be implemented. 
 * 
 */ 

double rootfind_and_mod_dose(RF_FUNC_PTR func, struct ProcessTrain *train, int dose_unit, int dose_location, 
                           char target_param, int target_unit, int target_location, double target, 
                           double x_lo, double x_up, FILE *fout) {

    /* Purpose: use Bisection or Secant method to find root of function. 
     * 
     * Inputs:
     *  func = function pointer
     *  x_0 = initial guess
     *  err_tol = error tolerance of search
     *  max_iter = maximum iterations allowed
     * 
     * Return:
     *  x = value of root found by search
     */
   
    double err_tol = ERROR_TOL;
    double max_iter = 25; 
    double x, x_0, x_1, f_0, f_1, x_mid, f_mid, error;
    int i, j;
    double x_lo_guess = x_lo;  // save initial lower guess
    double x_up_guess = x_up;  // save initial upper guess
    int debug = TRUE;  // (TRUE/FALSE): TRUE to output debugging print statements, FALSE for no output
    
    if (debug==TRUE) fprintf(fout,"=========== START rootfind_and_mod_dose() =========== \n\n");
    
    f_0 = func(train, dose_unit, x_lo, dose_location, target_param, 
                    target_unit, target, target_location, fout); 
    if (debug==TRUE) fprintf(fout, "rootfind_and_mod_dose: initial f(x_min) = %.3f \n", f_0);                
    
    f_1 = func(train, dose_unit, x_up, dose_location, target_param, 
                    target_unit, target, target_location, fout);    
    if (debug==TRUE) fprintf(fout, "rootfind_and_mod_dose: initial f(x_max) = %.3f \n", f_1); 
    
    if (f_0*f_1 < 0)  /* Check bracketing of lower and upper guesses */
    {
    /* ------------------ BISECTION METHOD ------------------ */ 
		/* Purpose: implement the Bisection method (https://en.wikipedia.org/wiki/Bisection_method) for root-finding.
		 * This bisection method is modified to solve for multiple roots, if they exist. 
		 * Of those roots, this code will return the minimum non-negative root. This method
		 * guarantees that a root is found given sufficient iterations. */
        
        if (debug==TRUE) fprintf(fout,"\nBracketing conditions met, using BISECTION METHOD... \n");
        
        /* Define array of roots to choose from if there are multiple roots */
        int max_starts = 3;  // maximum number of roots that can be found
        double roots[] = {DBL_MAX, DBL_MAX, DBL_MAX};  // initialize roots to maximum possible value
        double min_root; 
        
		/* Search for multiple roots */
        for (j=0; j<max_starts; j++) {
            
            fprintf(fout, "\n============= BISECTION METHOD: iteration %d ============== \n", j+1);
            
			/* After the first root is found (for j==0), look for other roots by updating lower and upper bounds */
            if (j==1) {  // once first root is found, check between the lower bound and just left of the root
                x_lo = x_lo_guess;
                x_up = roots[0] - 0.1;  // bracket the lower guess and just to the left of the root that was found
             
            } else if (j==2) { // then, check between just right of the root and the upper bound
                x_lo = roots[0];
                x_up = x_up_guess + 0.1; // bracket just to the right of the root that was found and the upper bound
            }
            
            if (debug==TRUE) {
                fprintf(fout, "\n");
                fprintf(fout, "low guess (x_lo) = %f\n", x_lo);
                fprintf(fout, "upper guess (x_up) = %f\n", x_up);
                fprintf(fout, "\n");
            }
            
			/* Perform Bisection Method */
            for (i=0; i<max_iter; i++) {
                x_mid = (x_lo + x_up)/2.0; 
                f_mid = func(train, dose_unit, x_mid, dose_location, target_param, 
                            target_unit, target, target_location, fout); 
                if (debug==TRUE) {
                    fprintf(fout, "rootfind_and_mod_dose: current x = %.8f \n", x_mid); 
                    fprintf(fout, "rootfind_and_mod_dose: current f(x) = %.8f \n", f_mid);  
                    fprintf(fout, "\n");
                }
                
                error = fabs(f_mid);
                
                if ((f_0 * f_mid) <= 0){
                    x_up = x_mid; 
                } else {
                    x_lo = x_mid; 
                }
                
                if (fabs(error) < err_tol) {
                    roots[j] = x_mid;
                    if (debug==TRUE) fprintf(fout, "Root found!\n");
                    break;
                }
            }
        }

		/* Choose minimum non-negative root to return from function */
        min_root = find_min_nonneg(roots, max_starts); 
        if (debug==TRUE) fprintf(fout,"rootfind_and_mod_dose: Minimum root found = %f\n\n", min_root); 
        
		/* Handle the case where no root has been found */
        if (min_root == DBL_MAX) { 
            fprintf(stderr, "No roots found with Bisection method. Increase maximum iterations! \n");
            exit(EXIT_FAILURE); 
        }
        
        /* Apply minimum, non-negative root found to change the chemical dose */
        double f_min_root = func(train, dose_unit, min_root, dose_location, target_param, 
                                target_unit, target, target_location, fout); 
        if (debug==TRUE) fprintf(fout, "rootfind_and_mod_dose: f(x) of minimum root = %f\n\n", f_min_root); 
        if (debug==TRUE) fprintf(fout, "============= END BISECTION METHOD: rootfind_and_mod_dose() ============\n\n");
        
        return  min_root;
    
    } else {    
    /* ------------------ SECANT METHOD ------------------ */ 
		/* Purpose: implement the Secant method (https://en.wikipedia.org/wiki/Secant_method) for root-finding.
		 * This Secant method is modified to solve for multiple roots, if they exist. 
		 * Of those roots, this code will return the minimum positive root. If no root exists, 
		 * the minimum |f(x)| which has a non-negative x-value is returned instead. */
		 
		fprintf(fout,"Bracketing conditions not met, using SECANT METHOD instead of Bisection... \n");
		        
        /* Keep track of minimum and maximum function evalutions during root-finding search */ 
        double bracket_array[] = {f_0, f_1};
        double f_min=DBL_MAX, f_max=-DBL_MAX;      // The minimum and maximum observed function values, respectively. Initialize to limits of values for doubles. 
        double x_min=DBL_MAX, x_max=-DBL_MAX;      // The x-values which corresponds to the minimum and maximum observed function values, repectively. Initialize to limits of values for doubles.
		int maximum;  // set maximum = 1 to get maximum, and maximum = 0 to get minimum 

        int min_index = find_extrema_index(bracket_array, sizeof(bracket_array)/sizeof(double), maximum=0);
        
		/* Keep track of the positive x-value which yields the minimum and maximum f(x) (i.e. x_min and x_max).
		 * In the case that there are no roots found, this value will be returned. */
        if (min_index == 0) {
            if (x_lo >= 0) x_min = x_lo;  // only set x_min value if x_lo is positive
            if (x_up >= 0) x_max = x_up;  // only set x_max value if x_up is positive
        } else if (min_index == 1) {
            if (x_lo >= 0) x_max = x_lo;  // only set x_min value if x_lo is positive
            if (x_up >= 0) x_min = x_up;  // only set x_max value if x_up is positive
        } else {
            fprintf(fout,"Error setting x_min and x_max. Exiting program...\n");
            exit(EXIT_FAILURE);
        }
        
        /* Define array of roots to choose from if there are multiple roots */
        int max_starts = 10; 
        double roots[max_starts];  
        double min_root; 

        for (j=0; j < max_starts; j++) { 
            
            fprintf(fout, "\n============= SECANT METHOD: iteration %d ==============\n", j+1);
            
            roots[j] = DBL_MAX;  // set default value of roots to NaN. If NaN remains after root-finding, no root was found. 
            x_0 = x_lo_guess + fabs(x_lo_guess-x_up_guess)*((double)(j)/(double)(max_starts)); 
         
            /* Quality check inputs */
            if (err_tol <= 0) {
                fprintf(stderr,"Error tolerance cannot be negative! Setting tolerance = 0.001 \n"); 
                err_tol = 0.001;
            }
               
            if (max_iter < 1) {
                fprintf(stderr,"Maximum iterations must be greater than 1! Setting maxiter = 30 \n");
                max_iter = 30; 
            }

            /* Change guesses for Secant Method */
            if (x_0 >= 0) {
                x_1 = x_0*(1 + 1E-4) + 1E-4;
            } else {
                x_1 = x_0*(1 + 1E-4) - 1E-4;
            }

            f_0 = func(train, dose_unit, x_0, dose_location, target_param, 
                      target_unit, target, target_location, fout); 
            
            f_1 = func(train, dose_unit, x_1, dose_location, target_param, 
                      target_unit, target, target_location, fout); 
            
            /* Determine x which corresponds to minimum and maximum observed function values 
			 * thus far which is non-negative*/ 
			x_min = update_xextrema_nonneg(f_min, f_0, f_1, x_min, x_0, x_1, maximum=0); 
            x_max = update_xextrema_nonneg(f_max, f_0, f_1, x_max, x_0, x_1, maximum=1);
            
            /* Update f_min and f_max based on x's found above */ 
            f_min = func(train, dose_unit, x_min, dose_location, target_param, 
                        target_unit, target, target_location, fout); 

            f_max = func(train, dose_unit, x_max, dose_location, target_param, 
                        target_unit, target, target_location, fout); 
            
			/* Perform Secant Method */
            for (i=1; i<= max_iter; i++) {
                if (f_0 == f_1) {
                    
                    if (x_0 != x_1) {
                        fprintf(fout,"Tolerance of %f reached. \n", (x_0 - x_1));
                    }
                    
                    x = (x_0 + x_1)/2.0; 
                    break; 
                    
                } else {
                        x = x_1 - f_1*(x_1 - x_0)/(f_1 - f_0);
                }
                
                if (fabs(x - x_1) < err_tol) {
                    roots[j] = x; 
                    break; 
                }
                
                x_0 = x_1;
                f_0 = f_1; 
                x_1 = x;
                f_1 = func(train, dose_unit, x_1, dose_location, target_param, 
                          target_unit, target, target_location, fout); 
                
                /* Determine x which corresponds to minimum and maximum observed function values thus far */ 
				x_min = update_xextrema_nonneg(f_min, f_0, f_1, x_min, x_0, x_1, maximum = 0); 
                x_max = update_xextrema_nonneg(f_max, f_0, f_1, x_max, x_0, x_1, maximum = 1);
            
                /* Update f_min and f_max based on x's found above */ 
                f_min = func(train, dose_unit, x_min, dose_location, target_param, 
                            target_unit, target, target_location, fout); 
                if (debug==TRUE) {
                    fprintf(fout, "rootfind_and_mod_dose: x_min = %.8f \n", x_min); 
                    fprintf(fout, "rootfind_and_mod_dose: f(x_min) = %.8f \n", f_min);  
                    fprintf(fout, "\n");
                }
                
                f_max = func(train, dose_unit, x_max, dose_location, target_param, 
                            target_unit, target, target_location, fout);
                            
                if (debug==TRUE) {
                    fprintf(fout, "rootfind_and_mod_dose: x_max = %.8f \n", x_max); 
                    fprintf(fout, "rootfind_and_mod_dose: f(x_max) = %.8f \n", f_max);  
                    fprintf(fout, "\n");
                }
			}
        }
        
        min_root = find_min_nonneg(roots, max_starts); 
        
        if (min_root == DBL_MAX) {  // this statement means that no root was found
            if (debug==TRUE) fprintf(fout,"No roots found: use x value which corresponds to function evalution closest to the root instead. \n");
            if ((fabs(f_min)) < fabs(f_max)) {
                /* Check to make sure that the x_min is not the default value */
                if (x_min == DBL_MAX) {  // this statement means that root was not found
                    fprintf(fout,"Error: no root found \n");
                    exit(EXIT_FAILURE); 
				} 
                
                /* Apply minimum, non-negative root found to change the chemical dose */
                func(train, dose_unit, x_min, dose_location, target_param, 
                        target_unit, target, target_location, fout); 
                if (debug==TRUE) {
                    fprintf(fout, "rootfind_and_mod_dose: x_min = %.8f \n", x_min); 
                    fprintf(fout, "rootfind_and_mod_dose: f(x_min) = %.8f \n", f_min);  
                    fprintf(fout, "\n");
                }
                
                return x_min;
                
            } else {
                /* Check to make sure that the x_max is not the default value */
                if (x_max == DBL_MAX) {  // this statement means that root was not found
                    fprintf(fout,"Error: no root found \n");
                    exit(EXIT_FAILURE); 
				} 
                
                /* Apply minimum, non-negative root found to change the chemical dose */
                func(train, dose_unit, x_max, dose_location, target_param, 
                        target_unit, target, target_location, fout); 
                if (debug==TRUE) {
                    fprintf(fout, "rootfind_and_mod_dose: x_max = %.8f \n", x_max); 
                    fprintf(fout, "rootfind_and_mod_dose: f(x_max) = %.8f \n", f_max);  
                    fprintf(fout, "\n");
                }
                
                return x_max;
            }
        }
        
        /* Apply minimum, non-negative root found to change the chemical dose */
        double f_min_root = func(train, dose_unit, min_root, dose_location, target_param, 
                                target_unit, target, target_location, fout); 
        if (debug==TRUE) {
            fprintf(fout, "rootfind_and_mod_dose: minimum root = %.8f \n", min_root); 
            fprintf(fout, "rootfind_and_mod_dose: f(x) of minimum root = %.8f \n", f_min_root);  
            fprintf(fout, "\n");
        }
        
        return min_root; 
    }
}

double mod_dose_check_target(struct ProcessTrain *train, int dose_unit,   double dose,   int dose_location, 
                           char target_param, int target_unit, double target, int target_location, FILE *fout)
     /* Purpose: 
     *  Modify the chemical dose at particular unit process, then calculate the difference between the unit process 
     *  effluent water quality and the target (or setpoint) value at a target unit process. 
     * 
     * Inputs:
     *  train           = Process train control structure.
     *  dose_unit       = Unit process at which the chemical dose is being added.
     *  dose            = Magnitude of dose (mg/L).  
     *  target_param    = Target water quality parameter.
     *  target_unit     = Unit process at which to check whether the target setpoint is met. 
     *  target_location = For the case that there are multiple of the same processes in a treatment train (e.g. 
     *                    two CO2 addition points), input the location for the particular instance of the unit process
     *                    to which the setpoint refers. 0 for the first unit process, 1 for the second, 2 for the third, etc. 
     *  fout            = Output file to store print statements (for debugging). 
     * 
     * Outputs: 
     *  eff - target = Difference between effluent and target water quality.
     * 
     * Note: 
     *  'unit process effluent' refers to the water quality leaving an individual unit process. This is not to be confused with the 
     *   effluent water quality of the treatment plant, overall. 
     */
{
    int    dose_cntr = 0; 
    int    target_cntr = 0; 
    double eff; 
    register struct UnitProcess *unit;
    struct Effluent             *influent;
    
    /* Debugging: print train output */ 
    int    debug = TRUE; 

    /* Set chemical dosing (supports lime, CO2, alum, and hypochlorite) */ 
    for( unit=FirstUnitProcess(train); unit; unit=NextUnitProcess(unit) )
          {
              switch (unit->type)
              {
                  case LIME:
                    if (dose_unit == LIME){
                        if (dose_cntr == dose_location) { 
                            unit->data.lime->dose = dose;  // change chemical dose for unit process 
                        }
                        dose_cntr++;
                    }
                    break; 

                  case CARBON_DIOXIDE:
                    if (dose_unit == CARBON_DIOXIDE){
                        if (dose_cntr == dose_location) { 
                            unit->data.chemical->co2 = dose;  // change chemical dose for unit process 
                        }
                        dose_cntr++;
                    }
                    break; 
                  
                  case ALUM: 
                    if (dose_unit == ALUM){
                        if (dose_cntr == dose_location) { 
                            unit->data.alum->dose = dose;  // change chemical dose for unit process 
                        }
                    dose_cntr++;
                    }
                    break; 
                    
                  case HYPOCHLORITE: 
                    if (dose_unit == HYPOCHLORITE){
                        if (dose_cntr == dose_location) {
                            unit->data.chemical->naocl = dose;  // change chemical dose for unit process 
                        }
                    dose_cntr++; 
                    }
                    break; 
                    
                  default:break; 
              }
          }
    
    /* Run model*/
    runmodel(train);

    /* Check to see if effluent water quality from unit process matches target setpoint */
    for( unit=FirstUnitProcess(train); unit; unit=NextUnitProcess(unit) )
        {
          switch (unit->type)
          {
              case INFLUENT:
                influent = &unit->eff; // get influent water quality to calculate % removal for TOC
                break;

              case CARBON_DIOXIDE:
                if (target_unit == CARBON_DIOXIDE){  
                    if (target_location == target_cntr){  // if there are multiple unit processes of this type
                        /* Get value for CO2 unit process effluent: either pH or alkalinity */
                        if (target_param == 'P'){ 
                            eff = unit->eff.pH;
                            if (debug == TRUE) 
                            {
                                fprintf(fout, "mod_dose_check_target(): pH at co2 = %f \n", eff);     
							    fprintf(fout, "mod_dose_check_target(): alk at co2 = %f \n", unit->eff.Alk * MW_CaCO3 / 2.0);
                            }   
                        } else if (target_param == 'A') {
                            eff = unit->eff.Alk * MW_CaCO3 / 2.0;
                            if (debug == TRUE) 
                            {
                                fprintf(fout, "mod_dose_check_target(): alk at co2 = %f \n", eff);
                            }
                        } else {
                            fprintf(stderr, "Incorrect target parameter input for CO2! \n"); 
                            exit(EXIT_FAILURE);
                        }
                    }
                    target_cntr++; 
                }
                break; 
              
              case RAPID_MIX:
                if (target_unit == RAPID_MIX){
                    
                    if (debug==TRUE) fprintf(fout, "pH at rapid mix = %f\n", unit->eff.pH);
                    
                    /* Get value for rapid mix effluent: pH */
                    if (target_param == 'P')       eff = unit->eff.pH;
                    else {
                        fprintf(stderr, "Incorrect target parameter input for Rapid Mix! \n"); 
                        exit(EXIT_FAILURE);
                    }
                }
                break;   
                
              case WTP_EFFLUENT:
                if (target_unit == WTP_EFFLUENT){
                    /* Get value for water treatment plant effluent: either pH or alkalinity */
                    if (target_param == 'P') {
                        eff = unit->eff.pH;
                        if (debug==TRUE) fprintf(fout, "pH at effluent = %f \n", eff);
                    } else if (target_param == 'A')  {
                        eff = unit->eff.Alk * MW_CaCO3 / 2.0;
                        if (debug==TRUE) fprintf(fout, "alk at effluent = %f \n", eff);
                    } else {
                        fprintf(stderr, "Incorrect target parameter input for Effluent! \n"); 
                        exit(EXIT_FAILURE);
                    }
                }
                break; 
                
              case END_OF_SYSTEM:
                if (target_unit == END_OF_SYSTEM){
                    /* Get value for end of distribution system effluent: either pH, alkalinity, TTHM, HAA5,
                        % TOC removal, or chlorine residual */
                    if (target_param == 'P') {
                        eff = unit->eff.pH;
                        if (debug==TRUE) fprintf(fout, "pH at EOS = %.8f \n", eff);
                    } else if (target_param == 'A') {
                        eff = unit->eff.Alk * MW_CaCO3 / 2.0;
                        if (debug==TRUE) fprintf(fout, "alk at EOS = %.8f \n", eff);
                    } else if (target_param == 'T') {
                        eff = unit->eff.TTHM;
                        if (debug==TRUE) fprintf(fout, "TTHM at EOS = %.8f \n", eff);
                    } else if (target_param == 'H') {
                        eff = unit->eff.HAA5; 
                        if (debug==TRUE) fprintf(fout, "HAA5 at EOS = %.8f \n", eff);
                    } else if (target_param == 'O') {
                        eff = ((influent->TOC - unit->eff.TOC) / influent->TOC * 100.0); 
                        if (debug==TRUE) fprintf(fout, "TOC rmv at EOS = %.8f \n", eff);
                    } else if (target_param == 'C') {
                        eff = unit->eff.FreeCl2 * MW_Cl2; 
                        if (debug==TRUE) fprintf(fout, "cl2 at EOS = %.8f \n", eff);
                    } else {
                        fprintf(stderr, "Incorrect target parameter input for End of System! \n"); 
                        exit(EXIT_FAILURE);
                    }
                }
                break; 
                    
                default:break; 
            }
        }
        
    return eff-target;  // difference between unit process effluent and target setpoint for a particular water quality parameter
}

/* WJR: moved this to extrema.c */ 

//double find_min_nonneg(double array[], int n) {
//    
//    /* Purpose: Finds minimum non-negative value of an array with n elements */
//  
//    int i; 
//    double min = DBL_MAX;
//     
//    for (i=0; i<n; i++) {
//        if (array[i] < min && array[i] >= 0.0) {  // set min if next number is lower or has not been set before
//           min = array[i];
//        }
//    }
//    
//    return min;
//}
//
//int find_extrema_index(double array[], int n, int maximum) {
//    
//    /* Purpose: Finds the index of the extrema (either maximum or minimum) of an array with n elements */
//	/* If maximum=1, it will find the index for the maximum. If maximum=0, it will find the index for the minimum */
//    /* Source: http://www.programmingsimplified.com/c/source-code/c-program-find-minimum-element-in-array */
//    /* Source: https://stackoverflow.com/questions/25500199/find-smallest-positive-number-in-an-array */
//  
//    int i, index; 
//	index = 0;  // index which results in minimum
//    
//	if (maximum == 0) {
//		double min = DBL_MAX;
//		 
//		for (i=0; i<n; i++) {
//			if (array[i] < min) {  // set min if next number is lower
//			   min = array[i];
//			   index=i;
//			}
//		}
//	} else if (maximum == 1) {
//		double max = -DBL_MAX;
//		 
//		for (i=0; i<n; i++) {
//			if (array[i] > max) { // set max if next number is higher
//			   max = array[i];
//			   index=i;
//			}
//		}
//	} else {
//		printf("Incorrect input for 'maximum' in find_extrema_index()! Input maximum=1 to get x_max or maximum=0 to get x_min.\n"); 
//		exit(EXIT_FAILURE);
//	}
//    
//    return index;
//}
//
//double update_xextrema_nonneg(double f_extrema, double f_0, double f_1, double x_extrema, double x_0, double x_1, int maximum) {
//    
//    /* Purpose: update x_min or x_max (i.e. the non-negative x-values which corresponds to the minimum or maximum observed 
//	 * function values, f(x)) 
//	 * 
//	 * If maximum = 1, find x_max. If maximum = 0, find x_min. */
//	
//	/* Remove negative values of x from consideration and update x based on desired extrema of f(x)*/
//	if (x_extrema < 0.0) {  // if x_extrema is negative
//	
//		if (x_0 < 0.0) {  // if x_extrema and x_0 are negative
//		
//			if (x_1 < 0.0) {  // if all three (x_extrema, x_0, and x_1) are negative
//				// do not update x_extrema
//				
//			} else {  // if x_extrema and x_0 are negative but x_1 is non-negative
//				x_extrema = x_1;
//			}
//			
//		} else { // if x_0 is non-negative
//		
//			if (x_1 < 0.0) {  // if x_0 is non-negative but both x_extrema and x_1 are negative
//				x_extrema = x_0; 
//					
//			} else { // if x_extrema is negative but both x_0 and x_1 are non-negative
//				double extrema_array[] = {f_0, f_1};
//				int index = find_extrema_index(extrema_array, sizeof(extrema_array)/sizeof(double), maximum);
//				if (index==0) {
//					x_extrema = x_0;  
//				} else if (index==1) {
//					x_extrema = x_1;
//				} else {
//					printf("Error: Improper index found in update_x_nonneg()! \n");
//					exit(EXIT_FAILURE);
//				}
//			}
//		}
//		
//	} else {  // if x_extrema is non-negative
//		
//		if (x_0 < 0.0) {  // if x_0 is negative
//		
//			if (x_1 < 0.0) {  // if x_extrema is non-negative but both x_0 and x_1 are negative
//				// do not update x_extrema
//				
//			} else {  // if x_0 is negative but both x_extrema and x_1 are non-negative
//				double extrema_array[] = {f_extrema, f_1};
//				int index = find_extrema_index(extrema_array, sizeof(extrema_array)/sizeof(double), maximum); 
//				if (index==0) {
//					// do not update x_extrema
//				} else if (index==1) {
//					x_extrema = x_1;
//				} else {
//					printf("Error: Improper index found in update_xextrema_nonneg()! \n");
//					exit(EXIT_FAILURE);
//				}
//			}
//			
//		} else { // if x_0 is non-negative
//		
//			if (x_1 < 0.0) {  // if x_extrema and x_0 are non-negative and x_1 is negative
//				double extrema_array[] = {f_extrema, f_0};
//				int index = find_extrema_index(extrema_array, sizeof(extrema_array)/sizeof(double), maximum);
//				if (index==0) {
//					// do not update x_extrema
//				} else if (index==1) {
//					x_extrema = x_0; 
//				} else {
//					printf("Error: Improper index found in update_xextrema_nonneg()! \n");
//					exit(EXIT_FAILURE);
//				}
//					
//			} else { // if all three x's (x_extrema, x_0, and x_1) are non-negative
//				double extrema_array[] = {f_extrema, f_0, f_1};
//				int index = find_extrema_index(extrema_array, sizeof(extrema_array)/sizeof(double), maximum);
//				if (index==0) {
//					// do not update x_extrema
//				} else if (index==1) {
//					x_extrema = x_0;  // update x_extrema to x_0 only if x_0 is positive
//				} else if (index==2) {
//					x_extrema = x_1;  // update x_extrema to x_1 only if x_1 is positive
//				} else {
//					printf("Error: Improper index found in update_xextrema_nonneg()! \n");
//					exit(EXIT_FAILURE);
//				}
//			}
//		}
//	}
//    
//    return x_extrema; 
//}