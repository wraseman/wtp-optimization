/* extrema.c */
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

/* Purpose: this file contains functions which find extremas (min or max) in an array */

double min(double array[], int n) {
	
	/* Purpose: Finds minimum value of an array of doubles of length n */
	
	int i; 
	double min = DBL_MAX;
	 
	for (i=0; i<n; i++) {
		if (array[i] < min) {  // set min if next number is lower or has not been set before
		   min = array[i];
		}
	}
	
	return min;
}

double max(double array[], int n) {
	
	/* Purpose: Finds maximum value of an array of doubles of length n */
	
	int i; 
	double max = -DBL_MAX;
	 
	for (i=0; i<n; i++) {
		if (array[i] > max) {  // set max if next number is higher or has not been set before
		   max = array[i];
		}
	}
	
	return max;
}

double find_min_nonneg(double array[], int n) {
    
    /* Purpose: Finds minimum non-negative value of an array with n elements */
  
    int i; 
    double min = DBL_MAX;
     
    for (i=0; i<n; i++) {
        if (array[i] < min && array[i] >= 0.0) {  // set min if next number is lower or has not been set before
           min = array[i];
        }
    }
    
    return min;
}

int find_extrema_index(double array[], int n, int maximum) {
    
    /* Purpose: Finds the index of the extrema (either maximum or minimum) of an array with n elements */
	/* If maximum=1, it will find the index for the maximum. If maximum=0, it will find the index for the minimum */
    /* Source: http://www.programmingsimplified.com/c/source-code/c-program-find-minimum-element-in-array */
    /* Source: https://stackoverflow.com/questions/25500199/find-smallest-positive-number-in-an-array */
  
    int i, index; 
	index = 0;  // index which results in minimum
    
	if (maximum == 0) {
		double min = DBL_MAX;
		 
		for (i=0; i<n; i++) {
			if (array[i] < min) {  // set min if next number is lower
			   min = array[i];
			   index=i;
			}
		}
	} else if (maximum == 1) {
		double max = -DBL_MAX;
		 
		for (i=0; i<n; i++) {
			if (array[i] > max) { // set max if next number is higher
			   max = array[i];
			   index=i;
			}
		}
	} else {
		printf("Incorrect input for 'maximum' in find_extrema_index()! Input maximum=1 to get x_max or maximum=0 to get x_min.\n"); 
		exit(EXIT_FAILURE);
	}
    
    return index;
}

double update_xextrema_nonneg(double f_extrema, double f_0, double f_1, double x_extrema, double x_0, double x_1, int maximum) {
    
    /* Purpose: update x_min or x_max (i.e. the non-negative x-values which corresponds to the minimum or maximum observed 
	 * function values, f(x)) 
	 * 
	 * If maximum = 1, find x_max. If maximum = 0, find x_min. */
	
	/* Remove negative values of x from consideration and update x based on desired extrema of f(x)*/
	if (x_extrema < 0.0) {  // if x_extrema is negative
	
		if (x_0 < 0.0) {  // if x_extrema and x_0 are negative
		
			if (x_1 < 0.0) {  // if all three (x_extrema, x_0, and x_1) are negative
				// do not update x_extrema
				
			} else {  // if x_extrema and x_0 are negative but x_1 is non-negative
				x_extrema = x_1;
			}
			
		} else { // if x_0 is non-negative
		
			if (x_1 < 0.0) {  // if x_0 is non-negative but both x_extrema and x_1 are negative
				x_extrema = x_0; 
					
			} else { // if x_extrema is negative but both x_0 and x_1 are non-negative
				double extrema_array[] = {f_0, f_1};
				int index = find_extrema_index(extrema_array, sizeof(extrema_array)/sizeof(double), maximum);
				if (index==0) {
					x_extrema = x_0;  
				} else if (index==1) {
					x_extrema = x_1;
				} else {
					printf("Error: Improper index found in update_x_nonneg()! \n");
					exit(EXIT_FAILURE);
				}
			}
		}
		
	} else {  // if x_extrema is non-negative
		
		if (x_0 < 0.0) {  // if x_0 is negative
		
			if (x_1 < 0.0) {  // if x_extrema is non-negative but both x_0 and x_1 are negative
				// do not update x_extrema
				
			} else {  // if x_0 is negative but both x_extrema and x_1 are non-negative
				double extrema_array[] = {f_extrema, f_1};
				int index = find_extrema_index(extrema_array, sizeof(extrema_array)/sizeof(double), maximum); 
				if (index==0) {
					// do not update x_extrema
				} else if (index==1) {
					x_extrema = x_1;
				} else {
					printf("Error: Improper index found in update_xextrema_nonneg()! \n");
					exit(EXIT_FAILURE);
				}
			}
			
		} else { // if x_0 is non-negative
		
			if (x_1 < 0.0) {  // if x_extrema and x_0 are non-negative and x_1 is negative
				double extrema_array[] = {f_extrema, f_0};
				int index = find_extrema_index(extrema_array, sizeof(extrema_array)/sizeof(double), maximum);
				if (index==0) {
					// do not update x_extrema
				} else if (index==1) {
					x_extrema = x_0; 
				} else {
					printf("Error: Improper index found in update_xextrema_nonneg()! \n");
					exit(EXIT_FAILURE);
				}
					
			} else { // if all three x's (x_extrema, x_0, and x_1) are non-negative
				double extrema_array[] = {f_extrema, f_0, f_1};
				int index = find_extrema_index(extrema_array, sizeof(extrema_array)/sizeof(double), maximum);
				if (index==0) {
					// do not update x_extrema
				} else if (index==1) {
					x_extrema = x_0;  // update x_extrema to x_0 only if x_0 is positive
				} else if (index==2) {
					x_extrema = x_1;  // update x_extrema to x_1 only if x_1 is positive
				} else {
					printf("Error: Improper index found in update_xextrema_nonneg()! \n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}
    
    return x_extrema; 
}