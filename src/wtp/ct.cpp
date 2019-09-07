/* CT.c */
#include "wtp.h"

void ct(struct UnitProcess *unit, double t_ten)

/*
*  Purpose: Use values of log_inactivation required (calculated once in influent
*           subroutine) to estimate inactivation ratios for Giardia, Virus, and
*           Cryptosporidium.  The overall inactivation ratio for a given pathogen
*           is composed of chlorine/chloramine, ozone, and chlorine dioxide
*           components.
*
*  ct() is called by basn_dbp() and filt_dbp()
*/

{

	//FILE *fptr;

	/* Inputs: */
	double pH;			 /* Disinfection pH                      (-) */
	double DegC;		 /* Disinfection temperature         (Deg C) */
	double free_cl2;	 /* Free chlorine residual     (mg/L) as Cl2 */
	double comb_cl2;	 /* Combined chlorine residual (mg/L) as Cl2 */
	double o3_res;		 /* Ozone Residual                    (mg/L) */
	double clo2_res;	 /* Chlorine dioxide residual         (mg/L) */
	double log_required; /* These  three are ultimately the log-inactivations required */
	double log_required_v;
	double log_required_c;

	/* Outputs: */
	double CT_Achieved;
	double CT_Required;
	double CT_Required_v;

	/* Internal: */
	double o3_CT_Achieved;
	double o3_CT_Required;
	double o3_CT_Required_v;
	double o3_CT_Required_c;
	double clo2_CT_Achieved;
	double clo2_CT_Required;
	double clo2_CT_Required_v;
	double clo2_CT_Required_c;

	double CT;
	double CT_v;

	/* Internal variables used in CT integration*/
	int n, i;					  /* Number of time increments in CT integration, counter var.*/
	double o3_conc[11], time[11]; /* Arrays storing calculated ozone residuals and times      */
	double start_time;			  /* Time since last ozone application at unit influent       */
	double decay_period;		  /* Mean residence time of unit over which ozone decay occurs*/
	double delta_t;				  /* Increment of time used in CT integration                 */
	double inf_o3;				  /* Influent ozone residual for unit process (mg/L)          */

	struct Effluent *eff;

	/* Get input parameters from UnitProcess data structure */
	eff = &unit->eff;
	pH = eff->pH;
	DegC = eff->DegK - 273.15;
	free_cl2 = eff->FreeCl2 * MW_Cl2; /* mg/L as Cl2 */
	comb_cl2 = eff->NH2Cl * MW_Cl2;   /* mg/L as Cl2 */
	o3_res = eff->o3_res;
	clo2_res = eff->clo2_res;

	log_required = eff->log_required_g;
	log_required_v = eff->log_required_v;
	log_required_c = eff->log_required_c;

	if (free_cl2 > 0.0 ||
		comb_cl2 > 0.0 || /* Self-Protect */
		(o3_res >= 0.05 && unit->eff.o3_minutes <= 120.0) ||
		(clo2_res >= 0.05 && unit->eff.clo2_minutes <= 120.0))
	{ /* We should only calculate CT values here if there are disinfectant residuals! */

		//  if( swflag == TRUE || gw_virus_flag == TRUE )   /* Surface Water System or GW requiring virus inactivation */
		//    {

		/******************* Determine CHLORINE CT_Achieved and CT_Required ***********************/

		if (free_cl2 > 0.0)
		{ /********************************FREE CHLORINE****************************/

			CT_Achieved = free_cl2 * t_ten;
			eff->ct_cl2 += CT_Achieved; //Updates effluent data structure cum. tracking variable

			if (log_required > 0.0)
			{ /***** Giardia - Free Chlorine **********/

				/* Estimate Free Chlorine CT Required for Giardia in SWTR guidance manual. */
				/* These two equations are regressions developed by Smith et al. (1995) */
				/* of SWTRGM table values */
				if (pH < 6.0)
					pH = 6.0;
				if (pH > 9.0)
					pH = 9.0;
				if (free_cl2 < 0.4)
					free_cl2 = 0.4;

				if (DegC < 12.5)
				{
					if (DegC < 1.0)
						DegC = 1.0;
					CT_Required = log_required / 2.832 * (12.006 + exp(2.46 - 0.073 * DegC + 0.125 * free_cl2 + 0.389 * pH));
				}
				else /*(DegC >= 12.5)*/
				{
					CT_Required = log_required / 2.770 * (-2.261 + exp(2.69 - 0.065 * DegC + 0.111 * free_cl2 + 0.361 * pH));
				}

				// Update cumulative giardia inactivation achieved (based on assumption of linearity
				// with CT_Achieved)
				if (CT_Required > 0.0)
					eff->log_inact_ach_g += log_required * CT_Achieved / CT_Required;

			} /*End "if(log_required > 0.0)"*/

			if (log_required_v > 0.0)
			{ /***** Virus - Free Chlorine **********/

				/* Estimate 2 log virus CT with chlorine in SWTR guidance manual. */
				if (pH <= 9.0)
				{
					if (DegC < 5.0)
						CT_v = 6.0;
					else if (DegC >= 5.0 && DegC < 10.0)
						CT_v = 4.0;
					else if (DegC >= 10.0 && DegC < 15.0)
						CT_v = 3.0;
					else if (DegC >= 15.0 && DegC < 20.0)
						CT_v = 2.0;
					else /*  DegC>=20.0 */
						CT_v = 1.0;
				}
				else if (pH > 9.0 && pH < 10.0)
				{ /* Use linear interpolation in this range */
					if (DegC < 5.0)
						CT_v = (pH - 9.0) * (45.0 - 6.0) + 6.0;
					else if (DegC >= 5.0 && DegC < 10.0)
						CT_v = (pH - 9.0) * (30.0 - 4.0) + 4.0;
					else if (DegC >= 10.0 && DegC < 15.0)
						CT_v = (pH - 9.0) * (22.0 - 3.0) + 3.0;
					else if (DegC >= 15.0 && DegC < 20.0)
						CT_v = (pH - 9.0) * (15.0 - 2.0) + 2.0;
					else if (DegC >= 20.0 && DegC < 25.0)
						CT_v = (pH - 9.0) * (11.0 - 1.0) + 1.0;
					else /*  DegC>=25.0 */
						CT_v = (pH - 9.0) * (7.0 - 1.0) + 1.0;
				}
				else /*pH >= 10.0*/
				{
					if (DegC < 5.0)
						CT_v = 45.0;
					else if (DegC >= 5.0 && DegC < 10.0)
						CT_v = 30.0;
					else if (DegC >= 10.0 && DegC < 15.0)
						CT_v = 22.0;
					else if (DegC >= 15.0 && DegC < 20.0)
						CT_v = 15.0;
					else if (DegC >= 20.0 && DegC < 25.0)
						CT_v = 11.0;
					else /*  DegC>=25.0 */
						CT_v = 7.0;
				}

				CT_Required_v = CT_v * (log_required_v / 2.0);

				// Update cumulative virus inactivation achieved
				if (CT_Required_v > 0.0)
					eff->log_inact_ach_v += log_required_v * CT_Achieved / CT_Required_v;

			} /*End "if(log_required_v > 0.0)"*/

		} /* End of "if(free_cl2 > 0.0)..."*/

		else if (comb_cl2 > 0.0)
		{ /******************************** CHLORAMINES *************************************/

			CT_Achieved = comb_cl2 * t_ten;

			eff->ct_nh2cl += CT_Achieved; //Updates effluent data structure cum. tracking variable

			if (log_required > 0.0)
			{ /***** Giardia - Chloramines **********/

				/* Calculate Giardia CT values required for 3 log inactivation */
				if (DegC < 5.0)
					CT = (5.0 - DegC) / 4.0 * (635 - 365) + 365;
				else if (DegC >= 5.0 && DegC < 10.0)
					CT = (10.0 - DegC) / 5.0 * (365 - 310) + 310;
				else if (DegC >= 10.0 && DegC < 15.0)
					CT = (15.0 - DegC) / 5.0 * (310 - 250) + 250;
				else if (DegC >= 15.0 && DegC < 20.0)
					CT = (20.0 - DegC) / 5.0 * (250 - 185) + 185;
				else if (DegC >= 20.0 && DegC < 25.0)
					CT = (25.0 - DegC) / 5.0 * (185 - 125) + 125;
				else /*  DegC>=25.0 */
					CT = 125;

				CT_Required = CT * (log_required / 0.5);

				// Update cumulative giardia inactivation achieved
				if (CT_Required > 0.0)
					eff->log_inact_ach_g += log_required * CT_Achieved / CT_Required;

			} /* End "if(log_required >0)"*/

			if (log_required_v > 0.0)
			{ /***** Viruses - Chloramines **********/

				/* Calculate Virus CT Values with Chloramine using linear interpolation (wrto 
	        temperature) of SWTR tabulated values of inactivation */
				if (DegC < 5.0)
					CT_v = (5.0 - DegC) / 4.0 * (1243 - 365) + 857;
				else if (DegC >= 5.0 && DegC < 10.0)
					CT_v = (10.0 - DegC) / 5.0 * (857 - 643) + 643;
				else if (DegC >= 10.0 && DegC < 15.0)
					CT_v = (15.0 - DegC) / 5.0 * (643 - 428) + 428;
				else if (DegC >= 15.0 && DegC < 20.0)
					CT_v = (20.0 - DegC) / 5.0 * (428 - 321) + 321;
				else if (DegC >= 20.0 && DegC < 25.0)
					CT_v = (25.0 - DegC) / 5.0 * (321 - 214) + 214;
				else /*  DegC>=25.0 */
					CT_v = 214;

				if (log_required_v > 2.0 && log_required_v <= 3.0)
				{
					CT_Required_v = CT_v * (log_required_v / 2.0) * 1.1065;
				}
				else if (log_required_v > 3.0)
				{
					CT_Required_v = CT_v * (log_required_v / 2.0) * 1.1600;
				}
				else /*log_required <=2 --> Use CT_v developed for 2 log inactivation */
				{
					CT_Required_v = CT_v;
				}

				// Update cumulative virus inactivation achieved
				if (CT_Required_v > 0.0)
					eff->log_inact_ach_v += log_required_v * CT_Achieved / CT_Required_v;

			} /* End "if(log_required_v >0)"*/
		}	 /* End of "if(comb_cl2 > 0.0)..."*/

		else
		{ /****** No free or combined chlorine ******/
			CT_Achieved = 0.0;
			CT_Required = 1.0;   /* Avoid a divide by zero. */
			CT_Required_v = 1.0; /* Avoid a divide by zero. */
		}						 /* End of free and combined chlorine if-else statement*/

		/************************ Determine OZONE CT_Achieved and CT_Required ************************/

		if (unit->type == O3_CONTACTOR && o3_res >= 0.1 &&
			unit->eff.first_o3_chamber == FALSE)
		{ /* Ozone CT only to be calculated in O3_CONTACTOR (but not the 1st
	      one) where there's a measurable residual (assumed to be 0.05 mg/L) */

			//Ozone CT Achieved is based on AVERAGE ozone residual in the reactor, not
			//simply the reactor effluent concentration

			/* First, set some basic integration parameters */
			n = 10.0;
			decay_period = unit->eff.o3_minutes - PrevUnitProcess(unit)->eff.o3_minutes;
			start_time = PrevUnitProcess(unit)->eff.o3_minutes;
			inf_o3 = PrevUnitProcess(unit)->eff.o3_res; //unit process influent ozone conc.

			/* Set integration increment size*/
			delta_t = decay_period / n;

			/* Perform integration */
			o3_CT_Achieved = 0.0;
			for (i = 0; i <= n; i++)
			{
				// First, calculate time since O3 addition at this point
				time[i] = start_time + (i * delta_t);

				//Next, calculate ozone concentration at this point and
				//CT for the increment
				if (i == 0)
					o3_conc[i] = inf_o3;
				else
				{
					o3_conc[i] = inf_o3 - unit->eff.K_o3decay *
											  (pow(time[i], 0.068) - pow(start_time, 0.068));

					if (o3_conc[i] > o3_conc[i - 1])
						o3_conc[i] = o3_conc[i - 1];

					o3_CT_Achieved += (o3_conc[i] + o3_conc[i - 1]) / 2 * delta_t;
				}

			} /* End "for(i=0..." */

			/* Adjust so CT is based on T10, not avg. res. time */
			o3_CT_Achieved *= (t_ten / decay_period);

			eff->ct_o3 += o3_CT_Achieved; //Updates effluent data structure cum. tracking variable

			if (log_required > 0.0)
			{ /********* Giardia - Ozone **********/

				/* Calculate CT Required for Giardia with Ozone by starting with 0.5-log numbers
	        and interpolating between temperatures based on SWTR Guidance Manual numbers */

				if (DegC < 5.0)
					CT = (5.0 - DegC) / 4.0 * (0.48 - 0.32) + 0.32;
				else if (DegC >= 5.0 && DegC < 10.0)
					CT = (10.0 - DegC) / 5.0 * (0.32 - 0.23) + 0.23;
				else if (DegC >= 10.0 && DegC < 15.0)
					CT = (15.0 - DegC) / 5.0 * (0.23 - 0.16) + 0.16;
				else if (DegC >= 15.0 && DegC < 20.0)
					CT = (20.0 - DegC) / 5.0 * (0.16 - 0.12) + 0.12;
				else if (DegC >= 20.0 && DegC < 25.0)
					CT = (25.0 - DegC) / 5.0 * (0.12 - 0.08) + 0.08;
				else /*  DegC>=25.0 */
					CT = 0.08;

				o3_CT_Required = CT * (log_required / 0.5);

				// Update cumulative Giardia inactivation achieved
				if (o3_CT_Required > 0.0)
					eff->log_inact_ach_g += log_required * o3_CT_Achieved / o3_CT_Required;

			} /* End "if(log_required > 0.0)"*/

			if (log_required_c > 0.0)
			{ /******** Crypto. - Ozone **********/

				/* Calculate CT Required for Crypto with Ozone based on LT2 Guidance Manuals'
	        CT Table; this table can be closely approxmiated with a regression eqn. */

				if (DegC > 0.5 && DegC <= 25.0)
					CT = log_required_c * 25.16 * exp(-0.0929 * DegC); //R2=0.99  (WJS, 7/2005)
				else if (DegC <= 0.5)
					CT = log_required_c * 25.16 * exp(-0.0929 * 0.5);
				else /* (DegC > 25.0) */
					CT = log_required_c * 25.16 * exp(-0.0929 * 25.0);
				o3_CT_Required_c = CT;

				/*Update Crypto log inactivation achieved based on linearity between CT ratio
	       and log ratio */
				if (o3_CT_Required_c > 0.0)
					eff->log_inact_ach_c += log_required_c * o3_CT_Achieved / o3_CT_Required_c;

			} /* End "if(log_required_c > 0.0)" */

			if (log_required_v > 0.0)
			{ /********* Viruses - Ozone **********/

				/* Calculate Virus CT Values with Ozone using linear interpolation (wrto
	        temperature) of SWTR tabulated values starting with 2-log numbers */
				if (DegC < 5.0)
					CT_v = (5.0 - DegC) / 4.0 * (0.90 - 0.60) + 0.60;
				else if (DegC >= 5.0 && DegC < 10.0)
					CT_v = (10.0 - DegC) / 5.0 * (0.60 - 0.50) + 0.50;
				else if (DegC >= 10.0 && DegC < 15.0)
					CT_v = (15.0 - DegC) / 5.0 * (0.50 - 0.30) + 0.30;
				else if (DegC >= 15.0 && DegC < 20.0)
					CT_v = (20.0 - DegC) / 5.0 * (0.30 - 0.25) + 0.25;
				else if (DegC >= 20.0 && DegC < 25.0)
					CT_v = (25.0 - DegC) / 5.0 * (0.25 - 0.15) + 0.15;
				else /*  DegC>=25.0 */
					CT_v = 0.15;

				o3_CT_Required_v = CT_v * (log_required_v / 2.0);

				/*Update virus log inactivation achieved*/
				if (o3_CT_Required_v > 0.0)
					eff->log_inact_ach_v += log_required_v * o3_CT_Achieved / o3_CT_Required_v;
			}

			/*
	  fprintf(fptr,"t_ten/decay_period         %f\n",t_ten/decay_period);
	  fprintf(fptr,"Corrected CT Achieved      %f\n",o3_CT_Achieved);
	  fprintf(fptr,"o3_CT_Required             %f\n",o3_CT_Required);   */

		} /* End 'if(o3_res > 0.05...)...' */
		else
		{ /****** No Ozone *******/
			o3_CT_Achieved = 0.0;
			o3_CT_Required = 1;
			o3_CT_Required_c = 1; /* Avoid zero division */
			o3_CT_Required_v = 1;
		}

		/****************** Determine CHLORINE DIOXIDE CT_Achieved and CT_Required ******************/

		if (clo2_res >= 0.10)
		{ /* Assumes that any process with an effluent clo2_res < 0.10 mg/L would
	      not be counted as having a measurable residual for CT calculation. */

			clo2_CT_Achieved = unit->eff.clo2_avg * t_ten;

			eff->ct_clo2 += clo2_CT_Achieved; //Updates effluent data structure cum. tracking variable

			if (log_required > 0.0)
			{ /******** Giardia - ClO2 **********/

				/* Calculate CT Required for Giardia with ClO2 by starting with 0.5-log numbers
             and interpolating between temperatures based on SWTR Guidance Manual numbers */
				if (DegC < 5.0)
					CT = (5.0 - DegC) / 4.0 * (10.0 - 4.3) + 4.3;
				else if (DegC >= 5.0 && DegC < 10.0)
					CT = (10.0 - DegC) / 5.0 * (4.3 - 4.0) + 4.0;
				else if (DegC >= 10.0 && DegC < 15.0)
					CT = (15.0 - DegC) / 5.0 * (4.0 - 3.2) + 3.2;
				else if (DegC >= 15.0 && DegC < 20.0)
					CT = (20.0 - DegC) / 5.0 * (3.2 - 2.5) + 2.5;
				else if (DegC >= 20.0 && DegC < 25.0)
					CT = (25.0 - DegC) / 5.0 * (2.5 - 2.0) + 2.0;
				else /*  DegC>=25.0 */
					CT = 2.0;

				clo2_CT_Required = CT * (log_required / 0.5);

				/*Update Giardia log inactivation achieved*/
				if (clo2_CT_Required > 0.0)
					eff->log_inact_ach_g += log_required * clo2_CT_Achieved / clo2_CT_Required;

			} /*End "if(log_required > 0.0)"*/

			if (log_required_v > 0.0)
			{ /******** Viruses - ClO2 **********/

				/* Calculate Virus CT Values with ClO2 using linear interpolation (wrto
             temperature) of SWTR tabulated values */
				if (log_required_v <= 2.0)
				{
					if (DegC < 5.0)
						CT_v = (5.0 - DegC) / 4.0 * (0.90 - 0.60) + 0.60;
					else if (DegC >= 5.0 && DegC < 10.0)
						CT_v = (10.0 - DegC) / 5.0 * (0.60 - 0.50) + 0.50;
					else if (DegC >= 10.0 && DegC < 15.0)
						CT_v = (15.0 - DegC) / 5.0 * (0.50 - 0.30) + 0.30;
					else if (DegC >= 15.0 && DegC < 20.0)
						CT_v = (20.0 - DegC) / 5.0 * (0.30 - 0.25) + 0.25;
					else if (DegC >= 20.0 && DegC < 25.0)
						CT_v = (25.0 - DegC) / 5.0 * (0.25 - 0.15) + 0.15;
					else /*  DegC>=25.0 */
						CT_v = 0.15;
				}
				else if (log_required_v <= 3.0 && log_required_v > 2.0)
				{
					if (DegC < 5.0)
						CT_v = (5.0 - DegC) / 4.0 * (25.6 - 17.1) + 17.1;
					else if (DegC >= 5.0 && DegC < 10.0)
						CT_v = (10.0 - DegC) / 5.0 * (17.1 - 12.8) + 12.8;
					else if (DegC >= 10.0 && DegC < 15.0)
						CT_v = (15.0 - DegC) / 5.0 * (12.8 - 8.6) + 8.6;
					else if (DegC >= 15.0 && DegC < 20.0)
						CT_v = (20.0 - DegC) / 5.0 * (8.6 - 6.4) + 6.4;
					else if (DegC >= 20.0 && DegC < 25.0)
						CT_v = (25.0 - DegC) / 5.0 * (6.4 - 4.3) + 4.3;
					else /*  DegC>=25.0 */
						CT_v = 4.3;
				}
				else /* log_required_v > 3.0 */
				{
					if (DegC < 5.0)
						CT_v = (5.0 - DegC) / 4.0 * (50.1 - 33.4) + 33.4;
					else if (DegC >= 5.0 && DegC < 10.0)
						CT_v = (10.0 - DegC) / 5.0 * (33.4 - 25.1) + 25.1;
					else if (DegC >= 10.0 && DegC < 15.0)
						CT_v = (15.0 - DegC) / 5.0 * (25.1 - 16.7) + 16.7;
					else if (DegC >= 15.0 && DegC < 20.0)
						CT_v = (20.0 - DegC) / 5.0 * (16.7 - 12.5) + 12.5;
					else if (DegC >= 20.0 && DegC < 25.0)
						CT_v = (25.0 - DegC) / 5.0 * (12.5 - 8.4) + 8.4;
					else /*  DegC>=25.0 */
						CT_v = 8.4;
				}
				clo2_CT_Required_v = CT_v;

				/*Update virus log inactivation achieved*/
				if (clo2_CT_Required_v > 0.0)
					eff->log_inact_ach_v += log_required_v * clo2_CT_Achieved / clo2_CT_Required_v;

			} /*End "if(log_required_v > 0.0)"*/

			if (log_required_c > 0.0)
			{ /******** Crypto. - ClO2 **********/

				/* Calculate CT Required for Crypto. with ClO2 based on LT2 Guidance Manuals table */

				if (DegC > 0.5 && DegC <= 25.0)
					CT = log_required_c * 664.34 * exp(-0.0873 * DegC);
				else if (DegC <= 0.5)
					CT = log_required_c * 664.34 * exp(-0.0873 * 0.5);
				else /*(DegC > 25.0)*/
					CT = log_required_c * 664.34 * exp(-0.0873 * 25.0);
				clo2_CT_Required_c = CT;

				/*Update Crypto. log inactivation achieved*/
				if (clo2_CT_Required_c > 0.0)
					eff->log_inact_ach_c += log_required_c * clo2_CT_Achieved / clo2_CT_Required_c;
			} /* End "if(log_required_c > 0.0)"*/

		} /* End of "if(clo2_res >= 0.10)..." */
		else
		{ /********* No ClO2 ********/
			clo2_CT_Achieved = 0.0;
			clo2_CT_Required = 1;
			clo2_CT_Required_v = 1;
		}

		/*************************** End of CT and Inactivation Calculations ************************/

		/* Update Inactivation Ratios in Unit Process Data Structure */
		/* Giardia */
		if (log_required > 0.0)
		{
			eff->ct_ratio += (CT_Achieved / CT_Required + o3_CT_Achieved / o3_CT_Required + clo2_CT_Achieved / clo2_CT_Required);
		}

		/* Virus   */
		if (log_required_v > 0.0)
		{
			eff->ct_ratio_v += (CT_Achieved / CT_Required_v + o3_CT_Achieved / o3_CT_Required_v + clo2_CT_Achieved / clo2_CT_Required_v);
		}

		/* Cryptosporidium   */
		if (log_required_c > 0.0)
		{
			eff->ct_ratio_c += o3_CT_Achieved / o3_CT_Required_c + clo2_CT_Achieved / clo2_CT_Required_c;
		}

		//    }	 /* End "if(swflag == TRUE)...." */

		//else
		/* groundwater  -  Do nothing, just remember that these variables = 999999.0
		eff->log_required_g;
  		eff->log_required_v;
  		eff->log_required_c; */

	} /* End " if( any residuals > 0)..." */

} /* End subroutine */
