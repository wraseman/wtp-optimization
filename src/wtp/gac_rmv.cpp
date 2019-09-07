/* GAC_Rmv.c   September 7, 1993 */
#include "wtp.h"

void gac_rmv(struct UnitProcess *unit)
/*
*  Purpose:
*   This subroutine calculates the average TOC and UV-254
*   removals by granular activated carbon adsorption for both single
*   and multiple contactor systems.  The type of system is selected by
*   the user at the input screen and stored in the variable "config".
*   If the GAC unit is a parallel system of contactors, the average, steady-
*   state effluent TOC and UV values are calculated. If the GAC unit is composed
*   of a single contactor, the program calculates the running average effluent
*   TOC and UV or the maximum expected TOC and UV (based on the values at
*   the end of the column regen period.  The decision of which value to
*   calculate is based on user input and stored in the variable "toc_calc"
*
*  Notes:
*   1. Equations developed from Summers and Solarik 1999.
*   2. Coded and annoted by WJS, 2/99
*/
{
  //  FILE *fptr2;

  /* Inputs: */
  double ebct;              /* Empty bed contact time (minutes) */
  double regen;             /* Reeneration frequency     (days) */
  char config;              /* GAC contactor config., 'S' for single
						    'B' for blended */
  char toc_calc;            /* Eff. TOC calc. basis for single contactor system
  			     'M' for max., 'A' for average */
  double inf_toc;           /* Influent TOC              (mg/L) */
  double inf_uv;            /* Influent UV (for internal calcs.) (1/cm) */
  double inf_uv_out;        /* Influent UV (for internal calcs.) (1/cm) */
  double peak_factor = 1.0; /* Set equal to eff->Peak/eff->Flow if coldflag == TRUE */
  double inf_pH;

  /* Outputs: */
  double eff_toc;
  double eff_uv;     /*UV for internal calcs.     (1/cm)*/
  double eff_uv_out; /*UV for outputting purposes (1/cm)*/

  /* Internal: */
  double coeff_a, coeff_b, coeff_d;

  struct Gac *gac;      /* Design and operating data */
  struct Effluent *eff; /* Water quality */

  if (unit == NULL || unit->data.ptr == NULL || unit->type != GAC)
    return;

  /* Get inputs from UnitProcess data structure */
  gac = unit->data.gac;
  eff = &unit->eff;
  if (coldflag == TRUE)
    peak_factor = eff->Peak / eff->influent->data.influent->avg_flow;
  ebct = gac->ebct;
  regen = gac->regen;
  config = gac->config;
  toc_calc = gac->toc_calc;
  inf_toc = eff->TOC;
  inf_uv = eff->UV;
  inf_uv_out = eff->UV_out;
  inf_pH = eff->pH;

  /* Adjust EBCT for peak flow scenario based on whether or not it is a single contactor
     or a system of parallel contactors */
  if (coldflag == TRUE && config == 'S')
    ebct /= peak_factor;
  else
  {
  } /* ebct stays same because we are at average flow or we are at peak flow, but
	     we are assuming that parallel contactors are brought on- and off-line to
	     maintain the EBCT at the user-specified value as flow changes*/

  /* Compute coefficients */
  coeff_a = 0.682 * inf_toc;
  coeff_b = 0.167 * pow(inf_pH, 2.0) - 0.808 * inf_pH + 19.086;
  coeff_d = inf_toc * (inf_pH * (-0.0000058 * pow(ebct, 2.0) + 0.000111 * ebct + 0.00125) + 0.0001444 * pow(ebct, 2.0) - 0.005486 * ebct + 0.06005);

  if (config == 'S' && toc_calc == 'M')
  /* Compute Eff. TOC based on breakthrough at end of regen. period
     this is a worst-case (Max. Eff. TOC scenario) */
  {
    eff_toc = coeff_a / (1 + coeff_b * exp(-coeff_d * regen));
  }
  else
  { /* We have either a blended effluent scenario or a time-averaged
	  single contactor scenario */
    eff_toc = coeff_a + coeff_a / coeff_d / regen * log(1.0 + coeff_b * exp(-coeff_d * regen)) - coeff_a / coeff_d / regen * log(1.0 + coeff_b);
  }

  /* Check for reasonableness - we are not creating TOC here */
  if (eff_toc > inf_toc)
    eff_toc = inf_toc;

  /* Calculate effluent UV based upon what processes preceded GAC
     updated from G. Solarik in May 2001 */
  if (pre_o3flag == TRUE || int_o3flag == TRUE)
  {
    eff_uv = 0.0114 * eff_toc - 0.0014;
    eff_uv_out = 0.0114 * eff_toc - 0.0014;
  }
  else
  {
    eff_uv = 0.0195 * eff_toc - 0.0077;
    eff_uv_out = 0.0195 * eff_toc - 0.0077;
  }

  /* Check for reasonableness - we are not creating UV here */
  if (eff_uv > inf_uv)
    eff_uv = inf_uv;
  if (eff_uv < 0.0)
    eff_uv = 0.0;
  if (eff_uv_out > inf_uv_out)
    eff_uv_out = inf_uv_out;
  if (eff_uv_out < 0.0)
    eff_uv_out = 0.0;

  /* Copy outputs to UnitProcess data structure */
  eff->TOC = eff_toc;
  eff->UV = eff_uv;
  eff->UV_out = eff_uv_out;

  /* Zero disinfect residuals in UnitProcess data structure */
  eff->FreeCl2 = 0.0;
  eff->NH2Cl = 0.0;
  eff->NHCl2 = 0.0;
  eff->cl2dose = 0.0;
  eff->o3_res = 0.0;
  eff->clo2_res = 0.0;
  eff->last_o3_inf = NULL;
  /* This line assures that ozonate is not called until ozone is
	    added again after this GAC unit */

  /*DEBUGGING CODE*/
  //      fptr2=fopen("debug.dat","a+");
  //      fprintf(fptr2,"Module: %d\n",unit->type);
  //      fprintf(fptr2,"coeff_a:              %f\n",coeff_a);
  //      fprintf(fptr2,"coeff_b:              %f\n",coeff_b);
  //      fprintf(fptr2,"coeff_d:              %f\n",coeff_d);
  //      fclose(fptr2);
  /*DEBUGGING CODE*/
}

/*************************************************************************************/
