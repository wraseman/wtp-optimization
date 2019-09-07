/* ClO2.c */
#include "wtp.h"

void clo2_decay(struct UnitProcess *unit, double cfstr_time)
/*
* Purpose:
*   Estimate ClO2 Residuals (for output)
*
* Input: "cfstr_time" is the residence time of the cfstr in which this algorithm is
*        calculating ClO2 decay (minutes)
*  
* Documention and Code by WJS, 5/2001
*/
{
  double clo2_minutes; /* Time since the last clo2 dose */

  /*Unit Process Influent Parameters*/
  double degC; /* Temperature (Degrees C)*/
  double pH;
  double uv;   /* UV-254 (1/cm)*/
  double toc;  /* TOC (mg/L) */
  double dose; /* ClO2 dose (mg/L)*/
  double k1;   /* First order rate constant */

  /*Internal Variables*/
  double term1, term2, term3, term4;
  double clo2_in;  /* ClO2 residual at influent */
  double clo2_res; /* intermediate */

  /*Outputs*/
  double clo2_out; /* clo2 residual at effluent*/

  register struct Effluent *eff;

  eff = &unit->eff;

  /*Get inputs from Effluent data structure*/
  degC = eff->DegK - 273.15;
  pH = eff->pH;
  uv = eff->UV;
  toc = eff->TOC;
  dose = eff->clo2dose;
  clo2_in = eff->clo2_res;
  clo2_minutes = eff->clo2_minutes;

  //Initialize this way in case clo2_in = 0.0
  clo2_out = clo2_in;

  /*These limits are established to keep predicted ClO2 decay from being
    unrealistically slow.  They are based on the model development
    database parameter value ranges. */
  if (degC < 4.3)
    degC = 4.3;
  if (uv < 0.027)
    uv = 0.027;
  if (toc < 1.3)
    toc = 1.3;

  if (clo2_minutes == 0.0)
  { //Take care of initial decay first, but only if this is right after a dose

    //Determine initial residual based upon raw or coagulated water equations
    if (swflag == TRUE && (rwdbpflag == TRUE || owdbpflag == TRUE))
    { // use raw water equation

      term1 = pow(dose, 1.802);
      term2 = pow(degC, -0.0475);
      term3 = pow(pH, 1.47);
      term4 = pow((toc * uv), -0.284);
      clo2_res = 0.0157 * term1 * term2 * term3 * term4;
    }
    else
    { //use coagulated water equation

      term1 = pow(dose, 1.415);
      term2 = pow(degC, -0.0395);
      term3 = pow(pH, 1.85);
      term4 = pow((toc * uv), -0.182);
      clo2_res = 0.0124 * term1 * term2 * term3 * term4;
    }

    //no increasing of ClO2 concentration here
    if (clo2_res > clo2_in)
      clo2_res = clo2_in;

  } //End"if(clo2_minutes = 0.0)"
  else
  {
    clo2_res = clo2_in;
  }

  //--------------------------------------------------------------------------------------

  if (clo2_res > 0.0 && dose > 0.0)
  { //Calculate decay with 1st order constant based upon raw or coagulated water equations
    //and fact that we are talking about a CFSTR where Ceff = Cin / (1-kt)

    if (swflag == TRUE && (rwdbpflag == TRUE || owdbpflag == TRUE))
    { // use raw water equation

      term1 = pow(dose, -0.673);
      term2 = pow(degC, 0.398);
      term3 = pow((toc * uv), 0.441);
      k1 = -0.014 * term1 * term2 * term3;
    }
    else
    { //use coagulated water equation

      term1 = pow(dose, -0.544);
      term2 = pow(degC, 0.237);
      term3 = pow((toc * uv), 0.764);
      k1 = -0.0282 * term1 * term2 * term3;
    }

    clo2_out = clo2_res / (1.0 - k1 * cfstr_time);

    //No residuals less than 0.1 mg/L
    if (clo2_out < 0.1)
      clo2_out = 0.0;

  } //End"if(clo2_res > 0.0 && dose > 0.0)"

  /*Update Unit Process Effluent data structure*/
  eff->clo2_res = clo2_out;
}