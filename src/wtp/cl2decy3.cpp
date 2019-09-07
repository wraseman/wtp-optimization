/* Cl2Decay3.c */
#include "wtp.h"

void new_cl2decay(struct UnitProcess *unit, double cfstr_time)
/*
* Purpose:
*   Estimate chlorine and chloramine decay in a reactor.
*   The calling routine is responsible for determining the mean residence time.
*   The calling routine is also responsible for updating unit.hours
*
* Inputs:
*   cfstr_time
*     Residence time of one of the CFSTRs in series for which this routine calculates
*     chlorine or chloramine decay.
*   hours
*     hours is the cumulative time from the point of chlorine addition
*     to the INFLUENT of the CFSTR.  The cumulative time to the effluent
*     is (hours+cfstr_time).  After calling cl2decay2(), the calling routine
*     must add 'cfstr_time' to unit->hours.
*   cl2dose
*     Cumulative free chlorine added to process train to INFLUENT of unit
*     process.
*   cl2res
*     On input, cl2res is free chlorine residual (mg/L as Cl2) in INFLUENT of
*     reactor.
*   comb_cl2
*     On input, comb_cl2 is monochloramine residual (mg/L as Cl2) in
*     INFLUENT of reactor.
*   cl2cnt
*     Number of times chlorine added to process train as of INFLUENT of
*     unit process.
*   pH, toc, uv
*     Water quality parameters in INFLUENT of unit process.
*
* Outputs: (Inputs updated)
*   cl2res
*     On return, cl2res is free chlorine residual (mg/L as Cl2) in EFFLUENT
*     of reactor.
*   comb_cl2
*     On return, comb_cl2 is monochloramine residual (mg/L as Cl2) in
*     EFFLUENT of reactor.
*
* Method:
*   The decay equations are those developed at the University of Colorado.
*   These are the saturation-type decay models of the following form:
*                   dC/dt = -alpha2*C / (alpha1 + C)
*   This rate expression is used in the equation for concentration in a
*   CFSTR:
*                    Ceff = Cin + dC/dt * dt
*   Ceff can be solved based on this these and use of the solution to a
*   quadratic equation to yield:
*       Ceff = -(alpha1-Cin+alpha2*dt)/2 - sqrt((alpha1-Cin+alpha2*dt)^2 -4(-alpha1*Cin))/2
*
* Documentation and Code by WJS - 05/2001
*/
{
  //FILE *fptr2;

  /* Inputs: ******************************************************/
  double toc;     /*  "      "    "                     (mg/L) */
  double uv;      /*  "      "    "                     (1/cm) */
                  //  double  hours;    /* Cumulative hours to CFSTR influent. (hrs) */
  double cl2dose; /* Total chlorine                     (mg/l) */

  /* Inputs changed (Outputs) *************************************/
  double cl2res;   /* Free chlorine residual     (mg/L) as Cl2 */
  double comb_cl2; /* Combined chlorine residual (mg/L) as Cl2 */

  /* Internal *****************************************************/
  double alpha1, alpha2, k2;
  double a, b, c;
  double Co;
  //  double  Ctin;
  double Ctout;
  //  double  DeltaCt;
  //  double  t;
  //  double guessCt,prev_guessCt,dose_fract,calcCt;

  struct Effluent *eff;

  /* Update [FreeCl] and [NH2Cl] using breakpoint estimation. */
  breakpt(unit);

  /* Get inputs from UnitProcess structure */
  eff = &unit->eff;
  toc = eff->TOC;
  uv = eff->UV;
  //  hours    = eff->hours;
  cl2res = eff->FreeCl2 * MW_Cl2;
  cl2dose = eff->cl2dose;
  comb_cl2 = eff->NH2Cl * MW_Cl2;

  /*Self-protection Statements*/
  if (uv < 0.001)
    uv = 0.001;
  if (toc < 0.1)
    toc = 0.1;
  if (cl2dose > 50.0)
    cl2dose = 50.0;
  if (cl2dose <= 0.0)
    cl2dose = 0.09;

  Co = cl2dose;

  /* Estimate free chlorine decay *********************************************/
  if (cl2res > 0.0)
  {

    //Parameters for Free Cl2 decay
    if (swflag == TRUE && (rwdbpflag == TRUE || owdbpflag == TRUE))
    { // use raw water decay equations
      alpha1 = -0.817 * Co;
      k2 = -2.2808 * pow((Co / uv), -1.2971);
    }
    else
    { //use coagulated water decay equations for coag. water, gac-treated or mem-treated water, or groundwater
      alpha1 = -0.8408 * Co;
      k2 = -0.404 * pow((Co / uv), -0.918);
    }

    alpha2 = k2 * toc;

    // Estimate free chlorine at CFSTR effluent ---------------

    //Quadratic equation parameters
    a = 1.0;
    b = alpha1 - cl2res + alpha2 * cfstr_time;
    c = -alpha1 * cl2res;

    Ctout = (-b - pow((b * b - 4.0 * a * c), 0.5)) / (2.0 * a);

    //Check that answer makes sense
    if (Ctout > cl2res)
      Ctout = cl2res;
    if (Ctout < 0.0)
      Ctout = 0.0;

    //Calculate new chlorine residual
    cl2res = Ctout;

    /* Limit chlorine residual to 0.1 mg/L */
    if (cl2res < 0.1)
    {
      cl2res = 0.0;
    }

  } //end   if( cl2res > 0.0 )

  /*****************************************************************************
  * Estimate chloramine decay...  Section A.6.3 and A-6.4  */

  if (comb_cl2 > 0.0)
  {
    //Parameters for chloramine decay
    alpha1 = -0.990 * Co;
    alpha2 = -0.015 * uv;

    // Estimate combined chlorine at CFSTR effluent---------------

    //Quadratic equation parameters
    a = 1.0;
    b = alpha1 - comb_cl2 + alpha2 * cfstr_time;
    c = -alpha1 * comb_cl2;

    Ctout = (-b - pow((b * b - 4.0 * a * c), 0.5)) / (2.0 * a);

    //Check that answer makes sense
    if (Ctout > comb_cl2)
      Ctout = comb_cl2;
    if (Ctout < 0.0)
      Ctout = 0.0;

    //Calculate new combined chlorine residual
    comb_cl2 = Ctout;

    /* Limit chlorine residual to 0.1 mg/L */
    if (comb_cl2 < 0.1)
    {
      comb_cl2 = 0.0;
    }

  } //end   if( comb_cl2>0.0 )

  /**************************************************************/

  /* Copy outputs to UnitProcess data structure */
  eff->FreeCl2 = cl2res / MW_Cl2;
  eff->NH2Cl = comb_cl2 / MW_Cl2;
}
