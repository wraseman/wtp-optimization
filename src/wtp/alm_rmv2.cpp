/* alm_rmv2.c */
#include "wtp.h"

void alum_rmv(struct UnitProcess *unit)
/*
*  Purpose: Estimate TOC removal (based on Marc Edwards model), UV removal,
*           and sludge production due to alum coagulation
*
* Notes:
*  1. Is there a real difference between 'raw_toc' and 'inf_toc'?
*     Same for raw_uv.
*  2. AlumDose>0 is used as a flag to indicate that this routine should
*     do its thing. AlumDose is set to zero at the conclusion of this 
*     routine.
*  3. TOC and UV removal equations are A-13 and A-14 in 1.21 version 
*     of WTP manual.  The equations predict finished water quality (plant
*     effluent) - NOT the effluent from the settling tank.
*  4. I found no reference for the sludge production equation.
*  5. The inputs are not range checked within this fuction.  Table A-11
*     lists the ranges of the inputs that were used in developing the
*     equations.
*
* Documentation and code by WJS, December 1998.
* Model parameters corrected to 'General' on Table 3 of Edwards (1997) by WJR, August 2017.
*/
{
  /* Inputs to model: */
  double pH;       /* In influent to rapid mix                         */
  double AlumDose; /* Cumulative             (meq/L)                   */

  /* Empirical Model Coefficients: */
  double K1 = -0.054; /* Used to calculate FNSDOC based on SUVA         */
  double K2 = 0.54;   /* Used to calculate FNSDOC based on SUVA         */
  double a;           /* Max. DOC sorption per mM of coagulant [mgDOC/(L-mM coag.)] */
  double x1 = 383.0;  /* Used to calculate parameter "a"                */
  double x2 = -98.6;  /* Used to calculate parameter "a"                */
  double x3 = 6.42;   /* Used to calculate parameter "a"                */
  double b = 0.145;   /* Sorption constant for sorbable DOC [L/mgDOC]   */

  /* Variables for model computations */
  double suva;         /* SUVA of uncoagulated water [L/mg-m]            */
  double doc_sorb_inf; /* Sorbable DOC in water before coagulation [mg/L]*/
  double doc_non_sorb; /* Non-sorbable DOC                         [mg/L]*/
  double fns_doc_inf;  /* Fraction of uncoagulated DOC that is non-sorbable */
  double term1, term1_out, term2, term3;
  double doc_sorb_eff1;
  double doc_sorb_eff2;
  double delta_uv;
  double delta_uv_out;

  /* The following two parameters are used to range check the output of 
  *  the model. ie the effluent should not be greater than the influent.
  */
  double inf_toc;    /* Influent TOC to rapid mix (mg/L) */
  double inf_uv;     /* Influent UV to rapid mix (for calcs.) (1/cm) */
  double inf_uv_out; /* Influent UV to rapid mix (for output) (1/cm) */

  /* Internal: */
  register struct Effluent *eff;

  if (unit != NULL)
  {
    /* Get inputs from UnitProcess data structure */
    eff = &unit->eff;
    pH = eff->pH;
    inf_toc = eff->TOC;
    inf_uv = eff->UV;
    inf_uv_out = eff->UV_out;

    /* turn alum+ferric dosage to milliequivalents of Al3+/Liter here */
    AlumDose = (eff->AlumDose + eff->FericDose / 270.0 * 297.0) / 297.0;

    /*For TOC removal model purposes, set 8.0 as max pH considered*/
    if (pH > 8.0)
      pH = 8.0;

    if (unit->type == RAPID_MIX && AlumDose > 0.0 && inf_toc > 0.0)
    {
      suva = 100 * inf_uv / inf_toc;
      fns_doc_inf = K1 * suva + K2;
      doc_non_sorb = fns_doc_inf * inf_toc;
      doc_sorb_inf = inf_toc - doc_non_sorb;
      a = x3 * pow(pH, 3.0) + x2 * pow(pH, 2.0) + x1 * pH;

      term1 = pow(((doc_sorb_inf * b) - (b * a * AlumDose) - 1), 2.0) + (4.0 * b * doc_sorb_inf);

      term2 = 1 + (a * b * AlumDose) - (doc_sorb_inf * b);

      /* Calculate both solutions of quadratic equation */
      doc_sorb_eff1 = (term2 - pow(term1, 0.5)) / (-2.0 * b);
      doc_sorb_eff2 = (term2 + pow(term1, 0.5)) / (-2.0 * b);

      /* If one of the answers is negative, choose the other as the solution */
      if (doc_sorb_eff1 >= 0.0)
        eff->TOC = doc_sorb_eff1 + doc_non_sorb;
      if (doc_sorb_eff2 >= 0.0)
        eff->TOC = doc_sorb_eff2 + doc_non_sorb;

      /* If both answers are positive, choose the one that is less than the influent*/
      if (doc_sorb_eff1 >= 0.0 && doc_sorb_eff2 >= 0.0)
      {
        if (doc_sorb_eff2 < doc_sorb_inf)
          eff->TOC = doc_sorb_eff2 + doc_non_sorb;
        if (doc_sorb_eff1 < doc_sorb_inf)
          eff->TOC = doc_sorb_eff1 + doc_non_sorb;
      } /* If both were pos. and less than the influent, the first will be selected */

      if (eff->TOC > inf_toc)
        eff->TOC = inf_toc;
      if (eff->TOC < 0.0)
        eff->TOC = 0.0;

      /* Update Unit Process Data Structure*/
      eff->toc_sludge += inf_toc - eff->TOC;
      eff->alum_to_sludge += eff->AlumDose;
      eff->iron_to_sludge += eff->FericDose;
      eff->AlumDose = 0.0;
      eff->FericDose = 0.0;

      /* UV-254 removal by Coagulation.  */

      /*Old equation*/
      /* eff->UV = exp( -4.639
                        +0.8793 * log(inf_uv)
                        -0.1846 * log(AlumDose*297)
			+0.5639 * pH 
		      );  */
      /*New equation*/
      if (pH < 3.0)
        pH = 3.0; /*self-protection*/
      term1 = pow(inf_uv, 1.0894);
      term1_out = pow(inf_uv_out, 1.0894);
      term2 = pow(AlumDose, 0.305);
      term3 = pow(pH, -0.9513);
      delta_uv = 5.7154 * term1 * term2 * term3;
      delta_uv_out = 5.7154 * term1_out * term2 * term3;
      eff->UV -= delta_uv;
      eff->UV_out -= delta_uv_out;
      if (eff->UV < 0.0)
        eff->UV = 0.0;
      if (eff->UV_out < 0.0)
        eff->UV_out = 0.0;
      if (eff->UV > inf_uv)
        eff->UV = inf_uv;
      if (eff->UV_out > inf_uv_out)
        eff->UV_out = inf_uv_out;

    } /* End if(unit->type == RAPID_MIX && AlumDose > 0.0 && inf_toc > 0.0)..*/

  } /* End if(unit != NULL)...*/
}
