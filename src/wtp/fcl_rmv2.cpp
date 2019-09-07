/* Fcl_Rmv2.c
*/
#include "wtp.h"
#define MW_Feric 270L /* (mg/mmol) */

void fecl_rmv(struct UnitProcess *unit)
/*
*  Purpose: Estimate TOC removal (based on Marc Edwards model) and UV removal,
*           due to ferric coagulation
*
* Notes:
*  1. Is there a real differance between 'raw_toc' and 'inf_toc'?
*     Same for raw_uv.
*  2. AlumDose>0 is used as a flag to indicate that this routine should
*     do its thing. AlumDose is set to zero at the conclusion of this 
*     routine.
*  3. TOC and UV removal equations are A-13 and A-14 in 1.21 version 
*     of WTP manual.  The equations predict finished water quality (plant
*     effluent) - NOT the effluent from the settling tank.
*  4. The inputs are not range checked within this function.  Table A-11
*     lists the ranges of the inputs that were used in developing the
*     equations.
*
* Documentation and code by WJS, December 1998.
*/
{
  /* Inputs to model: */
  double pH;         /* In influent to rapid mix                         */
  double FerricDose; /* Cumulative             (mg/L) as Al2(SO4)3-14H2O */

  /* Empirical Model Coefficients: */
  double K1 = -0.028; /* Used to calculate FNSDOC based on SUVA         */
  double K2 = 0.23;   /* Used to calculate FNSDOC based on SUVA         */
  double a;           /* Max. DOC sorption per mM of coagulant [mgDOC/(L-mM coag.)] */
  double x1 = 280.0;  /* Used to calculate parameter "a"                */
  double x2 = -73.9;  /* Used to calculate parameter "a"                */
  double x3 = 4.96;   /* Used to calculate parameter "a"                */
  double b = 0.068;   /* Sorption constant for sorbable DOC [L/mgDOC    */

  /* Variables for model computations */
  double suva;         /* SUVA of uncoagulated water [L/mg-m]            */
  double doc_sorb_inf; /* Sorbable DOC in water before coagulation [mg/L]*/
  double doc_non_sorb; /* Non-sorbable DOC                         [mg/L]*/
  double fns_doc_inf;  /* Fraction of uncoagulated DOC that is non-sorbable */
  double term1, term2, term1_out, term3;
  double doc_sorb_eff1;
  double doc_sorb_eff2;
  double delta_uv;
  double delta_uv_out;

  /* The following two parameters are used to range check the output of 
  *  the model. ie the effluent should not be greater than the influent.
  */
  double inf_toc;    /* Influent TOC to RM                        (mg/L) */
  double inf_uv;     /* Influent UV to RM (for calcs.)            (1/cm) */
  double inf_uv_out; /* Influent UV to RM (for outputting)        (1/cm) */

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

    /* turn coagulant dosage to milliequivalents/Liter here */
    FerricDose = eff->FericDose / 270.0;

    /*For TOC removal model purposes, set 8.0 as max pH considered*/
    if (pH > 8.0)
      pH = 8.0;

    if (unit->type == RAPID_MIX && FerricDose > 0.0 && inf_toc > 0.0)
    {
      suva = 100 * inf_uv / inf_toc;
      fns_doc_inf = K1 * suva + K2;
      if (fns_doc_inf < 0.0)
        fns_doc_inf = 0.0; /* Self-protection */
      if (fns_doc_inf > 1.0)
        fns_doc_inf = 1.0; /* Self-protection */
      doc_non_sorb = fns_doc_inf * inf_toc;
      doc_sorb_inf = inf_toc - doc_non_sorb;
      a = x3 * pow(pH, 3.0) + x2 * pow(pH, 2.0) + x1 * pH;
      term1 = pow(((doc_sorb_inf * b) - (b * a * FerricDose) - 1), 2.0) + (4.0 * b * doc_sorb_inf);

      term2 = 1 + (a * b * FerricDose) - (doc_sorb_inf * b);

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

      /* UV-254 removal by coagulation. */
      /* Old Equation for ferric*/
      /*eff->UV = pow(10.0,( 0.228
                              +1.025*log10(inf_uv)
                              -0.033*log10(FerricDose)
                              -2.367 / pH)
			); */
      /*New equation*/
      if (pH < 3.0)
        pH = 3.0; /*self-protection*/
      term1 = pow(inf_uv, 1.0894);
      term1_out = pow(inf_uv_out, 1.0894);
      term2 = pow(FerricDose, 0.305);
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

      /* Update UnitProcess data structure. */
      eff->toc_sludge += inf_toc - eff->TOC;
      eff->iron_to_sludge += eff->FericDose;
      eff->FericDose = 0.0;

    } /* End if(unit->type == RAPID_MIX ...)*/

  } /* End if(unit != NULL)...*/
}
