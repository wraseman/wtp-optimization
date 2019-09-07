/* mdrw2dbp.c */
#include "wtp.h"

void modrw2dbp(struct UnitProcess *unit)
/*
* Purpose: Estimate trihalomethane and haloacetic acid concentration in current
*	   unit process based upon the following combination of models and
*	   input conditions:
*   		1) Raw water DBP model from Jinsik Sohn, Univ. of Col., 1998
*		2) Ratio of raw water and Post-RM models from Gaby Solarik, 1998
*		3) Water quality input parameters from the RM influent
*		4) Time since influent of RM
*
* Notes:
*  1. Called from choose_model() whenever modrw2dbpflag == TRUE and current unit
*     is not RM.
*  2. The calling routine must update unit.effluent.hours to include the
*     contact time of this unit process.
*  3. inf_hours and eff_hours are cumulative from the last point of
*     chlorine addition.
*
* Documentation and code by WJS, 10/98
*/
{
   /* Inputs to model (from RM influent): */
   double doc;    /* Really, its TOC      (mg/L) */
   double uv;     /* UV254		       (1/cm) */
   double doc_uv; /* Product of DOC and UV-254 (mg/L/cm) */

   /* Inputs to model (from current unit process): */
   double pH;         /*                         (-) */
   double DegC;       /* Temperature         (Deg C) */
   double doc_out;    /* TOC at effluent of unit process (mg/L) */
   double br;         /* Bromide              (mg/L) */
   double br_ugl;     /* Bromide              (ug/L) */
   double cl2dose;    /* Free chlorine (mg/L) as Cl2 */
   double inf_time;   /* Values for these depend on whether or not there */
   double eff_time;   /*   was chlorination before the RM */
   double inf_hours;  /* Cumulative time since last chlorination
                          to current unit process influent (hours) */
   double inf_FreeCl; /* Process Influent Free Chlorine      (Mole/L) */
   double inf_nh2cl;  /* Process Influent Monochloramine     (Mole/L) */

   /* Intermediates: */
   double doc_reduction;

   /* The factors are applied to the Raw Water Models to simulate
   effect of pre-chlorination */
   double thm_factor;
   double haa_factor;
   double tox_factor;

   /* Outputs: */
   double delta_chcl3 = 0.0;   /* Chloroform           (ug/L) */
   double delta_chbrcl2 = 0.0; /* Bromodichloromethane (ug/L) */
   double delta_chbr2cl = 0.0; /* Dibromochloromethane (ug/L) */
   double delta_chbr3 = 0.0;   /* Bromoform            (ug/L) */
   double delta_tthm = 0.0;    /* Total THM            (ug/L) */

   /* Outputs: */
   double delta_mcaa = 0.0; /* Monochloroacetic acid (ug/L) */
   double delta_dcaa = 0.0; /* Dichloroacetic   acid (ug/L) */
   double delta_tcaa = 0.0; /* Trichloroacetic  acid (ug/L) */
   double delta_mbaa = 0.0; /* Monobromoacetic  acid (ug/L) */
   double delta_dbaa = 0.0; /* Dibromoacetic    acid (ug/L) */
   double delta_haa5 = 0.0; /* sum of above 5 HAAs   (ug/L) */
   double delta_bcaa = 0.0; /* Bromochloroacetic acid (ug/L)*/
   double delta_haa6 = 0.0; /* sum of above 6 HAAs   (ug/L) */
   double delta_tbaa = 0.0;
   double delta_dbcaa = 0.0;
   double delta_bdcaa = 0.0;
   double delta_haa9 = 0.0;
   double delta_haa6_for_haa9 = 0.0;
   double delta_haa9_haa6 = 0.0;
   double delta_tox = 0.0;

   /* Internal: */
   double factor, term1, term2, term3, term4, term5, term6;
   double eff_Br, moles_br_incorporated;
   register struct Effluent *eff;

   /* Self-protection: there should have been an RM encountered by now */
   if (unit->eff.last_rm_inf != NULL)
   {

      /* Get most water quality inputs from the current unit process */
      pH = unit->eff.pH;
      DegC = unit->eff.DegK - 273.15;
      br = unit->eff.Br_at_last_cl2 * MW_Br;
      br_ugl = br * 1000;
      // FreeCl    = unit->eff.FreeCl2;
      // nh2cl     = unit->eff.NH2Cl;
      doc_out = unit->eff.TOC;

      /* Get TOC and UV from effluent of process just preceding the Rapid Mix
     ince the pre-chlorination factors take organics removal during
     chlorination into account */
      eff = &unit->eff.last_rm_inf->eff;
      doc = eff->TOC;
      uv = eff->UV;

      /* The chlorine dose is retrieved from the global variable "modrw2dbpcl2" which stores
     the value of the cl2 residual at the RM effluent plus the cl2 dose just after the RM */
      cl2dose = modrw2dbpcl2;

      /*  Get Cumulative time to influent of UnitProcess.    */
      if (unit->eff.wtp_effluent != NULL)
      {
         /* 'unit' is a distribution sample -- get influent time from 
    *  WTP_EFFLUENT.
    */
         inf_hours = unit->eff.wtp_effluent->eff.hours;
         inf_FreeCl = unit->eff.wtp_effluent->eff.FreeCl2;
         inf_nh2cl = unit->eff.wtp_effluent->eff.NH2Cl;
      }
      else
      {
         if (PrevUnitProcess(unit) != NULL)
         {
            inf_hours = PrevUnitProcess(unit)->eff.hours;
            inf_FreeCl = PrevUnitProcess(unit)->eff.FreeCl2;
            inf_nh2cl = PrevUnitProcess(unit)->eff.NH2Cl;
         }
         else
         { /* 'unit' is first unit process in process train. */
            inf_hours = 0.0;
            inf_FreeCl = 0.0;
            inf_nh2cl = 0.0;
         }
      }

      /* This subroutine will calculate DBPs formed only in one process at a time
     and add them to those present in the influent of the unit process, so
     need to have the correct time term.  Also, if chlorination occurred before
     the RM, then (global) modrw2dbptime contains the res. time in the RM, so
     that we will not be adding the beginning part of the chlorination curve again */

      inf_time = modrw2dbptime + inf_hours;
      eff_time = modrw2dbptime + unit->eff.hours;

      /* Self-protection - Set minimum values */
      if (DegC <= 1.0)
         DegC = 1.0; /* Minimum Temperature */
      if (uv <= 0.001)
         uv = 0.001; /* Minimum uv          */
      if (br_ugl <= 1.0)
         br_ugl = 1.0; /* Min. br_ugl         */

      if (doc > 0.0 && doc < 0.1)
         doc = 0.1; /* Min. non-zero doc  */

      doc_uv = doc * uv;

      /* Calculate the Solarik factors to apply to the Sohn raw water model predictions */

      if ((inf_FreeCl > 0.0 || inf_nh2cl > 0.0) && doc > 0.0 && cl2dose > 0.0 && eff_time > 0.0)
      { /* THMs and/or HAAs should be forming */

         /* Calculate DBP reduction factors to apply to raw water formation model */

         doc_reduction = (doc - doc_out) / doc;
         if (doc_reduction < 0.0)
            doc_reduction = 0.0;

         thm_factor = 1.0 - (0.875 * doc_reduction);
         haa_factor = 1.0 - (0.776 * doc_reduction);
         tox_factor = 1.0 - (0.865 * doc_reduction);

         /* Now develop THM and HAA estimates based upon Jinsik Sohn's Raw Water Eq. */

         /* Estimate TTHM (ug/liter) ------------------------ Equation XX*/
         term1 = pow(doc, 1.098);
         term2 = pow(cl2dose, 0.152);
         term3 = pow(br_ugl, 0.068);
         term4 = pow(DegC, 0.609);
         term5 = pow(pH, 1.601);
         term6 = pow(eff_time, 0.263) - pow(inf_time, 0.263);
         delta_tthm = pow(10, -1.385) * term1 * term2 * term3 * term4 * term5 * term6;

         /* Estimate CHCl3 (ug/liter) ------------------------ Equation XX*/
         term1 = pow(doc, 1.617);
         term2 = pow(cl2dose, -0.094);
         term3 = pow(br_ugl, -0.175);
         term4 = pow(DegC, 0.607);
         term5 = pow(pH, 1.403);
         term6 = pow(eff_time, 0.306) - pow(inf_time, 0.306);
         delta_chcl3 = pow(10, -1.205) * term1 * term2 * term3 * term4 * term5 * term6;

         /* Estimate CHBrCl2 (ug/liter) ------------------------ Equation XX*/
         term1 = pow(doc, 0.901);
         term2 = pow(cl2dose, 0.017);
         term3 = pow(br_ugl, 0.733);
         term4 = pow(DegC, 0.498);
         term5 = pow(pH, 1.511);
         term6 = pow(eff_time, 0.199) - pow(inf_time, 0.199);
         delta_chbrcl2 = pow(10, -2.874) * term1 * term2 * term3 * term4 * term5 * term6;

         /* Estimate CHBr2Cl (ug/liter) ------------------------ Equation XX*/
         term1 = pow(doc, -0.226);
         term2 = pow(cl2dose, 0.108);
         term3 = pow(br_ugl, 1.81);
         term4 = pow(DegC, 0.512);
         term5 = pow(pH, 2.212);
         term6 = pow(eff_time, 0.146) - pow(inf_time, 0.146);
         delta_chbr2cl = pow(10, -5.649) * term1 * term2 * term3 * term4 * term5 * term6;

         /* Estimate CHBr3  (ug/liter) ------------------------ Equation XX*/
         term1 = pow(doc, -0.983);
         term2 = pow(cl2dose, 0.804);
         term3 = pow(br_ugl, 1.765);
         term4 = pow(DegC, 0.754);
         term5 = pow(pH, 2.139);
         term6 = pow(eff_time, 0.566) - pow(inf_time, 0.556);
         delta_chbr3 = pow(10, -7.83) * term1 * term2 * term3 * term4 * term5 * term6;

         /*  Estimate HAAs */

         /* Estimate HAA6  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc, 0.935);
         term2 = pow(cl2dose, 0.443);
         term3 = pow(br_ugl, -0.031);
         term4 = pow(DegC, 0.387);
         term5 = pow(pH, -0.655);
         term6 = pow(eff_time, 0.178) - pow(inf_time, 0.178);
         delta_haa6 = 9.98 * term1 * term2 * term3 * term4 * term5 * term6;

         /* Estimate HAA5  (ug/liter) ------------------------ Equation XX */
         /*    term1 = pow( doc      , 0.997  );
      term2 = pow( cl2dose  , 0.278  );
      term3 = pow( br_ugl   ,-0.138  );
      term4 = pow( DegC     , 0.341  );
      term5 = pow( pH       , -0.799 );
      term6 = pow( eff_time, 0.169 ) - pow( inf_time, 0.169 );
      delta_haa5 = 30.0*term1*term2*term3*term4*term5*term6; */

         /* Estimate MCAA  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc, 0.173);
         term2 = pow(cl2dose, 0.397);
         term3 = pow(br_ugl, 0.029);
         term4 = pow(DegC, 0.573);
         term5 = pow(pH, -0.279);
         if (inf_time > 0.0)
            term6 = pow(eff_time, -0.009) - pow(inf_time, -0.009);
         else
            term6 = pow(eff_time, -0.009);
         delta_mcaa = 0.45 * term1 * term2 * term3 * term4 * term5 * term6;
         if (delta_mcaa > 11.0)
            delta_mcaa = 11.0; //based on eqn. and model development input database

         /* Estimate MBAA  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc, -0.584);
         term2 = pow(cl2dose, 0.754);
         term3 = pow(br_ugl, 1.1);
         term4 = pow(DegC, 0.707);
         term5 = pow(pH, 0.604);
         term6 = pow(eff_time, 0.09) - pow(inf_time, 0.09);
         delta_mbaa = 6.21 * pow(10, -5.0) * term1 * term2 * term3 * term4 * term5 * term6;

         /* Estimate DCAA  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc, 1.396);
         term2 = pow(cl2dose, 0.379);
         term3 = pow(br_ugl, -0.149);
         term4 = pow(DegC, 0.465);
         term5 = pow(pH, 0.2);
         term6 = pow(eff_time, 0.218) - pow(inf_time, 0.218);
         delta_dcaa = 0.3 * term1 * term2 * term3 * term4 * term5 * term6;

         /* Estimate TCAA  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc, 1.152);
         term2 = pow(cl2dose, 0.331);
         term3 = pow(br_ugl, -0.229);
         term4 = pow(DegC, 0.299);
         term5 = pow(pH, -1.627);
         term6 = pow(eff_time, 0.18) - pow(inf_time, 0.18);
         delta_tcaa = 92.68 * term1 * term2 * term3 * term4 * term5 * term6;

         /* Estimate DBAA  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc, -1.086);
         term2 = pow(cl2dose, 0.673);
         term3 = pow(br_ugl, 2.052);
         term4 = pow(DegC, 0.38);
         term5 = pow(pH, -0.001);
         term6 = pow(eff_time, 0.095) - pow(inf_time, 0.095);
         delta_dbaa = 3.59 * pow(10, -5.0) * term1 * term2 * term3 * term4 * term5 * term6;

         /*  Estimate BCAA  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc, 0.463);
         term2 = pow(cl2dose, 0.522);
         term3 = pow(br_ugl, 0.667);
         term4 = pow(DegC, 0.379);
         term5 = pow(pH, 0.581);
         term6 = pow(eff_time, 0.22) - pow(inf_time, 0.22);
         delta_bcaa = 5.51 * pow(10, -3.0) * term1 * term2 * term3 * term4 * term5 * term6;

         //  The following are coag.-water HAA9 species and TOX models developed by MPI in 4/2001
         //  (used here in place of a suitable raw water model)

         /*  Estimate TBAA  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc_uv, 0.0657);
         term2 = pow(cl2dose, -2.51);
         term3 = pow(br_ugl, 2.32);
         term4 = pow(1.059, (DegC - 20.0));
         term5 = pow(0.555, (pH - 8.0));
         term6 = pow(eff_time, 1.26) - pow(inf_time, 1.26);
         delta_tbaa = 5.59 * pow(10, -6.0) * term1 * term2 * term3 * term4 * term5 * term6;

         /*  Estimate DBCAA  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc_uv, -0.0162);
         term2 = pow(cl2dose, -0.170);
         term3 = pow(br_ugl, 0.972);
         term4 = pow(1.054, (DegC - 20.0));
         term5 = pow(0.839, (pH - 8.0));
         term6 = pow(eff_time, 0.685) - pow(inf_time, 0.685);
         delta_dbcaa = 3.7 * pow(10, -3.0) * term1 * term2 * term3 * term4 * term5 * term6;

         /*  Estimate BDCAA  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc_uv, 0.230);
         term2 = pow(cl2dose, 0.140);
         term3 = pow(br_ugl, 0.301);
         term4 = pow(1.022, (DegC - 20.0));
         term5 = pow(0.700, (pH - 8.0));
         term6 = pow(eff_time, 0.422) - pow(inf_time, 0.422);
         delta_bdcaa = 0.589 * term1 * term2 * term3 * term4 * term5 * term6;

         /*  Estimate HAA9  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc_uv, 0.250);
         term2 = pow(cl2dose, 0.500);
         term3 = pow(br_ugl, 0.054);
         term4 = pow(1.015, (DegC - 20.0));
         term5 = pow(0.894, (pH - 8.0));
         term6 = pow(eff_time, 0.348) - pow(inf_time, 0.348);
         delta_haa9 = 10.783 * term1 * term2 * term3 * term4 * term5 * term6;

         /*  Estimate HAA6 for HAA9  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc_uv, 0.320);
         term2 = pow(cl2dose, 0.510);
         term3 = pow(br, -0.106);
         term4 = pow(1.017, (DegC - 20.0));
         term5 = pow(0.938, (pH - 8.0));
         term6 = pow(eff_time, 0.353) - pow(inf_time, 0.353);
         delta_haa6_for_haa9 = 18.6 * term1 * term2 * term3 * term4 * term5 * term6;

         /*  Estimate TOX  (ug/liter) ------------------------ Equation XX */
         term1 = pow(doc_uv, 0.362);
         term2 = pow(cl2dose, 0.129);
         term4 = pow(DegC, 0.211);
         term6 = pow(eff_time, 0.182) - pow(inf_time, 0.182);
         delta_tox = 109.0 * term1 * term2 * term4 * term6;

         /*******PRE-CL2 FACTORS************/
         /* Modify the delta_XXXXs (Jinsik Sohn Raw Water Model Predictions)
       by the appropriate XXX_factors (from G. Solarik Models) */
         delta_tthm *= thm_factor;
         delta_chcl3 *= thm_factor;
         delta_chbrcl2 *= thm_factor;
         delta_chbr2cl *= thm_factor;
         delta_chbr3 *= thm_factor;
         //    delta_haa5    *= haa_factor;
         delta_haa6 *= haa_factor;
         delta_mcaa *= haa_factor;
         delta_mbaa *= haa_factor;
         delta_dcaa *= haa_factor;
         delta_tcaa *= haa_factor;
         delta_dbaa *= haa_factor;
         delta_bcaa *= haa_factor;
         delta_tbaa *= haa_factor;
         delta_dbcaa *= haa_factor;
         delta_bdcaa *= haa_factor;
         delta_haa9 *= haa_factor;
         delta_tox *= tox_factor;

         /********CORRECTION FACTORS BASED ON ICR DATA CALIBRATION WORK (WJS, 5/2001)********/

         /* Correction factors to apply for free chlorine based on comparing
       predicted vs. observed in months 7-18 in AUX8 at the FINISH water
       location for plantmonths with Validation mode = VALID,
       WTP_DIS = CL2, predicted alk > 0 at FINISH, and predicted nh2cl = 0
       at FINISH */
         delta_tthm /= 0.77;
         delta_chcl3 /= 1.00;
         delta_chbr2cl /= 0.50;
         delta_chbrcl2 /= 0.86;
         delta_chbr3 /= 1.00;
         //      delta_haa5    /= 1.17;
         delta_haa6 /= 1.0;
         delta_mcaa /= 1.00;
         delta_dcaa /= 0.71;
         delta_tcaa /= 1.3;
         delta_mbaa /= 1.00;
         delta_dbaa /= 1.00;
         delta_bcaa /= 0.82;

         /* In Distribution system, need several different correction factors */
         if (unit->type == AVG_TAP || unit->type == END_OF_SYSTEM ||
             unit->type == LOCATION_1)
         {
            delta_dcaa *= (0.71 / 1.3);
            delta_bcaa *= (0.82 / 2.0);
            delta_chcl3 /= 1.1;
            delta_chbr2cl /= 0.80;
         }

         /* If chloramines exist, use appropriate factors for the ratio of formation
       with chloramines to that with free chlorine (20% is the default, the
       others were set based on calibration with AUX8 ICR plants with CLM as
       the DS disinfectant) */

         if (inf_FreeCl <= 0.0 && inf_nh2cl > 0.0)
         {
            delta_tthm *= 0.3;
            delta_chcl3 *= 0.3;
            delta_chbr2cl *= 0.3;
            delta_chbrcl2 *= 0.3;
            delta_chbr3 *= 0.3;
         }

         if (inf_nh2cl > 0.0)
         {
            //      delta_haa5 *= 0.05;
            delta_haa6 *= 0.35;
            delta_haa9 *= 0.2;
            delta_haa6_for_haa9 *= 0.2;
            delta_mcaa *= 0.2;
            delta_dcaa *= 0.35;
            delta_tcaa *= 0.05;
            delta_mbaa *= 0.2;
            delta_dbaa *= 0.2;
            delta_bcaa *= 0.3;
            delta_tbaa *= 0.2;
            delta_dbcaa *= 0.2;
            delta_bdcaa *= 0.2;
            delta_tox *= 0.325;
         }

         //MCAA should only be able to go to zero
         if ((delta_mcaa + unit->eff.MCAA) < 0.0)
            delta_mcaa = -unit->eff.MCAA;

         /******************************END OF CORRECTION FACTORS*******************************/

         /* Proportion the species based on the bulk TTHM */
         if (delta_chcl3 > 0.0 || delta_chbrcl2 > 0.0 || delta_chbr2cl > 0.0 || delta_chbr3 > 0.0)
         {
            factor = delta_tthm / (delta_chcl3 + delta_chbrcl2 +
                                   delta_chbr2cl + delta_chbr3);
            delta_chcl3 *= factor;
            delta_chbrcl2 *= factor;
            delta_chbr2cl *= factor;
            delta_chbr3 *= factor;
         }
         else
            delta_tthm = 0.0;

         /* Proportion-out 6 HAA species based on HAA6 equation */
         if (delta_mcaa > 0.0 || delta_mbaa > 0.0 || delta_bcaa > 0.0 ||
             delta_dcaa > 0.0 || delta_dbaa > 0.0 || delta_tcaa > 0.0)
         {
            factor = delta_haa6 / (delta_mcaa + delta_mbaa + delta_bcaa +
                                   delta_dcaa + delta_dbaa + delta_tcaa);
            delta_mcaa *= factor;
            delta_mbaa *= factor;
            delta_dcaa *= factor;
            delta_dbaa *= factor;
            delta_tcaa *= factor;
            delta_bcaa *= factor;
         }
         else
            delta_haa6 = 0.0;

         /* Calc. HAA5 based on HAA6 equation, minus proportioned BCAA */
         delta_haa5 = delta_haa6 - delta_bcaa;

         /* Proportion the HAA9 species based on HAA9-HAA6 */
         if (delta_haa9 > delta_haa6_for_haa9)
         {
            delta_haa9_haa6 = delta_haa9 - delta_haa6_for_haa9;
         }
         else
         {
            delta_haa9_haa6 = 0.0;
            delta_dbcaa = 0.0;
            delta_bdcaa = 0.0;
            delta_tbaa = 0.0;
         }

         if (delta_dbcaa > 0.0 || delta_bdcaa > 0.0 || delta_tbaa > 0.0)
         {
            factor = delta_haa9_haa6 / (delta_dbcaa + delta_bdcaa + delta_tbaa);
            delta_dbcaa *= factor;
            delta_bdcaa *= factor;
            delta_tbaa *= factor;
         }

         /* Calc. HAA9 based on HAA6, plus species */
         delta_haa9 = delta_haa6 + delta_tbaa + delta_dbcaa + delta_bdcaa;

         /* Add incremental outputs in UnitProcess data structure. */

         /* THMs */
         /* THMs */
         unit->eff.CHCl3 += delta_chcl3;
         unit->eff.CHBrCl2 += delta_chbrcl2;
         unit->eff.CHBr2Cl += delta_chbr2cl;
         unit->eff.CHBr3 += delta_chbr3;
         unit->eff.TTHM += delta_tthm;

         /* HAAs */
         unit->eff.MCAA += delta_mcaa;
         unit->eff.DCAA += delta_dcaa;
         unit->eff.TCAA += delta_tcaa;
         unit->eff.MBAA += delta_mbaa;
         unit->eff.DBAA += delta_dbaa;
         unit->eff.BCAA += delta_bcaa;
         unit->eff.TBAA += delta_tbaa;
         unit->eff.DBCAA += delta_dbcaa;
         unit->eff.BDCAA += delta_bdcaa;
         unit->eff.HAA5 += delta_haa5;
         unit->eff.HAA6 += delta_haa6;
         unit->eff.HAA9 += delta_haa9;

         /*TOX */
         unit->eff.tox += delta_tox;

         /*Stoichiometric Bromide incorporation (non-THM or HAA TOX not considered)*/
         moles_br_incorporated =
             (delta_chbrcl2 / 1000.0) / MW_BDCM * 1.0 +
             (delta_chbr2cl / 1000.0) / MW_DBCM * 2.0 +
             (delta_chbr3 / 1000.0) / MW_CHBR3 * 3.0 +
             (delta_mbaa / 1000.0) / MW_MBAA * 1.0 +
             (delta_bcaa / 1000.0) / MW_BCAA * 1.0 +
             (delta_dbaa / 1000.0) / MW_DBAA * 2.0 +
             (delta_bdcaa / 1000.0) / MW_BDCAA * 1.0 +
             (delta_dbcaa / 1000.0) / MW_DBCAA * 2.0 +
             (delta_tbaa / 1000.0) / MW_TBAA * 3.0;

         eff_Br = unit->eff.Br - moles_br_incorporated;

         if (eff_Br > unit->eff.Br)
            eff_Br = unit->eff.Br; /* No increasing Br allowed */
         if (eff_Br < 0.0)
            eff_Br = 0.0; /* No negative Br allowed   */

         unit->eff.Br = eff_Br; /*Update UnitProcess data structure */

         /***********************End of Bromide Incorporation**********************/

         /* Because of possibility of a negative delta_mcaa, the HAA5 species
       could all possibly be negative due to proportioning*/
         if (unit->eff.MCAA < 0.0)
            unit->eff.MCAA = 0.0;
         if (unit->eff.DCAA < 0.0)
            unit->eff.DCAA = 0.0;
         if (unit->eff.TCAA < 0.0)
            unit->eff.TCAA = 0.0;
         if (unit->eff.MBAA < 0.0)
            unit->eff.MBAA = 0.0;
         if (unit->eff.DBAA < 0.0)
            unit->eff.DBAA = 0.0;
         if (unit->eff.BCAA < 0.0)
            unit->eff.BCAA = 0.0;

      } /* End of "if(free chlorine residual OR combined chlorine residual exists)" */
   }    /* End of "if(last_rm_inf != NULL)" */
}
