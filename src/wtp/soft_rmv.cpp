/* Soft_Rmv.c -- September 8, 1993 */
#include "wtp.h"

void soft_rmv(struct UnitProcess *unit)
/*
* Purpose:
*   Estimate removal of TOC and UV during lime softening (with or without
*   coagulant addition.  Equations are based upon regressions applied to
*   combined ICR and WITAF database.  Equations were developed 1/99 by
*   WJS.  
*
*   Documentation and Code by WJS, 1/99 */
{
  /* Inputs: *****************************************/
  double raw_toc;   /* Raw water TOC           (mg/L)  */
  double raw_uv;    /* Raw water UV            (1/cm)  */
  double pH_soft;   /* pH of softening    (in the RM)  */
  double lime_dose; /* Cum. Lime Dose (mg/L as Ca(OH)2 */
  double coag_dose; /* Cum. Coag. Dose (meq/L as metal)*/
  double suva_raw;  /* Raw water SUVA         (L/mg-m) */

  /* Internals: ***************************************/
  double toc_removed; /* TOC removed by combined lime softening/coagulation (mg/L) */
  double uv_removed;  /* UV-254 removed by combined lime softening/coagulation (1/cm) */
  double term1, term2, term3, term4;

  struct Effluent *eff;
  struct Influent *raw;

  /* Self protection */
  if (unit == NULL)
    return;
  if (unit->eff.influent == NULL)
    return;

  /* Get inputs from UnitProcess data structure. */
  eff = &unit->eff;
  raw = eff->influent->data.influent;
  if (raw == NULL)
    return;

  raw_toc = raw->toc;
  raw_uv = raw->uv254;
  pH_soft = eff->pH;
  lime_dose = eff->LimeDose;
  coag_dose = eff->AlumDose / 297.0 + eff->FericDose / 270.0;
  suva_raw = raw_uv / raw_toc * 100.0;

  if (unit->type == RAPID_MIX)
  { /* This is the only place TOC/UV removal will take place for softening */

    /* Estimate TOC Removal*/
    term1 = pow(raw_toc, 1.3843);
    term2 = pow(pH_soft, 2.2387);
    term3 = pow(lime_dose, 0.1707);
    term4 = pow((1.0 + coag_dose), 2.4402);
    toc_removed = 0.0004657 * term1 * term2 * term3 * term4;

    /* Estimate UV-254 Removal*/
    term1 = pow(toc_removed, 0.8367);
    term2 = pow(suva_raw, 1.2501);
    uv_removed = 0.01685 * term1 * term2;

    /*Update Data Structure*/
    eff->TOC -= toc_removed;

    /****** Softening Plant TOC Adjustment Function  - WJS May, 2001***********/
    /*   (based on AUX8 periods 7-18, softening plants validation) */
    eff->TOC = f_adj_toc(eff->TOC, unit);
    /*************************************************************************/

    eff->UV -= uv_removed;
    eff->UV_out -= uv_removed;
    eff->toc_sludge += toc_removed;
    eff->alum_to_sludge += eff->AlumDose;
    eff->iron_to_sludge += eff->FericDose;
    eff->AlumDose = 0.0;
    eff->FericDose = 0.0;

    if (eff->TOC < 0.0)
      eff->TOC = 0.0;
    if (eff->UV < 0.0)
      eff->UV = 0.0;
    if (eff->UV_out < 0.0)
      eff->UV_out = 0.0;
  }
}
