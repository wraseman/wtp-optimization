// solids.c

#include "wtp.h"

void solids_rmv(struct UnitProcess *unit)
/*
* Purpose:
*   Estimate solids production of any type of filter or sed basin.
*   Sludge is generated in the sedimentation basin and floc carry-over
*   of 40 mg/L of CaCO3 is implemented to reduce pH swings following softening.
*   The floc carry-over is based upon recommendations of a Sep. '92 Tech. Memo.
*   by MPI which detailed the results of an analysis of WTP model results versus
*   actual behavior in softening plants.
* Documentation and Code by WJS, 7/05
*/

{
  double sludge = 0.0;
  double caco3_floc = 0.0;

  struct Effluent *eff;

  eff = &unit->eff;

  /* Estimate solids production based on type of sludge to produce */

  if (eff->floc_carry_cntr == FALSE && eff->Ca_solid > 0.0)
  { /* there's calcium carbonate flocc to remove and we're only going to have
	  flocc carry-over one time in the train */
    caco3_floc = (eff->Ca_solid - 0.0004) * MW_CaCO3;
    if (caco3_floc < 0.0)
      caco3_floc = 0.0;
    eff->Ca_solid = 0.0004; /*floc carry-over term = 40 mg/L of CaCO3*/
    eff->floc_carry_cntr = TRUE;
  }
  else
  {
    caco3_floc = eff->Ca_solid * MW_CaCO3;
    eff->Ca_solid = 0.0;
  }

  sludge += caco3_floc +
            eff->Mg_solid * MW_MgOH2 +
            eff->toc_sludge +
            eff->alum_to_sludge * 0.44 +
            eff->iron_to_sludge / 270 * 55.85 * 2.9 + /* to convert to mg of Fe/L*/
            eff->Turbidity * 1.4;

  eff->solids += sludge;
  eff->Mg_solid = 0.0;
  eff->iron_to_sludge = 0.0;
  eff->alum_to_sludge = 0.0;
  eff->toc_sludge = 0.0;
  eff->Turbidity = 0.0;

} /* End of solids_rmv() */