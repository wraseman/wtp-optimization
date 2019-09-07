/* BreakPt.c */
#include "wtp.h"

void breakpt(struct UnitProcess *unit)
/*
* Purpose: Breakpoint chlorination.  This function will convert
*   free chlorine and ammonia into monochloramine and dichloramine.
*
* Inputs:
*   Free Cl2, NH3, NH2Cl, NHCl2, CBminusCA, pH
*     On input, these are molar concentrations in the influent to the
*     unit process.  The units are (Mole/L) execpt CBminusCA which is
*     (equ/L) and pH.
*
* Outputs:
*   Free Cl2, NH3, NH2Cl, NHCl2, CBminusCA, pH
*     On output, these inputs are converted to effluent concentration
*     from the unit process.
*
* Method:
*   Breakpoint chemistry is assumed to be the following:
*      [HOCl] +  [NH3]   -> [H2O] + [NH2Cl]                      Eq 1.
*      [HOCl] + 2[NH2Cl] -> [H2O] + [N2]    + 3[H+] + 3[Cl-]     Eq 2.
*
* Notes:
*  1. CBminusCA is changed via Eq 2.  This function contains logic to call
*     phchange() when necessary.
*  2. This function assumes instantious reactions.  Reactions 1 and 2 are
*     known to have rates that are temperature and pH dependent however,
*     this function does not include any reaction rates.  
*  3. I could not find a discription of breakpoint in the 1.21 version of
*     WTP manual.  Eq 1 and Eq 2 were found in the 1.1 code from MPI.
*  4. Formation of dichloramine is not included in this version. The input
*     concentration of dichloramine is returned as the effluent
*     concentration.
*
* Documentation and code by M.Cummins  May 19, 1993
*/
{
  /* Input and outputs: ************************************/
  double freecl2;   /* Free chlorine            (Mole/L) */
  double nh3;       /* Ammonia                  (Mole/L) */
  double nh2cl;     /* Monochloramine           (Mole/L) */
  double nhcl2;     /* Dichloramine         (Mole/Liter) */
  double CBminusCA; /* Charge of untracked ions  (equ/L) */
  /***** pH          * passed to phchange()          (-) */

  /* Internal:*/
  int change_ph_flag = FALSE;
  struct Effluent *eff;

  /* Get inputs from UnitProcess data structure. */
  eff = &unit->eff;
  freecl2 = eff->FreeCl2;
  nh3 = eff->NH3;
  nh2cl = eff->NH2Cl;
  nhcl2 = eff->NHCl2;
  CBminusCA = eff->CBminusCA;

  /* Will reaction 1 (formation of monochloramine) occur? */
  if ((freecl2 > 0.0) && (nh3 > 0.0))
  {
    /* Reaction 1 will occour... Will Free Cl2 or NH3 limit reaction? */
    if (freecl2 < nh3)
    {
      /* Free chlorine limits amount of monochloramine formed. */
      nh2cl += freecl2;
      nh3 -= freecl2;
      freecl2 = 0.0;
    }
    else
    { /* NH3 limits amount of monochloramine formed. */
      nh2cl += nh3;
      freecl2 -= nh3;
      nh3 = 0.0;
    }
  }

  /* Will reaction 2 (breakpoint) occur? */
  if ((freecl2 > 0.0) && (nh2cl > 0.0))
  {
    /* Reaction 2 will occur... Will free cl2 or NH2Cl limit reaction? */
    if (freecl2 < 0.5 * nh2cl)
    /*
      * I have stoped every time I see the above to re-think the logic of
      * freecl2<0.5*nh2cl.  Reaction 2 is.. 1 mole of free chlorine will 
      * react with 2 mole of NH2Cl.  Therefore the balance point, all
      * reactants consumed, is [FreeCl2]=0.5*[NH2Cl]. Example, 
      * (1 mole)=0.5*(2 mole).  If free chlorine is less than 1 mole
      * say 0.9 mole then free chlorine will limit the reaction.
      */
    {
      /* Free chlorine limits breakpoint. */
      nh2cl -= 2.0 * freecl2;
      CBminusCA -= 3.0 * freecl2;
      freecl2 = 0.0;
    }
    else
    { /* Monochloramine limits breakpoint. */
      freecl2 -= 0.5 * nh2cl;
      CBminusCA -= 1.5 * nh2cl;
      nh2cl = 0.0;
    }
    change_ph_flag = TRUE;
  }

  /* Copy outputs to UnitProcess data structure */
  eff->FreeCl2 = freecl2;
  eff->NH3 = nh3;
  eff->NH2Cl = nh2cl;
  eff->NHCl2 = nhcl2;
  eff->CBminusCA = CBminusCA;

  if (change_ph_flag)
    phchange(unit, TRUE);
}
