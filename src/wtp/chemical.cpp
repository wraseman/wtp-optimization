/* Chemical.c -- Chemical addition functions
*
*  Purpose:
*    Calculate alkalinity change and solids production due to various
*    chemical additions.
*
*/
#include "wtp.h"
/*
* The following functions are in this file:
* void chem_reset ( struct ProcessTrain *train );
* void alumadd(  struct UnitProcess *unit );
* void fecladd(  struct UnitProcess *unit );
* void chloradd( struct UnitProcess *unit );
* void naohadd(  struct UnitProcess *unit );
* void h2so4add( struct UnitProcess *unit );
* void nh3add(   struct UnitProcess *unit );
* void limeadd(  struct UnitProcess *unit );
* void co2_add(  struct UnitProcess *unit );
* void kmno4add( struct UnitProcess *unit );
* void soda_add( struct UnitProcess *unit );
* void clo2add(  struct UnitProcess *unit );
* void so2add(   struct UnitProcess *unit );
*/

void chem_reset(struct ProcessTrain *train)

/* Purpose: Set chemical doses in treatment train to zero. Used for 
 * 	simulation and simulation-optimization mode to clear chemical doses
 * 	from previous model runs.
 * 
 * Input: 
 * 	train = treatment train structure
 *  unit = unit process structure
 * 
 * Return: void
 * 
 * Note: this function currently changes lime, carbon dioxide, alum, and 
 * 	sodium hypochlorite back to zero. Must change to include other
 * 	chemicals. 
 * 
 */

{

  struct UnitProcess *unit;

  for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
  {
    switch (unit->type)
    {
    case LIME:
      unit->data.lime->dose = 0.0;
      break;

    case CARBON_DIOXIDE:
      unit->data.chemical->co2 = 0.0;
      break;

    case ALUM:
      unit->data.alum->dose = 0.0;
      break;

    case HYPOCHLORITE:
      unit->data.chemical->naocl = 0.0;
      break;

    default:
      break;
    }
  }

  return;
}

void alumadd(struct UnitProcess *unit)
/*
*  Purpose:  Calculate change in alkalinity due to alum addition.
*
*
*  Input:
*    unit->data.alum->dose = alum dose (mg/L) as Al2(SO4)3-14H2O
*
*  Return:
*    AlumDose = Cumulative alum dose (mg/L) as Al2(SO4)3-14H2O
*    raw_toc  = TOC at the point where alum is first added.
*    raw_uv   = UV   "  "   "     "     "   "   "     "
*    CBminusCA
*
*  Notes:
*   1. Model currently assumes all aluminum ions form an insoluble
*      precipitate.  Therefore, actual alkalinity decreases should be
*      smaller than those predicted here.
*
*  Hummmm....
*         Al2(SO4)3*14H2O -> 2 Al(OH)3 + 3 SO4-- + 6 H+ + 8 H2O
*
*  Documentation and code by M.Cummins April 28, 1993
*/
{
  double dose;
  register struct Effluent *eff;

  if (unit == NULL || unit->type != ALUM)
    return;
  eff = &unit->eff;

  dose = unit->data.alum->dose; /* Alum dose: (mg/L) Al2(SO4)3*14H2O */

  if (dose > 0.0)
  {
    eff->AlumDose += dose;
    eff->CBminusCA -= (6 * dose / 594364L);
  }
}

void fecladd(struct UnitProcess *unit)
/*
*   Calculate change in alkalinity due to ferric addition.  Model
*   currently assumes all ferric ions form an insoluble
*   precipitate.  Therefore, actual alkalinity decreases should be
*   smaller than those predicted here.
*/
{
  double dose;
  register struct Effluent *eff;

  if (unit == NULL || unit->type != IRON)
    return;
  eff = &unit->eff;

  dose = unit->data.iron->dose;

  if (dose > 0.0)
  {
    /*  if( eff->FericDose==0.0 )
        {
          eff->raw_toc = eff->TOC;
          eff->raw_uv  = eff->UV;
        }  */
    eff->FericDose += dose;
    eff->CBminusCA -= (3 * dose / 270296L);
  }
}

void chloradd(struct UnitProcess *unit)
{
  double dose;
  register struct Effluent *eff;

  if (unit == NULL || (unit->type != CHLORINE && unit->type != HYPOCHLORITE))
    return;
  eff = &unit->eff;

  if (unit->type == CHLORINE)
  {
    dose = unit->data.chemical->chlor;
  }
  else
  { /*unit->type == HYPOCHLORITE*/
    dose = unit->data.chemical->naocl;
  }

  if (dose > 0.0)
  {
    eff->Br_at_last_cl2 = eff->Br;
    eff->cl2cnt += 1;
    eff->hours = 0.0;
    eff->pre_chlor_dose_track += dose;
    eff->FreeCl2 += (dose / MW_Cl2);
    breakpt(unit);
    eff->cl2dose = MW_Cl2 * (eff->FreeCl2 + eff->NH2Cl);
    if (unit->type == CHLORINE)
    {
      eff->CBminusCA -= dose / MW_Cl2;
    }
    else
    { /*unit->type == HYPOCHLORITE*/
      eff->CBminusCA += (dose * 2) / MW_Cl2;
    }
  }
}

void so2_add(struct UnitProcess *unit)
{
  double dose;
  double mols_reacted;
  register struct Effluent *eff;

  if (unit == NULL || unit->type != SULFUR_DIOXIDE)
    return;
  eff = &unit->eff;

  dose = unit->data.chemical->so2;

  if (dose > 0.0)
  {
    mols_reacted = dose / MW_SO2; /* Assume all SO2 will be consumed by Cl2*/
    eff->FreeCl2 -= (dose / MW_SO2);
    if (eff->FreeCl2 < 0.0)
    {
      mols_reacted = dose / MW_SO2 + eff->FreeCl2; /*Actual mols of SO2 consumed by Cl2*/
      eff->FreeCl2 = 0.0;                          /*All Cl2 used up*/
    }
    eff->CBminusCA -= 4 * mols_reacted; /*Its the mols consumed that produce acid*/
  }
}

void naohadd(struct UnitProcess *unit)
{
  unit->eff.CBminusCA += (unit->data.chemical->naoh / MW_NaOH);
}

void h2so4add(struct UnitProcess *unit)
{
  unit->eff.CBminusCA -= (2 * unit->data.chemical->h2so4 / MW_H2SO4);
}

void nh3add(struct UnitProcess *unit)
{
  unit->eff.NH3 += (unit->data.chemical->nh3 / MW_NH3);
  if (unit->type == AMMONIUM_SULFATE)
    unit->eff.CBminusCA -= unit->data.chemical->nh3 / MW_NH3;
}

void limeadd(struct UnitProcess *unit)
{
  unit->eff.Ca_aq += (unit->data.lime->dose / MW_LIME);
  unit->eff.LimeDose += unit->data.lime->dose;
  if (unit->data.lime->purpose == 'S')
    unit->eff.limesoftening = TRUE;
}

void co2_add(struct UnitProcess *unit)
{
  unit->eff.CO2_aq += (unit->data.chemical->co2 / MW_CO2);
}

void kmno4add(struct UnitProcess *unit)
{
  double dose;

  dose = unit->data.chemical->kmno4;

  unit->eff.CBminusCA -= (dose / 158034L);
  unit->eff.solids += (86937L * dose / 158034L);
}

void soda_add(struct UnitProcess *unit)
{
  double dose;

  dose = unit->data.chemical->soda;
  unit->eff.CBminusCA += (dose / (MW_SODA / 2));
  unit->eff.CO2_aq += (dose / MW_SODA);
}

void clo2add(struct UnitProcess *unit)
{

  /*FILE *fptr;*/
  double dose;
  double conversion;
  register struct Effluent *eff;

  if (unit == NULL || unit->type != CHLORINE_DIOXIDE)
    return;
  eff = &unit->eff;

  dose = unit->data.clo2->dose;
  conversion = unit->data.clo2->conversion;

  if (dose > 0.0)
  {
    eff->clo2_cntr++;
    eff->clo2_minutes = 0.0;
    eff->chlorite += dose * conversion / 100.0;
    eff->clo2dose = dose + eff->clo2_res;
    eff->clo2_res = eff->clo2dose;
  }
}

void o3add(struct UnitProcess *unit)
{

  double uv, uv_out, eff_uv, eff_uv_out, toc, o3_dose;
  register struct Effluent *eff;

  if (unit == NULL || unit->type != OZONE)
    return;
  eff = &unit->eff;

  eff->o3_minutes = 0.0; /* Set time since last O3 add. pt.= 0 */

  /* calculate the ozone residual at this point as the residual going
    in plus the dose applied here; this will be used as the dose at
    the last point of ozone addition in ozonate() */
  eff->o3_res += unit->data.chemical->o3;

  /* save WQ data at point of ozone application for use in ozonate() */
  eff->last_o3_inf = unit;

  /*Calculate UV oxidation as an instantaneous event at any
   point of ozonation.  Based on algorithm developed by G. Solarik
   from pilot/full-scale data in May, 2001. Note that in ozonate(),
   the UV is referenced as the UV in the unit process effluent
   upstream of this point.*/

  uv = eff->UV;
  eff_uv = eff->UV;
  uv_out = eff->UV_out;
  eff_uv_out = eff->UV_out;
  toc = eff->TOC;
  o3_dose = eff->o3_res;

  if (toc < 0.1)
    toc = 0.1;

  if (o3_dose > 0.0)
  {
    eff_uv = 0.622 * pow(uv, 0.931) * pow(o3_dose / toc, -0.252);
    eff_uv_out = 0.622 * pow(uv_out, 0.931) * pow(o3_dose / toc, -0.252);
    if (eff_uv > uv)
      eff_uv = uv;
    if (eff_uv_out > uv_out)
      eff_uv_out = uv_out;
  }

  //Update Unit Process Data Structure
  eff->UV = eff_uv;
  eff->UV_out = eff_uv_out;
}