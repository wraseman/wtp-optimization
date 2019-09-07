/* basin.c */
#include "wtp.h"

void basn_dbp(struct UnitProcess *unit)
/*
* Purpose:
*   Estimate chlorine decay, DBP formation, and CT ratio in a basin.
*
*
* Documention and code by WJS - May, 2001
*/
{
  //  FILE *fptr2;
  double theo_res_time;     /* This is the EBCT of the unit process basin (calc.)*/
  double mean_theo;         /* This is the ratio of the mean residence time to the
			   EBCT of the unit process basin (as entered by user) */
  double t_ten_theo;        /* This is the ratio of T10 to the EBCT of the unit
  			   process basin (as entered by the user */
  double t_ten;             /* This is T10 of the unit process basin (calc.) */
  double rxnhours = 0.00;   /* This is the mean res. time of each of the small CFSTRs
			   that the unit process basin is modeled as (calc.)*/
  double rxnminutes = 0.00; /* This is the mean res. time of each of the small CFSTRs
			   that the unit process basin is modeled as (calc.)*/
  double ncstrs;            /* Number of CSTRs that the unit process basin is
			   modeled as (calc.) */
  double clo2_avg = 0.0;    /* Average concentration for ClO2 for CT calcs */
  int i;

  register struct Effluent *eff;
  struct Basin *basin;
  struct Presed *presed;

  basin = unit->data.basin;
  presed = unit->data.presed;
  eff = &unit->eff;

  /*Zero ozone residual, if not an ozone contactor */
  if (unit->type != O3_CONTACTOR)
    eff->o3_res = 0.0;

  if (unit->type == PRESED_BASIN)
    theo_res_time = (presed->volume) / (eff->Flow / 1440); //minutes
  else
    theo_res_time = (basin->volume) / (eff->Flow / 1440); //minutes

  if (theo_res_time > 0.0)
  { // Disinfectant decay and DBP formation only if there is detention time

    if (unit->type == PRESED_BASIN)
      mean_theo = presed->sb_mean;
    else
      mean_theo = basin->sb_mean;

    if (unit->type == PRESED_BASIN)
      t_ten_theo = presed->sb_t_ten;
    else
      t_ten_theo = basin->sb_t_ten;

    t_ten = t_ten_theo * theo_res_time; //minutes

    // Using new_cl2decay() in the framework of CFSTRs in series
    //--------------------------------------------------------------------------------
    if (eff->cl2cnt > 0.0)
    {
      ncstrs = ncstr(t_ten_theo);
      rxnhours = mean_theo * theo_res_time / (60.0 * ncstrs); //hours
      for (i = 1; i <= ncstrs; i++)
      {
        new_cl2decay(unit, rxnhours);
        if ((eff->FreeCl2 + eff->NH2Cl) > 0.0)
          eff->hours += rxnhours;
      }

      //Calculate DBP formation
      choose_dbpmodel(unit); /* NEW CODE (WJS, 10/98) */
    }

    // Using clo2_decay() in the framework of CFSTRs in series
    //--------------------------------------------------------------------------------
    if (eff->clo2_cntr > 0)
    {
      ncstrs = ncstr(t_ten_theo);
      rxnminutes = mean_theo * theo_res_time / ncstrs; //minutes

      for (i = 1; i <= ncstrs; i++)
      {
        clo2_decay(unit, rxnminutes);
        clo2_avg += eff->clo2_res;
        eff->clo2_minutes += rxnminutes;
      }
      eff->clo2_avg = clo2_avg / ncstrs;
    }

    /*DEBUGGING CODE*/
    //      fptr2=fopen("debug.txt","a+");
    //      fprintf(fptr2,"Module: %d  after  'if'  \n",unit->type);
    //        fprintf(fptr2,"Free_Cl2:              %f\n",unit->eff.FreeCl2);
    //        fclose(fptr2);
    /*DEBUGGING CODE*/

    //Calculate CT
    if (eff->cl2cnt > 0.0 || eff->o3_res > 0.0 || eff->clo2_res > 0.0)
      ct(unit, t_ten);

  } //End "if(theo_res_time > 0.0)"
}
