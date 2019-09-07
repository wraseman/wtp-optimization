/* Dist_dbp.c */
#include "wtp.h"

void dist_dbp(struct UnitProcess *unit)
/*
*  Purpose: Estimate chlorine decay, THM, and HAA in distribution sample.
*/
{
  //FILE *fptr;
  /* Inputs: */
  double theo_res_time; /* (minutes)                          */
  double t_ten_theo;    /* Ratio of t10/theo_res_time (0-1.0) */
  double mean_theo;     /* Ratio                              */

  /* Outputs: */
  double eff_hours; /* Cumulative time to effluent.       */

  /* Internal: */
  double ncstrs;   /* Number of equ. constant stired tanks in series. */
  double rxnhours; /* Contact time for one constant stired tank.      */
  double rxnminutes;
  double peakfactor = 1.0; /* Adjust dist. sys. travel time for max. flow case*/
  int i;
  struct Effluent *eff;
  struct UnitProcess *end;

  /* Self protection: */
  if (unit == NULL || unit->data.ptr == NULL)
    return;

  /* Make sure we have no O3 residual*/
  unit->eff.o3_res = 0.0;

  /* Mark end of process train. */
  if (unit->eff.wtp_effluent == NULL)
  {
    /* This is the 1st time dist_dbp() has been called. */
    end = PrevUnitProcess(unit);
    if (end != NULL)
    {
      end->eff.wtp_effluent = end;
    }
    unit->eff.wtp_effluent = end;
  }
  if (unit->eff.wtp_effluent == NULL)
    return; /*unit is 1st in process train.*/

  /* Load unit.eff with data from end of process train.    */
  eff = &unit->eff;
  *eff = unit->eff.wtp_effluent->eff;

  if (coldflag == TRUE)
  {
    peakfactor = unit->eff.Peak / unit->eff.influent->data.influent->avg_flow;
  }

  /* Get inputs from UnitProcess data structure: */
  eff_hours = eff->hours;

  /*DEBUGGING CODE*/
  //      fptr=fopen("debug.dat","a+");
  //      fprintf(fptr,"Module: %d\n",unit->type);
  //      fprintf(fptr,"eff_hours(before):     %f\n",eff_hours);
  //      fclose(fptr);
  /*DEBUGGING CODE*/

  switch (unit->type)
  {

  case LOCATION_1:
    theo_res_time = unit->data.avg_tap->days * 1440 / peakfactor;
    mean_theo = 1.0;
    t_ten_theo = 1.0;
    break;

  case AVG_TAP:
    theo_res_time = unit->data.avg_tap->days * 1440 / peakfactor;
    mean_theo = 1.0;
    t_ten_theo = 1.0;
    break;

  case END_OF_SYSTEM:
    theo_res_time = unit->data.end_of_system->days * 1440 / peakfactor;
    mean_theo = 1.0;
    t_ten_theo = 1.0;
    break;

  default:
    theo_res_time = 0.0;
    mean_theo = 1.0;
    t_ten_theo = 1.0;
    break;
  }

  if (theo_res_time > 0.0)
  {
    //Cl2_decay based on series of CFSTRs
    if (eff->cl2cnt > 0)
    {
      phchange(unit, TRUE);

      /* Update eff_hours */
      ncstrs = ncstr(t_ten_theo);
      rxnhours = mean_theo * theo_res_time / (60.0 * ncstrs);
      for (i = 1; i <= ncstrs; i++)
      {
        new_cl2decay(unit, rxnhours);
        if ((eff->FreeCl2 + eff->NH2Cl) > 0.0)
          eff_hours += rxnhours;
        /*- Copy outputs to UnitProcess data structure - */
        eff->hours = eff_hours;
      }

      /* Call dbp formation equations: */
      /*DEBUGGING CODE*/
      //      fptr=fopen("debug.dat","a+");
      //      fprintf(fptr,"eff_hours(after):     %f\n\n",eff_hours);
      //      fclose(fptr);
      /*DEBUGGING CODE*/
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
        eff->clo2_minutes += rxnminutes;
      }
    }

  } //End"if(theo_res_time > 0.0)"
}
