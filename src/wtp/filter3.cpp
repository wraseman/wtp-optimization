/* Filter.c */

#include "wtp.h"

void filt_dbp(struct UnitProcess *unit)
/*
* Purpose:
*   Estimate chlorine decay, DBP formation and CT ratio.
*
* Documentation and code by WJS - May, 2001
*/

{
   double theo_res_time, mean_theo, t_ten_theo, t_ten, rxnhours, rxnminutes;
   double ncstrs;         /* Number of constant stired tank reactors */
   double clo2_avg = 0.0; /* Average concentration for ClO2 for CT calcs */
   int i;
   struct Filter *filt;  /* Design and Operating data for FILTERs */
   struct Ssf *ssf;      /* Design and Operating data for SLOW_FILTERs*/
   struct Def *def;      /* Design and Operating data for DE_FILTERs */
   struct Effluent *eff; /* Effluent data             */

   ssf = unit->data.ssf;
   def = unit->data.def;
   filt = unit->data.filter;
   eff = &unit->eff;

   /* No O3 residual here */
   eff->o3_res = 0.0;

   if (unit->type == FILTER)
      theo_res_time = (filt->volume) / (eff->Flow / 1440); //minutes
   else if (unit->type == SLOW_FILTER)
      theo_res_time = (ssf->volume) / (eff->Flow / 1440); //minutes
   else                                                   /*(unit->type == DE_FILTER)*/
      theo_res_time = (def->volume) / (eff->Flow / 1440); //minutes

   if (theo_res_time > 0.0)
   {
      if (unit->type == FILTER)
         mean_theo = filt->fi_mean;
      else if (unit->type == SLOW_FILTER)
         mean_theo = ssf->fi_mean;
      else /*(unit->type == DE_FILTER)*/
         mean_theo = def->fi_mean;

      if (unit->type == FILTER)
         t_ten_theo = filt->fi_t_ten;
      else if (unit->type == SLOW_FILTER)
         t_ten_theo = ssf->fi_t_ten;
      else /*(unit->type == DE_FILTER)*/
         t_ten_theo = def->fi_t_ten;

      phchange(unit, TRUE);
      t_ten = t_ten_theo * theo_res_time;

      // Apply new_cl2decay() in framework of CFSTRs in series
      //--------------------------------------------------------------------------------
      if (eff->cl2cnt > 0)
      {
         ncstrs = ncstr(t_ten_theo);
         rxnhours = mean_theo * theo_res_time / (60.0 * ncstrs);
         for (i = 1; i <= ncstrs; i++)
         {
            new_cl2decay(unit, rxnhours);
            if ((eff->FreeCl2 + eff->NH2Cl) > 0.0)
               eff->hours += rxnhours;
         }

         choose_dbpmodel(unit); /* NEW CODE (WJS, 10/98) */

      } // end if (eff->cl2cnt > 0)

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

      // Calculate CT
      if (eff->cl2cnt > 0.0 || eff->o3_res > 0.0 || eff->clo2_res > 0.0)
         ct(unit, t_ten);

   } //End"if(theo_res_time > 0.0)"

} /* End of filt_dbp() */

void biofilt_rmv(struct UnitProcess *unit)

/* Most of this subroutine is commented-out because of a major change to
   the TOC reduction routine through rapid rate biofilters */

/* For a rapid filter, this subroutine is called by runmodel only if FILTER is part
   of the process train and OZONE has appeared in the process
   train before FILTER and the free+combined chlorine residual
   at the filter influent is less than 0.1 mg/L.  It calculates
   TOC removal across the filter only.  Equations from Summers
   and Sohn based on literature data. Module coded by WJS 11/98.   */

// For a slow sand filter, this subroutine is called by runmodel to remove
// TOC based on temperature according to values given by Robin Collins (WJS 7/05)

{
   /* FILE *fptr2;*/

   /*double     temp; */       /* Temperature in deg. C of filter influent water       */
   /*double   o3_toc; */       /* Cumulative ozone:TOC ratio by the point of entry to
			   the filter                                           */
   /*double     ebct; */       /* Empty bed contact time (= theo. res. time) of filter (min.) */
   /*double  eff_toc; */       /* Holds influent TOC when receiving inputs and effluent
			   TOC at end of subroutine                             */
   /*double  removal = 0.0; */ /* The percent TOC removal across the biofilter         */

   struct Effluent *eff; /* Effluent data structure pointer                      */
   eff = &unit->eff;

   // For a rapid filter...
   if (unit->type == FILTER)
   {
      if (unit->data.filter->media == 'S')
         eff->TOC *= 0.85; //..with A/S media
      else                 /* unit->data.filter->media =='G'*/
         eff->TOC *= 0.80; //..with GAC media
   }

   // For a slow sand filter...
   if (unit->type == SLOW_FILTER)
   {
      if (eff->DegK >= 293.15)
         eff->TOC *= ((100.0 - unit->data.ssf->toc_rem_hight) / 100.0);
      else if (eff->DegK < 293.15 && eff->DegK >= 283.15)
         eff->TOC *= ((100.0 - unit->data.ssf->toc_rem_midt) / 100.0);
      else /* eff->DegK < 283.15 */
         eff->TOC *= ((100.0 - unit->data.ssf->toc_rem_lowt) / 100.0);
   }

   /* Zero disinfectant residuals through the biofilter */
   eff->FreeCl2 = 0.0;
   eff->NH2Cl = 0.0;
   eff->NHCl2 = 0.0;
   eff->cl2dose = 0.0;
   eff->o3_res = 0.0;

   // OLD STUFF...
   /* Get inputs from Unit Process Data Structure */
   /*temp    = eff->DegK - 273.15; */

   /*if(eff->toc_at_o3 > 0.0 ) */ /* self-protect */
                                  /*  {
   o3_toc  = eff->cum_o3dose / eff->toc_at_o3;*/
                                  /* O3:TOC Ratio will be defined as */
   /*  }	  */                     /* the cumulative O3 dosed divided */
   /*else*/ /* toc_at_o3 <= 0 */  /* by the TOC at first point of O3 */
                                  /*  {
   o3_toc  = 1.0;*/
                                  /* Max o3_TOC ratio used to develop equations was 0.77 */
                                  /*  }
if(o3_toc > 1.0) o3_toc = 1.0;

ebct    = (unit->data.filter->volume / eff->Flow) * 1440;   */
                                  /*if( coldflag==TRUE) ebct *= eff->Flow / eff->Peak;  */
                                  /*eff_toc = eff->TOC;  */

   /* Choose from one of four models based on Ozonation Location and Temperature */
   /*if (pre_o3flag == TRUE)
   {
    if (temp >= 15.0)
       {
	removal = 2.46 * ebct * o3_toc + 0.6875;   */
   /* n=8, R^2 = 0.80 */
   /*       }
    else
       {                                                
	removal = 2.27 * ebct - 11.273;            */
   /* n=4, R^2 = 0.63 */
   /*       }
    }
else if (int_o3flag == TRUE || post_o3flag == TRUE)
    {
     if (temp > 15.0)
	{
	 removal = 0.55679 * ebct * o3_toc + 13.58;  */
   /* n=34, R^2 = 0.06!! */
   /*	}
     else
	{
	 removal = 0.6442 * ebct * o3_toc + 2.925;  */
   /* n=6, R^2 = 0.62 */
   /*	}
     }
else     removal = 0.0;*/
   /* DO NOTHING, if no model is selected */;

   /*DEBUGGING CODE*/
   /*  fprintf(fptr2,"Module: %d\n",unit->type); */
   /*   fprintf(fptr2,"Bromide:  %f\n",unit->eff.Br * MW_Br);   */
   /*    fptr2=fopen("debug.dat","a+");
      fprintf(fptr2,"     WE ARE IN BIOFILTRMV()     \n");
      fprintf(fptr2,"removal:       %f\n",        removal);
      fprintf(fptr2,"ebct (min)     %f\n",           ebct);
      fprintf(fptr2,"temp           %f\n",           temp);
      fprintf(fptr2,"cum_o3_dose    %f\n",eff->cum_o3dose);
      fprintf(fptr2,"eff_toc_at o3  %f\n", eff->toc_at_o3);
      fprintf(fptr2,"O3:TOC ratio   %f\n",         o3_toc);
      fclose(fptr2);   */
   /*DEBUGGING CODE*/

   /* Self-protection */
   /*if (removal > 100.0) removal = 100.0; */

   /* Calculate filter effluent TOC */
   /*eff_toc *= (100.0 - removal)/100.0; */

   /* Update Unit Process Data Structure */
   /*eff->TOC = eff_toc;*/

   //eff->last_o3_inf = NULL;
   /* This line assures that ozonate is not called until ozone is
	    added again after this biofilter */

} /* End of biofilt_rmv() */