/* Ozone1.c */
#include "wtp.h"

void ozonate(struct UnitProcess *unit)
/*
* Purpose:
*   Estimate Bromate Formation, Bromide conversion, and Ozone Residual
*   in an ozone contactor. 
*
* Documention and Code by WJS, 98/99
*/
{
  //FILE *fptr;
  double theo_res_time;          /* This is the EBCT of the ozone contactor (calc.)*/
  double mean_theo;              /* This is the ratio of the mean residence time to the
  			   EBCT of the ozone contactor (as entered by user) */
  double t_mean;                 /* This is the mean res. time of the ozone contactor
			   [min.] (calc.) */
  double inf_time;               /* Minutes since last ozone dose at unit process influent */
  double eff_time;               /* Minutes since last ozone dose at unit process effluent */
  double pH;                     /* pH in unit process*/
  double br;                     /* bromide concentration at point of ozonation    (ug/L)  */
  double uv;                     /* UV-254 at point of ozonation (for internal calcs.) (1/cm)  */
  double nh3_n;                  /* Ammonia nitrogen in unit process  (mg/L)  */
  double temp;                   /* Temperature at point of ozonation   (deg C) */
  double alk;                    /* Total alkalinity in unit process (mg/L as CaCO3) */
  double toc;                    /* TOC at point of ozonation (mg/L) */
  double tsuva;                  /* TSUVA at point of ozonation = UV/TOC*100 (l/mg-m) */
  double o3_dose;                /* Ozone Dose               (mg/L)  */
  double o3_demand = 0.0;        /* Ozone demand             (mg/L)  */
  double o3_in;                  /* Ozone residual at process influent (mg/L) */
  double o3_eff;                 /* Ozone residual at process effluent (mg/L) */
  double bro3_in;                /* Bromate at process influent (ppb) */
  double bro3_out;               /* Bromate at process effluent (ppb) */
  double br_in;                  /* Bromide at process influent (ppb) */
  double br_out;                 /* Bromide at process effluent (ppb) */
  double K_o3_decay = 0.0;       /* Ozone decay constant */
  double delta_bro3_wnh3 = 0.0;  /* Bromate formed since last point of ozonation (ug/L)
  			   based on "with nh3" model  */
  double delta_bro3_wonh3 = 0.0; /* Bromate formed since last point of ozonation (ug/L)
			   based on "without nh3 model  */
  double delta_bro3 = 0.0;       /* Final bromate estimate (ug/L) */

  double term1, term2, term3, term4, term5, term6, term7, term8;

  register struct Effluent *eff;
  struct UnitProcess *last_ozone;

  eff = &unit->eff;

  /* Self-protection: there should have been ozone application before an
   ozone contactor, but... */
  if (eff->last_o3_inf != NULL && unit->type == O3_CONTACTOR)
  {

    last_ozone = eff->last_o3_inf; /* last_o3_inf is updated as the
						address of the last point of
                                                ozone application when an OZONE
						unit process is encountered in the
						main runmodel() loop */

    theo_res_time = (unit->data.basin->volume) / (eff->Flow / 1440);
    mean_theo = unit->data.basin->sb_mean;

    /* These model inputs are based on WQ parameters at point of ozonation */
    br = last_ozone->eff.Br * MW_Br * 1000;
    uv = PrevUnitProcess(last_ozone)->eff.UV; //because UV already oxidized in o3add()
    temp = last_ozone->eff.DegK - 273.15;
    o3_dose = last_ozone->eff.o3_res; /* eff.o3_res is the sum of the residual and the dose */
    toc = last_ozone->eff.TOC;
    //  pH		= last_ozone->eff.pH;
    //  nh3_n	= last_ozone->eff.NH3 * MW_NH3;
    //  alk		= last_ozone->eff.Alk * MW_CaCO3 / 2.0;

    /* These inputs can change through chemical addition between ozone points,
     and thus, should be based on the WQ in the chamber.  Or, they are
     quantities that need to be known at the influent to the unit process.*/
    pH = eff->pH;
    alk = eff->Alk * MW_CaCO3 / 2.0;
    nh3_n = eff->NH3 * MW_NH3;
    o3_in = eff->o3_res;
    bro3_in = eff->BrO3;
    br_in = eff->Br * MW_Br * 1000;

    /* Set Limits for parameters */
    if (uv < 0.001)
      uv = 0.001;
    if (alk < 1.0)
      alk = 1.0;
    if (nh3_n < 0.01)
      nh3_n = 0.01;
    if (toc < 0.1)
      toc = 0.1;

    /*Calculate TSUVA */
    tsuva = uv / toc * 100.0;

    /* Bromate and ozone residual's pH sensitivity --> need for a cut-off on this one */
    if (pH < 5.8)
      pH = 5.8; /* 5.8 is the minimum pH in the pilot database used to
			    develop ozone decay algorithm and correction factors for bromate
			    algorithm */
    if (pH > 8.8)
      pH = 8.8; /* 8.8 is the maximum pH in the pilot database used to
			    develop ozone decay algorithm and correction factors for bromate
			    algorithm */

    /* Calculate mean res. time as ratio of mean to theoretical times theoretical*/
    t_mean = mean_theo * theo_res_time;

    /* Time expressions */
    inf_time = eff->o3_minutes;
    eff_time = eff->o3_minutes + t_mean;
    if (inf_time > 120.0)
      inf_time = 120.0;
    if (eff_time > 120.0)
      eff_time = 120.0;

    /*********************** Begin O3 residual predictions ****************************/

    /* Estimate ozone residual at process effluent by subtracting predicted demand
     (only if time has elapsed) */

    if (o3_dose > 0.0 && eff_time > inf_time)
    {
      term1 = pow(o3_dose, 1.312);
      term2 = pow((o3_dose / uv), -0.386);
      term3 = pow(tsuva, -0.184);
      term4 = pow(alk, 0.023);
      term5 = pow(pH, 0.229);
      term6 = pow(temp, 0.087);
      term7 = pow(eff_time, 0.068) - pow(inf_time, 0.068);
      o3_demand = 0.996 * term1 * term2 * term3 * term4 * term5 * term6 * term7;
      K_o3_decay = o3_demand / term7; //to support CT calc.
    }

    o3_eff = o3_in - o3_demand;

    /*  
      fptr = fopen("debug.dat","a+");
      fprintf(fptr,"Module:            %d\n",unit->type);
      fprintf(fptr,"theo_res_time      %f\n",theo_res_time);
      fprintf(fptr,"t_mean             %f\n",t_mean);
      fprintf(fptr,"inf_time               %f\n",inf_time);
      fprintf(fptr,"eff_time               %f\n",eff_time);
      fprintf(fptr,"o3_dose            %f\n",o3_dose);
      fprintf(fptr,"inf_o3            %f\n",o3_in);
      fprintf(fptr,"eff_o3            %f\n\n",o3_eff);
      fclose(fptr);
*/

    //Make sure predictions are reasonable
    if (o3_eff > o3_dose)
      o3_eff = o3_dose;
    if (o3_eff > o3_in)
      o3_eff = o3_in;

    //Set limits for the amount of decay time and min. residual
    if (o3_eff < 0.10 || eff_time >= 120.0)
      o3_eff = 0.0;

    //Update UnitProcess data structure
    eff->o3_res = o3_eff;

    /************************ End O3 Residual predictions ********************************/

    /*************************BEGIN Bromate Calculations *********************************/

    /* Estimate effluent bromate concentration from this process based upon
     bromate and WQ at last point of ozonation and in unit process.  Model for
     bromate from Jinsik Sohn, Univ. of Colorado, 1998.  Model for ozone
     demand from G. Solarik, 2001 */

    if (br > 0.0 && t_mean > 0.0 && o3_dose > 0.0 && o3_eff > 0.0)
    { /* Bromate can form */

      /* First, calculate the estimate based on the model with ammonia */
      term1 = pow(uv, -0.593);
      term2 = pow(pH, 5.81);
      term3 = pow(o3_dose, 1.28);
      term4 = pow(br, 0.94);
      term5 = pow(alk, -0.167);
      term6 = pow(1.035, (temp - 20));
      term7 = pow(nh3_n, -0.051);
      term8 = pow(eff_time, 0.337) - pow(inf_time, 0.337);
      delta_bro3_wnh3 = 8.71 * pow(10, -8.0) * term1 * term2 * term3 * term4 * term5 * term6 * term7 * term8;

      /* Next, calculate the estimate based on the model without ammonia */
      term1 = pow(uv, -0.623);
      term2 = pow(pH, 5.68);
      term3 = pow(o3_dose, 1.31);
      term4 = pow(br, 0.96);
      term5 = pow(alk, -0.201);
      term6 = pow(1.035, (temp - 20));
      term7 = pow(eff_time, 0.336) - pow(inf_time, 0.336);
      delta_bro3_wonh3 = 1.19 * pow(10, -7.0) * term1 * term2 * term3 * term4 * term5 * term6 * term7;

      if (nh3_n > 0.1)
      { /* Use "Model with ammonia" only if it gives a lower bromate result than the
	     "Model without ammonia" */
        if (delta_bro3_wnh3 < delta_bro3_wonh3)
          delta_bro3 = delta_bro3_wnh3;
        else
          delta_bro3 = delta_bro3_wonh3;
      }
      else
      { /* Use "Model without ammonia" */
        delta_bro3 = delta_bro3_wonh3;
      } /* End "if/else (nh3_n > 0.1)" */

      //   Apply bromate correction factor based on Spyros' and Zaid's pilot data analysis
      delta_bro3 /= 2.52; /* R^2=0.80 */

      //Calculate first estimate of effluent bromate concentration
      bro3_out = bro3_in + delta_bro3;

      // Set stoichiometric max. based on bromide available to convert
      if ((bro3_out * 0.625) > br)
        bro3_out = (br / 0.625);

      //Make sure result is reasonable
      if (bro3_out < bro3_in)
        bro3_out = bro3_in;

    } /* End "if(br > 0.0 && time > 0.0 && o3_dose > 0.0)" */

    //Update UnitProcess data structure
    eff->BrO3 = bro3_out;

    /* Estimate effluent bromide concentration in process effluent and enforce limits*/
    br_out = br_in - 0.625 * delta_bro3;
    if (br_out < 0.0)
      br_out = 0.0;
    eff->Br = br_out / 1000.0 / MW_Br; /* Convert from ug/L to Moles/L */

    //***************************************End Bromate Calcs*********************************************

    /* Update Unit Process Data Structure */
    eff->K_o3decay = K_o3_decay;
    eff->o3_minutes = eff_time;

  } /* End of "if(last_o3_inf != NULL && unit->type == O3_CONTACTOR)" */

} //End of subroutine

/*DEBUGGING CODE*/
/*      
      fptr = fopen("debug.dat","a+");
      fprintf(fptr,"Module:            %d\n",unit->type);
      fprintf(fptr,"theo_res_time      %f\n",theo_res_time);
      fprintf(fptr,"t_mean             %f\n",t_mean);
      if(eff->wtp_effluent != NULL)
      {
      fprintf(fptr,"wtpeffo3minutes    %f\n",eff->wtp_effluent->eff.o3_minutes);
      }
      else
      {
      fprintf(fptr,"wtpeffo3minutes    NULL\n");
      }
      fprintf(fptr,"time               %f\n",time);
      fprintf(fptr,"o3_res (in ozonate)%f\n\n",eff->o3_res);

      fclose(fptr);
*/
/*DEBUGGING CODE*/

//      fptr = fopen("ozone.dat","a+");
//      fprintf(fptr,"made it!             \n");
//      fclose(fptr);