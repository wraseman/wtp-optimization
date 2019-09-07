/* Influent.c */
#include "wtp.h"

void influent(struct UnitProcess *unit)
/*
*  Purpose:
*    1. Copy Influent data to Effluent data structure.
*    2. Compute temperature dependent ionization & solubility coefficients.
*    3. Determine log_inactivation required for Giardia, Virus, and Cryptosporidium.
*/
{
  //FILE *fptr2;
  double kw;     /* Ionization coefficient of water                */
  double k1, k2; /* Ionization coefficients of carbonic acid.      */
  double H;      /* [H+]                              (Mole/Liter) */
  double OH;     /* [OH-]                             (Mole/Liter) */
  double denom, alpha_one, alpha_two;
  int bin34flag = FALSE; // Source water falls into LT2 bin 3 or 4

  double log_required_g;
  double log_required_v;
  double log_required_c;

  register struct Effluent *eff;
  register struct Influent *inf;

  if (unit == NULL || unit->type != INFLUENT)
    return;

  /* Move the pointers into registers. */
  inf = unit->data.influent;
  eff = &unit->eff;

  /* Copy input variables to effluent variables and convert units.   */
  if (coldflag == TRUE)
  {
    eff->DegK = inf->low_temp + 273.15;
    eff->Flow = inf->peak_flow; //Still store this in eff->Peak below for use
                                //in dist_dbp(), gac_rmv_dbp()
  }
  else
  {
    eff->DegK = inf->temp + 273.15;
    eff->Flow = inf->avg_flow;
  }

  eff->pH = inf->pH;
  eff->TOC = inf->toc;
  eff->UV = inf->uv254;
  eff->UV_out = inf->uv254;
  eff->Br = inf->bromide / MW_Br;
  eff->Br_at_last_cl2 = eff->Br;
  eff->Alk = inf->alkalinity / (MW_CaCO3 / 2);
  eff->NH3 = inf->nh3 / MW_NH3;
  eff->Turbidity = inf->ntu;

  eff->CO2_aq = 0.0;
  eff->Ca_aq = inf->calcium / MW_CaCO3;
  eff->Ca_solid = 0.0;
  eff->Ca_sludge = 0.0;
  eff->Mg_aq = (inf->hardness - inf->calcium) / MW_CaCO3;
  eff->Mg_solid = 0.0;
  eff->Mg_sludge = 0.0;

  eff->Peak = inf->peak_flow;

  /* Make sure swflag gets set somewhere */
  swflag = inf->swflag;

  /* Zero the following variables in the Effluent Data Structure */
  if (eff->TOC > 0.0)
    eff->SUVA = eff->UV_out / eff->TOC * 100.0;
  else
    eff->SUVA = 999999.9;

  eff->FreeCl2 = 0.0;
  eff->NH2Cl = 0.0;
  eff->NHCl2 = 0.0;
  eff->solids = 0.0;
  eff->cl2dose = 0.0;

  eff->CHCl3 = 0.0;
  eff->CHBrCl2 = 0.0;
  eff->CHBr2Cl = 0.0;
  eff->CHBr3 = 0.0;
  eff->TTHM = 0.0;

  eff->MCAA = 0.0;
  eff->DCAA = 0.0;
  eff->TCAA = 0.0;
  eff->MBAA = 0.0;
  eff->DBAA = 0.0;
  eff->BCAA = 0.0;
  eff->BDCAA = 0.0;
  eff->DBCAA = 0.0;
  eff->TBAA = 0.0;
  eff->HAA5 = 0.0;
  eff->HAA6 = 0.0;
  eff->HAA9 = 0.0;

  eff->ct_ratio = 0.0;
  eff->ct_ratio_v = 0.0;
  eff->ct_ratio_c = 0.0;

  eff->ct_cl2 = 0.0;
  eff->ct_nh2cl = 0.0;
  eff->ct_clo2 = 0.0;
  eff->ct_o3 = 0.0;

  eff->log_inact_ach_g = 0.0;
  eff->log_inact_ach_v = 0.0;
  eff->log_inact_ach_c = 0.0;

  eff->AlumDose = 0.0;
  eff->FericDose = 0.0;
  eff->LimeDose = 0.0;
  eff->floc_carry_cntr = FALSE;
  eff->EquivAlumDose = 0.0;
  eff->toc_sludge = 0.0;
  eff->alum_to_sludge = 0.0;
  eff->iron_to_sludge = 0.0;
  eff->wtp_effluent = NULL;
  eff->influent = unit;

  eff->cl2cnt = 0;
  eff->hours = 0.0;

  /* New Effluent data structure variables added by WJS, 10/98 that
     require initializing */
  eff->BrO3 = 0.0;
  eff->o3_res = 0.0;
  eff->K_o3decay = 0.0;
  eff->first_o3_chamber = FALSE;
  eff->o3_minutes = 0.0;
  eff->last_o3_inf = NULL;
  eff->last_rm_inf = NULL;
  eff->dbpmodel = "Raw Water Model";

  eff->processtime = 0.0;
  eff->traintime = 0.0;

  eff->clo2_res = 0.0;
  eff->clo2_minutes = 0.0;
  eff->clo2_avg = 0.0;
  eff->chlorite = 0.0;
  eff->clo2dose = 0.0;
  eff->clo2_cntr = 0;
  eff->clo2_init_decay = FALSE;

  eff->tox = 0.0;
  eff->limesoftening = FALSE;

  eff->ec_meeting_step1 = FALSE;
  eff->ec_exempt = FALSE;

  eff->pre_re_chlor_flag = FALSE;
  eff->pre_chlor_flag = FALSE;
  eff->pre_chlor_TTHM = 0.0;
  eff->pre_chlor_HAA6 = 0.0;
  eff->pre_chlor_dose_track = 0.0;
  eff->pre_chlor_dose_ratio = 0.0;

  /***********************************pH Stuff**********************************/

  /* PChem Coefficients: */
  kw = Kw(eff->DegK);
  k1 = K_HCO3(eff->DegK);
  k2 = K_CO3(eff->DegK);

  H = pow(10.0, -(eff->pH));
  OH = kw / H;

#if 1 == 0
  printf("Influent():\n");
  printf("DegC     :%f\n", eff->DegK - 273.15);
  printf("Kw       :%e\n", kw);
  printf("K1       :%e\n", k1);
  printf("K2       :%e\n", k2);
#endif

  /*
  * Calculate dissolved carbonate concentration [HCO3-]+[CO3--] via eq A-37
  * (mole/L).  Note: this does not include the [H2CO3] term.
  */
  denom = (H * H) + (k1 * H) + (k1 * k2);
  alpha_one = (k1 * H) / denom;
  alpha_two = (k1 * k2) / denom;
  eff->CO2_aq = (eff->Alk + H - OH) / (alpha_one + 2 * alpha_two);

  /*
  *  Call phchange() to set CBminusCA.
  */
  eff->CBminusCA = 0.0;
  phchange(unit, FALSE);

  /****** Calculate log_inactivation required for Giardia, Virus, and Crypto.********/

  if (swflag == TRUE)
  {
    /* From SWTR, for filtered or unfiltered systems, we have total disinfection
         requirements for Giardia and Viruses: */
    log_required_g = 3.0;
    log_required_v = 4.0;

    if (filtflag == TRUE)
    { //For filtered systems....

      //First, set total disinfection requirements for Crypto.
      if (dir_filtflag == TRUE) //Total Disinfection Requirements for LT2 for DIRECT FILTER
      {
        if (inf->crypto_conc < 0.075)
          log_required_c = 3.0; //Bin #1
        else if (inf->crypto_conc >= 0.075 &&
                 inf->crypto_conc < 1.00)
          log_required_c = 4.5; //Bin #2
        else if (inf->crypto_conc >= 1.00 &&
                 inf->crypto_conc < 3.00)
        {
          log_required_c = 5.5;
          bin34flag = TRUE;
        } //Bin #3
        else
        {
          log_required_c = 6.0;
          bin34flag = TRUE;
        } //Bin #4
      }
      else
      { //Total Disinfection Requirements for LT2 for all other filter techs.
        if (inf->crypto_conc < 0.075)
          log_required_c = 3.0; //Bin #1
        else if (inf->crypto_conc >= 0.075 &&
                 inf->crypto_conc < 1.00)
          log_required_c = 4.0; //Bin #2
        else if (inf->crypto_conc >= 1.00 &&
                 inf->crypto_conc < 3.00)
        {
          log_required_c = 5.0;
          bin34flag = TRUE;
        } //Bin #3
        else
        {
          log_required_c = 5.5;
          bin34flag = TRUE;
        } //Bin #4
      }   //end if-else (dir_filtflag == TRUE)

      /*Set these globals to use in echoing to run_wtp() output*/
      tot_dis_req_g = log_required_g;
      tot_dis_req_v = log_required_v;
      tot_dis_req_c = log_required_c;

      //Next, begin turning the "log_required_x" variables into log inactivations required

      //First, check on LT2 toolbox source water protection credit for Crypto.
      if (inf->lt2_wscp_flag == TRUE)
      {
        log_required_c -= 0.5;
        lt2_wscp_flag = TRUE;
      }

      //Now, apply the total log removals from treatment (includes, 1st & 2nd stage filtration,
      //bank filtration, optimized filter performance, UV, 2nd stage softening and pre-sed); these
      //were calculated at start of runmodel()
      log_required_g -= tot_giardia_lr;
      log_required_v -= tot_virus_lr;
      log_required_c -= tot_crypto_lr;

      //Need to make sure that if we are in Bins 3 and 4 that we will get enough CT credit
      // to fulfill the 1.0-log minimum credit via ozone, ClO2, membranes, bags/cartridges
      // bank filters or UV disinfection requirement
      if (bin34flag == TRUE && (log_required_c + nonconv_discredit) < 1.0)
      {
        log_required_c = 1.0 - nonconv_discredit;
        bin34_inactreqd = 1.0 - nonconv_discredit;
      }
    }
    else
    { //For unfiltered SW systems...

      //Crypto Inactivation Requirement for Unfiltered System = f(source Crypto)
      if (inf->crypto_conc <= 0.01)
        log_required_c = 2.0;
      else
        log_required_c = 3.0;

    } //end if-else (filtflag==TRUE)

    /* Then, check to see that log_required is >= 0*/
    if (log_required_g <= 0.0)
    {
      log_required_g = 0.0;
      eff->ct_ratio = 1.0;
    }
    if (log_required_v <= 0.0)
    {
      log_required_v = 0.0;
      eff->ct_ratio_v = 1.0;
    }
    if (log_required_c <= 0.0)
    {
      log_required_c = 0.0;
      eff->ct_ratio_c = 1.0;
    }

    /* Lastly, initialize UnitProcess data structure*/
    eff->log_required_g = log_required_g;
    eff->log_required_v = log_required_v;
    eff->log_required_c = log_required_c;
  }
  else /* Groundwater - no disinfection as yet required for cysts, but need to check
	      if based on inputs, disinfection is required for viruses, per proposed
	      GW Rule */
  {

    /* Six nines serves as indicator that no CT is required due to GW */
    eff->log_required_g = 999999.0;
    eff->log_required_c = 999999.0;

    /*Set these globals to use in echoing to run_wtp() output*/
    tot_dis_req_g = 0.0;
    tot_dis_req_c = 0.0;
    tot_dis_req_v = 0.0;

    if (inf->gw_virus_flag == TRUE)
    { //This is a groundwater needing virus disinfection
      gw_virus_flag = TRUE;
      log_required_v = inf->gw_tot_dis_req_v; //Sets total amount of disinfection from input value
      tot_dis_req_v = log_required_v;         //Saves this amount for use in output tables
      log_required_v -= tot_virus_lr;         //Now takes out removal to calc. inact. req'd

      if (log_required_v <= 0.0) //Check to see if we have more removal credit than
      {                          //  the total disinfection required (unlikely)
        log_required_v = 0.0;
        eff->ct_ratio_v = 1.0;
      }
      eff->log_required_v = log_required_v; //Saves the amount of inact. req'd for CT calcs.
    }                                       // end if(inf->gw_virus_flag == TRUE)

  } // if-else (swflag == TRUE)

} //End "influent()"

//Old stuff
#if 1 == 0
printf("CBminusCA:%e\n", eff->CBminusCA);
printf("Alk:     :%e\n", eff->Alk);
printf("CO2_aq   :%e\n", eff->CO2_aq);
#endif

/* Orginal 1.2 code to set CBminusCA. */
#if 1 == 0

ca_term = eff->Ca_aq * (2 + (k_caoh / H)) / (1 + (k_caoh / H) + (k_caoh2aq / (H * H)));
mg_term = eff->Mg_aq * (2 + (k_mgoh / H)) / (1 + (k_mgoh / H) + (k_mgoh2aq / (H * H)));
nh3_term = eff->NH3 / (1 + (k_nh3 / H));
ocl_term = eff->FreeCl2 / (1 + (H / k_hocl));
eff->CBminusCA = eff->Alk - ca_term - mg_term - nh3_term + ocl_term;
#endif
